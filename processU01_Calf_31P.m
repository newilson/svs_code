function processU01_Calf_31P(SessionDir, FinalOutDir, plt, extRef)
% PROCESSU01_CALF_31P  Absolute 31P quantification for one U01 calf session.
%
%   processU01_Calf_31P(SessionDir, FinalOutDir)
%   processU01_Calf_31P(SessionDir, FinalOutDir, plt)
%   processU01_Calf_31P(SessionDir, FinalOutDir, plt, extRef)
%
% Reads the highest-MID *31p*.dat from <SessionDir>\datfiles\, runs the
% basis-Varpro fit from demo_dat2mat_svs31p.m via dat2mat_svs31p.m, then
% converts each fitted metabolite area to absolute [mM] using PCr as an
% internal reference with per-metabolite T1/T2 relaxation correction:
%
%   [M] = (S_M / S_PCr) * (R_PCr / R_M) * (nP_PCr / nP_M) * [PCr]
%
%   R_x = (1 - exp(-TR/T1_x)) * exp(-TE/T2_x)
%
% Calf is a single muscle compartment, so (unlike the brain pipeline) there
% is no anatomic segmentation and no (1 - f_CSF) partial-volume term.
%
% [PCr] reference: 33 mM is the standard whole-muscle (wet-weight) value for
% skeletal-muscle 31P MRS (Kemp; Meyerspeer et al. consensus, NMR Biomed
% 2020).  To report per litre of intracellular water instead, use ~42.5 mM
% ( = 33 / 0.77 ) -- edit refConc_mM below.
%
% OPTIONAL external-reference (phantom) quant, mirroring processU01_Brain_31P:
% pass extRef.enable=true to add the ~8 ppm phantom singlet 'RefExt' to the
% basis and additionally quantify every metabolite against the phantom of
% known concentration.  The calf tissue segmentation is the leg (largest
% bright connected component in the B1 map) -- see extRefQuant31P, which also
% documents the fat/marrow caveat on V_M.  B1+ correction is DISABLED
% (extRef.applyB1=false); a missing B1 series is NOT fatal, it just leaves the
% external concentrations unset while the PCr-referenced quant proceeds.
%
% Output: <FinalOutDir>\<scanID>_calf_31P_abs.mat
%   (output, parsOut, hdr, absQ, pars, metProps, extRef)

arguments
    SessionDir  (1,:) char
    FinalOutDir (1,:) char
    plt         (1,1) logical = false
    extRef      struct = struct()   % OPTIONAL external-reference quant (off unless .enable=true)
end

% External-reference (phantom) absolute quant is entirely optional and OFF by
% default: with extRef.enable=false nothing below changes and the internal
% PCr reference path runs exactly as before.
extRef    = extRefDefaultsCalf(extRef, SessionDir);
useExtRef = extRef.enable;
if useExtRef, extRef.plt = extRef.plt || plt; end

assert(isfolder(SessionDir), 'SessionDir not found: %s', SessionDir);
if ~isfolder(FinalOutDir), mkdir(FinalOutDir); end

datDir = fullfile(SessionDir, 'datfiles');
assert(isfolder(datDir), 'datfiles\\ subfolder not found in %s', SessionDir);

datFile = pick_latest_31p(datDir);
fprintf('31P : %s\n', datFile);

% =====================================================================
% Processing/Fitting parameters (from demo_dat2mat_svs31p.m)
% =====================================================================
pars = struct();
pars.lb         = 3;
pars.blp.flag   = true;
pars.PC         = 0;
pars.plt        = plt;

pars.hsvd.flag  = false;
pars.hsvd.sv    = 30;

% peak-based phase + frequency correction: PCr / gATP / aATP
pars.peaks      = [0, -2.53, -7.56];
pars.peakRanges = [0.5 0.5;
                   0.4 0.4;
                   0.3 0.3];

pars.dofit            = true;
pars.fit.mode         = 'basisvarpro';
pars.block_align.size = [];

pars.fit.ppm_range    = [-19.5 19.5];
pars.fit.metabs       = { ...
    'PE','PC','DPG1','DPG2','Pi','GPE','GPC','BroadMP', ...
    'PCr','gATP','aATP','bATP', ...
    'NADH','NADplus', ...
    'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc', ...
    'Broad10_5'};
pars.fit.phase1Enable = true;

% Add the external-reference singlet to the basis subset only when enabled.
% Optionally add a broad background component (BroadRefBG) to absorb the
% positive-side envelope so RefExt stays a clean narrow peak.
if useExtRef
    pars.fit.metabs{end+1} = 'RefExt';
    pars.fit.refExtPpm     = extRef.refPpm;   % ~8 ppm for the D2O/PBS phantom
    if extRef.addBroadBG
        pars.fit.metabs{end+1} = 'BroadRefBG';
    end
end

nMet    = numel(pars.fit.metabs);
isUDP   = startsWith(pars.fit.metabs,'UDP');
isNAD   = ismember(pars.fit.metabs,{'NADH','NADplus'});
isBroad = startsWith(pars.fit.metabs,'Broad');

shB = repmat([-60 60], nMet, 1);
shB(isUDP,:)   = repmat([-6   6],   sum(isUDP),   1);
shB(isNAD,:)   = repmat([-12  12],  sum(isNAD),   1);
shB(strcmp(pars.fit.metabs,'Broad10_5'),:) = [-300 60];

lbB = repmat([0 100], nMet, 1);
lbB(isUDP,:)   = repmat([15 50],  sum(isUDP),   1);
lbB(isNAD,:)   = repmat([15 50],  sum(isNAD),   1);
lbB(isBroad,:) = repmat([20 200], sum(isBroad), 1);

lbInit = 5 * ones(nMet,1);
lbInit(isUDP)   = 25;
lbInit(isNAD)   = 25;
lbInit(isBroad) = 100;

% External-reference singlet bound overrides (kept near +8 ppm; phantom line
% may be broad, so allow a wider linewidth band).
if useExtRef
    ri = strcmp(pars.fit.metabs,'RefExt');
    shB(ri,:)  = extRef.refShiftBoundsHz;
    lbB(ri,:)  = extRef.refLbBoundsHz;
    lbInit(ri) = extRef.refLbInitHz;
end

pars.fit.fitOpt.shiftBounds = shB;
pars.fit.fitOpt.lbLBounds   = lbB;
pars.fit.fitOpt.lbLInit     = lbInit;
pars.fit.fitOpt.phaseBounds = [-20 20];

sc(1) = struct('idxA','UDPGal',    'idxB','UDPGlcNAc','ratio',0.53,'lambda',1e4);
sc(2) = struct('idxA','UDPGlc',    'idxB','UDPGlcNAc','ratio',1.00,'lambda',1e4);
sc(3) = struct('idxA','UDPGalNAc', 'idxB','UDPGlcNAc','ratio',0.80,'lambda',1e4);
pars.fit.softConstraints = sc;

pars.fit.tieGroups   = [];
pars.fit.lineShapeOpt = struct('enable',false,'nSide',3,'asymmetric',true,'maxSide',0.3);

% =====================================================================
% Per-metabolite T1 / T2 / number of equivalent 31P nuclei.
% PLACEHOLDER values at 7T compiled from Lu 2014, Bogner 2009, Ren 2015,
% Mariano 2014. EDIT BEFORE PUBLICATION USE -- and note these are brain
% relaxation values; muscle PCr/ATP/Pi relaxation differs.
% =====================================================================
metProps = struct( ...
    'name',  {'PE','PC','DPG1','DPG2','Pi','GPE','GPC','BroadMP', ...
              'PCr','gATP','aATP','bATP', ...
              'NADH','NADplus', ...
              'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc', ...
              'Broad10_5'}, ...
    'T1_ms', {6000, 6000, 3000, 3000, 3500, 6000, 6000, 3000, ...
              3800, 1200, 1200, 1500, ...
               400,  400, ...
               500,  500,  500,  500, ...
              3000}, ...
    'T2_ms', {200, 200, 100, 100, 150, 200, 200, 150, ...
              180,  40,  40,  30, ...
               10,  10, ...
               10,  10,  10,  10, ...
              50}, ...
    'nPhos', {1,1,1,1,1,1,1,1,  1,1,1,1,  1,1,  1,1,1,1,  1});

% External-reference relaxation/nuclei properties (placeholders; edit when
% the phantom compound's T1/T2 are measured).
if useExtRef
    metProps(end+1) = struct('name','RefExt', ...
        'T1_ms', extRef.T1_ms, 'T2_ms', extRef.T2_ms, 'nPhos', extRef.nP);
    if extRef.addBroadBG
        metProps(end+1) = struct('name','BroadRefBG','T1_ms',3000,'T2_ms',50,'nPhos',1);
    end
end

% =====================================================================
% Reference metabolites.  Every scan is reported against BOTH internal
% references (and, when enabled, the external phantom).  Changing reference
% rescales all concentrations by one common factor and never alters a ratio
% between metabolites within a scan -- the three sets exist to show how much
% the absolute scale depends on which reference you believe.
%   PRIMARY (the one written to absQ.metabolites) is PCr for calf.
% =====================================================================
refName    = 'PCr';        % primary reference for calf
refConc_mM = 33;           % [PCr] muscle, whole-muscle wet weight (placeholder)
refIdx     = find(strcmp({metProps.name}, refName), 1);
assert(~isempty(refIdx), 'Reference %s missing from metProps', refName);

% Secondary internal reference, reported alongside.
% [aATP]_muscle PLACEHOLDER -- whole-muscle ATP is conventionally 8.2 mM
% (Kemp; Meyerspeer et al. consensus, NMR Biomed 2020), i.e. PCr/ATP ~ 4.
% EDIT BEFORE PUBLICATION USE, exactly like the T1/T2 table above.
altRefName    = 'aATP';
altRefConc_mM = 8.2;

% =====================================================================
% Process / fit
% =====================================================================
[output, parsOut, hdr] = dat2mat_svs31p(pars, datFile);

% dat2mat_svs31p rmpaths utils at exit; re-add for any utils helpers
% used downstream.
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

% sequence timing (us -> ms)
TR_ms = hdr.MeasYaps.alTR{1} / 1000;
TE_ms = hdr.MeasYaps.alTE{1} / 1000;
fprintf('TR=%.0f ms  TE=%.2f ms\n', TR_ms, TE_ms);

% =====================================================================
% Absolute quant: PCr internal reference (single muscle compartment)
% =====================================================================
% dat2mat_svs31p stores results in a different field per fit mode:
%   basisvarpro -> output.basisVarproResults.{names, amplScaled}
%                  (amplScaled = ampl ./ basis nVec, so it's already
%                   per-spin-normalized -> keep metProps.nPhos = 1)
%   amares      -> output.fit.{names, areas}
if isfield(output,'basisVarproResults')
    fitNames = output.basisVarproResults.names;
    fitAreas = output.basisVarproResults.amplScaled;
elseif isfield(output,'fit')
    fitNames = output.fit.names;
    fitAreas = output.fit.areas;
else
    error('processU01_Calf_31P: no fit results found in dat2mat_svs31p output');
end

absQ = struct();
absQ.TR_ms     = TR_ms;
absQ.TE_ms     = TE_ms;

% Quantify against BOTH internal references.
absQ.byRef.(refName)    = refQuant31P(refName,    refConc_mM, ...
                                      fitNames, fitAreas, metProps, TR_ms, TE_ms);
absQ.byRef.(altRefName) = refQuant31P(altRefName, altRefConc_mM, ...
                                      fitNames, fitAreas, metProps, TR_ms, TE_ms);

% Primary reference (PCr for calf) stays in absQ.reference/absQ.metabolites
% so every existing consumer of this .mat keeps working unchanged.
absQ.reference   = absQ.byRef.(refName).reference;
absQ.metabolites = absQ.byRef.(refName).metabolites;

fprintf('\n--- 31P concentrations [mM] by reference ---\n');
printByRef(absQ.byRef, {refName, altRefName});

[~, scanID] = fileparts(SessionDir);

% =====================================================================
% Optional external-reference (phantom) absolute quant, for comparison
% against the internal PCr reference above.  Phantom/leg volumes come from
% the B1 map slab (same orientation as the ISIS slab).  B1+ correction is
% disabled; a missing or unreadable B1 series leaves absQ.external.ok=false
% and does NOT stop the pipeline.
% =====================================================================
if useExtRef
    % Prescribed excitation flip angle (deg).  Recorded for QC even though the
    % B1+ efficiency term is currently switched off in extRefQuant31P.
    if isempty(extRef.flipDeg)
        try
            extRef.flipDeg = hdr.MeasYaps.adFlipAngleDegree{1};
        catch
            extRef.flipDeg = 90;
            warning('extRef: could not read prescribed flip from hdr; assuming 90 deg.');
        end
    end
    fprintf('extRef: prescribed excitation flip = %.1f deg\n', extRef.flipDeg);
    try
        absQ.external = extRefQuant31P(fitNames, fitAreas, metProps, ...
            TR_ms, TE_ms, extRef, absQ, FinalOutDir, scanID);
    catch ME
        warning('External-reference quant failed: %s', ME.message);
        absQ.external = struct('ok', false, 'error', ME.message);
    end
end

% =====================================================================
% Save
% =====================================================================
outPath = fullfile(FinalOutDir, sprintf('%s_calf_31P_abs.mat', scanID));
save(outPath, 'output','parsOut','hdr','absQ','pars','metProps','extRef');
fprintf('Saved -> %s\n', outPath);

% remove path
rmpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

end

% =====================================================================
function printByRef(byRef, order)
% Side-by-side table of the internal-reference concentration sets.
names = {};
for i = 1:numel(order)
    s = byRef.(order{i});
    if s.ok, names = union(names, {s.metabolites.name}, 'stable'); end
end
if isempty(names), fprintf('  (no reference set available)\n'); return; end

fprintf('  %-12s', 'metabolite');
for i = 1:numel(order)
    s = byRef.(order{i});
    fprintf(' %14s', sprintf('%s=%gmM', order{i}, s.reference.conc_mM));
end
fprintf('\n');
for k = 1:numel(names)
    fprintf('  %-12s', names{k});
    for i = 1:numel(order)
        s = byRef.(order{i});
        v = NaN;
        if s.ok
            j = find(strcmp({s.metabolites.name}, names{k}), 1);
            if ~isempty(j), v = s.metabolites(j).conc_mM; end
        end
        fprintf(' %14.4g', v);
    end
    fprintf('\n');
end
end

% =====================================================================
function R = relax(TR_ms, TE_ms, T1_ms, T2_ms)
R = (1 - exp(-TR_ms / T1_ms)) * exp(-TE_ms / T2_ms);
end

% =====================================================================
function [area, T1_ms, T2_ms, nPhos] = lookup_fit(name, fitNames, fitAreas, metProps)
iF = find(strcmp(fitNames, name), 1);
iP = find(strcmp({metProps.name}, name), 1);
assert(~isempty(iF), 'Metabolite %s not in fit names', name);
assert(~isempty(iP), 'Metabolite %s not in metProps',  name);
area  = fitAreas(iF);
T1_ms = metProps(iP).T1_ms;
T2_ms = metProps(iP).T2_ms;
nPhos = metProps(iP).nPhos;
end

% =====================================================================
function datFile = pick_latest_31p(datDir)
% Highest-MID *31p*.dat in datDir.
allF = dir(fullfile(datDir, '*31p*.dat'));
assert(~isempty(allF), 'No *31p*.dat found in %s', datDir);

mids = nan(1, numel(allF));
for k = 1:numel(allF)
    tok = regexp(allF(k).name, 'MID(\d+)', 'tokens', 'once');
    if ~isempty(tok), mids(k) = str2double(tok{1}); end
end
assert(any(~isnan(mids)), 'No MID# in any *31p*.dat filename in %s', datDir);

[~, iMax] = max(mids);
datFile = fullfile(datDir, allF(iMax).name);
end

% =====================================================================
function e = extRefDefaultsCalf(e, SessionDir)
% Defaults for the optional external-reference quant, calf variant.  Mirrors
% extRefDefaults in processU01_Brain_31P but with anatomy='calf' and none of
% the brain-only fields (brainMaskNii / t1Nii / isisDcm), which only feed the
% HD-BET INT ISIS-slab volume path that has no calf analogue.
if ~isstruct(e), e = struct(); end
d = struct( ...
    'enable',           false, ...
    'anatomy',          'calf', ... % selects the leg segmentation in extRefQuant31P
    'applyB1',          false, ...  % B1+ / receive correction DISABLED (QC-reported only)
    'b1Dir',            '', ...     % auto-located under SessionDir\E* if empty; '' is OK
    'refPpm',           8.0, ...    % external-reference chemical shift (PCr=0)
    'refConc_mM',       500, ...    % known phantom concentration -- CHECK PER SESSION
    'nP',               1, ...      % equivalent 31P per molecule (sodium phosphate = 1)
    'flipDeg',          [], ...     % prescribed excitation flip (deg); auto-read from hdr if empty
    'T1_ms',            3000, ...   % PLACEHOLDER (unknown)
    'T2_ms',            200, ...    % PLACEHOLDER (unknown)
    'refShiftBoundsHz', [-40 40], ...% ~+-0.33 ppm at 7T, keeps RefExt on the 8 ppm peak
    'refLbBoundsHz',    [0 40], ... % NARROW phantom line (liquid phantom)
    'refLbInitHz',      15, ...
    'dropEndSlices',    1, ...      % B1 slab = all slices minus this many at EACH end
    'brainThr',         0.15, ...   % leg-mask threshold (fraction of max)
    'phantThr',         0.05, ...   % phantom-mask threshold (fraction of max)
    'phantMinVox',      5, ...      % min phantom size (voxels)
    'phantomMask',      [], ...     % optional explicit phantom mask override
    'addBroadBG',       true, ...   % add broad component under the +8 ppm envelope
    'phantomB1Method',  'extrap', ...% QC only while applyB1=false
    'b1ExtrapOrder',    2, ...
    'phantomB1Value',   [], ...
    'reciprocityReceive', true, ... % QC only while applyB1=false
    'brainVol_mL',      [], ...     % explicit leg-in-slab volume (mL) override
    'plt',              false);
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(e, fn{i}) || isempty(e.(fn{i}))
        e.(fn{i}) = d.(fn{i});
    end
end
if e.enable && isempty(e.b1Dir)
    e.b1Dir = findB1DirCalf(SessionDir);
end
end

% =====================================================================
function b1 = findB1DirCalf(SessionDir)
% Locate the B1 map DICOM folder under SessionDir\E* (E100/E200/E300).
% Returns '' when the session has no B1 scan -- NOT an error: the external
% quant then reports ok=false and the PCr-referenced quant proceeds.
% Calf sessions carry either PREP_TFL3D_B1MAP_10SLICES (3D, preferred) or
% PREP_TFL_2D_B1MAP (single slice); some also have a 3-file 31P MoCo map that
% must not win, hence the explicit ranking.
b1 = '';
exams = dir(fullfile(SessionDir,'E*'));
exams = exams([exams.isdir]);
cand  = [];
for k = 1:numel(exams)
    c = dir(fullfile(exams(k).folder, exams(k).name, '*B1MAP*'));
    if ~isempty(c), cand = [cand; c(:)]; end %#ok<AGROW>
end
if isempty(cand), return; end
cand = cand([cand.isdir]);
if isempty(cand), return; end

names = {cand.name};
score = zeros(numel(names),1);
score(contains(names,'TFL','IgnoreCase',true))   = 1;
score(contains(names,'TFL3D','IgnoreCase',true)) = 2;
nfile = zeros(numel(names),1);
for k = 1:numel(names)
    f = dir(fullfile(cand(k).folder, names{k}, '*'));
    nfile(k) = sum(~[f.isdir]);
end
[~, ord] = sortrows([score nfile], [-1 -2]);
b1 = fullfile(cand(ord(1)).folder, cand(ord(1)).name);
end
