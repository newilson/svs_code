function processU01_Brain_31P(SessionDir, FinalOutDir, segTsv, plt, extRef)
% PROCESSU01_BRAIN_31P  Absolute 31P quantification for one U01 brain session.
%
%   processU01_Brain_31P(SessionDir, FinalOutDir, segTsv, plt)
%
% segTsv is OPTIONAL (QC only for 31P): pass '' or omit it to skip the
% tissue-fraction read; the absolute quant does not depend on it.
%
% Reads the highest-MID *31p*.dat from <SessionDir>\datfiles\, runs the
% basis-Varpro fit from demo_dat2mat_svs31p.m via dat2mat_svs31p.m, then
% converts each fitted metabolite area to absolute [mM] using aATP as an
% internal reference ([aATP]_brain = 2.7 mM, single brain-parenchyma
% compartment) with per-metabolite T1/T2 relaxation correction:
%
%   [M] = (S_M / S_aATP) * (R_aATP / R_M) * (nP_aATP / nP_M) * [aATP]_brain
%
%   R_x = (1 - exp(-TR/T1_x)) * exp(-TE/T2_x)
%
% No CSF partial-volume term: 31P metabolites and the aATP reference are
% both CSF-absent, so the (1 - f_CSF) dilution is common to numerator and
% denominator and cancels in the S_M/S_aATP ratio.  f_CSF is read from the
% seg TSV (when provided) and stored for QC only.  (A (1 - f_CSF) correction would apply
% only to an external or tissue-water reference, where the reference is
% not diluted by CSF the same way the metabolites are.)
%
% Single-compartment brain (no GM/WM split): GM/WM 31P differences are
% documented (Hetherington MRM 2001;45:46; Ren NMR Biomed 2015;28:1455;
% Lu NMR Biomed 2014;27:1135) but cannot be separated within one SVS
% voxel.
%
% Output: <FinalOutDir>\<scanID>_brain_31P_abs.mat
%   (output, parsOut, hdr, absQ, segFracs, pars)

arguments
    SessionDir  (1,:) char
    FinalOutDir (1,:) char
    segTsv      char = ''     % OPTIONAL: seg is QC-only for 31P (see below)
    plt         (1,1) logical = false
    extRef      struct = struct()   % OPTIONAL external-reference quant (off unless .enable=true)
end

% External-reference (phantom) absolute quant is entirely optional and OFF by
% default: with extRef.enable=false nothing below changes and the internal
% aATP reference path runs exactly as before.  See extRefDefaults (bottom) for
% every tunable field.
extRef    = extRefDefaults(extRef, SessionDir);
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
    pars.fit.refExtPpm     = extRef.refPpm;   % ~8 ppm for the D2O phantom
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
% Mariano 2014. EDIT BEFORE PUBLICATION USE.
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

refName    = 'aATP';
refConc_mM = 2.7;          % [aATP]_brain (parenchyma), placeholder
refIdx     = find(strcmp({metProps.name}, refName), 1);
assert(~isempty(refIdx), 'Reference %s missing from metProps', refName);

% =====================================================================
% Tissue fractions from brainseg.sh TSV (stored for QC; NOT used in quant
% -- the CSF dilution cancels with the internal aATP reference).  The seg
% is OPTIONAL: if segTsv is empty, missing, or unreadable the fractions are
% left as NaN and the 31P quant proceeds unchanged.
%   columns: labelval, voxels, volume(mm^3), labelname (CSF/GM/WM)
% =====================================================================
segFracs = [NaN NaN NaN];   % [CSF GM WM]
fCSF     = NaN;
if isempty(segTsv)
    fprintf('No seg TSV provided; tissue fractions left as NaN (QC only).\n');
elseif ~isfile(segTsv)
    fprintf('Seg TSV not found (%s); tissue fractions left as NaN (QC only).\n', segTsv);
else
    try
        T = readtable(segTsv,'FileType','text','Delimiter','\t', ...
                      'VariableNamingRule','preserve');
        voxCol  = find(strcmpi(T.Properties.VariableNames,'voxels'),    1);
        nameCol = find(strcmpi(T.Properties.VariableNames,'labelname'), 1);
        assert(~isempty(voxCol) && ~isempty(nameCol), ...
            'unexpected TSV columns: %s', strjoin(T.Properties.VariableNames, ','));
        nm  = string(T{:,nameCol});
        vx  = double(T{:,voxCol});
        nCSF = sum(vx(nm=="CSF"));
        nGM  = sum(vx(nm=="GM"));
        nWM  = sum(vx(nm=="WM"));
        tot  = nCSF + nGM + nWM;
        assert(tot > 0, 'no CSF/GM/WM voxels in table');
        segFracs = [nCSF nGM nWM] / tot;   % [CSF GM WM] (GM/WM kept for reference)
        fCSF     = segFracs(1);
        fprintf('Tissue fractions: CSF=%.3f  GM=%.3f  WM=%.3f\n', segFracs);
    catch ME
        warning('Could not read seg TSV %s (%s); tissue fractions left as NaN.', ...
                segTsv, ME.message);
    end
end

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
% Absolute quant: aATP internal reference (no CSF term -- see header)
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
    error('processU01_Brain_31P: no fit results found in dat2mat_svs31p output');
end

[areaRef, T1ref, T2ref, nPref] = lookup_fit(refName, fitNames, fitAreas, metProps);
Rref = relax(TR_ms, TE_ms, T1ref, T2ref);

nFit = numel(fitNames);
absQ = struct();
absQ.reference = struct('name',refName,'conc_mM',refConc_mM, ...
                        'T1_ms',T1ref,'T2_ms',T2ref,'nPhos',nPref,'R',Rref);
absQ.fCSF      = fCSF;
absQ.TR_ms     = TR_ms;
absQ.TE_ms     = TE_ms;
absQ.metabolites = struct('name',{},'area',{},'conc_mM',{}, ...
                          'T1_ms',{},'T2_ms',{},'nPhos',{},'R',{},'ratio',{});

for k = 1:nFit
    nm_k = fitNames{k};
    if strcmp(nm_k,'RefExt'), continue; end   % external ref: quantified separately below
    pIdx = find(strcmp({metProps.name}, nm_k), 1);
    if isempty(pIdx)
        warning('No T1/T2/nPhos entry for %s -- skipping quant', nm_k);
        continue
    end
    T1_k = metProps(pIdx).T1_ms;
    T2_k = metProps(pIdx).T2_ms;
    nP_k = metProps(pIdx).nPhos;
    R_k  = relax(TR_ms, TE_ms, T1_k, T2_k);

    ratio   = fitAreas(k) / areaRef;
    % No (1 - fCSF) term: with an internal reference (aATP) that is also
    % CSF-absent, the CSF dilution cancels in the S_M/S_aATP ratio.
    conc_mM = ratio * (Rref / R_k) * (nPref / nP_k) * refConc_mM;

    absQ.metabolites(end+1) = struct( ...
        'name', nm_k, 'area', fitAreas(k), 'conc_mM', conc_mM, ...
        'T1_ms', T1_k, 'T2_ms', T2_k, 'nPhos', nP_k, ...
        'R', R_k, 'ratio', ratio); %#ok<AGROW>
end

for k = 1:numel(absQ.metabolites)
    fprintf('%-10s : %.4g mM\n', absQ.metabolites(k).name, absQ.metabolites(k).conc_mM);
end

[~, scanID] = fileparts(SessionDir);

% =====================================================================
% Optional external-reference (phantom) absolute quant, for comparison
% against the internal aATP reference above.  B1+ map + phantom/brain
% volumes read from the coregistered prep_tfl3d_b1map slab.
% =====================================================================
if useExtRef
    % Prescribed excitation flip angle (deg) drives the B1+ efficiency
    % eta = sin(relB1 * flip).  1D-ISIS excites with a hard pulse at the
    % prescribed flip (70 deg on S095, not 90), so read it from the header.
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
outPath = fullfile(FinalOutDir, sprintf('%s_brain_31P_abs.mat', scanID));
save(outPath, 'output','parsOut','hdr','absQ','segFracs','pars','metProps','extRef');
fprintf('Saved -> %s\n', outPath);

% remove path
rmpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

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
function e = extRefDefaults(e, SessionDir)
% Fill defaults for the optional external-reference quant.  Every field is
% overridable by the caller; the whole feature is off unless e.enable=true.
if ~isstruct(e), e = struct(); end
d = struct( ...
    'enable',           false, ...
    'b1Dir',            '', ...      % auto-located under SessionDir\E100 if empty
    'refPpm',           8.0, ...     % external-reference chemical shift (PCr=0)
    'refConc_mM',       500, ...     % known phantom concentration (sodium phosphate)
    'nP',               1, ...       % equivalent 31P per molecule (sodium phosphate = 1)
    'flipDeg',          [], ...      % prescribed excitation flip (deg); auto-read from hdr if empty
    'T1_ms',            3000, ...    % PLACEHOLDER (unknown)
    'T2_ms',            200, ...     % PLACEHOLDER (unknown)
    'refShiftBoundsHz', [-40 40], ...% ~+-0.33 ppm at 7T, keeps RefExt on the 8 ppm peak
    'refLbBoundsHz',    [0 40], ...  % NARROW phantom line (liquid phantom); envelope -> baseline
    'refLbInitHz',      15, ...
    'nSlicesCentral',   8, ...       % central B1 slices spanning the excited slab
    'brainThr',         0.15, ...    % head-mask threshold (fraction of max), used as fallback
    'brainMaskNii',     '', ...      % HD-BET brain mask NIfTI (skull-stripped T1); overrides threshold
    'phantThr',         0.05, ...    % phantom-mask threshold (fraction of max)
    'phantMinVox',      5, ...       % min phantom size (voxels)
    'phantomMask',      [], ...      % optional explicit phantom mask override
    'addBroadBG',       true, ...    % add broad component under the +8 ppm envelope
    'phantomB1Method',  'extrap', ...% 'extrap' (poly-fit brain B1 -> phantom), 'measure', 'value'
    'b1ExtrapOrder',    2, ...       % polynomial order for B1 extrapolation
    'phantomB1Value',   [], ...      % relB1 used when phantomB1Method='value'
    'reciprocityReceive', true, ...  % add receive factor relB1_ref/relB1_M (reciprocity: assume B1-=B1+)
    'brainVol_mL',      [], ...      % explicit brain-in-slab volume (mL); overrides everything below
    't1Nii',            '', ...      % T1 NIfTI (for HD-BET brain INT ISIS-slab volume)
    'isisDcm',          '', ...      % one ISIS_SVS DICOM (auto-located under E100 if empty)
    'plt',              false);
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(e, fn{i}) || isempty(e.(fn{i}))
        e.(fn{i}) = d.(fn{i});
    end
end
if e.enable && isempty(e.b1Dir)
    e.b1Dir = findB1Dir(SessionDir);
end
if e.enable && isempty(e.isisDcm)
    e.isisDcm = findIsisDcm(SessionDir);
end
end

% =====================================================================
function dcm = findIsisDcm(SessionDir)
% One ISIS_SVS DICOM under SessionDir\E100 (defines the 31P VOI geometry).
dcm = '';
cand = dir(fullfile(SessionDir,'E100','*ISIS*'));
cand = cand([cand.isdir]);
if isempty(cand), return; end
d = dir(fullfile(cand(1).folder, cand(1).name, '*'));
d = d(~[d.isdir]);
if ~isempty(d)
    dcm = fullfile(d(1).folder, d(1).name);
end
end

% =====================================================================
function b1 = findB1Dir(SessionDir)
% Locate the prep_tfl3d_b1map DICOM folder under SessionDir\E100.
b1 = '';
cand = dir(fullfile(SessionDir,'E100','*B1MAP*'));
if isempty(cand), cand = dir(fullfile(SessionDir,'E100','*b1map*')); end
cand = cand([cand.isdir]);
if ~isempty(cand)
    b1 = fullfile(cand(1).folder, cand(1).name);
end
end
