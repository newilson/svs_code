function processU01_Brain_31P(SessionDir, FinalOutDir, segTsv, plt, extRef, ref31P)
% PROCESSU01_BRAIN_31P  Absolute 31P quantification for one U01 brain session.
%
%   processU01_Brain_31P(SessionDir, FinalOutDir, segTsv, plt)
%   processU01_Brain_31P(SessionDir, FinalOutDir, segTsv, plt, extRef, ref31P)
%
% segTsv is OPTIONAL (QC only for 31P): pass '' or omit it to skip the
% tissue-fraction read; the absolute quant does not depend on it.
%
% ref31P selects which reference populates absQ.reference/absQ.metabolites:
% 'aATP' (default), 'PCr', or 'External'.  Every available set is computed and
% saved under absQ.byRef / absQ.external regardless of the choice; switching
% rescales all concentrations by one common factor and never changes a ratio
% between metabolites within a scan.  'External' requires extRef.enable=true.
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
    ref31P      (1,:) char = 'aATP' % primary reference: aATP | PCr | External
end

% External-reference (phantom) absolute quant is entirely optional and OFF by
% default: with extRef.enable=false nothing below changes and the internal
% aATP reference path runs exactly as before.  See extRefDefaults (bottom) for
% every tunable field.
extRef    = extRefDefaults(extRef, SessionDir, FinalOutDir, segTsv);
useExtRef = extRef.enable;
if useExtRef, extRef.plt = extRef.plt || plt; end

% Which reference populates absQ.reference / absQ.metabolites.  All available
% sets are computed and saved regardless; this only picks the primary.
if isempty(ref31P), ref31P = 'aATP'; end
assert(any(strcmpi(ref31P, {'aATP','PCr','External'})), ...
       'ref31P must be aATP, PCr or External (got ''%s'')', ref31P);
assert(~strcmpi(ref31P,'External') || useExtRef, ...
       'ref31P=External requires the external quant to be enabled (extRef.enable).');

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
pars.blp.flag   = false;
pars.dropPts    = 8;
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
    'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'};
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
pars.fit.fitOpt.nStart      = 8;

sc(1) = struct('idxA','UDPGal',    'idxB','UDPGlcNAc','ratio',0.53,'lambda',1e4);
sc(2) = struct('idxA','UDPGlc',    'idxB','UDPGlcNAc','ratio',1.00,'lambda',1e4);
sc(3) = struct('idxA','UDPGalNAc', 'idxB','UDPGlcNAc','ratio',0.80,'lambda',1e4);
pars.fit.softConstraints = sc;

pars.fit.tieGroups   = struct('members',{{'aATP','NADplus','NADH','UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'}, {'aATP','NADplus','NADH'}, {'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'}}, ...
                              'params',{{'lbL','lbG'}, {'phi0'}, {'phi0'}});
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
              'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'}, ...
    'T1_ms', {6000, 6000, 3000, 3000, 3500, 6000, 6000, 3000, ...
              3800, 1200, 1200, 1500, ...
               400,  400, ...
               500,  500,  500,  500}, ...
    'T2_ms', {200, 200, 100, 100, 150, 200, 200, 150, ...
              180,  40,  40,  30, ...
               10,  10, ...
               10,  10,  10,  10}, ...
    ... nPhos here is DELIBERATELY 1 for every species, including the NAD and
    ... UDP pyrophosphates (2 P each).  The spin count is already divided out
    ... upstream: dat2mat_svs31p computes amplScaled = ampl ./ n from
    ... basisInfo.n (2 for NADH/NADplus/UDP*), and the drivers pass amplScaled
    ... as fitAreas.  Setting 2 here would apply the same correction twice and
    ... halve NAD+/NADH a second time.
    'nPhos', {1,1,1,1,1,1,1,1,  1,1,1,1,  1,1,  1,1,1,1});

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
%   PRIMARY (the one written to absQ.metabolites) is aATP for brain.
% =====================================================================
refName    = 'aATP';       % primary reference for brain
refConc_mM = 2.7;          % [aATP]_brain (parenchyma), placeholder
refIdx     = find(strcmp({metProps.name}, refName), 1);
assert(~isempty(refIdx), 'Reference %s missing from metProps', refName);

% Secondary internal reference, reported alongside.
% [PCr]_brain PLACEHOLDER -- brain PCr is roughly 3.0-3.5 mM at 7T
% (Hetherington MRM 2001;45:46; Lu NMR Biomed 2014;27:1135).  EDIT BEFORE
% PUBLICATION USE, exactly like the T1/T2 table above.
altRefName    = 'PCr';
altRefConc_mM = 3.5;

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

absQ = struct();
absQ.fCSF      = fCSF;
absQ.TR_ms     = TR_ms;
absQ.TE_ms     = TE_ms;

% Quantify against BOTH internal references.  No (1 - fCSF) term in either:
% aATP and PCr are themselves CSF-absent, so the dilution cancels in the
% S_M/S_ref ratio.
absQ.byRef.(refName)    = refQuant31P(refName,    refConc_mM, ...
                                      fitNames, fitAreas, metProps, TR_ms, TE_ms);
absQ.byRef.(altRefName) = refQuant31P(altRefName, altRefConc_mM, ...
                                      fitNames, fitAreas, metProps, TR_ms, TE_ms);

% Primary reference populates absQ.reference/absQ.metabolites, which is what
% the Stage 5 summary and checkU01_2026_relaxation read.  With ref31P=External
% this internal set is the basis for extRefQuant31P's ext_over_int column and
% is replaced by the external set after that call.
if strcmpi(ref31P,'PCr'), priName = altRefName; else, priName = refName; end
absQ.reference   = absQ.byRef.(priName).reference;
absQ.metabolites = absQ.byRef.(priName).metabolites;
absQ.primaryRef  = priName;

fprintf('\n--- 31P concentrations [mM] by reference ---\n');
printByRef(absQ.byRef, {refName, altRefName});

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
    if isempty(extRef.slabThickness_mm)
        try
            extRef.slabThickness_mm = double(hdr.MeasYaps.sSpecPara.sVoI.dThickness);
        catch
            warning(['extRef: could not read ISIS slab thickness from hdr; V_ref ' ...
                     'falls back to the segmented phantom mask.']);
        end
    end

    % Acquisition excitation geometry, for the spectral profile term.  These
    % are THIS scan's values and differ from the B1 calibration's, which uses
    % a longer pulse to reach flips past 90 deg.
    [extRef.acqPulseDur, extRef.acqExcCtr, extRef.hzPerPpm] = acqExcGeom(hdr);
    extRef.acqRefVolt = double(hdr.MeasYaps.sTXSPEC.asNucleusInfo{1}.flReferenceAmplitude);
    extRef.basisInfo = basisInfoOf(output);
    extRef.line      = output.refExtLine;
    extRef.refPpm    = output.refExtLine.ppm;

    try
        absQ.external = extRefQuant31P(fitNames, fitAreas, metProps, ...
            TR_ms, TE_ms, extRef, absQ);
    catch ME
        warning('External-reference quant failed: %s', ME.message);
        absQ.external = struct('ok', false, 'error', ME.message);
    end
end

% Promote the external set to primary when asked.  Hard error rather than a
% silent fall back to an internal reference: the two differ by a scale factor,
% so quietly substituting one would mislabel every concentration.
if strcmpi(ref31P,'External')
    assert(isfield(absQ,'external') && isfield(absQ.external,'ok') && absQ.external.ok, ...
           'ref31P=External but the external quant did not succeed (see absQ.external).');
    absQ.reference   = struct('name','RefExt','conc_mM',extRef.refConc_mM, ...
                              'T1_ms',NaN,'T2_ms',NaN,'nPhos',NaN,'R',NaN);
    absQ.metabolites = absQ.external.metabolites;
    absQ.primaryRef  = 'External';
end
fprintf('31P primary reference: %s\n', absQ.primaryRef);

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
% Highest-MID(or FID) *31p*.dat in datDir.
allF = dir(fullfile(datDir, '*31p*.dat'));
assert(~isempty(allF), 'No *31p*.dat found in %s', datDir);

mids = nan(1, numel(allF));
for k = 1:numel(allF)
    tok = regexp(allF(k).name, '[MF]ID(\d+)', 'tokens', 'once');
    if ~isempty(tok), mids(k) = str2double(tok{1}); end
end
assert(any(~isnan(mids)), 'No MID#/FID# in any *31p*.dat filename in %s', datDir);

[~, iMax] = max(mids);
datFile = fullfile(datDir, allF(iMax).name);
end

% =====================================================================
function e = extRefDefaults(e, SessionDir, FinalOutDir, segTsv)
% Fill defaults for the optional external-reference quant.  Every field is
% overridable by the caller; the whole feature is off unless e.enable=true.
if ~isstruct(e), e = struct(); end
d = struct( ...
    'enable',           false, ...
    'anatomy',          'brain', ... % 'brain' | 'calf' -- selects tissue segmentation
    'applyB1',          false, ...   % B1+ / receive correction DISABLED (QC-reported only)
    'b1CalibFiles',     [], ...     % explicit B1CALIB DICOMs (cellstr); auto-located under SessionDir\E* if empty
    'refPpm',           8.0, ...     % nominal phantom shift (PCr=0); used only if no line is located
    'refConc_mM',       500, ...     % known phantom concentration (sodium phosphate)
    'nP',               1, ...       % equivalent 31P per molecule (sodium phosphate = 1)
    'flipDeg',          [], ...      % prescribed excitation flip (deg); auto-read from hdr if empty
    'T1_ms',            3000, ...    % PLACEHOLDER (unknown)
    'T2_ms',            200, ...     % PLACEHOLDER (unknown)
    'refShiftBoundsHz', [-40 40], ...% ~+-0.33 ppm at 7T, keeps RefExt on the 8 ppm peak
    'refLbBoundsHz',    [0 40], ...  % NARROW phantom line (liquid phantom); envelope -> baseline
    'refLbInitHz',      15, ...
    'refArea_mm2',      200, ...     % measured phantom cross-section (180-220 observed)
    'slabThickness_mm', [], ...      % ISIS slab (mm); auto-read from hdr if empty
    'brainMaskNii',     '', ...      % HD-BET brain mask NIfTI (skull-stripped T1); overrides threshold
    'addBroadBG',       true, ...    % add broad component under the +8 ppm envelope
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
% Session root, so extRefQuant31P can locate the *B1CALIB* series itself.
e.sessionDir = SessionDir;

if e.enable && isempty(e.isisDcm)
    e.isisDcm = findIsisDcm(SessionDir);
end
if e.enable && (isempty(e.brainMaskNii) || isempty(e.t1Nii))
    anatDir = '';
    if ~isempty(segTsv), anatDir = fileparts(segTsv); end
    [mk, t1] = findAnatNii({anatDir, FinalOutDir});
    if isempty(e.brainMaskNii), e.brainMaskNii = mk; end
    if isempty(e.t1Nii),        e.t1Nii        = t1; end
end
end

% =====================================================================
function [maskNii, t1Nii] = findAnatNii(dirs)
% HD-BET brain mask + the T1 it was derived from, as staged by Stage 1/1b of
% processU01_Brain.ps1.  Both land in the same folder: during a live run that
% is the staging dir segTsv points into; when this stage is re-run on its own
% it is the session's processed\ dir, which Stage final copied them into.
%
% Without these there is no V_M at all and the external quant stops: the old
% intensity-threshold head mask (skull and scalp included, ~790 mL against a
% real brain-in-slab ~560 mL, so ~1.4x low) was removed with the rest of the
% B1 mapping precisely because it failed silently.
maskNii = ''; t1Nii = '';
for k = 1:numel(dirs)
    d = dirs{k};
    if isempty(d) || ~isfolder(d), continue; end
    m = dir(fullfile(d, '*_hdbet.nii.gz'));
    if isempty(m), continue; end
    cand = fullfile(m(1).folder, m(1).name);
    % the T1 hd-bet ran on: same name minus the _hdbet suffix
    stem = erase(m(1).name, '_hdbet.nii.gz');
    t = dir(fullfile(d, [stem '.nii*']));
    if isempty(t), continue; end          % a mask without its T1 is unusable
    maskNii = cand;
    t1Nii   = fullfile(t(1).folder, t(1).name);
    return
end
end

% =====================================================================
function dcm = findIsisDcm(SessionDir)
% One ISIS_SVS DICOM under any SessionDir\E* (defines the 31P VOI geometry).
%
% Searches every exam folder and skips *ISIS* series that exported with no
% files -- that happens (S092_07082026), and taking the empty one would leave
% V_M on the head-mask fallback with nothing in the log to say why.
dcm = '';
exams = dir(fullfile(SessionDir,'E*'));
exams = exams([exams.isdir]);
for k = 1:numel(exams)
    cand = dir(fullfile(exams(k).folder, exams(k).name, '*ISIS*'));
    cand = cand([cand.isdir]);
    for j = 1:numel(cand)
        d = dir(fullfile(cand(j).folder, cand(j).name, '*'));
        d = d(~[d.isdir]);
        if ~isempty(d)
            dcm = fullfile(d(1).folder, d(1).name);
            return
        end
    end
end
warning(['processU01_Brain_31P: no ISIS DICOM with any files under %s\\E* -- ' ...
         'the external quant cannot intersect the brain mask with the slab and ' ...
         'V_M falls back to the intensity-threshold head (~1.4x too large). ' ...
         'Check absQ.external.V_M_source.'], SessionDir);
end


% =====================================================================
function [Tus, ctrPpm, hzPerPpm] = acqExcGeom(hdr)
% Excitation pulse duration (us), center (ppm) and Hz/ppm for THIS scan, from
% the same WIP slots isis_svs_refVolt reads.  Any of them empty makes
% extRefQuant31P fall back to the on-resonance sine (no profile term), so a
% missing field degrades the correction rather than corrupting it.
Tus = []; ctrPpm = []; hzPerPpm = [];
try Tus    = wipEntry(hdr,'alFree', 6);  catch, end
try ctrPpm = wipEntry(hdr,'adFree',11);  catch, end
try hzPerPpm = double(hdr.MeasYaps.sTXSPEC.asNucleusInfo{1}.lFrequency) * 1e-6; catch, end
if isempty(ctrPpm), ctrPpm = 0; end       % absent slot means centered, as on scanner
end

% =====================================================================
function v = wipEntry(hdr, fld, idx)
% One sWipMemBlock entry, empty cells read as 0 (matching dat2mat_svs31p).
v = [];
w = hdr.MeasYaps.sWipMemBlock.(fld);
w(cellfun('isempty',w)) = {0};
if numel(w) >= idx, v = double(w{idx}); end
end

% =====================================================================
function bi = basisInfoOf(output)
% Per-component chemical shifts from the basis fit, used to give every
% metabolite its own excitation-profile factor.  Empty when the fit did not
% carry basisInfo, which falls back to one tissue value for all of them.
bi = [];
try
    if isfield(output,'basisVarproResults') && isfield(output.basisVarproResults,'basisInfo')
        bi = output.basisVarproResults.basisInfo;
    end
catch
end
end
