function processU01_Brain_1H(SessionDir, FinalOutDir, segTsv, plt, outStem)
% PROCESSU01_BRAIN_1H  Absolute NAD+ quantification for one U01 brain session.
%
%   processU01_Brain_1H(SessionDir, FinalOutDir, segTsv)
%   processU01_Brain_1H(SessionDir, FinalOutDir, segTsv, plt, outStem)
%
% Reads the water-suppressed and water-reference .dat files from
% <SessionDir>\datfiles\ (or <SessionDir> itself), parses CSF/GM/WM voxel
% counts from the brainseg.sh TSV, runs the processing/fitting in
% dat2mat_svsNAD.m, then converts each peak to absolute [mM] via the
% 3-compartment Gasparovic model in absoluteQuant_svs.m.
%
% Output: <FinalOutDir>\<scanID>_brain_1H_abs.mat
%   (output, parsOut, hdr, absQ, voxel, segFracs, pars)
%
% When several candidate .dat files are present, the highest MID number is
% used.  See pick_latest_pair for the recognised naming conventions.
%
% Design rationale, parameter sweeps, rejected alternatives and open issues:
% see processU01_Brain_1H_NOTES.md (same folder).  Section numbers referenced
% below as "NOTES sN" refer to that file.

arguments
    SessionDir  (1,:) char
    FinalOutDir (1,:) char
    segTsv      (1,:) char
    plt         (1,1) logical = false
    outStem     (1,:) char = ''    % overrides <scanID> in the output filename
end

assert(isfolder(SessionDir), 'SessionDir not found: %s', SessionDir);
assert(isfile(segTsv),       'segTsv not found: %s',     segTsv);
if ~isfolder(FinalOutDir), mkdir(FinalOutDir); end

% Sessions are normally laid out as <SessionDir>\datfiles\*.dat, but some
% exported studies keep the .dat files directly in the session folder.
datDir = fullfile(SessionDir, 'datfiles');
if ~isfolder(datDir)
    datDir = SessionDir;
end

[wsDat, watDat] = pick_latest_pair(datDir);
fprintf('WS  : %s\n', wsDat);
fprintf('WAT : %s\n', watDat);

% =====================================================================
% Processing/Fitting parameters
% =====================================================================
pars.lb_met  = 3;
pars.lb_wat  = 3;
pars.eccopt  = 0;
pars.PC      = 0;
pars.base    = nan;
pars.plt     = plt;
pars.pltSpec = plt;

pars.block_align.size     = 64;
pars.block_align.method   = 'fd';
pars.block_align.Nfit     = 256;
pars.block_align.ppmRange = [8.5 10];
pars.ccopt.minsig_frac    = 0.05;

pars.WatSupPre.opts.type         = 'none';
pars.WatSupPost.opts.type        = 'none';
pars.WatSupPost.opts.filt.type   = 'average';
pars.WatSupPost.opts.filt.N      = 11;
pars.WatSupPost.opts.hsvd.bounds = [-1 1] * 300/8000;
pars.WatSupPost.opts.hsvd.nsin   = 25;

pars.dofit         = 'varpro';
pars.fit.mode      = 5;                 % Voigt
% Inert while pars.fit.regions is defined (only seeds focusRange, and focus
% weighting needs focusWeight > 1).  Kept consistent with the region edges.
pars.fit.ppm_range = [8.7 10.45];
% Wide by design.  Do NOT narrow before the coarse phase step is fixed -- the
% width is currently absorbing a first-order phase error, not real distortion,
% and the May 2026 test that justified it was confounded.  NOTES s4.
pars.fit.ph_range  = [-60 60];
pars.peaks.name    = {'NADH2','NADH6','NADH4','Trp'};
pars.peaks.range   = [9.25 9.45; 9.05 9.25; 8.8 9.05; 10.05 10.15];

% Two independent fit regions -- the NAD+ cluster and the narrow Trp singlet at
% ~10.1 ppm -- each running its own coarse->phase->fine fit.  A single wide
% region cannot serve both: a stiff baseline inflates Trp, a flexible one eats
% NADH4.  Splitting also keeps low-SNR Trp out of the NAD phase ramp.  NOTES s1.
pars.fit.coarseMode = 'perRegion';
% Region 1 lower edge 8.7, not 8.6: the extra 0.1 ppm turns sharply upward into
% an unmodelled ~8.5 ppm ATP2/adenine peak and needs a baseline the spline
% cannot follow, and that misfit lands on NADH4.  Cost: ~15% of NADH4's
% Lorentzian mass falls outside the window.  NOTES s1.
pars.fit.regions(1).fitRange = [8.7 10.0];
pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4'};
pars.fit.regions(1).name     = 'NAD';
% 3 intervals over 1.3 ppm -> 6 spline coefficients (ncoef = nIntervals + 3),
% the median flexibility the measured baseline needs here.  NOTES s2.
pars.fit.regions(1).baselineOpt.knotSpacing = 0.43;
% 0.65 ppm window keeps window:width ~3.2:1 for the intrinsically broad Trp.
pars.fit.regions(2).fitRange = [9.8 10.45];
pars.fit.regions(2).peaks    = {'Trp'};
pars.fit.regions(2).name     = 'Trp';
% Sentinel -> max(2, round(0.65/2)) = 2 intervals = 5 spline coefficients (the
% floor).  A 6-coef variant (knotSpacing 0.22) was tried and REJECTED: it wins
% on residual but puts an interior minimum in the baseline in 12/12 subjects,
% right under Trp -- peak/baseline degeneracy, not a better model.  Do NOT
% re-select this on Trp-window residual or on SpecTickle agreement.  NOTES s2.
pars.fit.regions(2).baselineOpt.knotSpacing = 2;
% Single-region alternative, retired Jul 2026 (H4/H2 = 2.44 vs 1.083 expected;
% Trp/H2 CV 45%).  Kept for A/B against the two-region config above.  NOTES s1.
% pars.fit.regions(1).fitRange = [8.7 10.3];
% pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4','Trp'};
% pars.fit.regions(1).name     = 'NAD/Trp';

pars.fit.baselineOpt.enable      = true;
pars.fit.baselineOpt.knotSpacing = 2;
pars.fit.baselineOpt.lambda      = 1;
pars.fit.baselineOpt.lambdaAmpl  = 0;
pars.fit.baselineOpt.TolFun      = 1e-6;

% HSVD-derived SOFT PRIORS on the fine fit (penalties, not constraints, so the
% data can still overrule them).  A hard width cap was tried first and rejected:
% a bound pins the width and dumps the discrepancy into the amplitude.
%   lambdaPrior: baseline penalty ||w.*(B*beta - hsvdBaseline)||^2.  OFF -- it
%     degrades every metric, because HSVD's lw_max truncates broad wings and so
%     leaves genuine peak signal in its "baseline".  Hook kept for a better
%     baseline estimate (untruncated HSVD, or a measured MM spectrum).
%   lambdaWidth: Voigt total FWHM vs the HSVD linewidth.  Does nearly all the
%     work and saturates -- anything in 0.05-2.0 is equivalent.
% NOTES s3 for the sweep table.
pars.fit.baselineOpt.lambdaPrior = 0.0;
pars.fit.baselineOpt.lambdaWidth = 0.3;

% Shared lineshape kernel DISABLED: it is collinear with the Voigt Gaussian
% width and was winning, leaving widthG railed at its lower bound.  Settings
% kept (disabled) for A/B tests.  NOTES s5.
pars.fit.lineShapeOpt.enable     = false;
pars.fit.lineShapeOpt.nSide      = 4;
pars.fit.lineShapeOpt.asymmetric = true;
pars.fit.lineShapeOpt.maxSide    = 0.50;

pars.fit.hsvdClean.enable = false;

% ---------------------------------------------------------------------
% Coarse seeding / phase correction by HSVD (replaces the coarse VARPRO).
% K/npts match SpecTickle's defaults; lw and phase acceptance windows follow
% its 7T 1H brain peak table.  phaseTolDeg matters -- loosening it lets
% inverted components seed a peak and wreck the phase ramp.  NOTES s6.
% Vectors are in pars.peaks.name order and get subset per fit region.
pars.fit.hsvdCoarse.enable      = true;
pars.fit.hsvdCoarse.K           = 60;
pars.fit.hsvdCoarse.npts        = 1024;
pars.fit.hsvdCoarse.lwMinHz     = [ 10  10  10  10];   % NADH2 NADH6 NADH4 Trp
pars.fit.hsvdCoarse.lwMaxHz     = [ 90  90  90 120];
pars.fit.hsvdCoarse.phaseTolDeg = [100 100 100 100];

pars.watfit.method    = 'fit';
pars.watfit.nComp     = 1;
pars.watfit.ppm_range = [4 5.5];

% =====================================================================
% Literature relaxation / water content (Swago 2025, NMR Biomed 5324)
% =====================================================================
WM  = struct('name','WM', 'T1_ms',1130, 'T2_ms', 40, 'waterConc_M',38.5, 'flag_metSig',true);
GM  = struct('name','GM', 'T1_ms',1939, 'T2_ms', 40, 'waterConc_M',44,   'flag_metSig',true);
CSF = struct('name','CSF','T1_ms',4400, 'T2_ms',300, 'waterConc_M',55,   'flag_metSig',false);

% Trp: indole NH proton (1H) at 10.1 ppm.  Trp T1/T2 (75/19 ms) are NOT
% FINALISED -- placeholders matching SpecTickle's constants for comparability,
% pending a measured Trp T1/T2.  Revisit before reporting absolute Trp.
% NAD+ values are literature and unchanged.  NOTES s7.
met_T1_ms = [205.6, 211.6, 237.3, 75];   % NADH2, NADH6, NADH4, Trp
met_T2_ms = [ 33.6,  29.1,  42.3, 19];
met_nProt = [    1,     1,     1,   1];
peakNames = {'NADH2','NADH6','NADH4','Trp'};
nPeaks    = numel(peakNames);

% =====================================================================
% Tissue fractions from brainseg.sh TSV
%   columns: labelval, voxels, volume(mm^3), labelname  (CSF/GM/WM)
% =====================================================================
T = readtable(segTsv,'FileType','text','Delimiter','\t', ...
              'VariableNamingRule','preserve');
voxCol  = find(strcmpi(T.Properties.VariableNames,'voxels'),    1);
nameCol = find(strcmpi(T.Properties.VariableNames,'labelname'), 1);
assert(~isempty(voxCol) && ~isempty(nameCol), ...
    'Unexpected TSV columns in %s: %s', segTsv, ...
    strjoin(T.Properties.VariableNames, ','));
nm  = string(T{:,nameCol});
vx  = double(T{:,voxCol});
nCSF = sum(vx(nm=="CSF"));
nGM  = sum(vx(nm=="GM"));
nWM  = sum(vx(nm=="WM"));
tot  = nCSF + nGM + nWM;
assert(tot > 0, 'No CSF/GM/WM voxels in %s', segTsv);
segFracs = [nCSF nGM nWM] / tot;     % [CSF GM WM]
fprintf('Tissue fractions: CSF=%.3f  GM=%.3f  WM=%.3f\n', segFracs);

% =====================================================================
% Process/Fit
% =====================================================================
[output, parsOut, hdr] = dat2mat_svsNAD(pars, wsDat, watDat);

% dat2mat_svsNAD rmpaths utils at exit; re-add for absoluteQuant_svs and
% any other utils helpers used downstream.
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

% =====================================================================
% Voxel + absolute quant (Gasparovic 3-compartment)
% =====================================================================
voxel = struct();
voxel.tissues = [CSF GM WM];
voxel.tissues(1).fraction = segFracs(1);
voxel.tissues(2).fraction = segFracs(2);
voxel.tissues(3).fraction = segFracs(3);

data = struct();
data.water.area = output.wat.fit.areaTotal;
% Use areasConv, not areas: it integrates the convolved component and is
% phase-invariant (area of the ZERO-PHASE model peak), so quantitation does not
% depend on the fitted phase.  `areas` is base-Voigt (ampl .* width) only.
for k = 1:nPeaks
    data.metabolites(k).name     = peakNames{k};
    data.metabolites(k).area     = output.met.fit.areasConv(k);
    data.metabolites(k).nProtons = met_nProt(k);
    data.metabolites(k).T1_ms    = met_T1_ms(k);
    data.metabolites(k).T2_ms    = met_T2_ms(k);
end
seq = struct('TR_ms',   output.met.TR_ms, 'TE_ms',   output.met.TE_ms, ...
             'TR_ws_ms',output.wat.TR_ms, 'TE_ws_ms',output.wat.TE_ms);

absQ = absoluteQuant_svs(data, seq, voxel);

for k = 1:nPeaks
    fprintf('%-6s : %.4g uM\n', peakNames{k}, absQ.metabolites(k).conc_mM * 1e3);
end

% =====================================================================
% Save
% =====================================================================
% Output stem defaults to the session folder name; pass outStem to override
% it (used when re-running a session in place and the existing
% <scanID>_brain_1H_abs.mat must be preserved as a baseline).
[~, scanID] = fileparts(SessionDir);
if ~isempty(outStem), scanID = outStem; end
outPath = fullfile(FinalOutDir, sprintf('%s_brain_1H_abs.mat', scanID));
save(outPath, 'output','parsOut','hdr','absQ','voxel','segFracs','pars');
fprintf('Saved -> %s\n', outPath);

% remove path
rmpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

end

% =====================================================================
function [wsDat, watDat] = pick_latest_pair(datDir)
% Highest-MID water-suppressed and water-reference .dat among files
% containing 'svssel'.
%
% Two naming conventions are recognised, because different studies label the
% same two scans differently:
%
%   water-suppressed  : '_NAD_'    or  'OSrem'
%   water reference   : '_Water_'  or  'waterref'
%
% 'waterref' deliberately does not match the '_Water_' test (no trailing
% underscore), so the conventions cannot cross-match; the water-suppressed
% test is applied first regardless.  NOTES s9.

allF = dir(fullfile(datDir, '*svssel*.dat'));

wsPat  = {'_NAD_',   'OSrem'   };   % water-suppressed (metabolite)
watPat = {'_Water_', 'waterref'};   % water reference

nNAD = 0; nWAT = 0;
nadFiles = strings(0); watFiles = strings(0);
nadMid   = [];        watMid   = [];

for k = 1:numel(allF)
    fn = allF(k).name;
    tok = regexp(fn, 'MID(\d+)', 'tokens', 'once');
    if isempty(tok), continue; end          % require MID# to disambiguate
    midNum = str2double(tok{1});

    if any(contains(fn, wsPat, 'IgnoreCase', true))
        nNAD = nNAD + 1;
        nadFiles(nNAD) = string(fn); %#ok<AGROW>
        nadMid(nNAD)   = midNum;     %#ok<AGROW>
    elseif any(contains(fn, watPat, 'IgnoreCase', true))
        nWAT = nWAT + 1;
        watFiles(nWAT) = string(fn); %#ok<AGROW>
        watMid(nWAT)   = midNum;     %#ok<AGROW>
    end
end

assert(nNAD > 0, 'No water-suppressed (%s) .dat in %s', ...
       strjoin(wsPat, ' / '),  datDir);
assert(nWAT > 0, 'No water-reference (%s) .dat in %s',  ...
       strjoin(watPat, ' / '), datDir);

[~, iWS]  = max(nadMid);
[~, iWAT] = max(watMid);
wsDat  = fullfile(datDir, char(nadFiles(iWS)));
watDat = fullfile(datDir, char(watFiles(iWAT)));
end
