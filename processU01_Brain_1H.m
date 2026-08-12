function processU01_Brain_1H(SessionDir, FinalOutDir, segTsv, plt, outStem)
% PROCESSU01_BRAIN_1H  Absolute NAD+ quantification for one U01 brain session.
%
%   processU01_Brain_1H(SessionDir, FinalOutDir, segTsv)
%   processU01_Brain_1H(SessionDir, FinalOutDir, segTsv, plt, outStem)
%
% Reads the water-suppressed and water-reference .dat files from
% <SessionDir>\datfiles\ (or <SessionDir> itself), parses CSF/GM/WM tissue
% fractions from the brainseg.sh TSV (either emitted layout -- see
% read_seg_fracs below), runs the processing/fitting in
% dat2mat_svsNAD.m, then converts each peak to absolute [mM] via absoluteQuant_svs.m.
%
% Output: <FinalOutDir>\<scanID>_brain_1H_abs.mat
%   (output, parsOut, hdr, absQ, voxel, segFracs, pars)
%
% When several candidate .dat files are present, the highest MID number is
% used.  See pick_latest_pair for the recognised naming conventions.
%


arguments
    SessionDir  (1,:) char
    FinalOutDir (1,:) char
    segTsv      (1,:) char
    plt         (1,1) logical = false
    outStem     (1,:) char = ''    % overrides <scanID> in the output filename
end

assert(isfolder(SessionDir), 'SessionDir not found: %s', SessionDir);
assert(isfile(segTsv),       'segTsv not found: %s',     segTsv);
if isempty(FinalOutDir), FinalOutDir = SessionDir; end
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
pars.fit.ppm_range = [8.7 10.6];
pars.fit.ph_range  = [-60 60];
pars.peaks.name    = {'NADH2','NADH6','NADH4','Trp'};
pars.peaks.range   = [9.25 9.45; 9.05 9.25; 8.8 9.05; 10.05 10.15];

% Fine-fit linewidth bounds in Hz, PER PEAK (rows follow pars.peaks.name; ppm
% conversion happens once f0 is read).  Without this the default [10 60] Hz
% applies to every peak, but Trp is genuinely broader: the coarse HSVD is
% allowed 120 Hz for it (hsvdCoarse.lwMaxHz below) and routinely returns ~64 Hz
% -- ABOVE the fine-fit cap -- so Trp's width pinned to its upper bound in most
% scans, and harder baseline priors only pushed it further onto the bound.
pars.fit.baselineOpt.widthBoundsHz = [10 60; 10 60; 10 60; 10 120];

% Two fit regions -- the NAD+ cluster and the narrow Trp singlet at ~10.1 ppm.
% ONE coarse HSVD + ONE phase covers both ('union'): the Trp region holds a
% single peak, and nadCoarsePhase needs two to fit a first-order phase, so
% per-region phasing left the Trp region with pc1 = 0 in every session and no
% phase at all where HSVD failed to seed Trp.  The baseline partition is still
% done per region inside the union branch.
pars.fit.coarseMode = 'union';
pars.fit.regions(1).fitRange = [8.7 10.0];
pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4'};
pars.fit.regions(1).name     = 'NAD';
pars.fit.regions(1).baselineOpt.knotSpacing = 0.43;
pars.fit.regions(2).fitRange = [9.8 10.6];
pars.fit.regions(2).peaks    = {'Trp'};
pars.fit.regions(2).name     = 'Trp';
pars.fit.regions(2).baselineOpt.knotSpacing = 2;

% Single-region alternative
% pars.fit.regions(1).fitRange = [8.7 10.3];
% pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4','Trp'};
% pars.fit.regions(1).name     = 'NAD/Trp';

pars.fit.baselineOpt.enable      = true;
pars.fit.baselineOpt.knotSpacing = 2;
pars.fit.baselineOpt.lambda      = 1;
pars.fit.baselineOpt.lambdaAmpl  = 0;
pars.fit.baselineOpt.TolFun      = 1e-6;

% HSVD-derived SOFT PRIORS on the fine fit
%   lambdaPrior: baseline penalty ||w.*(B*beta - hsvdBaseline)||^2. 
%   lambdaWidth: Voigt total FWHM vs the HSVD linewidth.
pars.fit.baselineOpt.lambdaPrior = 0;
pars.fit.baselineOpt.lambdaWidth = 0.005;
pars.fit.baselineOpt.widthPriorPeaks = [true true true false];  % NADH2 NADH6 NADH4 Trp
pars.fit.baselineOpt.etaInit         = [0 0 0 1];               % NADH2 NADH6 NADH4 Trp

% Lineshape kernel DISABLED: it is collinear with the Voigt Gaussian width
pars.fit.lineShapeOpt.enable     = false;
pars.fit.lineShapeOpt.nSide      = 4;
pars.fit.lineShapeOpt.asymmetric = true;
pars.fit.lineShapeOpt.maxSide    = 0.50;

pars.fit.hsvdClean.enable = false;

% ---------------------------------------------------------------------
% Coarse seeding / phase correction by HSVD
pars.fit.coarseFirstOrder       = false;
pars.fit.hsvdCoarse.enable      = true;
pars.fit.hsvdCoarse.K           = 60;
pars.fit.hsvdCoarse.npts        = 1024;
pars.fit.hsvdCoarse.lwMinHz     = [ 10  10  10  10];   % NADH2 NADH6 NADH4 Trp
pars.fit.hsvdCoarse.lwMaxHz     = [ 90  90  90 120];
pars.fit.hsvdCoarse.phaseTolDeg = [100 100 100 100];

pars.watfit.method    = 'fit';
pars.watfit.nComp     = 1;
pars.watfit.ppm_range = [4 5.5];
pars.watfit.etaInit   = 0;      % fit is insensitive to this seed (tested 0-1)

% =====================================================================
% Literature relaxation / water content (Swago 2025, NMR Biomed 5324)
% =====================================================================
WM  = struct('name','WM', 'T1_ms',1130, 'T2_ms', 40, 'waterConc_M',38.5, 'flag_metSig',true);
GM  = struct('name','GM', 'T1_ms',1939, 'T2_ms', 40, 'waterConc_M',44,   'flag_metSig',true);
CSF = struct('name','CSF','T1_ms',4400, 'T2_ms',300, 'waterConc_M',55,   'flag_metSig',false);

met_T1_ms = [205.6, 211.6, 237.3, 75];   % NADH2, NADH6, NADH4, Trp
met_T2_ms = [ 33.6,  29.1,  42.3, 19];
met_nProt = [    1,     1,     1,   1];
peakNames = {'NADH2','NADH6','NADH4','Trp'};
nPeaks    = numel(peakNames);

% =====================================================================
% Tissue fractions from brainseg.sh TSV
% =====================================================================
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));
segFracs = read_seg_fracs(segTsv);   % [CSF GM WM], sums to 1
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
% Two naming conventions are recognized, because different studies label the
% same two scans differently:
%
%   water-suppressed  : '_NAD_'    or  'ppm'
%   water reference   : '_Water_'  or  'waterref'
%


allF = dir(fullfile(datDir, '*svssel*.dat'));

wsPat  = {'_NAD_',   'ppm'   };   % water-suppressed (metabolite)
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
