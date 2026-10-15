function processNAD_Chan_1H(SessionDir, FinalOutDir, segTsv, scan, plt, outStem, combineRuns)
% PROCESSNAD_CHAN_1H  Absolute NAD+ quantification for one NAD_Chan brain scan.
%
%   processNAD_Chan_1H(SessionDir, FinalOutDir, segTsv, scan)
%   processNAD_Chan_1H(SessionDir, FinalOutDir, segTsv, scan, plt, outStem, combineRuns)
%
%   scan        : 'slab'  whole slab       (*svssel*_sl_NAD.dat  / _sl_WAT.dat)
%                 'vox'   occipital voxel  (*svssel*_vox_NAD.dat / _vox_WAT.dat)
%   combineRuns : false (default) uses the highest-MID water-suppressed run;
%                 true concatenates the averages of every run of that scan
%
% Finds the water-suppressed / water-reference .dat pair for the requested
% scan anywhere under SessionDir (highest MID wins), parses CSF/GM/WM
% fractions from the brainseg.sh TSV for that scan's mask, runs
% dat2mat_svsNAD.m, then converts each peak to absolute [mM] via
% absoluteQuant_svs.m.  Design notes: processNAD_Chan_NOTES.md.
%
% Output: <FinalOutDir>\<scanID>_<scan>_1H_abs.mat
%   (output, parsOut, hdr, absQ, voxel, segFracs, pars)

arguments
    SessionDir  (1,:) char
    FinalOutDir (1,:) char
    segTsv      (1,:) char
    scan        (1,:) char {mustBeMember(scan,{'slab','vox'})}
    plt         (1,1) logical = false
    outStem     (1,:) char = ''    % overrides <scanID> in the output filename
    combineRuns (1,1) logical = false
end

assert(isfolder(SessionDir), 'SessionDir not found: %s', SessionDir);
assert(isfile(segTsv),       'segTsv not found: %s',     segTsv);
if isempty(FinalOutDir), FinalOutDir = SessionDir; end
if ~isfolder(FinalOutDir), mkdir(FinalOutDir); end

[wsDat, watDat] = pick_scan_pair(SessionDir, scan, combineRuns);
fprintf('WS  : %s\n', wsDat{:});
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
pars.fit.ppm_range = [8.7 10.6];
pars.fit.ph_range  = [-60 60];
pars.peaks.name    = {'NADH2','NADH6','NADH4','Trp'};
pars.peaks.range   = [9.25 9.45; 9.05 9.25; 8.8 9.05; 10.05 10.15];
pars.fit.baselineOpt.widthBoundsHz = [10 60; 10 60; 10 60; 10 120];

pars.fit.coarseMode = 'union';
pars.fit.regions(1).fitRange = [8.7 10.0];
pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4'};
pars.fit.regions(1).name     = 'NAD';
pars.fit.regions(1).baselineOpt.knotSpacing = 0.43;
pars.fit.regions(2).fitRange = [9.8 10.6];
pars.fit.regions(2).peaks    = {'Trp'};
pars.fit.regions(2).name     = 'Trp';
pars.fit.regions(2).baselineOpt.knotSpacing = 2;

pars.fit.baselineOpt.enable      = true;
pars.fit.baselineOpt.knotSpacing = 2;
pars.fit.baselineOpt.lambda      = 1;
pars.fit.baselineOpt.lambdaAmpl  = 0;
pars.fit.baselineOpt.TolFun      = 1e-6;

pars.fit.baselineOpt.lambdaPrior = 0;
pars.fit.baselineOpt.lambdaWidth = 0.005;
pars.fit.baselineOpt.widthPriorPeaks = [true true true false];  % NADH2 NADH6 NADH4 Trp
pars.fit.baselineOpt.etaInit         = [0 0 0 1];               % NADH2 NADH6 NADH4 Trp

pars.fit.lineShapeOpt.enable     = false;
pars.fit.lineShapeOpt.nSide      = 4;
pars.fit.lineShapeOpt.asymmetric = true;
pars.fit.lineShapeOpt.maxSide    = 0.50;

pars.fit.hsvdClean.enable = false;

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
pars.watfit.etaInit   = 0;

% ---- occipital voxel (processNAD_Chan_NOTES.md) ---------------------
if strcmp(scan, 'vox')
    pars.fit.mode        = 3;           % complex Lorentzian
    pars.fit.ph_range    = [-20 20];
    pars.watfit.mode     = 5;
    pars.watfit.ph_range = [-60 60];
    pars.fit.regions = [];
    pars.fit.regions(1).fitRange = [8.7 10.6];
    pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4','Trp'};
    pars.fit.regions(1).name     = 'NAD/Trp';
    pars.fit.regions(1).baselineOpt.knotSpacing = 0.43;
end

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

% dat2mat_svsNAD rmpaths utils at exit
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
outPath = fullfile(FinalOutDir, sprintf('%s_%s_1H_abs.mat', scanID, scan));
save(outPath, 'output','parsOut','hdr','absQ','voxel','segFracs','pars');
fprintf('Saved -> %s\n', outPath);

rmpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));
end

% =====================================================================
function [wsDat, watDat] = pick_scan_pair(SessionDir, scan, combineRuns)
% Highest-MID *svssel*_<tag>_NAD.dat (or, with combineRuns, every one in MID
% order) as a cell, and the highest-MID *svssel*_<tag>_WAT.dat, anywhere
% under SessionDir, with tag = 'sl' (slab) or 'vox' (voxel).

tag  = struct('slab','sl', 'vox','vox');
tag  = tag.(scan);
allF = dir(fullfile(SessionDir, '**', '*svssel*.dat'));

wsF = strings(0);  wsMid  = [];
watF = strings(0); watMid = [];
for k = 1:numel(allF)
    fn  = allF(k).name;
    tok = regexp(fn, '[MF]ID(\d+)', 'tokens', 'once');
    if isempty(tok), continue; end
    mid = str2double(tok{1});
    full = string(fullfile(allF(k).folder, fn));
    if ~isempty(regexpi(fn, ['_' tag '_NAD\.dat$'], 'once'))
        wsF(end+1)  = full; wsMid(end+1)  = mid; %#ok<AGROW>
    elseif ~isempty(regexpi(fn, ['_' tag '_WAT\.dat$'], 'once'))
        watF(end+1) = full; watMid(end+1) = mid; %#ok<AGROW>
    end
end

assert(~isempty(wsF),  'No *svssel*_%s_NAD.dat under %s', tag, SessionDir);
assert(~isempty(watF), 'No *svssel*_%s_WAT.dat under %s', tag, SessionDir);

[~, iWS]  = sort(wsMid);
if ~combineRuns, iWS = iWS(end); end
[~, iWAT] = max(watMid);
wsDat  = cellstr(wsF(iWS));
watDat = char(watF(iWAT));
end
