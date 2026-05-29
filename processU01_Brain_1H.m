function processU01_Brain_1H(SessionDir, FinalOutDir, segTsv, plt)
% PROCESSU01_BRAIN_1H  Absolute NAD+ quantification for one U01 brain session.
%
%   processU01_Brain_1H(SessionDir, FinalOutDir, segTsv)
%
% Reads the water-suppressed (*_NAD_*) and water-reference (*_Water_*)
% .dat files from <SessionDir>\datfiles\, parses CSF/GM/WM voxel counts
% from the brainseg.sh TSV, runs the processing/fitting in dat2mat_svsNAD.m, then converts each
% peak to absolute [mM] via the 3-compartment Gasparovic model in
% absoluteQuant_svs.m.
%
% Output: <FinalOutDir>\<scanID>_brain_1H_abs.mat
%   (output, parsOut, hdr, absQ, voxel, segFracs, pars)
%
% When multiple *_NAD_* or *_Water_* .dat files are present, the one with
% the highest MID number is used.

arguments
    SessionDir  (1,:) char
    FinalOutDir (1,:) char
    segTsv      (1,:) char
    plt         (1,1) logical = false
end

assert(isfolder(SessionDir), 'SessionDir not found: %s', SessionDir);
assert(isfile(segTsv),       'segTsv not found: %s',     segTsv);
if ~isfolder(FinalOutDir), mkdir(FinalOutDir); end

datDir = fullfile(SessionDir, 'datfiles');
assert(isfolder(datDir), 'datfiles\\ subfolder not found in %s', SessionDir);

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
% pars.fit.ppm_range = [8.7 10.0];
pars.fit.ppm_range = [8.7 10.3];
% NOTE: wide phase range is intentional.  Tightening ph_range to [-15 15]
% (and/or making lineShapeOpt symmetric / lowering its maxSide; see below)
% makes the per-component figures LOOK cleaner -- pure-absorptive peaks,
% baseline no longer floats above the data -- but tested (May 2026) it makes
% NADH2:NADH6:NADH4 LESS self-consistent across averages/sessions.  The wide
% phase + asymmetric lineShape kernel are absorbing real systematic
% distortion (eddy currents, B0 asymmetry, residual phase drift); take that
% freedom away and the distortion is forced into the peak amplitudes, where
% it pulls each NAD proton unequally.  Do not narrow without re-verifying
% NAD self-consistency.  See feedback-nad-1h-phase-and-lineshape-stay-wide.
pars.fit.ph_range  = [-60 60];
pars.peaks.name    = {'NADH2','NADH6','NADH4','Trp'};
pars.peaks.range   = [9.25 9.45; 9.05 9.25; 8.8 9.05; 10.05 10.15];

% Two independent fit regions: the NAD+ cluster and a narrow Trp singlet at
% ~10.1 ppm, each running its own coarse->phase->fine fit (coarseMode
% 'perRegion').  Trp's narrow window lets the local spline baseline flex enough
% to absorb the broad macromolecular hump under Trp without inflating Trp, while
% the NAD window keeps a stiff baseline that does not nibble the broad NADH4
% peak.  A single wide region cannot do both at once: a stiff baseline inflates
% Trp ~2x, a flexible one eats NADH4, and an explicit broad component
% (BroadMM, name-detected by dat2mat_svsNAD) is multimodal/fragile across
% sessions (verified on S007/S076).
pars.fit.coarseMode = 'perRegion';
% pars.fit.regions(1).fitRange = [8.7 10.0];
% pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4'};
% pars.fit.regions(1).name     = 'NAD';
% pars.fit.regions(2).fitRange = [9.95 10.3];   % fairly narrow, around 10.1
% pars.fit.regions(2).peaks    = {'Trp'};
% pars.fit.regions(2).name     = 'Trp';
pars.fit.regions(1).fitRange = [8.7 10.3];
pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4','Trp'};
pars.fit.regions(1).name = 'NAD/Trp';

pars.fit.baselineOpt.enable      = true;
pars.fit.baselineOpt.knotSpacing = 2;
pars.fit.baselineOpt.lambda      = 1;
pars.fit.baselineOpt.lambdaAmpl  = 0;
pars.fit.baselineOpt.TolFun      = 1e-6;

% NOTE: asymmetric shared lineshape kernel is intentional -- see the wide
% ph_range comment above.  Together they absorb real systematic line
% distortion; constraining either degrades NAD self-consistency.
pars.fit.lineShapeOpt.enable     = true;
pars.fit.lineShapeOpt.nSide      = 4;
pars.fit.lineShapeOpt.asymmetric = true;
pars.fit.lineShapeOpt.maxSide    = 0.50;

pars.fit.hsvdClean.enable = false;

pars.watfit.method    = 'fit';
pars.watfit.nComp     = 1;
pars.watfit.ppm_range = [4 5.5];

% =====================================================================
% Literature relaxation / water content (Swago 2025, NMR Biomed 5324)
% =====================================================================
WM  = struct('name','WM', 'T1_ms',1130, 'T2_ms', 40, 'waterConc_M',38.5, 'flag_metSig',true);
GM  = struct('name','GM', 'T1_ms',1939, 'T2_ms', 40, 'waterConc_M',44,   'flag_metSig',true);
CSF = struct('name','CSF','T1_ms',4400, 'T2_ms',300, 'waterConc_M',55,   'flag_metSig',false);

% Trp: indole NH proton (1H) at 10.1 ppm; T1/T2 set to 100/10 ms.
met_T1_ms = [205.6, 211.6, 237.3, 100];   % NADH2, NADH6, NADH4, Trp
met_T2_ms = [ 33.6,  29.1,  42.3,  10];
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
% Use kernel-aware areas: `areas` is base-Voigt (ampl .* width) only and
% under-counts whenever pars.fit.lineShapeOpt is enabled; `areasConv`
% integrates the convolved component (see dat2mat_svsNAD / curvefitAuto_
% varproBaseline).  When the water fit has no kernel and the metabolite
% fit does, swapping `areas`->`areasConv` raises the reported metabolite
% concentrations by the kernel factor (~10-40% with the current settings).
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
outPath = fullfile(FinalOutDir, sprintf('%s_brain_1H_abs.mat', scanID));
save(outPath, 'output','parsOut','hdr','absQ','voxel','segFracs','pars');
fprintf('Saved -> %s\n', outPath);

% remove path
rmpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

end

% =====================================================================
function [wsDat, watDat] = pick_latest_pair(datDir)
% Highest-MID *_NAD_* (water-suppressed) and *_Water_* (water-ref) .dat
% among files containing 'svssel'.

allF = dir(fullfile(datDir, '*svssel*.dat'));

nNAD = 0; nWAT = 0;
nadFiles = strings(0); watFiles = strings(0);
nadMid   = [];        watMid   = [];

for k = 1:numel(allF)
    fn = allF(k).name;
    tok = regexp(fn, 'MID(\d+)', 'tokens', 'once');
    if isempty(tok), continue; end          % require MID# to disambiguate
    midNum = str2double(tok{1});

    if contains(fn, '_NAD_', 'IgnoreCase', true)
        nNAD = nNAD + 1;
        nadFiles(nNAD) = string(fn); %#ok<AGROW>
        nadMid(nNAD)   = midNum;     %#ok<AGROW>
    elseif contains(fn, '_Water_', 'IgnoreCase', true)
        nWAT = nWAT + 1;
        watFiles(nWAT) = string(fn); %#ok<AGROW>
        watMid(nWAT)   = midNum;     %#ok<AGROW>
    end
end

assert(nNAD > 0, 'No water-suppressed (*_NAD_*) .dat in %s',   datDir);
assert(nWAT > 0, 'No water-reference (*_Water_*) .dat in %s',  datDir);

[~, iWS]  = max(nadMid);
[~, iWAT] = max(watMid);
wsDat  = fullfile(datDir, char(nadFiles(iWS)));
watDat = fullfile(datDir, char(watFiles(iWAT)));
end
