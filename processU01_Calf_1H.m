function processU01_Calf_1H(SessionDir, FinalOutDir, plt)
% PROCESSU01_CALF_1H  Absolute NAD+ quantification for one U01 calf session.
%
%   processU01_Calf_1H(SessionDir, FinalOutDir)
%   processU01_Calf_1H(SessionDir, FinalOutDir, plt)
%
% Reads the water-suppressed (*_NAD_*) and water-reference (*_Water_*)
% .dat files from <SessionDir>\datfiles\, runs the processing/fitting in
% dat2mat_svsNAD.m, then converts each peak to absolute [mM] via the
% water-reference model in absoluteQuant_svs.m.
%
% Unlike the brain pipeline, calf is treated as a SINGLE muscle compartment:
% there is no anatomic segmentation, so absoluteQuant_svs is called with one
% tissue (fraction = 1, flag_metSig = true), which recovers the un-segmented
% Gasparovic case (f_CSF = 0).  
%
% Output: <FinalOutDir>\<scanID>_calf_1H_abs.mat
%   (output, parsOut, hdr, absQ, voxel, pars)
%
% When multiple *_NAD_* or *_Water_* .dat files are present, the one with
% the highest MID number is used.

arguments
    SessionDir  (1,:) char
    FinalOutDir (1,:) char
    plt         (1,1) logical = false
end

assert(isfolder(SessionDir), 'SessionDir not found: %s', SessionDir);
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
pars.block_align.ppmRange = [8.95 10.7];
pars.ccopt.minsig_frac    = 0.05;

pars.WatSupPre.opts.type         = 'none';
pars.WatSupPost.opts.type        = 'none';
pars.WatSupPost.opts.filt.type   = 'average';
pars.WatSupPost.opts.filt.N      = 11;
pars.WatSupPost.opts.hsvd.bounds = [-1 1] * 300/8000;
pars.WatSupPost.opts.hsvd.nsin   = 25;

pars.dofit         = 'varpro';
pars.fit.mode      = 5;                 % Voigt
pars.fit.ppm_range = [8.95 10.7];
% NOTE: wide phase range is intentional.  Tightening ph_range to [-15 15]
% (and/or making lineShapeOpt symmetric / lowering its maxSide; see below)
% makes the per-component figures LOOK cleaner -- pure-absorptive peaks,
% baseline no longer floats above the data -- but tested (May 2026 on S095)
% it makes NADH2:NADH6:NADH4 LESS self-consistent across averages/sessions.
% The wide phase + asymmetric lineShape kernel are absorbing real systematic
% distortion (eddy currents, B0 asymmetry, residual phase drift); take that
% freedom away and the distortion is forced into the peak amplitudes, where
% it pulls each NAD proton unequally.  Do not narrow without re-verifying
% NAD self-consistency.  See feedback-nad-1h-phase-and-lineshape-stay-wide.
pars.fit.ph_range  = [-60 60];

pars.peaks.name    = {'NADH2','NADH6','NADH4','NR','Trp','Peak10_3'};
pars.peaks.range   = [9.25 9.45;     % NADH2     ~9.35
                      9.05 9.25;     % NADH6     ~9.15
                      8.80 9.05;     % NADH4     ~8.93
                      9.6 9.75;     % NR        ~9.70
                      10.05 10.15;   % Trp       ~10.1  (small)
                      10.25 10.35];  % Peak10_3  ~10.3  (larger, unidentified)

% Two independent fit groups, each running its own coarse->phase->fine fit
% (coarseMode 'perRegion') with a local baseline:
%   1) the NAD+ cluster (H2/H6/H4)
%   2) the downfield group: NR (~9.7), a small Trp singlet (~10.1), and a
%      larger unidentified peak (~10.3)
% Fitting them in separate windows keeps each group's baseline and phase
% local, so the weak Trp peak is not perturbed by the NAD cluster or the
% larger 10.3 peak.
% pars.fit.coarseMode = 'perRegion';
% pars.fit.regions(1).fitRange = [8.7 9.5];
% pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4'};
% pars.fit.regions(1).name     = 'NAD';
% pars.fit.regions(2).fitRange = [9.5 10.4];
% pars.fit.regions(2).peaks    = {'NR','Trp','Peak10_3'};
% pars.fit.regions(2).name     = 'NR/Trp/10.3';
% --- Alternative: single region over the whole downfield range -------------
pars.fit.regions(1).fitRange = [8.95 10.7];
pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4','NR','Trp','Peak10_3'};
pars.fit.regions(1).name     = 'all';

pars.fit.baselineOpt.enable      = true;
pars.fit.baselineOpt.knotSpacing = 2;
pars.fit.baselineOpt.lambda      = 1;
pars.fit.baselineOpt.lambdaAmpl  = 0;
pars.fit.baselineOpt.TolFun      = 1e-6;

% NOTE: asymmetric shared lineshape kernel is intentional -- see the wide
% ph_range comment above.  Together they absorb real systematic line
% distortion; constraining either degrades NAD self-consistency.  maxSide
% is 0.20 here (calf) vs 0.50 in the brain driver -- calf data appears to
% tolerate a slightly tamer kernel; do not lower further without re-checking
% NAD ratios.
pars.fit.lineShapeOpt.enable     = true;
pars.fit.lineShapeOpt.nSide      = 4;
pars.fit.lineShapeOpt.asymmetric = true;
pars.fit.lineShapeOpt.maxSide    = 0.20;

pars.fit.hsvdClean.enable = false;

pars.watfit.method    = 'fit';
pars.watfit.nComp     = 1;
pars.watfit.ppm_range = [4 5.5];

% =====================================================================
% Single muscle compartment: water content + relaxation.
%   waterConc_M : ~0.77 g water/g tissue * 55.5 M ~= 42.7 mol/L
%   T1/T2       : muscle water at 7T (rough literature placeholders)
% =====================================================================
Muscle = struct('name','Muscle', 'T1_ms', 2000, 'T2_ms', 25, ...
                'waterConc_M', 0.78 * 55, 'flag_metSig',true);

% Downfield-proton relaxation / proton counts.
% PLACEHOLDER values for H4, Trp, and Peak10_3.
met_T1_ms = [50, 50, 50, 100, 100, 100]; % NADH2,NADH6,NADH4,NR,Trp,Peak10_3
met_T2_ms = [ 14,  14,  14,  11,  11,  11];
met_nProt = [    1,     1,     1,   1,   1,   1];
peakNames = {'NADH2','NADH6','NADH4','NR','Trp','Peak10_3'};
nPeaks    = numel(peakNames);

% =====================================================================
% Process/Fit
% =====================================================================
[output, parsOut, hdr] = dat2mat_svsNAD(pars, wsDat, watDat);

% dat2mat_svsNAD rmpaths utils at exit; re-add for absoluteQuant_svs and
% any other utils helpers used downstream.
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

% =====================================================================
% Voxel + absolute quant (single-compartment water reference)
% =====================================================================
voxel = struct();
voxel.tissues = Muscle;
voxel.tissues(1).fraction = 1;

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
outPath = fullfile(FinalOutDir, sprintf('%s_calf_1H_abs.mat', scanID));
save(outPath, 'output','parsOut','hdr','absQ','voxel','pars');
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
