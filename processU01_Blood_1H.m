function processU01_Blood_1H(SessionDir, FinalOutDir, wsExp, nwsExp, plt)
% PROCESSU01_BLOOD_1H  Absolute NAD+ quantification for one U01 blood session
%                      acquired on the Bruker high-field NMR.
%
%   processU01_Blood_1H(SessionDir, FinalOutDir, wsExp, nwsExp)
%   processU01_Blood_1H(SessionDir, FinalOutDir, wsExp, nwsExp, plt)
%
% Reads the water-suppressed and water-reference acquisitions from
% <SessionDir>\<wsExp>\ and <SessionDir>\<nwsExp>\ (each Bruker experiment
% is its own numbered subfolder containing fid/ser + acqus), runs the
% processing/fitting in bru2mat_svsNAD.m, then converts each peak to
% absolute [mM] via the water-reference model in absoluteQuant_svs.m.
%
% Whole blood is treated as a single compartment: absoluteQuant_svs is
% called with one tissue (fraction = 1, flag_metSig = true), which recovers
% the un-segmented Gasparovic case (f_CSF = 0).  The fit recipe mirrors
% processU01_Calf_1H.m (Voigt, asymmetric lineshape, wide phase).
%
% Output: <FinalOutDir>\<scanID>_blood_1H_abs.mat
%   (output, parsOut, hdr, absQ, voxel, pars)

arguments
    SessionDir  (1,:) char
    FinalOutDir (1,:) char
    wsExp       (1,1) double {mustBeInteger,mustBePositive}
    nwsExp      (1,1) double {mustBeInteger,mustBePositive}
    plt         (1,1) logical = false
end

assert(isfolder(SessionDir), 'SessionDir not found: %s', SessionDir);
if ~isfolder(FinalOutDir), mkdir(FinalOutDir); end

wsDir  = fullfile(SessionDir, num2str(wsExp));
nwsDir = fullfile(SessionDir, num2str(nwsExp));
assert(isfolder(wsDir),  'WS  exp folder not found: %s', wsDir);
assert(isfolder(nwsDir), 'NWS exp folder not found: %s', nwsDir);
fprintf('WS  : %s\n', wsDir);
fprintf('NWS : %s\n', nwsDir);

% =====================================================================
% Processing/Fitting parameters (same recipe as processU01_Calf_1H.m)
% =====================================================================
pars = struct();
pars.lb_met  = 3;
pars.lb_wat  = 3;
pars.eccopt  = 0;
pars.PC      = 0;
pars.base    = nan;
pars.plt     = plt;
pars.pltSpec = plt;

pars.center_freq  = 4.72;
pars.initPhDeg    = 0;       % initial zero-order phase guess (degrees)
pars.initPPMShift = 0;       % initial chemical-shift offset (ppm)

pars.WatSupPre.opts.type         = 'none';
pars.WatSupPost.opts.type        = 'none';
pars.WatSupPost.opts.filt.type   = 'average';
pars.WatSupPost.opts.filt.N      = 11;
pars.WatSupPost.opts.hsvd.bounds = [-1 1] * 300/8000;
pars.WatSupPost.opts.hsvd.nsin   = 25;

pars.dofit         = 'varpro';
pars.fit.mode      = 5;                 % Voigt
pars.fit.ppm_range = [8.9 10.3];
% NOTE: wide phase range is intentional.  Tightening ph_range to [-15 15]
% (and/or making lineShapeOpt symmetric / lowering its maxSide) makes the
% per-component figures LOOK cleaner -- pure-absorptive peaks, baseline no
% longer floats above the data -- but degrades NADH2:NADH6:NADH4 self-
% consistency across averages/sessions.  The wide phase + asymmetric
% lineShape kernel are absorbing real systematic distortion (eddy currents,
% B0 asymmetry, residual phase drift); take that freedom away and the
% distortion is forced into the peak amplitudes, where it pulls each NAD
% proton unequally.  Do not narrow without re-verifying NAD self-
% consistency.  See feedback-nad-1h-phase-and-lineshape-stay-wide.
pars.fit.ph_range  = [-60 60];
% NAD peak windows shifted +0.05-0.1 ppm relative to brain/calf to track
% the observed positions in whole blood (the brain windows leave each NAD
% peak hugging the upper-ppm boundary, with amplitudes pegged tiny).
% Revise once Bruker blood NAD chemical shifts are nailed down.
pars.peaks.name    = {'NADH2','NADH6','NADH4','Trp'};
pars.peaks.range   = [9.35 9.55; 9.20 9.35; 8.90 9.20; 9.90 10.20];
% Per-peak FINE-fit width seed (ppm).  Default (NAN / <=0) keeps the coarse-
% fitted width; explicit value forces the fine optimizer to start there.
% Trp on blood parks in a narrow local minimum during the NAD-dominant
% coarse fit; seeding it wider escapes that basin.  0.25 ppm overshot,
% 0.15 ppm lands closer to the visible Trp FWHM at 500 MHz.
pars.peaks.widthInit = [NaN; NaN; NaN; 0.15];

% Two-region alternative (NAD cluster vs narrow Trp window) is documented
% in processU01_Brain_1H.m if local-baseline isolation is needed.  Default
% is a single wide region covering NAD + Trp, like the brain wrapper.
pars.fit.coarseMode = 'perRegion';
% pars.fit.regions(1).fitRange = [8.9 10.0];
% pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4'};
% pars.fit.regions(1).name     = 'NAD';
% pars.fit.regions(2).fitRange = [9.85 10.3];   % narrow, around 10.05
% pars.fit.regions(2).peaks    = {'Trp'};
% pars.fit.regions(2).name     = 'Trp';
pars.fit.regions(1).fitRange = [8.9 10.3];
pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4','Trp'};
pars.fit.regions(1).name     = 'NAD/Trp';

pars.fit.baselineOpt.enable      = true;
pars.fit.baselineOpt.knotSpacing = 2;
pars.fit.baselineOpt.lambda      = 1;
pars.fit.baselineOpt.lambdaAmpl  = 0;
pars.fit.baselineOpt.TolFun      = 1e-6;
% Per-peak linewidth bounds (ppm).  NAD trio caps at 0.25 ppm (~125 Hz at
% 500 MHz); Trp gets a wider 0.6 ppm cap because the indole NH at ~10.1 ppm
% is in fast exchange with water in whole blood and can carry a few-hundred-
% Hz line.  Row order matches pars.peaks.name = {NADH2, NADH6, NADH4, Trp}.
pars.fit.baselineOpt.widthBounds = [0.02 0.25; ...
                                    0.02 0.25; ...
                                    0.02 0.25; ...
                                    0.02 0.60];

% NOTE: asymmetric shared lineshape kernel is intentional -- see the wide
% ph_range comment above.  Together they absorb real systematic line
% distortion; constraining either degrades NAD self-consistency.
pars.fit.lineShapeOpt.enable     = true;
pars.fit.lineShapeOpt.nSide      = 4;
pars.fit.lineShapeOpt.asymmetric = true;
pars.fit.lineShapeOpt.maxSide    = 0.50;

pars.fit.hsvdClean.enable = false;

pars.watfit.method    = 'fit';
pars.watfit.nComp     = 1;                          % single peak; initial width set inside bru2mat_svsNAD
pars.watfit.ppm_range = [4 5.5];
% Water base is pure Lorentzian (mode 3): 1/x^2 wings reach the broad
% pedestal the Voigt (mode 5) Gaussian component pulls inward.  Verified
% against demo_bru2mat_svsNAD on Human_Blood_1; the wide symmetric
% lineshape kernel below handles the residual shoulder structure.
pars.watfit.mode      = 3;                          % 3 = Lorentzian (complex)
% Inherit baselineOpt from the metabolite fit but allow water to be much
% wider than NAD peaks (radiation damping / exchange give blood water a
% substantial broad pedestal that the metabolite cap can't reach).  The
% initial width is set inside bru2mat_svsNAD (default 0.1 ppm, matched to
% the visible water FWHM at 500 MHz); override via pars.watfit.widthInit.
pars.watfit.baselineOpt = pars.fit.baselineOpt;
pars.watfit.baselineOpt.widthBounds = [0.02 1.0];   % ~10-500 Hz at 500 MHz

% Symmetric lineshape kernel models the broad water pedestal on top of the
% narrow Lorentzian core.  See demo_bru2mat_svsNAD.m for tuning notes.
pars.watfit.widthInit = 0.05;                       % narrow base Voigt; kernel adds shoulders
pars.watfit.lineShapeOpt.enable     = true;
pars.watfit.lineShapeOpt.nSide      = 15;
pars.watfit.lineShapeOpt.asymmetric = false;
pars.watfit.lineShapeOpt.maxSide    = 1.0;
% Wider area integration window than the default 10x: Lorentzian 1/x^2
% wings converge slowly, so a pad of ~100x captures >99.7% of the peak.
% Only affects output.wat.fit.areaTotal (and downstream absolute quant),
% not the fit itself.
pars.watfit.lineShapeOpt.areaPadFactor = 100;

% =====================================================================
% Single blood compartment: water content + relaxation.
% PLACEHOLDER values for whole blood (high-field). EDIT BEFORE PUBLICATION.
%   waterConc_M : whole blood is ~83% water by mass -> ~0.83 * 55.5 M ~= 46 M
%   T1/T2       : blood water at high field (rough literature placeholders)
% =====================================================================
Blood = struct('name','Blood', 'T1_ms', 2000, 'T2_ms', 50, ...
               'waterConc_M', 0.8 * 55, 'flag_metSig',true);

% Downfield-proton relaxation / proton counts. We don't think these values
% matter much based on sequence parameters
met_T1_ms = [300, 300, 300, 300];   % NADH2, NADH6, NADH4, Trp
met_T2_ms = [ 10,  10,  10,  10];
met_nProt = [    1,     1,     1,   1];
peakNames = {'NADH2','NADH6','NADH4','Trp'};
nPeaks    = numel(peakNames);

% =====================================================================
% Process/Fit
% =====================================================================
[output, parsOut, hdr] = bru2mat_svsNAD(pars, wsDir, nwsDir);

% bru2mat_svsNAD rmpaths utils at exit; re-add for absoluteQuant_svs and
% any other utils helpers used downstream.
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

% =====================================================================
% Voxel + absolute quant (single-compartment water reference)
% =====================================================================
voxel = struct();
voxel.tissues = Blood;
voxel.tissues(1).fraction = 1;

data = struct();
data.water.area = output.wat.fit.areaTotal;
% Use kernel-aware areas: `areas` is base-Voigt (ampl .* width) only and
% under-counts whenever pars.fit.lineShapeOpt is enabled; `areasConv`
% integrates the convolved component and matches the kernel-aware water
% areaTotal so the NAD/water ratio is internally consistent.
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
outPath = fullfile(FinalOutDir, sprintf('%s_blood_1H_abs.mat', scanID));
save(outPath, 'output','parsOut','hdr','absQ','voxel','pars');
fprintf('Saved -> %s\n', outPath);

% remove path
rmpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

end
