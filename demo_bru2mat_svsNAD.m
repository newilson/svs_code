%% demo bru2mat_svsNAD - 1H NAD downfield on Bruker high-field data
clear,clc

%% Bruker data - 500 MHz Chemistry scanner

% Each acquisition lives in its own numbered Bruker experiment subfolder.
% Update wsExp / nwsExp once the 1H NAD experiment numbers are known --
% these are placeholders next to the 31P demo's bru_dir convention.
scan = 0;
switch scan
    case 0
        sessionDir = 'C:\Users\CAMIPM-NW\Downloads\rbnmr\rbnmr\data\Human_Blood_1';
        wsExp      = 4;   % water-suppressed 1H NAD 
        nwsExp     = 6;   % water reference  
    case 1
        sessionDir = 'C:\Users\CAMIPM-NW\Downloads\rbnmr\rbnmr\data\Human_Blood_2';
        wsExp      = 4;
        nwsExp     = 6;
end

ws_dir  = fullfile(sessionDir, num2str(wsExp));
nws_dir = fullfile(sessionDir, num2str(nwsExp));

use_nws = true;

%% Processing parameters
pars = struct();
pars.lb_met  = 3;
pars.lb_wat  = 3;
pars.eccopt  = 0;            % 0: Klose ECC (single-coil Bruker)
pars.PC      = [];            % skip manual phase GUI
pars.base    = nan;          % no semi-manual baseline
pars.plt     = true;
pars.pltSpec = true;

pars.center_freq  = 4.72;    % carrier on water for 1H NAD acquisitions
pars.initPhDeg    = 0;       % initial zero-order phase guess (degrees)
pars.initPPMShift = 0;       % initial chemical-shift offset (ppm)

pars.hsvd.flag  = false;
pars.hsvd.sv    = 30;

% Bruker is single coil + already-summed fid: no coil combination, no block
% averaging, no block alignment.  Leave water-suppression filters off (the
% sequence's CHESS handles it); switch on Post if residual water bleeds in.
pars.WatSupPre.opts.type         = 'none';
pars.WatSupPost.opts.type        = 'none';
pars.WatSupPost.opts.filt.type   = 'average';
pars.WatSupPost.opts.filt.N      = 11;
pars.WatSupPost.opts.hsvd.bounds = [-1 1] * 300/8000;
pars.WatSupPost.opts.hsvd.nsin   = 25;

%% Fitting parameters (multi-region VARPRO, same recipe as processU01_Brain_1H)
pars.dofit         = 'varpro';
pars.fit.mode      = 5;                 % Voigt
pars.fit.ppm_range = [8.9 10.3];

% Wide phase range + asymmetric lineshape kernel are intentional: they
% absorb real systematic distortion (eddy currents, B0 asymmetry, residual
% phase) without tying it into peak amplitudes.  See
% feedback-nad-1h-phase-and-lineshape-stay-wide.md before narrowing.
pars.fit.ph_range  = [-60 60];

pars.peaks.name    = {'NADH2','NADH6','NADH4','Trp'};
pars.peaks.range   = [9.35 9.55; 9.20 9.35; 8.9 9.20; 9.9 10.2]; % make these bigger
% Per-peak FINE-fit width seed (ppm).  Default (NAN / <=0) keeps the coarse-
% fitted width; explicit value forces the fine optimizer to start there.
% Trp on blood parks in a narrow local minimum during the NAD-dominant
% coarse fit; seeding it wider escapes that basin.  0.25 ppm overshot,
% 0.15 ppm lands closer to the visible Trp FWHM at 500 MHz.
pars.peaks.widthInit = [NaN; NaN; NaN; 0.15];

% Single wide NAD+Trp region (brain default).  Two-region alternative
% (NAD cluster vs narrow Trp window) is documented in processU01_Brain_1H.m
% if local-baseline isolation is needed.
pars.fit.coarseMode = 'perRegion';
% pars.fit.regions(1).fitRange = [8.7 10.0];
% pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4'};
% pars.fit.regions(1).name     = 'NAD';
% pars.fit.regions(2).fitRange = [9.95 10.3];   % narrow, around 10.1
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

pars.fit.lineShapeOpt.enable     = true;
pars.fit.lineShapeOpt.nSide      = 4;
pars.fit.lineShapeOpt.asymmetric = true;
pars.fit.lineShapeOpt.maxSide    = 0.50;   

pars.fit.hsvdClean.enable = false;

pars.watfit.method    = 'fit';
pars.watfit.nComp     = 1;                          % single peak; initial width set inside bru2mat_svsNAD
pars.watfit.ppm_range = [4 5.5];
% Test: force water to a pure Lorentzian (mode 3 = lorentzian_complex_basis).
% Hypothesis is that the observed broad water tails are 1/x^2 Lorentzian
% wings the Voigt (mode 5) cannot reach because its Gaussian component pulls
% the wings inward.  Mode 3 drops the wG / wT parameters; the optimizer then
% leans on wL + the lineshape kernel.  Comment out to revert to inherited
% Voigt (pars.fit.mode = 5).
pars.watfit.mode      = 3;                          % 3 = Lorentzian (complex)
% Inherit baselineOpt from the metabolite fit but allow water to be much
% wider than NAD peaks (radiation damping / exchange give blood water a
% substantial broad pedestal that the metabolite cap can't reach).  The
% initial width is set inside bru2mat_svsNAD (default 0.1 ppm, matched to
% the visible water FWHM at 500 MHz); override via pars.watfit.widthInit.
pars.watfit.baselineOpt = pars.fit.baselineOpt;
pars.watfit.baselineOpt.widthBounds = [0.02 1.0];   % ~10-500 Hz at 500 MHz

% Lineshape kernel for water: SYMMETRIC (water peak should be symmetric
% about the carrier) with enough side points to model the broad pedestal.
% Each side coef in [0, maxSide]; the optimizer learns positive shoulder
% weights that the convolution then turns into broad tails on top of the
% narrow Voigt core.  nSide=15 covers ~±0.06 ppm at 500 MHz, enough to
% model the inner shoulders; bump if visible tails are wider.
pars.watfit.widthInit = 0.05;                       % narrow base Voigt; kernel adds shoulders
pars.watfit.lineShapeOpt.enable     = true;
pars.watfit.lineShapeOpt.nSide      = 15;
pars.watfit.lineShapeOpt.asymmetric = false;
pars.watfit.lineShapeOpt.maxSide    = 1.0;
pars.watfit.lineShapeOpt.areaPadFactor = 100;       % wide trapz for Lorentzian 1/x^2 tails


%% process
tic
if use_nws
    [output,pars,fullHdr] = bru2mat_svsNAD(pars,ws_dir,nws_dir);
else
    [output,pars,fullHdr] = bru2mat_svsNAD(pars,ws_dir,[]);
end
toc

%% Quick report of fitted areas
if isfield(output,'met') && isfield(output.met,'fit') && isfield(output.met.fit,'areas')
    fprintf('\n*** 1H NAD areas (Bruker) ***\n');
    for k = 1:numel(output.met.fit.names)
        fprintf('  %-10s : %.4g\n', output.met.fit.names{k}, output.met.fit.areas(k));
    end
    if isfield(output,'wat') && isfield(output.wat,'fit') && isfield(output.wat.fit,'areaTotal')
        fprintf('  -----\n  water (%s) : %.4g\n', output.wat.fit.method, output.wat.fit.areaTotal);
    end
end
