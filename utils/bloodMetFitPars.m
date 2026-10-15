function pars = bloodMetFitPars()
% BLOODMETFITPARS  Processing + fitting recipe for the blood 1H pipeline:
% downfield NAD/Trp fit, water reference (bloodWatfitPars) and TSP reference.
% Shared by processU01_Blood_1H and the relaxometry analyses so they all fit
% the same way.  Callers set pars.plt / pars.pltSpec, and pars.serIndex to
% fit one row of a 2D ser instead of the average.

pars = struct();
pars.lb_met  = 6;
pars.lb_wat  = 6;
pars.eccopt  = 0;  % 0 = Klose (water-phase divide); restored to original f_water in bru2mat_svsNAD
pars.PC      = 0;
pars.base    = nan;

pars.initPhDeg    = 0;       % initial zero-order phase guess (degrees)
pars.initPPMShift = 0;       % initial chemical-shift offset (ppm)

pars.WatSupPre.opts.type         = 'none';
pars.WatSupPost.opts.type        = 'none';
pars.WatSupPost.opts.filt.type   = 'average';
pars.WatSupPost.opts.filt.N      = 11;
pars.WatSupPost.opts.hsvd.bounds = [-1 1] * 300/8000;
pars.WatSupPost.opts.hsvd.nsin   = 25;

pars.dofit         = 'varpro';
pars.fit.mode      = 5;                 % Lorentzian
pars.fit.ppm_range = [8.9 10.3];        % stops ABOVE the steep 8.55-9.0 aromatic envelope (and NAD+ H4 ~8.74) so the stiff baseline below stays well-behaved (validated via relaxometry run_t2relax_narrow)

% NOTE: wide phase range is intentional.  Tightening ph_range to [-15 15]
% (and/or making lineShapeOpt symmetric / lowering its maxSide) makes the
% per-component figures LOOK cleaner -- pure-absorptive peaks, baseline no
% longer floats above the data -- but degrades NADH2:NADH6 self-
% consistency across averages/sessions.  The wide phase + asymmetric
% lineShape kernel are absorbing real systematic distortion (eddy currents,
% B0 asymmetry, residual phase drift); take that freedom away and the
% distortion is forced into the peak amplitudes, where it pulls each NAD
% proton unequally.  Do not narrow without re-verifying NAD self-
% consistency.  See feedback-nad-1h-phase-and-lineshape-stay-wide.
pars.fit.ph_range  = [-60 60];

% ---------------------------------------------------------------------
% Metabolite model.  Optimized across the rbnmr blood datasets
% (Human_Blood_{3,5,6,8,11,12}); harness + diagnostics in _blood_opt/.
%
%  * SINGLE BROAD Trp.  The blood Trp feature is one broad ~0.12-0.18 ppm
%    blob (sometimes with a faint doublet shoulder); a single Voigt captures
%    its AREA robustly.  The earlier narrow TrpA/TrpB doublet was abandoned:
%    on these datasets the blob is broad, so the narrow pair pinned at its
%    width/range walls and undershot.
%  * NADH2 (~9.2-9.45) and NADH6 (~9.0-9.2) are the two PROMINENT NAD+ peaks.
%    Broad9_4 is a separate MINOR shoulder near ~9.5-9.6 (it is never one of
%    the 9.2/9.4 peaks).  The 'Broad' prefix routes Broad9_4 through the
%    broad-component path: excluded from the coarse phase fit, seeded broad,
%    re-added as a free component in the FINE fit, still quantified per-peak.
% ---------------------------------------------------------------------
pars.peaks.name    = {'Broad9_4','NADH2','NADH6','Trp'};
pars.peaks.range   = [9.50 9.70; 9.22 9.48; 9.00 9.22; 9.90 10.20];

% Per-peak FINE-fit width seed (ppm), aligned with pars.peaks.name
% {Broad9_4,NADH2,NADH6,Trp}.  NaN keeps the coarse-fitted width; Trp is
% seeded broad (0.10 ppm) so the fine fit starts inside the broad blob rather
% than a narrow local minimum.  Length MUST equal numel(pars.peaks.name)
% (bru2mat_svsNAD asserts this).
pars.peaks.widthInit = [NaN; NaN; NaN; 0.10];

% ---------------------------------------------------------------------
% Coarse fit: HSVD decomposition (same method as the Siemens brain/calf path
% in dat2mat_svsNAD).  Components landing inside a peak's range that pass the
% linewidth and phase filters seed that peak's center/width/phase; everything
% else is treated as baseline.  Set .enable = false to fall back to the legacy
% coarse VARPRO fit.
%
% The linewidth bounds are in Hz and are SPECIFIC TO THIS EXPERIMENT (Bruker
% 14.1 T, f0 = 600.13 MHz, selgpse_cwl).  They are NOT field-scaled -- a new
% experiment gets its own bounds, measured the same way: run once, dump
% output.met.fit.region(r).coarse.hsvd, and read off the linewidths of the
% components that land in each peak's range.  Measured across the six rbnmr
% use-this datasets (phase-qualifying components only):
%     NADH2  34.3 - 64.6 Hz   (one 190 Hz outlier in Blood_5 = background,
%                              deliberately excluded by the 90 Hz ceiling)
%     NADH6  27.8 - 61.0 Hz
%     Trp    65.9 - 118.7 Hz  (one 22.9 Hz sliver in Blood_6, amplitude 6 vs
%                              100 for the real peak, excluded by the floor)
% Bounds below take those ranges with margin.  Broad9_4 is excluded from the
% coarse phase fit entirely (the 'Broad' prefix), so its row is unused and is
% only present to keep the vectors aligned with pars.peaks.name.
%
% phaseTolDeg = 60: at the 180 default nothing is ever rejected on phase, so
% two near-opposed components in one peak range get averaged into a single
% seed and partially cancel (Blood_11: +88 and -54 deg both landing on Trp).
% zeroOrderOnly = true: see _blood_opt/coarse_phase_notes.md.
pars.fit.hsvdCoarse.enable        = true;
pars.fit.hsvdCoarse.K             = 60;
pars.fit.hsvdCoarse.npts          = 1024;
%                                Broad9_4  NADH2  NADH6   Trp
pars.fit.hsvdCoarse.lwMinHz       = [   40;    20;    20;    40];
pars.fit.hsvdCoarse.lwMaxHz       = [  300;    90;    90;   150];
pars.fit.hsvdCoarse.phaseTolDeg   = 60;
pars.fit.hsvdCoarse.zeroOrderOnly = true;

% DATA-DRIVEN per-dataset re-seeding of the NAD/Trp center ranges.  The
% downfield referencing drifts ~0.15 ppm between sessions; that drift parks a
% NAD component in a valley and the >=0 amplitude solve in the fitter then
% NULLS it (the "NADH6 collapse").  For each listed peak, bru2mat_svsNAD finds
% the strongest peak of the background-subtracted spectrum within 'search' and
% tightens that peak's range to detection +/- 'halfWidth'.  The NADH6/NADH2
% split at 9.23 separates the two NAD peaks in every dataset (NADH6 <=9.20,
% NADH2 >=9.25).  Broad9_4 is NOT auto-seeded (kept at its fixed ~9.5-9.7
% window so it cannot wander onto a NAD peak).
pars.peaks.autoSeed.enable      = true;
pars.peaks.autoSeed.bgSmoothPPM = 0.30;
pars.peaks.autoSeed.peaks = struct( ...
    'name',      {'NADH6',      'NADH2',      'Trp'        }, ...
    'search',    {[9.00 9.23],  [9.23 9.50],  [9.86 10.20] }, ...
    'halfWidth', {0.07,         0.07,         0.13         });

% SINGLE region spanning NAD + Trp.  One region matters: the strong Trp blob
% anchors the coarse phase line, so the derived phase correction is near
% identity and the NAD region stays absorptive.  A separate NAD-only region
% phases off just 2 peaks, mis-estimates the phase, and rotates the NAD peaks
% dispersive -- whose >=0 amplitudes then null.
pars.fit.coarseMode = 'perRegion';
pars.fit.regions(1).fitRange = [8.9 10.3];
pars.fit.regions(1).peaks    = {'NADH2','NADH6','Trp','Broad9_4'};
pars.fit.regions(1).name     = 'NAD/Trp';

pars.fit.baselineOpt.enable      = true;
% knotSpacing 0.5 (moderately stiff).  The [8.9 10.3] window stops ABOVE the
% steep 8.55-9.0 aromatic envelope, so the in-band background is gentle and a
% stiff spline tracks it cleanly -- no wiggle for peak tails to trade area with.
% Validated on the relaxometry data (run_t2relax_narrow): this range+stiffness
% gives a smooth baseline and a stable NADH6.  Do NOT widen the window below
% ~8.9 at this stiffness -- the steep envelope then needs a flexible
% (knot <= 0.25) spline or it leaks into the peaks (and would re-require H4).
pars.fit.baselineOpt.knotSpacing = 0.5;
pars.fit.baselineOpt.lambda      = 1;
pars.fit.baselineOpt.lambdaAmpl  = 0;
pars.fit.baselineOpt.TolFun      = 1e-6;
% Per-peak linewidth bounds (ppm).  Row order matches pars.peaks.name =
% {Broad9_4, NADH2, NADH6, Trp}.  NAD peaks cap at 0.20 ppm; Trp is allowed
% broad (up to 0.25 ppm) to fill the blob.  The Broad9_4 row is a placeholder:
% the 'Broad' prefix makes bru2mat_svsNAD override it with the broad-component
% width bounds set in pars.fit.broadWidthBounds below.
pars.fit.baselineOpt.widthBounds = [0.02  0.25; ...   % Broad9_4 (overridden by broad bounds)
                                    0.02  0.20; ...   % NADH2
                                    0.02  0.20; ...   % NADH6
                                    0.03  0.25];      % Trp (broad)

% Broad-component (Broad9_4) width control.  bru2mat_svsNAD's defaults
% ([65/f0 180/f0] = 65-180 Hz, init 150 Hz) were tuned for the 7T in-vivo
% macromolecular baseline and are far too broad for high-field blood: a
% ~150 Hz line at 9.5 ppm leaks into NADH2 and steals NAD area.  Constrain it
% to "a bit broader than NAD, not macromolecular" -- at 600 MHz these are
% ~24-90 Hz bounds, ~42 Hz init.  Both are in ppm.
pars.fit.broadWidthBounds = [0.04 0.15];   % ppm  (~24-90 Hz @ 600 MHz)
pars.fit.broadWidthInit   = 0.07;          % ppm  (~42 Hz @ 600 MHz)

% Shared lineshape kernel DISABLED for blood.  The asymmetric kernel is
% intentional for the in-vivo (brain/calf) pipeline -- it absorbs gradient
% eddy-current / B0 line distortion -- but the high-field NMR blood lines are
% clean and well-phased, so the kernel adds no benefit and slightly raises the
% residual (sweep: kernel-off mean score 0.147 vs 0.159 kernel-on across the 6
% datasets).  The wide ph_range above is kept.  See
% feedback-nad-1h-phase-and-lineshape-stay-wide (which is in-vivo-specific).
pars.fit.lineShapeOpt.enable     = false;
pars.fit.lineShapeOpt.nSide      = 4;
pars.fit.lineShapeOpt.asymmetric = true;
pars.fit.lineShapeOpt.maxSide    = 0.50;

pars.fit.hsvdClean.enable = false;

pars.watfit = bloodWatfitPars();   % 2 co-located Voigts, shared with processU01_Blood_waterT2

pars.tspfit.lb          = 6;
pars.tspfit.search      = [-0.5 0.35];   % raw-axis ppm window to locate TSP
pars.tspfit.window      = [-0.5 0.3];    % fit window relative to located TSP
pars.tspfit.mode        = 5;
pars.tspfit.centerHalf  = 0.05;
pars.tspfit.widthInit   = 0.06;
pars.tspfit.widthBounds = [0.02 0.30];
pars.tspfit.ph_range    = [-60 60];
pars.tspfit.baseline    = true;
end
