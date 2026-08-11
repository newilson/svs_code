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
pars.lb_met  = 6;
pars.lb_wat  = 6;
pars.eccopt  = 0;  % 0 = Klose (water-phase divide); restored to original f_water in bru2mat_svsNAD
pars.PC      = 0;
pars.base    = nan;
pars.plt     = plt;
pars.pltSpec = plt;

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
% linewidth and phase filters seed that peak's centre/width/phase; everything
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

pars.watfit.method    = 'fit';
pars.watfit.ppm_range = [4 5.5];

% Two-component Lorentzian water model: blood water is super-Lorentzian
% (a sharp ~0.029 ppm core sitting on heavier-than-Lorentzian wings).
% Physically this is a DISTRIBUTION of Lorentzians -- intra/extracellular
% (plasma vs RBC) compartments with different susceptibility/T2*, plus
% voxel B0 spread -- NOT the dipolar semisolid super-Lorentzian (water is
% mobile, so residual dipolar coupling averages to zero).  A single
% component (Voigt or Lorentzian) structurally cannot match both the sharp
% cusp and the heavy wings: one amplitude locks the core:wing ratio, so it
% undershoots the apex 30-40% and leaves a bipolar derivative-shaped
% residual at the apex.  TWO co-located Lorentzians (narrow core + broad
% pedestal, both pinned ON the water line) resolve it: peak residual drops
% from ~0.40 to ~0.018 of apex.  Components stay separated and stable via
% three controls below: (1) co-located init, (2) NON-OVERLAPPING per-
% component width bounds, (3) tight per-component center bounds.  No
% lineshape kernel -- the broad component already supplies the pedestal;
% adding a kernel on top is redundant and drives the broad Lorentzian to
% its width bound while its area collapses (verified).
pars.watfit.mode      = 3;                          % 3 = Lorentzian (complex)
pars.watfit.nComp     = 2;                          % narrow core + broad pedestal

pars.watfit.twoStagePhase = false;                  % no kernel => no warm-start needed

% Co-located seeds + tight, per-component bounds.  Both components start ON
% the water line (4.70); centerBounds keep them there (narrow tighter than
% broad).  widthInit seeds one narrow, one broad; the per-component
% widthBounds are DISJOINT ([0.015 0.05] vs [0.05 0.45]) so the optimizer
% physically cannot merge the two into a single width.
pars.watfit.centerInit   = [4.70 4.70];
pars.watfit.centerBounds = [4.64 4.76;              % narrow core
                            4.58 4.82];             % broad pedestal
pars.watfit.widthInit    = [0.025 0.10];            % [narrow broad] ppm
pars.watfit.baselineOpt  = pars.fit.baselineOpt;
pars.watfit.baselineOpt.enable      = false;        % flat baseline; 2 comps carry all area
pars.watfit.baselineOpt.widthBounds = [0.015 0.050; % narrow core band
                                       0.050 0.450];% broad pedestal band
pars.watfit.lineShapeOpt.enable     = false;        % no kernel (broad comp is the pedestal)
% Wider area integration window than the default 10x: Lorentzian 1/x^2
% wings converge slowly, so a pad of ~100x captures >99.7% of the peak.
% Only affects output.wat.fit.areaTotal (and downstream absolute quant),
% not the fit itself.
pars.watfit.lineShapeOpt.areaPadFactor = 100;

% =====================================================================
% Single blood compartment: water content + relaxation.
%   waterConc_M : whole blood is ~83% water by mass -> ~0.83 * 55.5 M ~= 46 M
%                 (0.8*55 = 44 M used below; still a rough placeholder)
%   T2_ms       : MEASURED on the water arm of the blood relaxometry series
%                 (_t2relax, selective spin-echo, TE 20-82 ms): 29.3 +/- 3.4 ms,
%                 R^2 = 0.996.  Unaffected by the TE=10 ms question that set the
%                 metabolite T2 above -- the water series has no 10 ms scan.
%   T1_ms       : still a PLACEHOLDER (rough literature value) and the dominant
%                 systematic on every absolute concentration here: conc scales
%                 with (1-exp(-TR/T1_water)), which is ~0.18 at TR = 500 ms.
%                 EDIT BEFORE PUBLICATION.
% =====================================================================
Blood = struct('name','Blood', 'T1_ms', 2500, 'T2_ms', 29.3, ...
               'waterConc_M', 0.8 * 55, 'flag_metSig',true);

% Downfield-proton relaxation / proton counts.  Placeholders for now (same
% for every peak); kept as explicit per-peak arrays because they will be
% replaced with measured per-peak T1/T2 later.  Trp is now a single broad
% component, quantified as one peak.
%
% ORDER MUST MATCH pars.peaks.name = {Broad9_4, NADH2, NADH6, Trp}, because
% the quant loop indexes output.met.fit.areasConv positionally and areasConv
% is returned in pars.peaks.name order (an assert below enforces this).
met_T1_ms = [300, 300, 300, 300];   % Broad9_4, NADH2, NADH6, Trp
% T2 (ms) measured on the blood relaxometry series (_t2relax\run_t2relax_narrow
% -> refit_t2_dropTE10), narrow [8.9 10.3] fit, TE = 20..82 ms.  The TE=10 ms
% point is EXCLUDED: that scan carries extra broad short-T2 background (gone by
% TE=20) which the fixed-lineshape model has no component for, so the
% amplitude-only solve dumps it into the broadest peak and starves NADH6 --
% biasing NADH6 long (53 ms) and NADH4 short (20 ms).  With it dropped the four
% window/reference-scan variants agree to 1-4% and both R^2 and CI improve.
% Broad9_4 is a fitting-only shoulder that is never quantified -> its T2 is
% irrelevant; 30 ms is an arbitrary placeholder.
met_T2_ms = [ 30.0,  40.6,  35.8,  18.9];   % Broad9_4 (placeholder), NADH2, NADH6, Trp
met_nProt = [   1,    1,    1,    1];
peakNames = {'Broad9_4','NADH2','NADH6','Trp'};
nPeaks    = numel(peakNames);

% =====================================================================
% Process/Fit
% =====================================================================
[output, parsOut, hdr] = bru2mat_svsNAD(pars, wsDir, nwsDir);

% bru2mat_svsNAD rmpaths utils at exit; re-add for absoluteQuant_svs and
% any other utils helpers used downstream.
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

% Guard: the quant loop below indexes output.met.fit.areasConv positionally,
% so the relaxation/proton arrays + peakNames must be in the SAME order as the
% fitted components (pars.peaks.name).  Fail loudly if they ever drift apart.
assert(isequal(output.met.fit.names(:).', peakNames), ...
    ['peakNames/met_T1_ms/met_T2_ms order must match the fitted components ' ...
     '(output.met.fit.names = {%s})'], strjoin(output.met.fit.names, ', '));

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

% Write the same metabolite concentrations to a TSV alongside the .mat file.
tsvPath = fullfile(FinalOutDir, sprintf('%s_blood_1H_abs.tsv', scanID));
fid = fopen(tsvPath, 'w');
if ~isequal(fid,-1)
    fprintf(fid, 'metabolite\tconc_uM\n');
    for k = 1:nPeaks
        fprintf(fid, '%s\t%.4g\n', peakNames{k}, absQ.metabolites(k).conc_mM * 1e3);
    end
    fclose(fid);
    fprintf('Saved -> %s\n', tsvPath);
else
    fprintf('Could not save TSV');
end

% remove path
rmpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

end
