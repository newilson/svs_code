function processU01_Blood_1H(SessionDir, FinalOutDir, plt)
% PROCESSU01_BLOOD_1H  Absolute NAD+/Trp quantification for one U01 blood
%                      session acquired on the Bruker high-field NMR.
%
%   processU01_Blood_1H(SessionDir, FinalOutDir)
%   processU01_Blood_1H(SessionDir, FinalOutDir, plt)
%
% Each scan is a labeled Bruker experiment folder under SessionDir
% (fid + acqus):
%   selse_9p5ppm  downfield NAD/Trp  (carrier 9.5 ppm)
%   selse_4p7ppm  water reference    (carrier 4.7 ppm)
%   selse_0ppm    TSP reference      (carrier 0 ppm)
% 31p and selse_4p7ppm_T2relax are not used here.
%
% The downfield scan is processed/fit by bru2mat_svsNAD.m with the water
% scan as the Klose ECC + water reference.  The TSP scan goes through the
% same reader (Klose ECC against the water scan) and its singlet is fit with
% one Voigt + linear baseline.  Each downfield peak is converted to absolute
% concentration against both internal references:
%   water : absoluteQuant_svs.m (single whole-blood compartment)
%   TSP   : [M] = (S_M/S_TSP) * (N_TSP/N_M) * [TSP] * R_TSP/R_M
% with R = (1 - exp(-TR/T1)) * exp(-TE/T2).
%
% <stem> = bloodOutStem(SessionDir): the session folder name, with '_blood'
% appended unless it already ends in '_blood'.
%
% Requires <FinalOutDir>\<stem>_waterT2.mat from processU01_Blood_waterT2
% (session water T2, used for the water reference).
% The water model is bloodWatfitPars (shared with the T2 series).
%
% Output: <FinalOutDir>\<stem>_1H_abs.mat
%   (output, parsOut, hdr, absQ, absQtsp, tsp, waterT2, waterCheck, voxel, pars)
%         <FinalOutDir>\<stem>_1H_abs.tsv, one row:
%   dataset, WaterT2_ms, WaterT2_CI95_ms, WaterT2_R2, WaterRef_over_T2pred,
%   WaterT1_implied_ms, <peak>_uM_water..., <peak>_uM_TSP...

arguments
    SessionDir  (1,:) char
    FinalOutDir (1,:) char
    plt         (1,1) logical = false
end

assert(isfolder(SessionDir), 'SessionDir not found: %s', SessionDir);
if ~isfolder(FinalOutDir), mkdir(FinalOutDir); end

metDir = fullfile(SessionDir, 'selse_9p5ppm');
watDir = fullfile(SessionDir, 'selse_4p7ppm');
tspDir = fullfile(SessionDir, 'selse_0ppm');
assert(isfolder(metDir), 'Downfield scan folder not found: %s', metDir);
assert(isfolder(watDir), 'Water scan folder not found: %s', watDir);
assert(isfolder(tspDir), 'TSP scan folder not found: %s', tspDir);
fprintf('MET : %s\n', metDir);
fprintf('WAT : %s\n', watDir);
fprintf('TSP : %s\n', tspDir);

addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

[~, scanID] = fileparts(SessionDir);
stem   = bloodOutStem(SessionDir);
t2Path = fullfile(FinalOutDir, [stem '_waterT2.mat']);
assert(isfile(t2Path), ...
    'Water T2 result not found: %s (run processU01_Blood_waterT2 first)', t2Path);
load(t2Path, 'waterT2');
fprintf('Water T2 = %.2f +/- %.2f ms (%s)\n', waterT2.T2_ms, waterT2.T2_ci_ms, t2Path);

% =====================================================================
% Processing/Fitting parameters (shared recipe: utils/bloodMetFitPars.m)
% =====================================================================
pars         = bloodMetFitPars();
pars.plt     = plt;
pars.pltSpec = plt;

% =====================================================================
% Single blood compartment: water content + relaxation.
%   waterConc_M : whole blood is ~83% water by mass -> ~0.83 * 55.5 M ~= 46 M
%                 (0.8*55 = 44 M used below; still a rough placeholder)
%   T2_ms       : measured per session from selse_4p7ppm_T2relax
%                 (processU01_Blood_waterT2).
%   T1_ms       : selective saturation recovery, exp 3 (4.7 ppm) of
%                 FromFlip\aswain\rbnmr\data\Human_Test_Relaxometry, fit by
%                 _t2relax\t1relax_satrec.m -> t1_fit_results\: 2207 +/- 138 ms
%                 (R^2 0.9998).  DIFFERENT SAMPLE from any U01 session, and it
%                 disagrees with the ~1.0 s implied by waterCheck below (RG
%                 linearity suspect).  conc scales with (1-exp(-TR/T1_water)).
% =====================================================================
Blood = struct('name','Blood', 'T1_ms', 2207, 'T2_ms', waterT2.T2_ms, ...
               'waterConc_M', 0.8 * 55, 'flag_metSig',true);

% Downfield-proton relaxation / proton counts.  Placeholders for now (same
% for every peak); kept as explicit per-peak arrays because they will be
% replaced with measured per-peak T1/T2 later.  Trp is now a single broad
% component, quantified as one peak.
%
% ORDER MUST MATCH pars.peaks.name = {Broad9_4, NADH2, NADH6, Trp}, because
% the quant loop indexes output.met.fit.areasConv positionally and areasConv
% is returned in pars.peaks.name order (an assert below enforces this).
% T1 (ms) by selective saturation recovery, exp 5 (9.5 ppm) of
% FromFlip\aswain\rbnmr\data\Human_Test_Relaxometry, fit by
% _t2relax\t1relax_satrec.m -> t1_fit_results\ (fixed-lineshape arm, all
% R^2 >= 0.99).  Different sample from any U01 session.
met_T1_ms = [291, 402, 465, 386];   % Broad9_4, NADH2, NADH6, Trp
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
% TSP internal reference (selse_0ppm)
% =====================================================================
% T1_ms: selective saturation recovery, exp 4 (0 ppm) of
% FromFlip\aswain\rbnmr\data\Human_Test_Relaxometry, fit by
% _t2relax\t1relax_satrec.m -> t1_fit_results\: 1295 +/- 610 ms.
% conc_mM and T2_ms still to be supplied (T2 ~23.6 ms was measured on the
% rbnmr 0 ppm T2 series, _t2relax\t2relax_tspT2.m, but on yet another sample).
TSP = struct('conc_mM', NaN, 'nProt', 9, 'T1_ms', 1295, 'T2_ms', NaN);

% =====================================================================
% Process/Fit
% =====================================================================
[output, parsOut, hdr] = bru2mat_svsNAD(pars, metDir, watDir);

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

% =====================================================================
% TSP fit + TSP-referenced quant
% =====================================================================
tsp = fitTSP(tspDir, watDir, pars.tspfit, plt);
absQtsp = tspQuant(data.metabolites, seq, tsp, TSP);

fprintf('TSP at %.4f ppm (raw axis), area %.5g, water area %.5g\n', ...
    tsp.center, tsp.area, data.water.area);
fprintf('%-8s %12s %12s\n', 'peak', 'water [uM]', 'TSP [uM]');
for k = 1:nPeaks
    fprintf('%-8s %12.4g %12.4g\n', peakNames{k}, ...
        absQ.metabolites(k).conc_mM * 1e3, absQtsp.metabolites(k).conc_mM * 1e3);
end
if any(isnan([TSP.conc_mM TSP.T1_ms TSP.T2_ms]))
    warning('processU01_Blood_1H:tspConst', ...
        'TSP conc_mM/T1_ms/T2_ms not set; TSP-referenced concentrations are NaN.');
end

% =====================================================================
% Water reference scan vs water T2 series
% =====================================================================
waterCheck = struct();
waterCheck.T2_ms         = waterT2.T2_ms;
waterCheck.T2_ci_ms      = waterT2.T2_ci_ms;
waterCheck.R2            = waterT2.R2;
waterCheck.refArea       = data.water.area;
waterCheck.predArea      = waterT2.S0 * exp(-output.wat.TE_ms / waterT2.T2_ms);
waterCheck.refOverPred   = waterCheck.refArea / waterCheck.predArea;
waterCheck.T1_implied_ms = NaN;
if waterCheck.refOverPred > 0 && waterCheck.refOverPred < 1
    waterCheck.T1_implied_ms = -output.wat.TR_ms / log(1 - waterCheck.refOverPred);
end
fprintf(['Water ref area %.5g vs T2-series prediction %.5g at TE %.2f ms: ratio %.3f ' ...
         '(implied T1 %.0f ms at TR %.0f ms)\n'], waterCheck.refArea, waterCheck.predArea, ...
         output.wat.TE_ms, waterCheck.refOverPred, waterCheck.T1_implied_ms, output.wat.TR_ms);

% =====================================================================
% Save
% =====================================================================
outPath = fullfile(FinalOutDir, [stem '_1H_abs.mat']);
save(outPath, 'output','parsOut','hdr','absQ','absQtsp','tsp','waterT2','waterCheck','voxel','pars');
fprintf('Saved -> %s\n', outPath);

cols = {'WaterT2_ms','WaterT2_CI95_ms','WaterT2_R2','WaterRef_over_T2pred','WaterT1_implied_ms'};
vals = [waterCheck.T2_ms, waterCheck.T2_ci_ms, waterCheck.R2, ...
        waterCheck.refOverPred, waterCheck.T1_implied_ms];
for k = 1:nPeaks
    cols{end+1} = sprintf('%s_uM_water', peakNames{k}); %#ok<AGROW>
    vals(end+1) = absQ.metabolites(k).conc_mM * 1e3;    %#ok<AGROW>
end
for k = 1:nPeaks
    cols{end+1} = sprintf('%s_uM_TSP', peakNames{k});   %#ok<AGROW>
    vals(end+1) = absQtsp.metabolites(k).conc_mM * 1e3; %#ok<AGROW>
end

tsvPath = fullfile(FinalOutDir, [stem '_1H_abs.tsv']);
fid = fopen(tsvPath, 'w');
if ~isequal(fid,-1)
    fprintf(fid, 'dataset');
    fprintf(fid, '\t%s', cols{:});
    fprintf(fid, '\n%s', scanID);
    fprintf(fid, '\t%.4g', vals);
    fprintf(fid, '\n');
    fclose(fid);
    fprintf('Saved -> %s\n', tsvPath);
else
    fprintf('Could not save TSV');
end

% remove path
rmpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

end

% =====================================================================
function tsp = fitTSP(tspDir, watDir, opt, plt)
% Process the selse_0ppm scan through bru2mat_svsNAD (Klose ECC against the
% water scan, no metabolite fit) and fit the TSP singlet.

p = struct('lb_met',opt.lb, 'lb_wat',opt.lb, 'eccopt',0, 'PC',0, 'base',nan, ...
           'plt',false, 'pltSpec',false, 'dofit',false);
p.WatSupPre.opts.type  = 'none';
p.WatSupPost.opts.type = 'none';
[o, ~, hdrT] = bru2mat_svsNAD(p, tspDir, watDir);
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

f0  = hdrT.acqus.BF1;
out = fitBloodTSP(o.met.ppm, o.met.spec, opt, f0);

tsp = struct();
tsp.area     = out.area;
tsp.areaWin  = out.areaWin;
tsp.pars     = out.pars;
tsp.parnames = out.parnames{1};
tsp.center   = out.center;
tsp.located  = out.located;
tsp.window   = out.window;
tsp.ppm      = out.ppm;
tsp.spec     = out.spec;
tsp.fit      = out.fit;
tsp.baseline = out.baseline;
tsp.TR_ms    = o.met.TR_ms;
tsp.TE_ms    = o.met.TE_ms;
tsp.nAvg     = o.met.nAvg;
tsp.opt      = opt;

if plt
    figure('Name','1H TSP fit');
    plot(out.ppm, out.spec, 'k', out.ppm, out.fit, 'r', out.ppm, out.baseline, 'b--', ...
         out.ppm, out.spec - out.fit, 'g');
    legend({'data','fit','baseline','residual'});
    set(gca, 'xdir', 'reverse');
    title(sprintf('TSP fit: %.4f ppm, area %.4g', tsp.center, tsp.area));
end
end

% =====================================================================
function absQ = tspQuant(mets, seq, tsp, TSP)
R_T = (1 - exp(-tsp.TR_ms/TSP.T1_ms)) * exp(-tsp.TE_ms/TSP.T2_ms);
absQ.metabolites = repmat(struct('name','','conc_mM',NaN,'R_M',NaN,'ratio',NaN), numel(mets), 1);
for k = 1:numel(mets)
    mt = mets(k);
    Rm = (1 - exp(-seq.TR_ms/mt.T1_ms)) * exp(-seq.TE_ms/mt.T2_ms);
    ratio = mt.area / tsp.area;
    absQ.metabolites(k).name    = mt.name;
    absQ.metabolites(k).conc_mM = ratio * (TSP.nProt/mt.nProtons) * TSP.conc_mM * R_T / Rm;
    absQ.metabolites(k).R_M     = Rm;
    absQ.metabolites(k).ratio   = ratio;
end
absQ.ref = struct('name','TSP', 'area',tsp.area, 'R',R_T, 'conc_mM',TSP.conc_mM, ...
                  'nProt',TSP.nProt, 'T1_ms',TSP.T1_ms, 'T2_ms',TSP.T2_ms, ...
                  'TR_ms',tsp.TR_ms, 'TE_ms',tsp.TE_ms);
end
