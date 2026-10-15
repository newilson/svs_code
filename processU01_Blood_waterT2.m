function waterT2 = processU01_Blood_waterT2(SessionDir, FinalOutDir, plt)
% PROCESSU01_BLOOD_WATERT2  Water T2 for one U01 blood session from the
%                           selse_4p7ppm_T2relax spin-echo series.
%
%   waterT2 = processU01_Blood_waterT2(SessionDir, FinalOutDir)
%   waterT2 = processU01_Blood_waterT2(SessionDir, FinalOutDir, plt)
%
% Each row of the 2D ser is one selective spin echo (selgpse_t2_cwl) with
% both echo delays taken from the vd list, so TE = 2*vd + 2*p16 + p12 (the
% same TE definition as the single-TE scans in bru2mat_svsNAD).  Every echo
% is processed like the water reference in bru2mat_svsNAD (RG/NS
% normalization, Klose ECC against itself, water restored to 4.7 ppm, lb)
% and fit with the shared blood water model (bloodWatfitPars).  The summed
% component areas are fit to S0*exp(-TE/T2).
%
% Output: <FinalOutDir>\<stem>_waterT2.mat  (waterT2)
%         <FinalOutDir>\<stem>_waterT2.tsv  (per-echo areas)
% with <stem> = bloodOutStem(SessionDir).
% processU01_Blood_1H reads the .mat for the session water T2.

arguments
    SessionDir  (1,:) char
    FinalOutDir (1,:) char
    plt         (1,1) logical = false
end

assert(isfolder(SessionDir), 'SessionDir not found: %s', SessionDir);
if ~isfolder(FinalOutDir), mkdir(FinalOutDir); end
t2Dir = fullfile(SessionDir, 'selse_4p7ppm_T2relax');
assert(isfolder(t2Dir), 'Water T2 series folder not found: %s', t2Dir);
fprintf('T2  : %s\n', t2Dir);

thisPath = fileparts(mfilename('fullpath'));
addpath(fullfile(thisPath, 'utils'));
addpath(fullfile(thisPath, 'bruker'));

lb     = 6;
watfit = bloodWatfitPars();

% =====================================================================
% Read the echo series
% =====================================================================
[P, ser] = bruker_read(t2Dir);
ser = bruker_remove_digital_filter(P, ser);
nTE = P.acqu2s.TD;
if size(ser,1) == nTE && size(ser,2) ~= nTE
    ser = ser.';
end
assert(size(ser,2) == nTE, 'ser has %d traces, expected TD1 = %d', size(ser,2), nTE);

vd = readVdList(t2Dir, P.acqus.VDLIST);
assert(numel(vd) >= nTE, 'vd list has %d entries, need %d', numel(vd), nTE);
vd    = vd(1:nTE);
TE_ms = 2*vd*1e3 + 2*P.acqus.P(17)/1e3 + P.acqus.P(13)/1e3;

bw  = P.acqus.SW_h;
f0  = P.acqus.BF1;
O1  = P.acqus.O1;
ser = ser / P.acqus.RG / P.acqus.NS;
N   = size(ser,1);
t   = (0:N-1)/bw;
hz  = (-1/2:1/N:1/2-1/N)*bw;
ppm = (O1 + hz)/f0;
f_water_Hz = 4.7*f0 - O1;

% =====================================================================
% Per-echo water fit
% =====================================================================
area = zeros(nTE,1);
fits = struct('spec',{},'fit',{},'comp',{},'pars',{},'parnames',{},'areasConv',{});
xw   = [];
for k = 1:nTE
    spec = bloodProcWaterFID(ser(:,k), t, lb, f_water_Hz);
    out  = fitBloodWater(ppm, spec, watfit, f0);
    xw  = out.ppm;

    area(k) = out.area;
    fits(k) = struct('spec',out.spec, 'fit',out.fit, 'comp',out.comp, 'pars',out.pars, ...
                     'parnames',out.parnames, 'areasConv',out.areasConv);
end

% =====================================================================
% Mono-exponential T2
% =====================================================================
T = fitMonoExpT2(TE_ms, area);

waterT2 = struct();
waterT2.T2_ms    = T.T2;
waterT2.T2_ci_ms = T.T2_ci;
waterT2.S0       = T.S0;
waterT2.S0_ci    = T.S0_ci;
waterT2.R2       = T.R2;
waterT2.TE_ms    = TE_ms(:);
waterT2.area     = area;
waterT2.areaFit  = T.Sfit;
waterT2.ppm      = xw;
waterT2.fits     = fits;
waterT2.lb       = lb;
waterT2.watfit   = watfit;
waterT2.d1_ms    = P.acqus.D(2)*1e3;
waterT2.RG       = P.acqus.RG;
waterT2.NS       = P.acqus.NS;
waterT2.srcDir   = t2Dir;

fprintf('%8s %12s %12s %8s\n', 'TE [ms]', 'area', 'fit', 'resid %');
for k = 1:nTE
    fprintf('%8.2f %12.5g %12.5g %8.2f\n', TE_ms(k), area(k), T.Sfit(k), ...
        100*(area(k) - T.Sfit(k))/T.Sfit(k));
end
fprintf('Water T2 = %.2f +/- %.2f ms (95%% CI), R^2 = %.5f, n = %d\n', ...
    T.T2, T.T2_ci, T.R2, nTE);

% =====================================================================
% Save
% =====================================================================
stem    = bloodOutStem(SessionDir);
outPath = fullfile(FinalOutDir, [stem '_waterT2.mat']);
save(outPath, 'waterT2');
fprintf('Saved -> %s\n', outPath);

tsvPath = fullfile(FinalOutDir, [stem '_waterT2.tsv']);
fid = fopen(tsvPath, 'w');
if ~isequal(fid,-1)
    fprintf(fid, 'TE_ms\tarea\tarea_fit\tresid_pct\n');
    for k = 1:nTE
        fprintf(fid, '%.4g\t%.6g\t%.6g\t%.4g\n', TE_ms(k), area(k), T.Sfit(k), ...
            100*(area(k) - T.Sfit(k))/T.Sfit(k));
    end
    fclose(fid);
    fprintf('Saved -> %s\n', tsvPath);
else
    fprintf('Could not save TSV');
end

if plt
    figure('Name','Water T2 spectra');
    tiledlayout('flow');
    for k = 1:nTE
        nexttile;
        plot(xw, fits(k).spec, 'k', xw, fits(k).fit, 'r', xw, fits(k).spec - fits(k).fit, 'g');
        set(gca, 'xdir', 'reverse'); xlim([4.4 5.0]);
        title(sprintf('TE %.1f ms', TE_ms(k)));
    end

    figure('Name','Water T2 decay');
    te = linspace(0, max(TE_ms)*1.05, 200);
    subplot(2,1,1);
    semilogy(TE_ms, area, 'ko', te, T.S0*exp(-te/T.T2), 'r-');
    xlabel('TE [ms]'); ylabel('water area');
    title(sprintf('Water T2 = %.1f \\pm %.1f ms, R^2 = %.4f', T.T2, T.T2_ci, T.R2));
    subplot(2,1,2);
    plot(TE_ms, 100*(area - T.Sfit)./T.Sfit, 'ko-'); yline(0, ':');
    xlabel('TE [ms]'); ylabel('residual [%]');
end

rmpath(fullfile(thisPath, 'bruker'));
rmpath(fullfile(thisPath, 'utils'));

end

% =====================================================================
function vd = readVdList(expDir, name)
name = regexprep(char(name), '[<>]', '');
cand = {fullfile(expDir, 'lists', 'vd', name), fullfile(expDir, 'vdlist')};
cand = cand(cellfun(@isfile, cand));
assert(~isempty(cand), 'vd list %s not found under %s', name, expDir);
lines = strtrim(splitlines(fileread(cand{1})));
vd = [];
for i = 1:numel(lines)
    tok = regexp(lines{i}, '^([-+]?[\d.]+(?:[eE][-+]?\d+)?)\s*([munsMUNS]?)$', 'tokens', 'once');
    if isempty(tok), continue; end
    v = str2double(tok{1});
    switch lower(tok{2})
        case 'm', v = v*1e-3;
        case 'u', v = v*1e-6;
        case 'n', v = v*1e-9;
    end
    vd(end+1,1) = v; %#ok<AGROW>
end
end
