function compareSubjectResultsU01(subject, dataRoot, outDir)
% COMPARESUBJECTRESULTSU01  Overlay the two U01 visits of one subject.
%
%   compareSubjectResultsU01(1)
%   compareSubjectResultsU01('S01')
%   compareSubjectResultsU01(1, dataRoot, outDir)
%
% Finds the subject's sessions under <dataRoot>\Brain, \Calf and \Blood
% (folders named <study>_S01_<mmddyyyy>_<anatomy>), loads their processed
% *_abs.mat outputs and, per anatomy, writes to outDir (default
% <dataRoot>\SummaryResults\S01):
%
%   S01_<anatomy>_table.png/.fig  both visits and the percent change
%   S01_<anatomy>_v<k>_<date>_1H_fine_fit_region<n>_<name>.png/.fig
%   S01_<anatomy>_v<k>_<date>_31P_fit.png/.fig
%                                 copies of visit k's newest 1H fine-fit figure
%                                 for each fit region and newest 31P fit
%                                 figure, from its processed folder
%   S01_<anatomy>.png/.fig        water spectrum, 1H downfield spectrum
%                                 (divided by the water integral), full 31P
%                                 spectrum and 31P zoomed on aATP/NAD/UDP(G), written
%                                 by the Save button of
%                                 compareSubjectResultsU01_GUI, where each
%                                 spectrum's phase and scale can be adjusted
%
% appendPngsToPptx(fullfile(outDir, 'S01 Results.pptx'), fullfile(outDir, 'S01_*.png'))
% then puts them on slides for brain, calf, then blood, each as table, fits
% (1H regions in order, then 31P; both visits of a fit on one slide), then
% overlay.
%
% The three GUIs open together; press Save in each when ready.  First
% visit blue, second red.  With more than two sessions the earliest and
% latest are compared.  A missing output leaves its panel empty and its table
% entries blank.
%
% Table rows
%   brain 1H : NADH2, NADH4, NADH6, NAD+ mean(H2,H4,H6), NAD+ mean(H2,H6), Trp
%   calf 1H  : NADH2, NADH6, NAD+ mean(H2,H6), NR, Trp, Peak10_3
%   blood 1H : NADH2, NADH6, NAD+ mean(H2,H6), Trp
%   all 1H   : water integral, water FWHM (total Voigt width of the water fit)
%   brain/calf 31P : NAD+, NADH, PCr, aATP against the external phantom, the
%              phantom area, the phantom B1+ relative to PCr from the B1
%              calibration, NAD+ and PCr against aATP, and NAD+ and aATP
%              against PCr (each internal reference labeled with its assumed
%              concentration)
%   blood 31P: NAD+, NADH against aATP (only if a blood 31P output exists)

sid = subjectID(subject);
if nargin < 2 || isempty(dataRoot)
    dataRoot = 'C:\Users\CAMIPM-NW\Documents\859875_U01';
end
if nargin < 3 || isempty(outDir)
    outDir = fullfile(dataRoot, 'SummaryResults', sid);
end
assert(isfolder(dataRoot), 'dataRoot not found: %s', dataRoot);
if ~isfolder(outDir), mkdir(outDir); end
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

anatomies = {'Brain', 'Calf', 'Blood'};
for a = 1:numel(anatomies)
    ses = findSessions(fullfile(dataRoot, anatomies{a}), sid);
    if isempty(ses)
        fprintf('%s %s: no sessions found\n', sid, anatomies{a});
        continue
    end
    if numel(ses) > 2
        warning('%s %s: %d sessions, comparing the earliest and latest', sid, anatomies{a}, numel(ses));
        ses = ses([1 end]);
    end
    anat = lower(anatomies{a});
    D = cellfun(@(s) loadSession(s, anat), {ses.dir});
    for k = 1:numel(D)
        fprintf('%s %-5s %s  1H:%d 31P:%d\n', sid, anat, ses(k).name, ~isempty(D(k).h1), ~isempty(D(k).p31));
    end
    labels = arrayfun(@(s) sprintf('Visit %d (%s)', s.visit, char(s.date, 'yyyy-MM-dd')), ...
        ses(:)', 'UniformOutput', false);

    rows = tableRows(D, anat);
    fig = drawTable(rows, labels, sprintf('%s %s', sid, anat));
    saveFig(fig, fullfile(outDir, sprintf('%s_%s_table', sid, anat)));
    copyFitFigures(ses, sid, anat, outDir);

    compareSubjectResultsU01_GUI(displaySpectra(D), labels, sprintf('%s %s', sid, anat), ...
        fullfile(outDir, sprintf('%s_%s', sid, anat)));
end
end

% =====================================================================
function sid = subjectID(subject)
if isnumeric(subject)
    sid = sprintf('S%02d', subject);
else
    tok = regexp(char(subject), '^[Ss]?(\d+)$', 'tokens', 'once');
    assert(~isempty(tok), 'subject must be a number or ''S<number>'': %s', char(subject));
    sid = sprintf('S%02d', str2double(tok{1}));
end
end

function ses = findSessions(anatDir, sid)
ses = struct('dir', {}, 'name', {}, 'date', {}, 'visit', {});
if ~isfolder(anatDir), return; end
d = dir(anatDir);
d = d([d.isdir] & ~startsWith({d.name}, '.'));
for k = 1:numel(d)
    tok = regexp(d(k).name, ['_' sid '_(\d{8})_'], 'tokens', 'once');
    if isempty(tok), continue; end
    ses(end+1) = struct('dir', fullfile(anatDir, d(k).name), 'name', d(k).name, ...
        'date', datetime(tok{1}, 'InputFormat', 'MMddyyyy'), 'visit', 0); %#ok<AGROW>
end
[~, i] = sort([ses.date]);
ses = ses(i);
for k = 1:numel(ses), ses(k).visit = k; end
end

% =====================================================================
function D = loadSession(sesDir, anat)
[~, scanID] = fileparts(sesDir);
stem = scanID;
if ~endsWith(scanID, ['_' anat]), stem = [scanID '_' anat]; end
procDir = fullfile(sesDir, 'processed');
D = struct('h1', [], 'p31', []);
D.h1  = loadFirst(procDir, {[scanID '_' anat '_1H_abs.mat'],  [stem '_1H_abs.mat']},  {'output', 'absQ'});
D.p31 = loadFirst(procDir, {[scanID '_' anat '_31P_abs.mat'], [stem '_31P_abs.mat']}, {'output', 'absQ'});
end

function copyFitFigures(ses, sid, anat, outDir)
% newest 1H fine-fit figure per fit region and newest 31P fit figure
% (.png and .fig) of each session, from its processed folder
kinds = {'_1H_fine_fit_-_region_*.png', '_1H_fine_fit_-_region_(\d+)_(.+)_\d+\.png$', ...
            @(t) sprintf('1H_fine_fit_region%s_%s', t{1}, t{2}); ...
         '_31P_basisVarpro_fit_*.png', '_31P_basisVarpro_fit_(\d+)\.png$', @(t) '31P_fit'};
for k = 1:numel(ses)
    procDir = fullfile(ses(k).dir, 'processed');
    pre = fullfile(outDir, sprintf('%s_%s_v%d_%s_', sid, anat, ses(k).visit, char(ses(k).date, 'yyyy-MM-dd')));
    nCopied = 0;
    for q = 1:size(kinds, 1)
        d = dir(fullfile(procDir, [ses(k).name kinds{q,1}]));
        tok = regexp({d.name}, kinds{q,2}, 'tokens', 'once');
        keep = ~cellfun(@isempty, tok);
        d = d(keep);
        names = cellfun(kinds{q,3}, tok(keep), 'UniformOutput', false);
        u = unique(names);
        for r = 1:numel(u)
            i = find(strcmp(names, u{r}));
            [~, j] = max([d(i).datenum]);
            src = fullfile(procDir, d(i(j)).name);
            copyOne(src, [pre u{r} '.png']);
            if isfile([src(1:end-4) '.fig']), copyOne([src(1:end-4) '.fig'], [pre u{r} '.fig']); end
        end
        nCopied = nCopied + numel(u);
    end
    fprintf('%s %-5s v%d: copied %d fit figure(s)\n', sid, anat, ses(k).visit, nCopied);
end
end

function copyOne(src, dst)
[ok, msg] = copyfile(src, dst, 'f');
if ~ok, warning('Could not copy %s to %s: %s', src, dst, msg); end
end

function m = loadFirst(procDir, names, vars)
m = [];
for k = 1:numel(names)
    f = fullfile(procDir, names{k});
    if isfile(f)
        m = load(f, vars{:});
        return
    end
end
end

% =====================================================================
function S = displaySpectra(D)
% complex spectra for the GUI, each cut to its display range: 1H downfield
% (fit phase, divided by the water integral, zero-filled), 1H water
% (zero-filled, with the water fit area) and 31P (acquisition-delay phase
% removed as in the fit, with the phantom area), each with a default PC1
% pivot (downfield: the fit's coarse-phase pivot, water 4.7 ppm, 31P 0 ppm)
S = repmat(struct('df', [], 'wat', [], 'p31', []), 1, numel(D));
for k = 1:numel(D)
    if ~isempty(D(k).h1)
        M = D(k).h1.output.met;
        [rng, band] = downfieldRange(M);
        piv = mean(band);
        if isfield(M.fit, 'phaseCorr'), piv = M.fit.phaseCorr.pivot; end
        [s, ppm] = zeroFill(phasedDownfield(M), M.ppm);
        i = ppm >= rng(1) & ppm <= rng(2);
        S(k).df = struct('ppm', ppm(i)', 'spec', s(i) / D(k).h1.output.wat.fit.areaTotal, 'band', band, 'range', rng, ...
            'pivot', piv);
        [s, ppm] = zeroFill(D(k).h1.output.wat.spec, D(k).h1.output.wat.ppm);
        i = ppm >= 4.2 & ppm <= 5.2;
        S(k).wat = struct('ppm', ppm(i)', 'spec', s(i), 'area', waterArea(D(k)), 'pivot', 4.7);
    end
    if ~isempty(D(k).p31)
        O = D(k).p31.output;
        s = O.spec(:, 1);
        if isfield(O, 'basisVarproResults') && isfield(O.basisVarproResults.timeInfo, 'delayPhase')
            s = s .* exp(-1i*O.basisVarproResults.timeInfo.delayPhase(O.ppm(:)));
        end
        i = O.ppm >= -20 & O.ppm <= 15;
        S(k).p31 = struct('ppm', O.ppm(i)', 'spec', s(i), 'area', phantomArea(D(k)), 'pivot', 0);
    end
end
end

function [rng, band] = downfieldRange(M)
% band: the excitation band, else the fit window; rng: band plus 0.3 ppm
if isfield(M, 'exc_freq') && isfield(M, 'exc_bw')
    band = M.exc_freq + [-1 1]*M.exc_bw/2;
else
    band = [min(M.fit.ppm) max(M.fit.ppm)];
end
rng = band + [-0.3 0.3];
end

function [s, ppm] = zeroFill(s, ppm)
% sinc-interpolate a 1H spectrum to at least 8192 points for display
N = numel(s);
nZF = max(N, 8192);
if nZF == N, return; end
fid = ifft(ifftshift(s(:)));
s = fftshift(fft(fid, nZF));
c0 = floor(N/2) + 1;
ppm = ppm(c0) + ((0:nZF-1) - floor(nZF/2)) * (ppm(2) - ppm(1)) * N/nZF;
end

function s = phasedDownfield(M)
% The fitted spectrum with the coarse phase the fit used.  The sign
% convention differs between coarse estimators and spec_full is stored
% before or after the phase depending on the region mode, so pick whichever
% candidate reproduces the stored fit data output.met.fit.spec.
F = M.fit;
s = F.spec_full(:);
if ~isfield(F, 'phaseCorr'), return; end
pc = [F.phaseCorr.pc0 F.phaseCorr.pc1];
[tf, loc] = ismember(F.ppm, M.ppm);
y = F.spec(:);
ok = tf(:) & y ~= 0;
cand = {s, shiftSpectrumPhase(s, pc, M.ppm, F.phaseCorr.pivot), shiftSpectrumPhase(s, -pc, M.ppm, F.phaseCorr.pivot)};
err = cellfun(@(c) norm(real(c(loc(ok))) - y(ok)), cand);
[e, i] = min(err);
if e > 1e-6*norm(y(ok))
    warning('downfield phase not reproduced (relative error %.2g)', e/norm(y(ok)));
end
s = cand{i};
end

% =====================================================================
function rows = tableRows(D, anat)
rows = struct('label', {}, 'v', {}, 'fmt', {}, 'header', {});
rows = addHeader(rows, '1H (water reference)');
switch anat
    case 'brain'
        rows = addRow(rows, 'NADH2 [\muM]',           arrayfun(@(d) h1Conc(d, {'NADH2'}), D), '%.0f');
        rows = addRow(rows, 'NADH4 [\muM]',           arrayfun(@(d) h1Conc(d, {'NADH4'}), D), '%.0f');
        rows = addRow(rows, 'NADH6 [\muM]',           arrayfun(@(d) h1Conc(d, {'NADH6'}), D), '%.0f');
        rows = addRow(rows, 'NAD^+ H2/H4/H6 mean [\muM]', arrayfun(@(d) h1Conc(d, {'NADH2','NADH4','NADH6'}), D), '%.0f');
        rows = addRow(rows, 'NAD^+ H2/H6 mean [\muM]',    arrayfun(@(d) h1Conc(d, {'NADH2','NADH6'}), D), '%.0f');
        rows = addRow(rows, 'Trp [\muM]',             arrayfun(@(d) h1Conc(d, {'Trp'}), D), '%.0f');
    case 'calf'
        rows = addRow(rows, 'NADH2 [\muM]',           arrayfun(@(d) h1Conc(d, {'NADH2'}), D), '%.0f');
        rows = addRow(rows, 'NADH6 [\muM]',           arrayfun(@(d) h1Conc(d, {'NADH6'}), D), '%.0f');
        rows = addRow(rows, 'NAD^+ H2/H6 mean [\muM]', arrayfun(@(d) h1Conc(d, {'NADH2','NADH6'}), D), '%.0f');
        rows = addRow(rows, 'NR [\muM]',              arrayfun(@(d) h1Conc(d, {'NR'}), D), '%.0f');
        rows = addRow(rows, 'Trp [\muM]',             arrayfun(@(d) h1Conc(d, {'Trp'}), D), '%.0f');
        rows = addRow(rows, '10.3 ppm peak [\muM]',   arrayfun(@(d) h1Conc(d, {'Peak10_3'}), D), '%.0f');
    case 'blood'
        rows = addRow(rows, 'NADH2 [\muM]',           arrayfun(@(d) h1Conc(d, {'NADH2'}), D), '%.0f');
        rows = addRow(rows, 'NADH6 [\muM]',           arrayfun(@(d) h1Conc(d, {'NADH6'}), D), '%.0f');
        rows = addRow(rows, 'NAD^+ H2/H6 mean [\muM]', arrayfun(@(d) h1Conc(d, {'NADH2','NADH6'}), D), '%.0f');
        rows = addRow(rows, 'Trp [\muM]',             arrayfun(@(d) h1Conc(d, {'Trp'}), D), '%.0f');
end
rows = addRow(rows, 'Water integral [a.u.]', arrayfun(@waterArea, D), '%.4g');
rows = addRow(rows, 'Water FWHM [Hz]',       arrayfun(@waterFWHM, D), '%.1f');

if strcmp(anat, 'blood')
    if all(arrayfun(@(d) isempty(d.p31), D)), return; end
    rows = addHeader(rows, ['31P (\alpha-ATP reference' refConcStr(D, 'primary', ', ') ')']);
    rows = addRow(rows, 'NAD^+ [mM]', arrayfun(@(d) p31Conc(d, 'primary', 'NADplus'), D), '%.3g');
    rows = addRow(rows, 'NADH [mM]',  arrayfun(@(d) p31Conc(d, 'primary', 'NADH'), D), '%.3g');
    return
end
rows = addHeader(rows, '31P (external phantom reference)');
rows = addRow(rows, 'NAD^+ [mM]',   arrayfun(@(d) p31Conc(d, 'external', 'NADplus'), D), '%.3g');
rows = addRow(rows, 'NADH [mM]',    arrayfun(@(d) p31Conc(d, 'external', 'NADH'), D), '%.3g');
rows = addRow(rows, 'PCr [mM]',     arrayfun(@(d) p31Conc(d, 'external', 'PCr'), D), '%.3g');
rows = addRow(rows, '\alpha-ATP [mM]', arrayfun(@(d) p31Conc(d, 'external', 'aATP'), D), '%.3g');
rows = addRow(rows, 'External phantom area [a.u.]', arrayfun(@phantomArea, D), '%.4g');
rows = addRow(rows, 'Phantom B1^+ / PCr B1^+ (B1 calibration)', arrayfun(@phantomRelB1, D), '%.3f');
rows = addHeader(rows, '31P (internal references)');
aRef = ['\alpha-ATP ref' refConcStr(D, 'aATP', ' ')];
pRef = ['PCr ref' refConcStr(D, 'PCr', ' ')];
rows = addRow(rows, ['NAD^+ [mM], ' aRef],     arrayfun(@(d) p31Conc(d, 'aATP', 'NADplus'), D), '%.3g');
rows = addRow(rows, ['PCr [mM], ' aRef],       arrayfun(@(d) p31Conc(d, 'aATP', 'PCr'), D), '%.3g');
rows = addRow(rows, ['NAD^+ [mM], ' pRef],     arrayfun(@(d) p31Conc(d, 'PCr', 'NADplus'), D), '%.3g');
rows = addRow(rows, ['\alpha-ATP [mM], ' pRef], arrayfun(@(d) p31Conc(d, 'PCr', 'aATP'), D), '%.3g');
end

function s = refConcStr(D, ref, sep)
% assumed concentration of an internal 31P reference ('aATP', 'PCr' or
% 'primary'), e.g. '2.7 mM', prefixed by sep; '' when no session has it
c = [];
for k = 1:numel(D)
    if isempty(D(k).p31), continue; end
    Q = D(k).p31.absQ;
    if strcmp(ref, 'primary')
        if isfield(Q, 'reference'), c(end+1) = Q.reference.conc_mM; end %#ok<AGROW>
    elseif isfield(Q, 'byRef') && isfield(Q.byRef, ref)
        c(end+1) = Q.byRef.(ref).reference.conc_mM; %#ok<AGROW>
    end
end
c = unique(c(isfinite(c)));
s = '';
if isempty(c), return; end
s = [sep strjoin(arrayfun(@(x) sprintf('%g', x), c, 'UniformOutput', false), '/') ' mM'];
end

function rows = addHeader(rows, label)
rows(end+1) = struct('label', label, 'v', [], 'fmt', '', 'header', true);
end

function rows = addRow(rows, label, v, fmt)
rows(end+1) = struct('label', label, 'v', v(:)', 'fmt', fmt, 'header', false);
end

function v = h1Conc(d, names)
% mean over the named peaks, in uM
v = NaN;
if isempty(d.h1), return; end
M = d.h1.absQ.metabolites;
c = arrayfun(@(n) concOf(M, n{1}), names);
v = mean(c) * 1e3;
end

function v = p31Conc(d, ref, name)
% ref: 'external', 'aATP', 'PCr' or 'primary' (absQ.metabolites)
v = NaN;
if isempty(d.p31), return; end
Q = d.p31.absQ;
switch ref
    case 'primary'
        M = Q.metabolites;
    case 'external'
        if ~isfield(Q, 'external') || ~isfield(Q.external, 'ok') || ~Q.external.ok, return; end
        M = Q.external.metabolites;
    otherwise
        if ~isfield(Q, 'byRef') || ~isfield(Q.byRef, ref), return; end
        M = Q.byRef.(ref).metabolites;
end
v = concOf(M, name);
end

function c = concOf(M, name)
k = find(strcmp({M.name}, name), 1);
if isempty(k), c = NaN; else, c = M(k).conc_mM; end
end

function v = phantomArea(d)
v = NaN;
if isempty(d.p31) || ~isfield(d.p31.absQ, 'external') || ~isfield(d.p31.absQ.external, 'areaRef'), return; end
v = d.p31.absQ.external.areaRef;
end

function v = phantomRelB1(d)
% external phantom relB1 over tissue (PCr) relB1, both from the 31P B1CALIB
% flip series
v = NaN;
if isempty(d.p31) || ~isfield(d.p31.absQ, 'external'), return; end
E = d.p31.absQ.external;
if isfield(E, 'relB1_ref_calib') && isfield(E, 'relB1_M_calib')
    v = E.relB1_ref_calib / E.relB1_M_calib;
end
end

function v = waterArea(d)
v = NaN;
if isempty(d.h1), return; end
v = d.h1.output.wat.fit.areaTotal;
end

function v = waterFWHM(d)
% total Voigt FWHM of a one-component water fit; with more components, the
% FWHM of their sum read off the fit grid
v = NaN;
if isempty(d.h1), return; end
W = d.h1.output.wat;
F = W.fit;
f0 = abs((W.hz(2) - W.hz(1)) / (W.ppm(2) - W.ppm(1)));
n = numel(F.names);
k = find(strcmp(F.parnames, 'width'), 1);
if n == 1 && ~isempty(k)
    v = F.pars(k) * f0;
    return
end
y = sum(F.comp, 2);
x = F.ppm(:);
[ym, im] = max(y);
half = ym / 2;
iL = find(y(1:im) < half, 1, 'last');
iR = im - 1 + find(y(im:end) < half, 1, 'first');
if isempty(iL) || isempty(iR), return; end
xL = interp1(y([iL iL+1]), x([iL iL+1]), half);
xR = interp1(y([iR-1 iR]), x([iR-1 iR]), half);
v = abs(xR - xL) * f0;
end

% =====================================================================
function fig = drawTable(rows, labels, ttl)
nR = numel(rows);
rowH = 0.32;
figH = 70 + 30*(nR + 1);
fig = figure('Name', [ttl ' table'], 'Color', 'w', 'Position', [100 100 760 figH]);
ax = axes(fig, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 10], 'YLim', [-(nR + 1)*rowH - 0.25 0.6]);
hold(ax, 'on');
xC = [0.15 5.6 7.5 9.85];
text(ax, 5, 0.3, ttl, 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'none');
hdr = [{'Measure'}, labels, {'Change [%]'}];
if numel(labels) == 1, hdr = [{'Measure'}, labels, {''}, {'Change [%]'}]; end
al = {'left', 'right', 'right', 'right'};
y = -0.5*rowH;
for c = 1:4
    text(ax, xC(c), y, hdr{c}, 'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', al{c});
end
plot(ax, [0.05 9.95], [y y] - 0.5*rowH, 'k', 'LineWidth', 1);
for r = 1:nR
    y = -(r + 0.5)*rowH;
    R = rows(r);
    if R.header
        patch(ax, [0.05 9.95 9.95 0.05], y + 0.5*rowH*[-1 -1 1 1], [0.9 0.92 0.96], 'EdgeColor', 'none');
        text(ax, xC(1), y, R.label, 'FontWeight', 'bold', 'FontSize', 10);
        continue
    end
    text(ax, xC(1) + 0.15, y, R.label, 'FontSize', 10);
    v = [R.v NaN(1, 2 - numel(R.v))];
    for c = 1:2
        text(ax, xC(c + 1), y, numStr(v(c), R.fmt), 'FontSize', 10, 'HorizontalAlignment', 'right');
    end
    text(ax, xC(4), y, numStr(100*(v(2) - v(1))/v(1), '%+.1f'), 'FontSize', 10, 'HorizontalAlignment', 'right');
end
plot(ax, [0.05 9.95], -(nR + 1)*rowH*[1 1], 'k', 'LineWidth', 1);
end

function s = numStr(v, fmt)
if ~isfinite(v), s = '–'; else, s = sprintf(fmt, v); end
end

function saveFig(fig, base)
try
    savefig(fig, [base '.fig']);
    exportgraphics(fig, [base '.png'], 'Resolution', 300);
catch ME
    warning('Could not save %s (is it open in another program?): %s', base, ME.message);
end
end
