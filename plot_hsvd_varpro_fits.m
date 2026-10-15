function plot_hsvd_varpro_fits(opt)
%PLOT_HSVD_VARPRO_FITS  Overlay the HSVD and VARPRO fits, metabolite and water.
%
%   plot_hsvd_varpro_fits                          % calf
%   plot_hsvd_varpro_fits('cohort','u01')
%   plot_hsvd_varpro_fits('cohort','calf','sessions',{'S004_04162026'})
%
% Reads the outputs written by compare_hsvd_varpro and draws, for every
% session, the metabolite fit region(s) and the water reference with:
%
%   data       phased spectrum the fine fit saw (metabolite) / real(nws) (water)
%   VARPRO     region fit (peaks + spline baseline), baseline dashed
%   HSVD       all components, peak-assigned + unassigned; unassigned dashed
%
% The HSVD metabolite traces are the ones the coarse step stored
% (region.coarse.fitPhased / .baselinePhased, same phase frame as the data);
% the water HSVD model is rebuilt from the saved component list with
% utils/hsvdSpectra.
%
% Writes <outDir>\fits\<cohort>\<session>.png plus two overview grids,
% <cohort>_metabolite.png and <cohort>_water.png.

arguments
    opt.cohort   (1,:) char    = 'calf'
    opt.outDir   (1,:) char    = fullfile(fileparts(mfilename('fullpath')), '_hsvd_varpro')
    opt.sessions (1,:) cell    = {}
    opt.overview (1,1) logical = true
end

REPO = fileparts(mfilename('fullpath'));
addpath(fullfile(REPO,'utils'));

srcDir = fullfile(opt.outDir, lower(opt.cohort));
assert(isfolder(srcDir), 'no processed outputs in %s (run compare_hsvd_varpro first)', srcDir);
f = dir(fullfile(srcDir,'*_1H_abs.mat'));
assert(~isempty(f), 'no *_1H_abs.mat in %s', srcDir);

figDir = fullfile(opt.outDir,'fits', lower(opt.cohort));
if ~isfolder(figDir), mkdir(figDir); end

S = {};
for i = 1:numel(f)
    tag = regexprep(f(i).name, '_(brain|calf)_1H_abs\.mat$', '');
    if ~isempty(opt.sessions) && ~any(strcmp(tag, opt.sessions)), continue; end
    T = traces(fullfile(f(i).folder, f(i).name), tag);
    S{end+1} = T; %#ok<AGROW>
    png = fullfile(figDir, [tag '.png']);
    sessionFigure(T, png);
    fprintf('%-46s %s\n', tag, png);
end
assert(~isempty(S), 'no sessions matched');

if opt.overview
    p1 = fullfile(opt.outDir, sprintf('fits_%s_metabolite.png', lower(opt.cohort)));
    p2 = fullfile(opt.outDir, sprintf('fits_%s_water.png',      lower(opt.cohort)));
    overviewFigure(S, 'met', p1);
    overviewFigure(S, 'wat', p2);
    fprintf('\noverview -> %s\n         -> %s\n', p1, p2);
end
end

% =====================================================================
function T = traces(mat, tag)
L = load(mat,'output');
o = L.output;  F = o.met.fit;
T.tag = tag;
T.names = F.hsvd.names;

% ---- metabolite region(s) -------------------------------------------
ppmAll = o.met.ppm(:);
for r = 1:numel(F.region)
    Rg = F.region(r);
    x  = ppmAll(ppmAll > min(Rg.fitRange) & ppmAll < max(Rg.fitRange));
    [tf, loc] = ismember(x, F.ppm(:));
    assert(all(tf) && numel(x) == numel(Rg.coarse.fit), ...
        '%s: region %d axis does not line up with the stored traces', tag, r);
    T.met(r).name   = Rg.name;
    T.met(r).x      = x(:);
    T.met(r).y      = real(F.spec(loc(:)));
    T.met(r).y      = T.met(r).y(:);
    T.met(r).varpro = Rg.fit.fit(:);
    T.met(r).vbase  = Rg.fit.baseline(:);
    % A coarse-phase QC fallback (pars.fit.coarseQC) replaces the HSVD coarse
    % fit with the VARPRO one, which carries no component traces to draw.
    T.met(r).qcFallback = ~isfield(Rg.coarse,'baselinePhased');
    if T.met(r).qcFallback
        T.met(r).hbase = nan(size(x(:)));
        T.met(r).hsvd  = nan(size(x(:)));
    else
        T.met(r).hbase = Rg.coarse.baselinePhased(:);
        T.met(r).hsvd  = Rg.coarse.fitPhased(:) + T.met(r).hbase;
    end
    sel = Rg.peakSelect;
    if all(isnan(F.hsvd.area(sel)))
        T.met(r).ratio = NaN;          % HSVD assigned nothing in this region
    else
        T.met(r).ratio = sum(F.hsvd.area(sel),'omitnan') / sum(F.areasConv(sel),'omitnan');
    end
end

% ---- water ----------------------------------------------------------
W  = o.wat;
h  = W.fit.hsvd;
xw = W.fit.ppm(:);
[tf, loc] = ismember(xw, W.ppm(:));
assert(all(tf), '%s: water axis does not line up', tag);

N    = numel(W.ppm);
dppm = W.ppm(2) - W.ppm(1);
f0w  = abs((W.hz(2) - W.hz(1)) / dppm);
ppm0 = W.ppm(floor(N/2) + 1);
Sc   = hsvdSpectra((h.compPpm - ppm0)*f0w, h.compLwHz*pi, h.compAmp, N, 1/(N*f0w*abs(dppm)));
asg  = logical(h.assigned(:)');

T.wat.x      = xw;
T.wat.y      = real(W.fit.spec(:));
T.wat.varpro = W.fit.spec_fit(:);
T.wat.vbase  = W.fit.baseline(:);
T.wat.hsvd   = real(sum(Sc(loc,:), 2));
T.wat.hbase  = real(sum(Sc(loc,~asg), 2));
T.wat.ratio  = h.area / W.fit.areaTotal;
T.wat.nComp  = h.nComp;
end

% =====================================================================
function sessionFigure(T, png)
nR  = numel(T.met);
fig = figure('Visible','off','Position',[100 100 520*(nR+1) 380]);
for r = 1:nR
    subplot(1, nR+1, r);
    if T.met(r).qcFallback
        ttl = sprintf('%s  (coarse-phase QC fallback: no HSVD fit)', T.met(r).name);
    else
        ttl = sprintf('%s  (HSVD/VARPRO area %.2f)', T.met(r).name, T.met(r).ratio);
    end
    panel(T.met(r), ttl);
    if r == 1, ylabel('signal'); end
end
subplot(1, nR+1, nR+1);
panel(T.wat, sprintf('water, %d comps  (HSVD/VARPRO area %.2f)', T.wat.nComp, T.wat.ratio));
lg = legend({'data','VARPRO','VARPRO baseline','HSVD','HSVD unassigned'}, ...
            'Location','northeast');
lg.FontSize = 7;  lg.Box = 'off';
sgtitle(strrep(T.tag,'_','\_'));
print(fig, png, '-dpng', '-r150');
close(fig);
end

function overviewFigure(S, which, png)
n  = numel(S);
nc = min(4, n);
nr = ceil(n/nc);
fig = figure('Visible','off','Position',[50 50 460*nc 300*nr]);
for i = 1:n
    subplot(nr, nc, i);
    if strcmp(which,'met')
        P  = S{i}.met;                      % every fit region in one panel
        tt = strtrim(sprintf('%.2f ', [P.ratio]));
        if any([P.qcFallback]), tt = [tt ' (QC fallback)']; end
    else
        P  = S{i}.wat;
        tt = sprintf('%.2f, %d comps', P.ratio, P.nComp);
    end
    panel(P, '');
    title(sprintf('%s\nHSVD/VARPRO %s', strrep(shortTag(S{i}.tag),'_','\_'), tt), ...
          'FontSize', 8, 'FontWeight','normal');
end
sgtitle(sprintf('%s: data (black), VARPRO (red), HSVD (blue); residuals below', ...
        upper(which)));
print(fig, png, '-dpng', '-r150');
close(fig);
end

function panel(P, ttl)
% P is one fit region, or a struct array of them drawn in the same axes.
yAll = vertcat(P.y);
xAll = vertcat(P.x);
span = max(yAll) - min(yAll);
if ~isfinite(span) || span <= 0, span = 1; end
off = min(yAll) - 0.30*span;                   % residuals sit under the data
hold on
for q = 1:numel(P)
    y = P(q).y;
    plot(P(q).x, y, 'k', 'LineWidth', 1);
    plot(P(q).x, P(q).varpro, 'r', 'LineWidth', 1);
    plot(P(q).x, P(q).vbase,  'r--', 'LineWidth', 0.5);
    plot(P(q).x, P(q).hsvd,   'b',  'LineWidth', 1);
    plot(P(q).x, P(q).hbase,  'b--', 'LineWidth', 0.5);
    plot(P(q).x, off + (y - P(q).varpro), 'r', 'LineWidth', 0.5);
    plot(P(q).x, off + (y - P(q).hsvd),   'b', 'LineWidth', 0.5);
end
plot([min(xAll) max(xAll)], [off off], 'k:', 'LineWidth', 0.5);
hold off
set(gca,'xdir','reverse');
xlim([min(xAll) max(xAll)]);
ylim([off - 0.35*span, max(yAll) + 0.15*span]);
xlabel('ppm');
grid on
if ~isempty(ttl), title(ttl, 'FontSize', 9); end
end

function s = shortTag(tag)
s = regexprep(tag, '_31P-1H-HeadCoil.*$', '');
s = regexprep(s, '_brain$|_calf$', '');
end
