function outFiles = plotU01_2026_fits(nucleus, anat, dataRoot, onlyLatest)
% PLOTU01_2026_FITS  Save a data/fit/residual PNG next to every processed scan.
%
%   plotU01_2026_fits('1H','brain')
%   plotU01_2026_fits('31P','calf')
%   plotU01_2026_fits(nucleus, anat, dataRoot, onlyLatest)
%
% Writes <session>\processed\<scanID>_<anat>_<nucleus>_fit.png.
%
% Two rows: data with the total fit overlaid, and the residual beneath.  No
% baseline or per-component traces -- those live in the fit structs if needed.
% Nothing is refitted; everything is read back out of the stored _abs.mat.
%
%   1H  -- one COLUMN per fitted region.  Brain is fitted in two regions
%          (NAD 8.7-10.0 ppm, Trp 9.8-10.45) that OVERLAP, and the stored
%          union-window model arrays overwrite each other where they do, so
%          each region has to be drawn on its own axis.  Calf is a single
%          region (8.95-10.7) and gets one column.
%   31P -- a single column zoomed on the aATP / NAD / UDP region.  The fit
%          window runs +-19 ppm, which compresses these lines to nothing.
%
%   onlyLatest  true (default) keeps latest-protocol sessions only
%   dataRoot    a folder of sessions, or a single session (has processed\)

if nargin < 1 || isempty(nucleus),   nucleus = '1H';   end
if nargin < 2 || isempty(anat),      anat    = 'brain'; end
if nargin < 3 || isempty(dataRoot)
    dataRoot = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026\2026';
end
if nargin < 4 || isempty(onlyLatest), onlyLatest = true; end
nucleus = upper(char(nucleus));  anat = lower(char(anat));

ZOOM31P = [-11.5 -5.5];    % ppm, PCr-referenced: covers aATP, NAD, all UDP

% ---- palette ----------------------------------------------------------
C = @(h) reshape(sscanf(h(2:end),'%2x'),1,3)/255;
surface = C('#fcfcfb');
gridClr = C('#e1e0d9');
muted   = C('#898781');
inkSec  = C('#52514e');
inkPri  = C('#0b0b0b');
BLUE    = C('#2a78d6');

sessions = findSessions(dataRoot);
outFiles = {};

for k = 1:numel(sessions)
    s = sessions{k};
    [~, scanID] = fileparts(s);
    f = fullfile(s, 'processed', sprintf('%s_%s_%s_abs.mat', scanID, anat, nucleus));
    if ~isfile(f), continue; end

    try
        [panels, ok] = readPanels(f, nucleus, onlyLatest, ZOOM31P);
    catch ME
        fprintf('  skip %s: %s\n', scanID, ME.message);  continue
    end
    if ~ok, fprintf('  skip (previous protocol): %s\n', scanID); continue; end

    nC = numel(panels);
    fig = figure('Color', surface, 'Units','pixels', ...
                 'Position',[80 80 460*nC 560], 'InvertHardcopy','off', 'Visible','off');
    tl = tiledlayout(fig, 3, nC, 'TileSpacing','compact', 'Padding','compact');

    for r = 1:nC
        P = panels(r);

        % ---- data + fit (2 of the 3 tile rows, so it dominates)
        ax = nexttile(tl, r, [2 1]);
        prep(ax, surface, gridClr, inkSec);
        hd = plot(ax, P.x, P.y,   '-', 'Color', inkPri, 'LineWidth', 1.2);
        hf = plot(ax, P.x, P.fit, '-', 'Color', BLUE,   'LineWidth', 2.0);
        title(ax, P.title, 'Color', inkPri, 'FontSize', 11, 'FontWeight','bold');
        if r == 1
            lg = legend(ax, [hd hf], {'data','fit'}, 'Box','off', ...
                        'TextColor', inkSec, 'FontSize', 9, 'Location','northwest');
            lg.ItemTokenSize = [14 8];
        end

        % ---- residual
        ax = nexttile(tl, 2*nC + r);
        prep(ax, surface, gridClr, inkSec);
        res = P.y - P.fit;
        plot(ax, [min(P.x) max(P.x)], [0 0], '-', 'Color', gridClr, 'LineWidth', 0.75);
        plot(ax, P.x, res, '-', 'Color', muted, 'LineWidth', 1.2);
        xlabel(ax, 'ppm', 'Color', inkSec, 'FontSize', 10);
        text(ax, 0.02, 0.86, sprintf('residual  rms %.3g', sqrt(mean(res.^2))), ...
             'Units','normalized', 'Color', inkSec, 'FontSize', 9);
    end

    ylabel(tl, 'signal (a.u.)', 'Color', inkSec, 'FontSize', 10);

    outPng = fullfile(s, 'processed', sprintf('%s_%s_%s_fit.png', scanID, anat, nucleus));
    exportgraphics(fig, outPng, 'Resolution', 180, 'BackgroundColor', surface);
    close(fig);
    outFiles{end+1} = outPng; %#ok<AGROW>
    fprintf('  wrote %s\n', outPng);
end

fprintf('\n%d %s %s fit figure(s) written.\n', numel(outFiles), anat, nucleus);
end

% =====================================================================
function [panels, ok] = readPanels(f, nucleus, onlyLatest, zoom31P)
panels = struct('x',{},'y',{},'fit',{},'title',{});
ok = true;

if strcmp(nucleus,'1H')
    m   = load(f, 'output');
    met = m.output.met;
    fit = met.fit;
    if onlyLatest && ~contains(lower(char(getf(met,'seqname',''))), 'svssel_20260101')
        ok = false; return
    end
    ppmU = fit.ppm(:);
    for r = 1:numel(fit.region)
        R  = fit.region(r);
        fr = R.fitRange;
        % Same predicate that built the region axis during fitting, applied to
        % the union window -- verified to reproduce the model length exactly.
        sel = ppmU > min(fr) & ppmU < max(fr);
        yf  = real(R.fit.fit(:));
        assert(sum(sel) == numel(yf), 'region %d: %d axis vs %d model points', ...
               r, sum(sel), numel(yf));
        panels(end+1) = struct('x', ppmU(sel), 'y', real(fit.spec(sel)), ...
            'fit', yf, 'title', sprintf('%s  [%.2f-%.2f ppm]', R.name, min(fr), max(fr))); %#ok<AGROW>
    end
else
    m = load(f, 'output', 'absQ');
    r = m.output.basisVarproResults;
    if onlyLatest && isfield(m,'absQ') && isfield(m.absQ,'TE_ms') ...
            && abs(m.absQ.TE_ms - 0.1) > 1e-9
        ok = false; return
    end
    x   = r.xFit(:);
    sel = x >= min(zoom31P) & x <= max(zoom31P);
    assert(any(sel), 'zoom window %g..%g contains no fitted points', zoom31P);
    panels(end+1) = struct('x', x(sel), 'y', real(r.yFit(sel)), ...
        'fit', real(r.fit(sel)), 'title', 'aATP / NAD / UDP');
end
end

% =====================================================================
function prep(ax, surface, gridClr, inkSec)
set(ax, 'Color', surface, 'Box','off', 'XDir','reverse', 'TickDir','out', ...
        'XColor', inkSec, 'YColor', inkSec, 'FontSize', 9, ...
        'LineWidth', 0.75, 'TickLength',[0.008 0.008]);
ax.YGrid = 'on';  ax.GridColor = gridClr;  ax.GridAlpha = 1;  ax.Layer = 'bottom';
hold(ax, 'on');
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% =====================================================================
function sessions = findSessions(dataRoot)
if isfolder(fullfile(dataRoot, 'processed')), sessions = {dataRoot}; return; end
sessions = {};
d = dir(dataRoot); d = d([d.isdir]);
d = d(~ismember({d.name}, {'.','..'}));
for k = 1:numel(d)
    p = fullfile(d(k).folder, d(k).name);
    if strcmpi(d(k).name, 'Repeatability')
        dd = dir(p); dd = dd([dd.isdir]);
        dd = dd(~ismember({dd.name}, {'.','..'}));
        for j = 1:numel(dd)
            sessions{end+1} = fullfile(dd(j).folder, dd(j).name); %#ok<AGROW>
        end
    else
        sessions{end+1} = p; %#ok<AGROW>
    end
end
end
