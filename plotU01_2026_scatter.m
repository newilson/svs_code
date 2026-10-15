function outPng = plotU01_2026_scatter(anat, nucleus, inDir, outPng)
% PLOTU01_2026_SCATTER  Per-scan concentration scatter, latest protocol.
%
%   plotU01_2026_scatter('brain','1H')
%   plotU01_2026_scatter('calf','31P')
%   plotU01_2026_scatter(anat, nucleus, inDir, outPng)
%
% One x position per scan, ordered earliest to most recent; one series per
% reported metabolite.  Subject identity is deliberately not encoded -- every
% scan is just the next point along.
%
% Series per panel:
%   brain 1H   NAD+ H2, H6, Trp        [uM]
%   calf  1H   NAD+ H2, H6             [uM]
%   * 31P      NAD+, NADH              [mM]
%
% Each legend entry carries that metabolite's mean +- SD over the scans shown.
%
% 31P uses INTERNAL aATP for both anatomies -- 2.7 mM in brain, 8.2 mM in
% calf.  No external referencing.
%
% Same protocol filter as statsU01_2026_latest: 31P TE 0.1 ms (flip 70, rect
% NOE) and 1H svssel_20260101.

if nargin < 1 || isempty(anat),    anat    = 'brain'; end
if nargin < 2 || isempty(nucleus), nucleus = '1H';    end
if nargin < 3 || isempty(inDir)
    inDir = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026';
end
anat = lower(char(anat));  nucleus = upper(char(nucleus));
if nargin < 4 || isempty(outPng)
    outPng = fullfile(inDir, sprintf('U01_2026_%s_%s_scatter.png', anat, nucleus));
end

% ---- palette ----------------------------------------------------------
% Slots 1-3 (blue/orange/aqua) validate all-pairs, which is what a scatter
% needs.  If a 4th series is ever added back, use VIOLET #4a3aa7, not the
% documented slot 4 (yellow): yellow-vs-orange fails the normal-vision floor
% all-pairs (dE 13.7 < 15), while blue/orange/violet/aqua passes (worst CVD
% dE 9.2, worst normal-vision dE 16.3, light mode).
C = @(h) reshape(sscanf(h(2:end),'%2x'),1,3)/255;
surface = C('#fcfcfb');
gridClr = C('#e1e0d9');
muted   = C('#898781');
inkSec  = C('#52514e');
inkPri  = C('#0b0b0b');

% ---- what to plot -----------------------------------------------------
% Shared with plotU01_2026_repeat so the two figure families cannot drift.
% Both anatomies reference internal aATP; collectU01_2026_TSV writes every
% reference set for every scan, so the aATP set comes straight out of the TSV
% -- no reprocessing needed to change reference.
[cols, labels, series, shapes, unit] = u01SeriesSpec(anat, nucleus);

% ---- data -------------------------------------------------------------
T = readtable(fullfile(inDir, sprintf('U01_2026_%s.tsv', anat)), 'FileType','text', ...
              'Delimiter','\t', 'VariableNamingRule','preserve');
T = T(abs(double(T.TE_31P_ms)-0.1) < 1e-9 & contains(string(T.seq_1H),'svssel_20260101'), :);

% Excluded by inspection, not by any rule the filter above can see: known
% protocol problems on this session (NW, 2026-08-07).  Applied to BOTH nuclei
% of that visit -- drop the entry if the 31P of that session is in fact sound.
EXCLUDE = "S007_07222026_31P-1H-HeadCoil_external_ref";
T = T(~ismember(string(T.session), EXCLUDE), :);

T = sortrows(T, {'scandate','session'});          % earliest -> most recent
n = height(T);
assert(n > 0, 'no %s %s scans survive the protocol filter', anat, nucleus);

V = nan(n, numel(cols));
for q = 1:numel(cols)
    if ismember(cols{q}, T.Properties.VariableNames), V(:,q) = double(T.(cols{q})); end
end
assert(any(isfinite(V(:))), 'no finite values for %s %s -- check column names', anat, nucleus);

% ---- axis range: nice round ticks with headroom -----------------------
[yMax, yStep] = u01NiceAxis(max(V(:), [], 'omitnan'));

% ---- figure -----------------------------------------------------------
% Legend entries carry mean +- SD, so they are long; give the legend its own
% full-width row under the title rather than squeezing it into the corner.
statLabels = labels;
for q = 1:numel(cols)
    v = V(isfinite(V(:,q)), q);
    if isempty(v), continue; end
    statLabels{q} = sprintf('%s  %s', labels{q}, fmtPair(mean(v), std(v)));
end

fig = figure('Color', surface, 'Units','pixels', 'Position',[100 100 780 520], ...
             'InvertHardcopy','off');
ax = axes(fig, 'Color', surface, 'Units','normalized', 'Position',[0.10 0.16 0.87 0.67]);
hold(ax,'on');

xr = [0.55 n+0.45];
yt = 0:yStep:yMax;
for y = yt
    plot(ax, xr, [y y], '-', 'Color', gridClr, 'LineWidth', 0.75);
end

h = gobjects(numel(cols),1);
for q = 1:numel(cols)
    h(q) = plot(ax, 1:n, V(:,q), shapes{q}, 'MarkerSize', 9, ...
                'MarkerFaceColor', series(q,:), 'MarkerEdgeColor', surface, ...
                'LineWidth', 2, 'LineStyle','none');
end
% No right-edge direct labels: they land on top of the last markers.  The
% legend carries identity; aqua's sub-3:1 contrast relief is the table view.

% ---- axes chrome ------------------------------------------------------
% XColor 'none' rather than a hidden axle: hiding the axle still leaves a stub
% painted past the last tick.  Scan numbers are drawn as text below instead.
set(ax, 'XLim',xr, 'YLim',[0 yMax], 'XTick',[], 'XColor','none', ...
        'YTick',yt, 'Box','off', 'YColor',inkSec, ...
        'FontSize',11, 'LineWidth',0.75, 'TickLength',[0 0]);
for k = 1:n
    text(ax, k, -0.033*yMax, sprintf('%d',k), 'Color', inkSec, 'FontSize', 11, ...
         'HorizontalAlignment','center', 'VerticalAlignment','top');
end
text(ax, mean(xr), -0.139*yMax, 'Scan No.', 'Color', inkSec, 'FontSize', 11, ...
     'HorizontalAlignment','center');
ylabel(ax, sprintf('Concentration (%s)', unit), 'Color', inkSec, 'FontSize', 11);

title(ax, sprintf('%s %s NAD metabolites', upperFirst(anat), nucleus), ...
      'Color', inkPri, 'FontSize', 13, 'FontWeight','bold', ...
      'HorizontalAlignment','left', 'Units','normalized', 'Position',[0 1.16 0]);

lg = legend(ax, h, statLabels, 'Orientation','horizontal', 'Box','off', ...
            'TextColor', inkSec, 'FontSize', 11);
lg.ItemTokenSize = [14 8];
lg.Units = 'normalized';
lg.Position = [0.10 0.855 0.87 0.05];

ax.YAxis.Axle.Visible = 'off';
ax.Clipping = 'off';

exportgraphics(fig, outPng, 'Resolution', 200, 'BackgroundColor', surface);
fprintf('Wrote %s\n', outPng);

for k = 1:n
    fprintf('  %d  %s  %-48s', k, datestr(T.scandate(k),'yyyy-mm-dd'), T.session{k}); %#ok<DATST>
    fprintf(' %9.4g', V(k,:));  fprintf('\n');
end
end

% =====================================================================
function s = upperFirst(s)
s(1) = upper(s(1));
end

function s = fmtPair(m, sd)
% Both numbers share the decimal count, chosen from the mean -- formatting
% them independently gives things like "1.00\pm0.094".  uM values run in the
% hundreds, mM values down to ~0.03, hence the spread of cases.
if ~isfinite(m), s = '-'; return; end
if     abs(m) >= 100, d = 0;
elseif abs(m) >= 10,  d = 1;
elseif abs(m) >= 1,   d = 2;
else,                 d = 3;
end
s = sprintf('%.*f\\pm%.*f', d, m, d, sd);
end
