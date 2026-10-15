function outPng = plotU01_2026_repeat(anat, nucleus, inDir, outPng)
% PLOTU01_2026_REPEAT  Test-retest scatter for the Repeatability\ cohort.
%
%   plotU01_2026_repeat('calf','1H')
%   plotU01_2026_repeat(anat, nucleus, inDir, outPng)
%
% One x position per SUBJECT.  Both visits are drawn at the same x in the same
% color, so the vertical gap between a pair IS that subject's test-retest
% spread -- the visits are deliberately not distinguished from each other.
%
% Legend reports, per metabolite, the MEAN OVER SUBJECTS of the per-subject
% CV between the two visits:
%       CV_i    = 100 * std([v1_i v2_i]) / mean([v1_i v2_i])
%       CV_mean = mean_i(CV_i)
% Note this is NOT the Bland-Altman within-subject CV that
% statsU01_2026_latest reports (wsSD = sqrt(sum(d^2)/2n) against the grand
% mean).  For a pair, std([v1 v2]) = |d|/sqrt(2), so a per-subject CV runs
% about 1/sqrt(2) of the equivalent Bland-Altman figure before the averaging
% is taken into account.  Quote whichever you name explicitly.
%
% NO PROTOCOL FILTER is applied: the whole Repeatability cohort was acquired
% on the PREVIOUS 31P protocol (TE 0.3 ms, flip 90, hann NOE), so filtering to
% the latest protocol would leave nothing.  These figures are therefore NOT
% comparable in absolute terms to the plotU01_2026_scatter ones.

if nargin < 1 || isempty(anat),    anat    = 'calf'; end
if nargin < 2 || isempty(nucleus), nucleus = '1H';   end
if nargin < 3 || isempty(inDir)
    inDir = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026';
end
anat = lower(char(anat));  nucleus = upper(char(nucleus));
if nargin < 4 || isempty(outPng)
    outPng = fullfile(inDir, sprintf('U01_2026_%s_%s_repeat.png', anat, nucleus));
end

[cols, labels, series, shapes, unit] = u01SeriesSpec(anat, nucleus);

C = @(h) reshape(sscanf(h(2:end),'%2x'),1,3)/255;
surface = C('#fcfcfb');
gridClr = C('#e1e0d9');
muted   = C('#898781');
inkSec  = C('#52514e');
inkPri  = C('#0b0b0b');

% ---- data -------------------------------------------------------------
T = readtable(fullfile(inDir, sprintf('U01_2026_%s.tsv', anat)), 'FileType','text', ...
              'Delimiter','\t', 'VariableNamingRule','preserve');
assert(ismember('cohort', T.Properties.VariableNames), ...
       'no cohort column -- re-run collectU01_2026_TSV');
T = T(strcmp(string(T.cohort),'repeat'), :);
assert(~isempty(T), ['no %s scans in the Repeatability cohort. The four brain ' ...
       'repeatability sessions have never been processed.'], anat);
T = sortrows(T, {'subject','scandate'});

subs = unique(string(T.subject), 'stable');
nS   = numel(subs);
nQ   = numel(cols);

% V(subject, visit, metabolite); only the first two visits are used
V = nan(nS, 2, nQ);
for s = 1:nS
    idx = find(string(T.subject) == subs(s));
    idx = idx(1:min(2,end));
    for q = 1:nQ
        if ~ismember(cols{q}, T.Properties.VariableNames), continue; end
        v = double(T.(cols{q})(idx));
        V(s, 1:numel(v), q) = v(:).';
    end
end
assert(any(isfinite(V(:))), 'no finite values for %s %s', anat, nucleus);

% ---- per-subject CV, then mean over subjects ---------------------------
cvMean = nan(1,nQ);  nPair = zeros(1,nQ);
for q = 1:nQ
    cv = nan(nS,1);
    for s = 1:nS
        pair = squeeze(V(s,:,q));
        pair = pair(isfinite(pair));
        if numel(pair) < 2, continue; end
        cv(s) = 100 * std(pair) / abs(mean(pair));
    end
    nPair(q)  = sum(isfinite(cv));
    cvMean(q) = mean(cv, 'omitnan');
end

statLabels = labels;
for q = 1:nQ
    if isfinite(cvMean(q))
        statLabels{q} = sprintf('%s  CV %.1f%%', labels{q}, cvMean(q));
    end
end

[yMax, yStep] = u01NiceAxis(max(V(:), [], 'omitnan'));

% ---- figure -----------------------------------------------------------
fig = figure('Color', surface, 'Units','pixels', 'Position',[100 100 780 520], ...
             'InvertHardcopy','off');
ax = axes(fig, 'Color', surface, 'Units','normalized', 'Position',[0.10 0.16 0.87 0.67]);
hold(ax,'on');

xr = [0.55 nS+0.45];
yt = 0:yStep:yMax;
for y = yt
    plot(ax, xr, [y y], '-', 'Color', gridClr, 'LineWidth', 0.75);
end

% A hairline joining each pair makes it obvious the two points belong
% together where they sit far apart; it carries no extra information.
h = gobjects(nQ,1);
for q = 1:nQ
    for s = 1:nS
        pair = squeeze(V(s,:,q));
        if all(isfinite(pair))
            plot(ax, [s s], pair, '-', 'Color', [series(q,:) 0.35], 'LineWidth', 1.25);
        end
    end
    p = plot(ax, repmat((1:nS).',1,2), squeeze(V(:,:,q)), shapes{q}, 'MarkerSize', 9, ...
             'MarkerFaceColor', series(q,:), 'MarkerEdgeColor', surface, ...
             'LineWidth', 2, 'LineStyle','none');
    h(q) = p(1);
end

% ---- axes chrome ------------------------------------------------------
set(ax, 'XLim',xr, 'YLim',[0 yMax], 'XTick',[], 'XColor','none', ...
        'YTick',yt, 'Box','off', 'YColor',inkSec, ...
        'FontSize',11, 'LineWidth',0.75, 'TickLength',[0 0]);
for s = 1:nS
    text(ax, s, -0.033*yMax, char(subs(s)), 'Color', inkSec, 'FontSize', 11, ...
         'HorizontalAlignment','center', 'VerticalAlignment','top');
end
text(ax, mean(xr), -0.139*yMax, 'Subject', 'Color', inkSec, 'FontSize', 11, ...
     'HorizontalAlignment','center');
ylabel(ax, sprintf('Concentration (%s)', unit), 'Color', inkSec, 'FontSize', 11);

title(ax, sprintf('%s %s repeatability', upperFirst(anat), nucleus), ...
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

fprintf('  CV is mean over subjects of per-subject CV between visits\n');
for q = 1:nQ
    fprintf('  %-10s CV %5.1f%%  (%d subject pairs)\n', labels{q}, cvMean(q), nPair(q));
end
for s = 1:nS
    fprintf('  %-6s', subs(s));
    for q = 1:nQ, fprintf('  %9.4g %9.4g |', V(s,1,q), V(s,2,q)); end
    fprintf('\n');
end
end

% =====================================================================
function s = upperFirst(s)
s(1) = upper(s(1));
end
