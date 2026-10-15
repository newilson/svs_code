function plotExtRefFits31P(dataRoot, outDir)
% PLOTEXTREFFITS31P  Overall 31P fit and the reference-peak component, all scans.
%
%   plotExtRefFits31P()
%   plotExtRefFits31P(dataRoot, outDir)
%
% One figure covering every session whose 31P fit included the RefExt singlet.
% Each scan gets a vertical pair of panels:
%
%   top     the whole fitted range -- data, total model, residual
%   bottom  zoomed on the phantom peak -- the RefExt component and, when
%           present, BroadRefBG, both drawn on the baseline so they sit where
%           they actually contribute
%
% The point of the pairing is that RefExt looks fine in isolation but is only
% interpretable against the broad envelope it is separated from: BroadRefBG
% underlies the whole +3..13 ppm region and is far larger than RefExt itself,
% so how the two split the area sets the external concentration scale.  The
% fitted RefExt area is printed on each panel because that is the number the
% external quant actually uses, and it spans ~59x across the calf scans.
%
% Spectra are normalized per panel to their own maximum, so shapes are
% comparable across scans; absolute scale lives in the annotated area.
%
% Output: <outDir>\U01_2026_extref_fits.png and .fig

if nargin < 1 || isempty(dataRoot)
    dataRoot = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026\2026';
end
if nargin < 2 || isempty(outDir)
    outDir = fileparts(dataRoot);
end

ZOOM = [2 16];        % ppm window around the +8 ppm phantom line

% ---- gather ------------------------------------------------------------
S = [];
sess = findSessions(dataRoot);
for k = 1:numel(sess)
    procDir = fullfile(sess{k}, 'processed');
    if ~isfolder(procDir), continue; end
    f = dir(fullfile(procDir, '*_31P_abs.mat'));
    for j = 1:numel(f)
        d = readOne(fullfile(f(j).folder, f(j).name));
        if isempty(d), continue; end
        if isempty(S), S = d; else, S(end+1) = d; end %#ok<AGROW>
    end
end
assert(~isempty(S), 'No external-reference fits found under %s', dataRoot);

[~, ord] = sortrows([string({S.anat})', string({S.subject})', string({S.scandate})']);
S = S(ord);
for a = {'brain','calf'}
    idx = find(strcmp({S.anat}, a{1}));
    for i = 1:numel(idx)
        S(idx(i)).id = sprintf('%s%d', upper(a{1}(1)), i);
    end
end

% ---- figure ------------------------------------------------------------
n    = numel(S);
ncol = 3;
npr  = ceil(n / ncol);                 % pair-rows
fig  = figure('Name','U01_2026 extref fits', 'Color','w', ...
              'Position', [40 40 1500 320*npr]);
t = tiledlayout(fig, 2*npr, ncol, 'TileSpacing','compact', 'Padding','compact');

for i = 1:n
    s   = S(i);
    col = mod(i-1, ncol) + 1;
    pr  = floor((i-1) / ncol) + 1;

    % ---------- top: whole fitted range ----------
    nexttile(t, (2*pr-2)*ncol + col);
    sc = max(abs(s.y));
    plot(s.x, s.y/sc, 'Color',[.65 .65 .65], 'LineWidth',1.6); hold on
    plot(s.x, s.fit/sc, 'k', 'LineWidth',1.0);
    plot(s.x, (s.y - s.fit)/sc - 0.35, 'Color',[.85 .33 .10], 'LineWidth',0.7);
    set(gca,'XDir','reverse','YTick',[],'Box','off','FontSize',8);
    xlim([min(s.x) max(s.x)]);
    yline(-0.35, ':', 'Color',[.7 .7 .7]);
    title(sprintf('%s  %s %s  (%s)', s.id, s.subject, s.scandate, s.anat), ...
          'FontSize',9, 'FontWeight','bold');
    if i == 1, legend({'data','fit','resid'},'FontSize',7,'Box','off', ...
                      'Location','northwest'); end

    % ---------- bottom: the phantom peak ----------
    nexttile(t, (2*pr-1)*ncol + col);
    m  = s.x >= ZOOM(1) & s.x <= ZOOM(2);
    sc = max(abs(s.y(m)));
    if ~(sc > 0), sc = 1; end
    plot(s.x(m), s.y(m)/sc,   'Color',[.65 .65 .65], 'LineWidth',1.6); hold on
    plot(s.x(m), s.fit(m)/sc, 'k', 'LineWidth',1.0);
    plot(s.x(m), s.base(m)/sc, ':', 'Color',[.4 .4 .4], 'LineWidth',0.9);
    if ~isempty(s.bg)
        plot(s.x(m), (s.bg(m)+s.base(m))/sc, 'Color',[.93 .69 .13], 'LineWidth',1.1);
    end
    if ~isempty(s.ref)
        area(s.x(m), (s.ref(m)+s.base(m))/sc, 'FaceColor',[.85 .10 .10], ...
             'FaceAlpha',0.30, 'EdgeColor',[.75 0 0], 'LineWidth',1.3, ...
             'BaseValue', min(s.base(m)/sc));
    end
    set(gca,'XDir','reverse','YTick',[],'Box','off','FontSize',8);
    xlim(ZOOM);
    xlabel('ppm (PCr = 0)','FontSize',8);
    txt = sprintf('RefExt area %.3g', s.areaRef);
    if ~s.ok, txt = [txt '  [not quantified]']; end %#ok<AGROW>
    text(0.03, 0.92, txt, 'Units','normalized', 'FontSize',7.5, ...
         'Color',[.6 0 0], 'FontWeight','bold');
    if i == 1
        lg = {'data','fit','baseline'};
        if ~isempty(s.bg), lg{end+1} = 'BroadRefBG'; end
        if ~isempty(s.ref), lg{end+1} = 'RefExt'; end
        legend(lg, 'FontSize',7, 'Box','off', 'Location','northeast');
    end
end

title(t, 'U01\_2026 external-reference scans: overall 31P fit (top) and the RefExt phantom peak (bottom)', ...
      'FontWeight','bold');

png = fullfile(outDir, 'U01_2026_extref_fits.png');
fgf = fullfile(outDir, 'U01_2026_extref_fits.fig');
exportgraphics(fig, png, 'Resolution', 150);
savefig(fig, fgf);
fprintf('Wrote %s\n       %s\n', png, fgf);
end

% =====================================================================
function d = readOne(matFile)
d = [];
try
    m = load(matFile);
catch ME
    warning('plotExtRefFits31P: load failed %s (%s)', matFile, ME.message); return
end
if ~isfield(m,'absQ') || ~isfield(m.absQ,'external'), return; end
if ~isfield(m,'output') || ~isfield(m.output,'basisVarproResults'), return; end
b = m.output.basisVarproResults;

[~, base] = fileparts(matFile);
tok = regexp(base, '^(.*)_(brain|calf)_31P_abs$', 'tokens', 'once');
if isempty(tok), return; end

d.id = ''; d.session = tok{1}; d.anat = tok{2};
[d.subject, d.scandate] = parseScanID(d.session);
d.x    = b.xFit(:);
d.y    = b.yFit(:);
d.fit  = b.fit(:);
d.base = b.baseline(:);
d.ref  = pickComp(b, 'RefExt');
d.bg   = pickComp(b, 'BroadRefBG');
e = m.absQ.external;
d.ok      = isfield(e,'ok') && logical(e.ok);
d.areaRef = NaN;
if isfield(e,'areaRef'), d.areaRef = e.areaRef; end
end

% =====================================================================
function v = pickComp(b, name)
v = [];
if ~isfield(b,'comp') || ~isfield(b,'names'), return; end
i = find(strcmp(b.names(:), name), 1);
if isempty(i) || i > size(b.comp,2), return; end
v = b.comp(:, i);
end

% =====================================================================
function sessions = findSessions(dataRoot)
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

% =====================================================================
function [subj, dstr] = parseScanID(scanID)
subj = ''; dstr = '';
tok = regexp(scanID, '^(S\d+)_(\d{2})(\d{2})(\d{4})', 'tokens', 'once');
if ~isempty(tok)
    subj = tok{1};
    dstr = sprintf('%s-%s-%s', tok{4}, tok{2}, tok{3});
else
    p = strsplit(scanID, '_'); subj = p{1};
end
end
