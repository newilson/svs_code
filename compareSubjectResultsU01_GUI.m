function f = compareSubjectResultsU01_GUI(S, labels, ttl, figName)
% COMPARESUBJECTRESULTSU01_GUI  Adjust and save the two-visit spectra overlay.
%
%   f = compareSubjectResultsU01_GUI(S, labels, ttl, figName)
%
% S(k).df, S(k).wat and S(k).p31 hold visit k's complex 1H downfield, 1H
% water and 31P spectra (.ppm, .spec, .pivot = default PC1 pivot [ppm];
% .df also .band, the ppm range its y-axis
% is scaled to, and .range, its x-axis; .wat and .p31 also .area, the water
% fit area and the external phantom area), as prepared by
% compareSubjectResultsU01.  An empty field is a missing output.  labels{k}
% names visit k, ttl titles the figure.  The water and full 31P legends show
% .area times the current scale factor.
%
% The overlay is drawn in the middle; visit 1 controls are on the left and
% visit 2 controls on the right.  Each spectrum gets a zero- and first-order
% phase and a scale factor, applied on top of the phase and scaling it
% arrives with (all zero and scale 1 leave it unchanged).  The phase is
% PC0 + PC1*(ppm - pivot) [deg], with PC1 in deg/ppm.  The two 31P panels
% share one setting.  Save writes <figName>.fig and
% <figName>.png (300 dpi) with the current settings, stored in the saved
% figure's UserData.

col   = {[0 0.35 0.85], [0.85 0.1 0.1]};
keys  = {'wat', 'df', 'p31'};
names = {'1H water', '1H downfield', '31P'};
R = struct('name', {'PC0', 'PC1', 'Pivot', 'Scale'}, 'field', {'pc0', 'pc1', 'piv', 'scl'}, ...
    'slider', {[-180 180 1/360 10/360], [-20 20 0.25/40 2.5/40], [], [0.1 3 0.01/2.9 0.1/2.9]}, ...
    'tip', {'zero-order phase [deg]', 'first-order phase [deg/ppm]', 'first-order pivot [ppm]', 'scale factor'}, ...
    'fmt', {'%.1f', '%.2f', '%.2f', '%.3g'});
nV = numel(S);
piv0 = NaN(nV, 3);
for k = 1:nV
    for g = 1:3
        if ~isempty(S(k).(keys{g})), piv0(k,g) = S(k).(keys{g}).pivot; end
    end
end
A = struct('pc0', zeros(nV, 3), 'pc1', zeros(nV, 3), 'piv', piv0, 'scl', ones(nV, 3));

f = figure('Name', [ttl ': adjust, Save, then close'], 'NumberTitle', 'off', 'Visible', 'off', ...
    'Position', [60 60 1600 900], 'Tag', 'compareSubjectResultsU01_GUI');
pnl = uipanel('Parent', f, 'Position', [.17 .075 .66 .92], 'BackgroundColor', 'w', 'BorderType', 'none');
G = drawPanels(pnl, S, A, labels, ttl);

ctl  = gobjects(nV, 3, numel(R), 2);
xCol = [.005 .835];
for k = 1:nV
    cp = uipanel('Parent', f, 'Title', labels{k}, 'ForegroundColor', col{k}, 'FontWeight', 'bold', ...
        'Position', [xCol(k) .075 .16 .92]);
    for g = 1:3
        y0 = 0.98 - (g-1)*0.325;
        on = 'on';
        if isempty(S(k).(keys{g})), on = 'off'; end
        uicontrol('Parent', cp, 'style', 'text', 'units', 'normalized', 'position', [.03 y0-.045 .94 .035], ...
            'string', names{g}, 'HorizontalAlignment', 'left', 'fontweight', 'bold');
        for r = 1:numel(R)
            yr = y0 - .045 - r*.05;
            uicontrol('Parent', cp, 'style', 'text', 'units', 'normalized', 'position', [.03 yr .19 .035], ...
                'string', R(r).name, 'HorizontalAlignment', 'left', 'Tooltip', R(r).tip);
            if ~isempty(R(r).slider)
                sl = R(r).slider;
                ctl(k,g,r,1) = uicontrol('Parent', cp, 'Style', 'slider', 'units', 'normalized', ...
                    'Position', [.22 yr .53 .035], 'min', sl(1), 'max', sl(2), 'sliderstep', sl(3:4), ...
                    'Enable', on, 'Tooltip', R(r).tip, 'Tag', sprintf('%sslider_%d_%s', lower(R(r).name), k, keys{g}), ...
                    'callback', @(src, ~) setVal(k, g, r, src.Value));
                addlistener(ctl(k,g,r,1), 'ContinuousValueChange', @(src, ~) setVal(k, g, r, src.Value));
            end
            ctl(k,g,r,2) = uicontrol('Parent', cp, 'style', 'edit', 'units', 'normalized', ...
                'position', [.77 yr .2 .035], 'Enable', on, 'Tooltip', R(r).tip, ...
                'Tag', sprintf('%sedit_%d_%s', lower(R(r).name), k, keys{g}), ...
                'callback', @(src, ~) setVal(k, g, r, str2double(src.String)));
        end
        syncControls(k, g);
    end
end

uicontrol('Parent', f, 'style', 'pushbutton', 'units', 'normalized', 'position', [.43 .015 .065 .045], ...
    'callback', @saveFigure, 'string', 'Save', 'fontweight', 'bold', 'Tag', 'save');
uicontrol('Parent', f, 'style', 'pushbutton', 'units', 'normalized', 'position', [.505 .015 .065 .045], ...
    'callback', @resetAll, 'string', 'Reset', 'Tag', 'reset');

set(f, 'visible', 'on', 'toolbar', 'figure')

    function setVal(k, g, r, v)
        if isfinite(v)
            switch R(r).field
                case 'pc0', A.pc0(k,g) = mod(v + 180, 360) - 180;
                case 'scl', if v > 0, A.scl(k,g) = v; end
                otherwise,  A.(R(r).field)(k,g) = v;
            end
        end
        syncControls(k, g);
        updatePanels(G, S, A, g);
    end

    function syncControls(k, g)
        for rr = 1:numel(R)
            v = A.(R(rr).field)(k,g);
            if isgraphics(ctl(k,g,rr,1))
                set(ctl(k,g,rr,1), 'value', min(max(v, R(rr).slider(1)), R(rr).slider(2)));
            end
            set(ctl(k,g,rr,2), 'string', sprintf(R(rr).fmt, v));
        end
    end

    function resetAll(~, ~)
        A = struct('pc0', zeros(nV, 3), 'pc1', zeros(nV, 3), 'piv', piv0, 'scl', ones(nV, 3));
        for kk = 1:nV
            for gg = 1:3, syncControls(kk, gg); end
        end
        for gg = 1:3, updatePanels(G, S, A, gg); end
    end

    function saveFigure(~, ~)
        fs = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1400 900], 'Name', ttl);
        drawPanels(fs, S, A, labels, ttl);
        fs.UserData = struct('spectra', {names}, 'visits', {labels(:)'}, 'pc0_deg', A.pc0, ...
            'pc1_degPerPpm', A.pc1, 'pivot_ppm', A.piv, 'scale', A.scl);
        set(fs, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')');
        savefig(fs, [figName '.fig']);
        exportgraphics(fs, [figName '.png'], 'Resolution', 300);
        close(fs);
        fprintf('Saved %s.fig/.png\n', figName);
        for kk = 1:nV
            for gg = 1:3
                if A.pc0(kk,gg) ~= 0 || A.pc1(kk,gg) ~= 0 || A.scl(kk,gg) ~= 1
                    fprintf('  %s %s: PC0 %+.1f deg, PC1 %+.2f deg/ppm (pivot %.2f ppm), scale %.3g\n', labels{kk}, ...
                        names{gg}, A.pc0(kk,gg), A.pc1(kk,gg), A.piv(kk,gg), A.scl(kk,gg));
                end
            end
        end
    end
end

% =====================================================================
function G = drawPanels(parent, S, A, labels, ttl)
col  = {[0 0.35 0.85], [0.85 0.1 0.1]};
keys = {'wat', 'df', 'p31'};
P = struct('key', {1, 2, 3, 3}, ...
    'title', {'1H water', '1H downfield', '31P', '31P: \alpha-ATP / NAD / UDP(G)'}, ...
    'ylabel', {'signal [a.u.]', 'signal / water integral', 'signal [a.u.]', 'signal [a.u.]'}, ...
    'xr', {[-Inf Inf], [-Inf Inf], [-20 15], [-11.5 -6.5]}, ...
    'top', {0.08, 0.12, 0.08, 0.08}, ...
    'area', {'area', '', 'phantom area', ''});
nV = numel(S);
G = struct('P', P, 'labels', {labels}, 'ax', gobjects(1, 4), 'line', gobjects(4, nV), ...
    'idx', {cell(4, nV)}, 'yidx', {cell(4, nV)});

tl = tiledlayout(parent, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, ttl, 'FontWeight', 'bold', 'Interpreter', 'none');
for p = 1:4
    ax = nexttile(tl);
    hold(ax, 'on');
    G.ax(p) = ax;
    set(ax, 'XDir', 'reverse', 'Box', 'on', 'FontSize', 10);
    xlabel(ax, 'ppm'); ylabel(ax, P(p).ylabel); title(ax, P(p).title);
    g = P(p).key;
    h = gobjects(0); xl = [Inf -Inf];
    for k = 1:nV
        X = S(k).(keys{g});
        if isempty(X), continue; end
        i = X.ppm >= P(p).xr(1) & X.ppm <= P(p).xr(2);
        G.idx{p,k} = i;
        if isfield(X, 'band')
            G.yidx{p,k} = X.ppm(i) >= X.band(1) & X.ppm(i) <= X.band(2);
        else
            G.yidx{p,k} = true(nnz(i), 1);
        end
        G.line(p,k) = plot(ax, X.ppm(i), adjusted(X, i, A, k, g), 'Color', col{k}, 'LineWidth', 1, ...
            'DisplayName', lineName(P(p), labels{k}, X, A.scl(k,g)));
        h(end+1) = G.line(p,k); %#ok<AGROW>
        if isfield(X, 'range'), xr = X.range; else, xr = [min(X.ppm(i)) max(X.ppm(i))]; end
        xl = [min(xl(1), xr(1)) max(xl(2), xr(2))];
    end
    if isempty(h)
        text(ax, 0.5, 0.5, 'no processed output', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5], 'FontSize', 12);
        set(ax, 'XTick', [], 'YTick', []);
        continue
    end
    xlim(ax, xl);
    legend(ax, h, 'Location', 'northeast', 'Box', 'off');
    if p == 4
        xline(ax, [-7.56 -8.30 -9.80], ':', {'\alpha-ATP', 'NAD', 'UDP(G)'}, 'Color', [0.5 0.5 0.5], ...
            'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', ...
            'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
    end
end
setYLim(G, 1:4);
end

function updatePanels(G, S, A, g)
keys = {'wat', 'df', 'p31'};
ps = find([G.P.key] == g);
for p = ps
    for k = 1:numel(S)
        if ~isgraphics(G.line(p,k)), continue; end
        X = S(k).(keys{g});
        set(G.line(p,k), 'YData', adjusted(X, G.idx{p,k}, A, k, g), ...
            'DisplayName', lineName(G.P(p), G.labels{k}, X, A.scl(k,g)));
    end
end
setYLim(G, ps);
end

function s = lineName(Pp, label, X, scl)
s = label;
if ~isempty(Pp.area) && isfield(X, 'area') && isfinite(X.area)
    s = sprintf('%s, %s %.4g', label, Pp.area, X.area * scl);
end
end

function setYLim(G, ps)
for p = ps
    ys = [];
    for k = 1:size(G.line, 2)
        if ~isgraphics(G.line(p,k)), continue; end
        y = G.line(p,k).YData(:);
        ys = [ys; y(G.yidx{p,k})]; %#ok<AGROW>
    end
    if isempty(ys) || max(ys) == min(ys), continue; end
    ylim(G.ax(p), [min(ys) max(ys)] + [-0.05 G.P(p).top]*(max(ys) - min(ys)));
end
end

function y = adjusted(X, idx, A, k, g)
ppm = X.ppm(idx);
ph  = A.pc0(k,g) + A.pc1(k,g) * (ppm(:) - A.piv(k,g));
y   = real(X.spec(idx) .* exp(-1i*pi/180*ph)) * A.scl(k,g);
end
