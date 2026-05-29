function fig = plot_varpro_1h(x, y, out, names, titleStr)
% PLOT_VARPRO_1H  Diagnostic plot for curvefitAuto_varproBaseline output (1H).
%
%   fig = plot_varpro_1h(x, y, out, names)
%   fig = plot_varpro_1h(x, y, out, names, titleStr)
%
%   Three stacked panels (data + total fit + spline baseline; per-component
%   contributions overlaid on the data; residual), all in the fit window.
%   x-axes are linked and the per-component legend floats in the right-margin
%   gap so it doesn't shrink the middle panel.
%
%   Inputs:
%       x        - frequency axis (ppm) over the fit window [np x 1]
%       y        - data (real) over the fit window         [np x 1]
%       out      - struct returned by curvefitAuto_varproBaseline; uses
%                  .fit, .baseline, .comp, .areas
%       names    - cell array of peak names {n x 1}
%       titleStr - optional super-title string

if nargin < 5 || isempty(titleStr), titleStr = ''; end

x        = x(:);
y        = real(y(:));
yfit     = out.fit(:);
baseline = out.baseline(:);
comp     = out.comp;
ampl     = out.areas(:);
resid    = y - yfit;

if nargin < 4 || isempty(names)
    names = arrayfun(@(k) sprintf('p%d',k), 1:numel(ampl), 'UniformOutput', false);
end

fig = figure('Position',[100 100 1100 800], 'Name','1H varpro fit');

% ---- Top: data + fit + baseline
ax1 = subplot(3,1,1);
plot(x, y, 'k', 'LineWidth', 1.0); hold on;
plot(x, real(yfit), 'r', 'LineWidth', 1.2);
if ~isempty(baseline) && any(baseline ~= 0)
    plot(x, real(baseline), 'b--', 'LineWidth', 1.0);
    legend('data','fit','baseline','Location','best');
else
    legend('data','fit','Location','best');
end
set(ax1,'XDir','reverse');
xlabel('ppm');
ylabel('a.u.');
if ~isempty(titleStr)
    title(ax1, titleStr, 'Interpreter','none');
else
    title(ax1, 'varpro fit');
end
grid on;
yl = ylim;

% ---- Middle: per-component contributions (comp is already amp-scaled)
ax2 = subplot(3,1,2);
hold on;
plot(x, y, 'Color',[0.6 0.6 0.6]);
nMet = size(comp,2);
cmap = lines(nMet);
labels = cell(nMet,1);
for k = 1:nMet
    plot(x, real(comp(:,k)), 'Color', cmap(k,:), 'LineWidth', 1.0);
    labels{k} = sprintf('%s (ampl=%.3g)', names{k}, ampl(k));
end
set(ax2,'XDir','reverse');
xlabel('ppm');
ylabel('a.u.');
title('per-component fits (ampl in legend; raw ampl in out.areas)');
grid on;
ylim(yl);

% ---- Bottom: residual
ax3 = subplot(3,1,3);
plot(x, real(resid), 'Color', [0.2 0.2 0.2]);
set(ax3,'XDir','reverse');
xlabel('ppm');
ylabel('a.u.');
title('residual (data - fit)');
grid on;
ylim(yl);

linkaxes([ax1 ax2 ax3], 'x');

% ---- Float the per-component legend in a reserved right-margin gap so
%      the middle axes width matches the top/bottom axes (vertical alignment).
pos1 = get(ax1,'Position');
pos2 = get(ax2,'Position');
pos3 = get(ax3,'Position');

legGap   = 0.18;
axRight  = 1 - legGap - 0.02;
axWidth  = axRight - pos1(1);
set(ax1,'Position',[pos1(1) pos1(2) axWidth pos1(4)]);
set(ax2,'Position',[pos2(1) pos2(2) axWidth pos2(4)]);
set(ax3,'Position',[pos3(1) pos3(2) axWidth pos3(4)]);

hLeg = legend(ax2, [{'data'}, labels(:).'], 'FontSize', 8, 'Interpreter','none');
set(hLeg, 'Location', 'none');
legX = pos1(1) + axWidth + 0.01;
legY = pos2(2);
legW = legGap;
legH = pos2(4);
set(hLeg, 'Position', [legX legY legW legH]);

end
