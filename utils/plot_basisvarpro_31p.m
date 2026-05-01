function fig = plot_basisvarpro_31p(out, titleStr)
% PLOT_BASISVARPRO_31P  Diagnostic plot for curvefitAuto_basisVarpro output.
%
%   fig = plot_basisvarpro_31p(out)
%   fig = plot_basisvarpro_31p(out, titleStr)
%
%   Plots data + total fit + spline baseline (top), per-component
%   contributions (middle), and residual (bottom). All in the fit window.

if nargin < 2 || isempty(titleStr), titleStr = ''; end

x        = out.xFit;
y        = out.yFit;
yfit     = out.fit;
baseline = out.baseline;
comp     = out.comp;
ampl     = out.ampl;
names    = out.names;
resid    = out.residualFit;

% Per-component nP for legend display (ampl/n is concentration-proportional).
nVec = ones(numel(names),1);
if isfield(out,'basisInfo') && ~isempty(out.basisInfo)
    for kk = 1:min(numel(out.basisInfo), numel(names))
        if isfield(out.basisInfo(kk),'n') && ~isempty(out.basisInfo(kk).n) ...
                && ~isnan(out.basisInfo(kk).n) && out.basisInfo(kk).n > 0
            nVec(kk) = out.basisInfo(kk).n;
        end
    end
end
amplDisp = ampl ./ nVec;

fig = figure('Position',[100 100 1100 800]);

% ---- Top: data + fit + baseline
ax1 = subplot(3,1,1);
plot(x, real(y), 'k', 'LineWidth', 1.0); hold on;
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
    title(ax1, 'basisVarpro fit');
end
grid on;
yl = ylim;

% ---- Middle: per-component contributions (comp is already amp-scaled)
ax2 = subplot(3,1,2);
hold on;
plot(x, real(y), 'Color',[0.6 0.6 0.6]);
nMet = size(comp,2);
cmap = lines(nMet);
labels = cell(nMet,1);
for k = 1:nMet
    plot(x, real(comp(:,k)), 'Color', cmap(k,:), 'LineWidth', 1.0);
    labels{k} = sprintf('%s (ampl/n=%.3g)', names{k}, amplDisp(k));
end
set(ax2,'XDir','reverse');
xlabel('ppm');
ylabel('a.u.');
title('per-component fits (ampl/n in legend; raw ampl in out.ampl)');
grid on;
ylim(yl);

% ---- Bottom: residual
ax3 = subplot(3,1,3);
plot(x, real(resid), 'Color', [0.2 0.2 0.2]);
set(ax3,'XDir','reverse');
xlabel('ppm');
ylabel('a.u.');
title('residual (fit - data)');
grid on;
ylim(yl);

linkaxes([ax1 ax2 ax3], 'x');

% ---- Force vertical x-axis alignment, then float the per-component
%      legend in the gap to the right of the middle panel without
%      letting it shrink that panel's axes.
pos1 = get(ax1,'Position');
pos2 = get(ax2,'Position');
pos3 = get(ax3,'Position');

% pull the figure's right margin back so we have room for the legend
legGap   = 0.18;     % fraction of figure width reserved for legend
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
