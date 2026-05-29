function fig = plot_componentSpectra_1h(x, y, out, names, titleStr)
% PLOT_COMPONENTSPECTRA_1H  Per-peak gallery of 1H component fits.
%
%   fig = plot_componentSpectra_1h(x, y, out, names)
%   fig = plot_componentSpectra_1h(x, y, out, names, titleStr)
%
%   One subplot per peak. Each shows:
%     - data on the fit window (gray)
%     - amplitude-scaled fitted component (red)
%   Title for each panel shows peak name and fitted amplitude.
%
%   Useful for asking "is this component doing real work, or has the spline
%   baseline absorbed it?".
%
%   Inputs:
%       x        - frequency axis (ppm) over the fit window [np x 1]
%       y        - data (real) over the fit window          [np x 1]
%       out      - struct returned by curvefitAuto_varproBaseline; uses
%                  .comp and .areas
%       names    - cell array of peak names {n x 1}
%       titleStr - optional super-title string

if nargin < 5 || isempty(titleStr), titleStr = ''; end

x    = x(:);
y    = real(y(:));
comp = out.comp;
ampl = out.areas(:);

if nargin < 4 || isempty(names)
    names = arrayfun(@(k) sprintf('p%d',k), 1:numel(ampl), 'UniformOutput', false);
end

nMet  = numel(names);
ncols = min(nMet, 4);
nrows = ceil(nMet/ncols);

fig = figure('Position',[10 10 320*ncols 260*nrows], 'Name','1H component spectra gallery');

ymax = max(abs(y));
if ymax == 0, ymax = 1; end
yLim = [-0.4*ymax, 1.1*ymax];

axList = gobjects(nMet,1);
for k = 1:nMet
    ax = subplot(nrows, ncols, k);
    axList(k) = ax;

    ck = real(comp(:,k));

    plot(x, y,  'Color',[0.6 0.6 0.6], 'LineWidth', 1.0); hold on;
    plot(x, ck, 'r', 'LineWidth', 1.4);

    set(ax, 'XDir','reverse');
    ylim(yLim);
    grid on;
    title(sprintf('%s  (ampl=%.3g)', names{k}, ampl(k)), 'Interpreter','none');

    if k == 1
        legend({'data','ampl \times peak'}, 'Location','northwest', 'FontSize', 8);
    end
    if k > (nrows-1)*ncols, xlabel('ppm'); end
end

linkaxes(axList, 'x');

if ~isempty(titleStr)
    try, sgtitle(titleStr, 'Interpreter','none'); catch, end
else
    try, sgtitle('1H component fits: per-peak contribution'); catch, end
end

end
