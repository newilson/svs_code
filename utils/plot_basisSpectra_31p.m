function fig = plot_basisSpectra_31p(out, titleStr)
% PLOT_BASISSPECTRA_31P  Per-metabolite gallery of basis spectra.
%
%   fig = plot_basisSpectra_31p(out, titleStr)
%
%   One subplot per metabolite. Each shows:
%     - data on the fit window (gray)
%     - modified basis spectrum (shift/lb/phase applied) at unit amplitude,
%       rescaled so its peak matches the data peak (blue dashed) - shows
%       the basis SHAPE
%     - amplitude-scaled fitted component (red) - shows the basis
%       CONTRIBUTION to the fit
%   Title for each panel shows metabolite name and fitted amplitude.
%
%   Useful for asking "is this component doing real work, or has the spline
%   baseline absorbed it?".

if nargin < 2 || isempty(titleStr), titleStr = ''; end

x         = out.xFit;
y         = out.yFit;
comp      = out.comp;
basisSpec = out.basisSpec;
ampl      = out.ampl;
names     = out.names;

% Per-component nP for title display.
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

nMet  = numel(names);
ncols = min(nMet, 4);
nrows = ceil(nMet/ncols);

fig = figure('Position',[10 10 320*ncols 260*nrows], 'Name','31P basis spectra gallery');

ydat   = real(y);
ymax   = max(abs(ydat));
yLim   = [-0.4*ymax, 1.1*ymax];

axList = gobjects(nMet,1);
for k = 1:nMet
    ax = subplot(nrows, ncols, k);
    axList(k) = ax;

    bk = real(basisSpec(:,k));
    ck = real(comp(:,k));

    % rescale basis shape to peak of data
    pkBasis = max(abs(bk));
    if pkBasis > 0
        bkPlot = bk * (ymax / pkBasis);
    else
        bkPlot = bk;
    end

    plot(x, ydat,    'Color',[0.6 0.6 0.6], 'LineWidth', 1.0); hold on;
    % plot(x, bkPlot,  'Color',[0.1 0.4 0.8], 'LineStyle','--', 'LineWidth', 0.9);
    plot(x, ck,      'r', 'LineWidth', 1.4);

    set(ax, 'XDir','reverse');
    ylim(yLim);
    grid on;
    title(sprintf('%s  (ampl/n=%.3g)', names{k}, amplDisp(k)), 'Interpreter','none');

    if k == 1
        % legend({'data','basis shape (peak-norm)','ampl \times basis'}, 'Location','northwest', 'FontSize', 8);
        legend({'data','ampl \times basis'},'Location','northwest', 'FontSize', 8);
    end
    if k > (nrows-1)*ncols, xlabel('ppm'); end
end

linkaxes(axList, 'x');

if ~isempty(titleStr)
    try, sgtitle(titleStr, 'Interpreter','none'); catch, end
else
    try, sgtitle('Basis spectra: shape vs fitted contribution'); catch, end
end

end
