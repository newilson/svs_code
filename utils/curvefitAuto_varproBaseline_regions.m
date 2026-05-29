function output = curvefitAuto_varproBaseline_regions(x, y, peakTable, regions, shared)
% curvefitAuto_varproBaseline_regions - fit several independent spectral
% regions, each as a separate VARPRO peak + spline-baseline fit.
%
% This is a thin multi-region wrapper around curvefitAuto_varproBaseline.
% Each region is fit completely independently: its own ppm window, its own
% subset of peaks, its own lineshape mode, phase bounds, and baseline /
% lineshape options.  No parameters are shared across regions and there is
% no joint solve.  Per-region settings are formed by merging a common
% SHARED block with each region's overrides, so options common to every
% region are written once and only the differences are listed.
%
% Because curvefitAuto_varproBaseline has no fit-range argument (it fits all
% of x), this wrapper restricts x and y to each region's window before
% calling it - which also gives every region an independent baseline.
%
% Inputs:
%   x         - frequency axis (ppm or Hz), full length [N x 1]
%   y         - spectral data, full length [N x 1] (real part is fit)
%   peakTable - struct holding the global peak list (may be [] if every
%               region supplies its own peaks):
%                  .center0 - peak positions [n x 1]
%                  .width0  - peak FWHM      [n x 1]
%                  .amp0    - optional initial amplitudes [n x 1] (only a
%                             size/upper-bound helper; defaults to ones)
%                  .names   - optional names {n x 1} (for peakSelect by name
%                             and for labeling)
%   regions   - [R x 1] struct array, one element per region:
%                  .fitRange    - [min max] window (REQUIRED; defines region)
%                  .peakSelect  - which global peaks are active here: numeric
%                                 indices or names into peakTable.  Omit to
%                                 auto-select peaks whose center0 falls in
%                                 fitRange.
%                  .center0/.width0/.amp0/.names - optional region-local peak
%                                 list; if .center0 is given it REPLACES the
%                                 global table / peakSelect for this region.
%                  .name        - optional region label (char)
%                  .mode        - lineshape mode (overrides shared.mode)
%                  .ph_range    - [min max] phase bounds (overrides shared)
%                  .minw        - minimum linewidth (overrides shared)
%                  .baselineOpt - overrides merged onto shared.baselineOpt
%                  .lineShapeOpt- overrides merged onto shared.lineShapeOpt
%   shared    - optional struct of defaults applied to every region:
%                  .mode, .ph_range, .minw, .baselineOpt, .lineShapeOpt
%               Omit or pass [] for no shared defaults.
%
% Output (fields):
%   .regions       - [R x 1] struct array; regions(r) is the full output of
%                    curvefitAuto_varproBaseline for region r (unchanged
%                    single-region contract)
%   .fit           - stitched fitted spectrum, full length [N x 1], NaN
%                    outside every region window
%   .fit_peaks     - stitched peak-only fit [N x 1], NaN outside regions
%   .baseline      - stitched spline baseline [N x 1], NaN outside regions
%   .residualFit   - stitched (fit - real(data)) [N x 1], NaN outside regions
%   .fitIdx        - union logical mask of points used by any region
%   .x, .y         - echoes of the full input axis / spectrum
%   .areas         - peak areas from all regions, concatenated [M x 1]
%   .center        - fitted peak positions, concatenated [M x 1]
%   .names         - matching peak names [M x 1 cell]
%   .region        - region index for each peak row [M x 1]
%   .regionName    - region label for each peak row [M x 1 cell]
%   .regionRanges  - [R x 2] fit ranges
%   .peakTable, .sharedOpt, .regionsSpec - echoes of inputs
%
% Single-region behaviour is unchanged: to fit one window, call
% curvefitAuto_varproBaseline directly.  This wrapper only loops it.
%
% NW - 05/2026

% ------------------------------------------------------------------- %
% input checks
% ------------------------------------------------------------------- %
x = x(:);
y = y(:);
N = length(x);
if length(y) ~= N
    error('x and y length mismatch');
end
if nargin < 4 || isempty(regions)
    error('regions struct array is required (use curvefitAuto_varproBaseline for a single region)');
end
if ~isstruct(regions)
    error('regions must be a struct array');
end
regions = regions(:);
R = numel(regions);

if nargin < 3 || isempty(peakTable)
    peakTable = struct();
end
if nargin < 5 || isempty(shared)
    shared = struct();
end
sharedBase = getOptStruct(shared, 'baselineOpt');
sharedLine = getOptStruct(shared, 'lineShapeOpt');

% ------------------------------------------------------------------- %
% per-region fits
% ------------------------------------------------------------------- %
regOut   = cell(R, 1);
areasAll = [];
areasConvAll = [];
centerAll = [];
namesAll = {};
regIdx   = [];
regNameAll = {};
fitFull      = nan(N, 1);
peaksFull    = nan(N, 1);
baseFull     = nan(N, 1);
residFull    = nan(N, 1);
unionMask    = false(N, 1);
regionRanges = nan(R, 2);

for r = 1:R
    reg = regions(r);
    if ~isfield(reg, 'fitRange') || isempty(reg.fitRange) || numel(reg.fitRange) ~= 2
        error('regions(%d).fitRange must be a [min max] pair', r);
    end
    lo = min(reg.fitRange);
    hi = max(reg.fitRange);
    regionRanges(r, :) = [lo hi];

    % restrict data to this region's window
    idx = (x >= lo) & (x <= hi);
    if nnz(idx) < 10
        error('regions(%d).fitRange produced <10 points; check axis units', r);
    end
    xr = x(idx);
    yr = y(idx);

    % resolve this region's peaks (region-local table, peakSelect, or auto)
    [center0, width0, amp0, pkNames] = resolveRegionPeaks(reg, peakTable, lo, hi, r);

    % per-region scalar options (shared default, region override)
    mode     = getField(reg, 'mode',     getField(shared, 'mode',     5));
    ph_range = getField(reg, 'ph_range', getField(shared, 'ph_range', []));
    minw     = getField(reg, 'minw',     getField(shared, 'minw',     0));

    % per-region struct options (merge shared with region overrides)
    bo  = mergeStruct(sharedBase, getOptStruct(reg, 'baselineOpt'));
    lso = mergeStruct(sharedLine, getOptStruct(reg, 'lineShapeOpt'));

    % run the single-region fit (existing, untouched function)
    out = curvefitAuto_varproBaseline(xr, yr, mode, amp0, center0, width0, ...
        ph_range, minw, bo, lso);
    regOut{r} = out;

    % place this region's curves into the stitched full-length arrays
    if any(unionMask & idx)
        warning('curvefitAuto_varproBaseline_regions:overlap', ...
            'Region %d overlaps a previous region; later region wins in stitched output.', r);
    end
    fitFull(idx)   = out.fit;
    peaksFull(idx) = out.fit_peaks;
    baseFull(idx)  = out.baseline;
    residFull(idx) = out.fit - real(yr);
    unionMask = unionMask | idx;

    % concatenate area / position / name tables, tagged by region.  In every
    % mode, pars layout begins [amp(1:n) pos(n+1:2n) ...], so positions are
    % the second block; areas come straight from out.areas.
    nLoc = numel(center0);
    pos  = out.pars(nLoc + (1:nLoc));
    areasAll  = [areasAll;  out.areas(:)];               %#ok<AGROW>
    if isfield(out, 'areasConv')
        areasConvAll = [areasConvAll; out.areasConv(:)]; %#ok<AGROW>
    else
        areasConvAll = [areasConvAll; out.areas(:)];     %#ok<AGROW>
    end
    centerAll = [centerAll; pos(:)];                     %#ok<AGROW>
    namesAll  = [namesAll;  pkNames(:)];                 %#ok<AGROW>
    regIdx    = [regIdx;    r * ones(nLoc, 1)];          %#ok<AGROW>
    rname     = getField(reg, 'name', sprintf('region%02d', r));
    regNameAll = [regNameAll; repmat({rname}, nLoc, 1)]; %#ok<AGROW>
end

% ------------------------------------------------------------------- %
% assemble output
% ------------------------------------------------------------------- %
output.regions      = [regOut{:}];
output.fit          = fitFull;
output.fit_peaks    = peaksFull;
output.baseline     = baseFull;
output.residualFit  = residFull;
output.fitIdx       = unionMask;
output.x            = x;
output.y            = y;
output.areas        = areasAll;
output.areasConv    = areasConvAll;   % kernel-aware (see curvefitAuto_varproBaseline)
output.center       = centerAll;
output.names        = namesAll;
output.region       = regIdx;
output.regionName   = regNameAll;
output.regionRanges = regionRanges;
output.peakTable    = peakTable;
output.sharedOpt    = shared;
output.regionsSpec  = regions;

end % ---- end of main function ------------------------------------- %


% =================================================================== %
% local helpers
% =================================================================== %

function [center0, width0, amp0, names] = resolveRegionPeaks(reg, peakTable, lo, hi, r)
% Determine the peaks active in a region.  A region-local .center0 takes
% precedence; otherwise peaks are pulled from the global peakTable using
% .peakSelect (indices or names) or, if absent, auto-selected by center.
    if isfield(reg, 'center0') && ~isempty(reg.center0)
        center0 = reg.center0(:);
        if ~isfield(reg, 'width0') || numel(reg.width0) ~= numel(center0)
            error('regions(%d): region-local width0 must match center0 length', r);
        end
        width0 = reg.width0(:);
        amp0   = getField(reg, 'amp0', []);
        names  = getField(reg, 'names', {});
        src    = 'region';
        sel    = (1:numel(center0)).';
    else
        center0g = getPT(peakTable, 'center0');
        width0g  = getPT(peakTable, 'width0');
        if isempty(center0g)
            error('regions(%d): no region-local peaks and peakTable.center0 is empty', r);
        end
        sel = resolvePeakSelect(getField(reg, 'peakSelect', []), peakTable, lo, hi, r);
        center0 = center0g(sel);
        width0  = width0g(sel);
        amp0    = getPT(peakTable, 'amp0');
        if ~isempty(amp0), amp0 = amp0(sel); end
        names   = getPT(peakTable, 'names');
        if ~isempty(names), names = names(sel); end
        src     = 'table';
    end

    n = numel(center0);
    if numel(width0) ~= n
        error('regions(%d): width0 length does not match center0', r);
    end

    % amp0 is only a size/upper-bound helper inside the fitter; default ones.
    if isempty(amp0)
        amp0 = ones(n, 1);
    else
        amp0 = amp0(:);
        if numel(amp0) ~= n
            error('regions(%d): amp0 length does not match center0', r);
        end
    end

    % names: default to source-aware labels.
    if isempty(names)
        names = cell(n, 1);
        for k = 1:n
            if strcmp(src, 'table')
                names{k} = sprintf('peak%02d', sel(k));
            else
                names{k} = sprintf('r%d_p%02d', r, k);
            end
        end
    else
        names = names(:);
        if numel(names) ~= n
            error('regions(%d): names length does not match center0', r);
        end
    end
end

function sel = resolvePeakSelect(peakSelect, peakTable, lo, hi, r)
% Resolve a region's peakSelect into row indices into the global peakTable.
    centers = getPT(peakTable, 'center0');
    nAll = numel(centers);
    if isempty(peakSelect)
        sel = find(centers >= lo & centers <= hi);   % auto by center-in-range
        if isempty(sel)
            error(['regions(%d): no peakTable centers fall within fitRange; ' ...
                   'set regions(%d).peakSelect explicitly'], r, r);
        end
        sel = sel(:);
        return
    end
    if isnumeric(peakSelect)
        sel = peakSelect(:);
        if any(sel < 1) || any(sel > nAll) || any(sel ~= round(sel))
            error('regions(%d): peakSelect indices must be integers in 1..%d', r, nAll);
        end
        return
    end
    if ischar(peakSelect) || isstring(peakSelect)
        peakSelect = cellstr(peakSelect);
    end
    if ~iscell(peakSelect)
        error('regions(%d): peakSelect must be names, indices, or empty', r);
    end
    names = getPT(peakTable, 'names');
    if isempty(names)
        error('regions(%d): peakSelect given by name, but peakTable has no .names', r);
    end
    sel = zeros(numel(peakSelect), 1);
    for k = 1:numel(peakSelect)
        ii = find(strcmpi(names, char(peakSelect{k})), 1);
        if isempty(ii)
            error('regions(%d): peakSelect name "%s" not found in peakTable.names', ...
                r, char(peakSelect{k}));
        end
        sel(k) = ii;
    end
end

function v = getPT(peakTable, name)
    if isstruct(peakTable) && isfield(peakTable, name) && ~isempty(peakTable.(name))
        v = peakTable.(name);
        if ~iscell(v), v = v(:); end
    else
        v = [];
    end
end

function s = getOptStruct(parent, name)
    if isstruct(parent) && isfield(parent, name) && isstruct(parent.(name))
        s = parent.(name);
    else
        s = struct();
    end
end

function v = getField(s, name, default)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end

function out = mergeStruct(base, override)
% Copy every field of OVERRIDE onto BASE (override wins).
    out = base;
    if ~isstruct(override)
        return
    end
    fn = fieldnames(override);
    for k = 1:numel(fn)
        out.(fn{k}) = override.(fn{k});
    end
end
