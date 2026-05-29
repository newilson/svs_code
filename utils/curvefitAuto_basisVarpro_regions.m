function output = curvefitAuto_basisVarpro_regions(x, y, basisFIDs, basisInfo, timeInfo, regions, shared)
% curvefitAuto_basisVarpro_regions - fit several independent spectral
% regions, each as a separate LCModel/Osprey-style basis VARPRO fit.
%
% This is a thin multi-region wrapper around curvefitAuto_basisVarpro.  Each
% region is fit completely independently (its own ppm window, its own basis
% subset, and its own fit / baseline / lineshape / parametric options); no
% parameters are shared across regions and there is no joint solve.  The
% per-region option sets are formed by merging a common SHARED block with
% each region's overrides, so settings common to every region are written
% once and only the per-region differences are listed.
%
% Inputs:
%   x          - full ppm axis [N x 1] (same as curvefitAuto_basisVarpro)
%   y          - full measured spectrum [N x 1]
%   basisFIDs  - simulated basis FIDs [Nb x nMet] (pass [] if every region
%                uses only parametric peaks)
%   basisInfo  - struct array [nMet x 1] describing each basis (.name etc.)
%   timeInfo   - acquisition metadata struct (.dwellTime, .f0, ...); shared
%                by all regions
%   regions    - [R x 1] struct array, one element per region:
%                  .fitRange     - [min max] ppm window (REQUIRED; this is
%                                  the authoritative region definition)
%                  .basisSelect  - which basis columns are active in this
%                                  region: cell of basisInfo names, or
%                                  numeric indices into basisFIDs columns.
%                                  Default = all basis columns.  (Parametric
%                                  peaks are region-local via .parametricOpt
%                                  and are not affected by basisSelect.)
%                  .name         - optional region label (char)
%                  .fitOpt       - per-region overrides merged onto
%                                  shared.fitOpt
%                  .baselineOpt  - per-region overrides merged onto
%                                  shared.baselineOpt
%                  .lineShapeOpt - per-region overrides merged onto
%                                  shared.lineShapeOpt
%                  .parametricOpt- per-region overrides merged onto
%                                  shared.parametricOpt
%   shared     - optional struct of defaults applied to every region, with
%                fields .fitOpt, .baselineOpt, .lineShapeOpt, .parametricOpt
%                (each optional).  Omit or pass [] for no shared defaults.
%
% Output (fields):
%   .regions       - [R x 1] struct array; regions(r) is the full output of
%                    curvefitAuto_basisVarpro for region r (unchanged
%                    single-region contract)
%   .fit           - stitched fitted spectrum, full length [N x 1], NaN
%                    outside every region window
%   .fit_peaks     - stitched peak-only fit [N x 1], NaN outside regions
%   .baseline      - stitched spline baseline [N x 1], NaN outside regions
%   .residualFit   - stitched (fit - data) [N x 1], NaN outside regions
%   .fitIdx        - union logical mask of all points used by any region
%   .x, .y         - echoes of the full input axis / spectrum
%   .ampl          - amplitudes from all regions, concatenated [M x 1]
%   .names         - matching component names [M x 1 cell]
%   .region        - region index for each amplitude row [M x 1]
%   .regionName    - region label for each amplitude row [M x 1 cell]
%   .regionRanges  - [R x 2] fit ranges
%   .timeInfo, .sharedOpt, .regionsSpec - echoes of inputs
%
% Single-region behaviour is unchanged: to fit one window, call
% curvefitAuto_basisVarpro directly.  This wrapper never edits or bypasses
% it - it only loops it.
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
if nargin < 5 || isempty(timeInfo)
    error('timeInfo (with .dwellTime and .f0) is required');
end
if nargin < 6 || isempty(regions)
    error('regions struct array is required (use curvefitAuto_basisVarpro for a single region)');
end
if ~isstruct(regions)
    error('regions must be a struct array');
end
regions = regions(:);
R = numel(regions);

if nargin < 7 || isempty(shared)
    shared = struct();
end
sharedFit   = getOptStruct(shared, 'fitOpt');
sharedBase  = getOptStruct(shared, 'baselineOpt');
sharedLine  = getOptStruct(shared, 'lineShapeOpt');
sharedParam = getOptStruct(shared, 'parametricOpt');

nB = size(basisFIDs, 2);

% ------------------------------------------------------------------- %
% per-region fits
% ------------------------------------------------------------------- %
regOut   = cell(R, 1);
amplAll  = [];
namesAll = {};
regIdx   = [];
regNameAll = {};
fitFull       = nan(N, 1);
peaksFull     = nan(N, 1);
baseFull      = nan(N, 1);
residFull     = nan(N, 1);
unionMask     = false(N, 1);
regionRanges  = nan(R, 2);

for r = 1:R
    reg = regions(r);
    if ~isfield(reg, 'fitRange') || isempty(reg.fitRange) || numel(reg.fitRange) ~= 2
        error('regions(%d).fitRange must be a [min max] pair', r);
    end
    regionRanges(r, :) = [min(reg.fitRange) max(reg.fitRange)];

    % which basis columns are active in this region
    sel = resolveBasisSelect(getField(reg, 'basisSelect', []), basisInfo, nB);
    if isempty(basisFIDs)
        bf = [];
        bi = basisInfo;
    else
        bf = basisFIDs(:, sel);
        if isempty(basisInfo) || numel(basisInfo) < max([sel(:); 0])
            bi = basisInfo;            % let the fitter assign default names
        else
            bi = basisInfo(sel);
        end
    end

    % merge shared defaults with this region's overrides
    fo  = mergeStruct(sharedFit,   getOptStruct(reg, 'fitOpt'));
    bo  = mergeStruct(sharedBase,  getOptStruct(reg, 'baselineOpt'));
    lso = mergeStruct(sharedLine,  getOptStruct(reg, 'lineShapeOpt'));
    po  = mergeStruct(sharedParam, getOptStruct(reg, 'parametricOpt'));
    fo.fitRange = reg.fitRange;        % top-level fitRange is authoritative

    % run the single-region fit (existing, untouched function)
    out = curvefitAuto_basisVarpro(x, y, bf, bi, timeInfo, fo, bo, lso, po);
    regOut{r} = out;

    % place this region's curves into the stitched full-length arrays
    idx = out.fitIdx(:);
    if any(unionMask & idx)
        warning('curvefitAuto_basisVarpro_regions:overlap', ...
            'Region %d overlaps a previous region; later region wins in stitched output.', r);
    end
    fitFull(idx)   = out.fit;
    peaksFull(idx) = out.fit_peaks;
    baseFull(idx)  = out.baseline;
    residFull(idx) = out.residualFit;
    unionMask = unionMask | idx;

    % concatenate amplitude / name tables, tagged by region
    nLoc = numel(out.ampl);
    amplAll  = [amplAll;  out.ampl(:)];                 %#ok<AGROW>
    namesAll = [namesAll; out.names(:)];                %#ok<AGROW>
    regIdx   = [regIdx;   r * ones(nLoc, 1)];           %#ok<AGROW>
    rname    = getField(reg, 'name', sprintf('region%02d', r));
    regNameAll = [regNameAll; repmat({rname}, nLoc, 1)];%#ok<AGROW>
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
output.ampl         = amplAll;
output.names        = namesAll;
output.region       = regIdx;
output.regionName   = regNameAll;
output.regionRanges = regionRanges;
output.timeInfo     = timeInfo;
output.sharedOpt    = shared;
output.regionsSpec  = regions;

end % ---- end of main function ------------------------------------- %


% =================================================================== %
% local helpers
% =================================================================== %

function s = getOptStruct(parent, name)
% Return parent.(name) if it is a struct, else an empty struct().
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
% Copy every field of OVERRIDE onto BASE (override wins).  BASE is the
% shared default; OVERRIDE is the per-region customization.
    out = base;
    if ~isstruct(override)
        return
    end
    fn = fieldnames(override);
    for k = 1:numel(fn)
        out.(fn{k}) = override.(fn{k});
    end
end

function sel = resolveBasisSelect(basisSelect, basisInfo, nB)
% Resolve a region's basisSelect (names, indices, or empty) into a vector
% of column indices into basisFIDs.
    if isempty(basisSelect)
        sel = 1:nB;
        return
    end
    if isnumeric(basisSelect)
        sel = basisSelect(:).';
        if any(sel < 1) || any(sel > nB) || any(sel ~= round(sel))
            error('basisSelect indices must be integers in 1..%d', nB);
        end
        return
    end
    if ischar(basisSelect) || isstring(basisSelect)
        basisSelect = cellstr(basisSelect);
    end
    if ~iscell(basisSelect)
        error('basisSelect must be names (cell/char), numeric indices, or empty');
    end
    if isempty(basisInfo) || ~isfield(basisInfo, 'name')
        error('basisSelect given by name, but basisInfo has no .name field');
    end
    allNames = {basisInfo.name};
    sel = zeros(1, numel(basisSelect));
    for k = 1:numel(basisSelect)
        idx = find(strcmpi(allNames, char(basisSelect{k})), 1);
        if isempty(idx)
            error('basisSelect name "%s" not found in basisInfo', char(basisSelect{k}));
        end
        sel(k) = idx;
    end
end
