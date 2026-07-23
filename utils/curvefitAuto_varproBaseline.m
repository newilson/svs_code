
function output = curvefitAuto_varproBaseline(x,y,mode,amp0,center0,width0,ph_range,minw,baselineOpt,lineShapeOpt)
% curvefitAuto_varproBaseline - nonlinear curve fitting for
% multiple peaks with variable-projection spline baseline fitting and an
% optional shared extra lineshape distortion kernel (LCModel/Osprey-like).
%
% Inputs:
%   x        - frequency axis (ppm or Hz)
%   y        - spectral data to fit (currently fits real(y))
%   mode     - lineshape model:
%                1 = Gaussian
%                2 = real Lorentzian
%                3 = complex Lorentzian
%                4 = complex Gaussian
%                5 = complex Voigt (default)
%                6 = magnitude Lorentzian
%                7 = magnitude Voigt
%   amp0     - initial amplitudes [n x 1]
%              (used only as starting guess / upper bound helper)
%   center0  - initial peak positions [n x 1]
%   width0   - initial linewidths FWHM [n x 1]
%   ph_range - [min max] phase bounds in degrees (modes 3,4,5 only)
%   minw     - minimum linewidth (optional, default 0)
%   baselineOpt - struct with spline-baseline options:
%                .enable      = true/false     (default true)
%                .style       = 'osprey' or 'lcmodel' (default 'osprey')
%                .knotSpacing = spline knot spacing in x units
%                .lambda      = spline roughness penalty
%                .normalize   = normalize spline columns (default true)
%                .amplUB      = amplitude upper bound (scalar or [n x 1])
%                .linearTol   = lsqlin tolerance
%                .verbose     = display iterations
%                .weights     = weight vector [np x 1] for weighted
%                               least squares (default ones = uniform)
%                .widthBounds = [min max] linewidth bounds in x-units
%                               (scalar row or [n x 2] per-peak)
%                .priorBaseline / .lambdaPrior = soft prior pulling the spline
%                               toward a supplied baseline estimate
%                .priorWidth  / .lambdaWidth  = soft prior pulling each peak's
%                               total linewidth toward a supplied estimate
%                               (modes 5/7 only).  See parseBaselineOpt.
%   lineShapeOpt - struct with shared additional lineshape options:
%                .enable      = true/false     (default false)
%                .nSide       = number of kernel points on one side
%                               of the center point (default 2)
%                .asymmetric  = true/false     (default false)
%                .maxSide     = upper bound for each side coefficient
%                               (default 0.50)
%                .normalize   = normalize kernel area to 1 (default true)
%
% Output:
%   output.pars              - optimized parameters (same style as
%                              curvefitAuto + appended lineshape pars)
%   output.parnames          - parameter names (describes order of pars)
%   output.fit               - total fitted spectrum [npts x 1]
%   output.fit_peaks         - peak-only fit [npts x 1]
%   output.baseline          - fitted spline baseline [npts x 1]
%   output.comp              - individual peak components [npts x n]
%   output.areas             - peak areas [1 x n]
%   output.beta              - spline coefficients
%   output.B                 - spline basis matrix
%   output.lineShapePars     - nonlinear side coefficients of extra kernel
%   output.lineShapeKernel   - shared extra kernel actually used
%   output.lb              	- lower bounds (same order as pars)
%   output.ub              	- upper bounds (same order as pars)
%   output.mode              - fitting mode
%
% Notes:
%   1) This is a variable-projection version: peak amplitudes and spline
%      coefficients are solved linearly at each nonlinear iteration.
%   2) The nonlinear optimizer fits only position / width / phase terms and,
%      optionally, the side coefficients of a shared extra lineshape kernel.
%   3) The extra lineshape kernel is convolved with metabolite peak basis
%      functions only, not with the spline baseline, similar in spirit to
%      LCModel/Osprey supplementary lineshape distortion handling.
%   4) For now, the fitting is performed on real(y), matching the original
%      curvefitAuto style that returns real-valued spectra in all modes.
%
% NW - 03/2026

x = x(:); y = y(:);
np = length(x);
if (length(y) ~= np)
    error('x and y array size mismatch')
end

if ~isreal(y)
    warning('curvefitAuto_varproBaseline currently fits real(y) only.')
    y = real(y);
end

if nargin<8 || isempty(minw)
    minw = 0;
end

if (nargin < 3) || isempty(mode)
    mode = 5; % complex Voigt
end

if nargin<9 || isempty(baselineOpt)
    baselineOpt = struct();
end
baselineOpt = parseBaselineOpt(x, baselineOpt);

if nargin<10 || isempty(lineShapeOpt)
    lineShapeOpt = struct();
end
lineShapeOpt = parseLineShapeOpt(lineShapeOpt);

% Weight vector for weighted least squares
if ~isempty(baselineOpt.weights)
    w = baselineOpt.weights(:);
    if length(w) ~= np
        error('baselineOpt.weights must have the same length as x (%d)', np)
    end
else
    w = ones(np, 1);
end

% rms weighted signal -- the scale the soft priors are normalised against, so
% that lambdaPrior / lambdaWidth are dimensionless and portable between
% datasets with different absolute scaling.
rmsWY = sqrt(mean((w .* y).^2));
if ~isfinite(rmsWY) || rmsWY <= 0
    rmsWY = 1;
end

% Initial guesses
if nargin<6 || isempty(amp0) || isempty(center0) || isempty(width0)
    [~,center,width,amp] = findpeaks(y,x);
    figure('Name','Auto peak detection (findpeaks)'), plot(x,y,'-k',center,amp,'*r')
    n = length(center);
else
    amp = amp0(:)';
    center = center0(:)';
    width = width0(:)';
    n = length(amp);
end

if length(center) ~= n || length(width) ~= n
    error('amp0, center0, and width0 array size mismatch')
end

if isfield(baselineOpt,'centerBounds') && ~isempty(baselineOpt.centerBounds)
    centermin = baselineOpt.centerBounds(:,1)';
    centermax = baselineOpt.centerBounds(:,2)';
else
    centermax = center+width/2;
    centermin = center-width/2;
end

% Optional per-peak or global linewidth bounds [min max] in x-units (ppm)
% Overrides the default width bounds derived from initial width guess.
if isfield(baselineOpt,'widthBounds') && ~isempty(baselineOpt.widthBounds)
    wb = baselineOpt.widthBounds;
    if size(wb,1) == 1  % scalar [min max] applied to all peaks
        widthBoundsMin = wb(1) * ones(1, n);
        widthBoundsMax = wb(2) * ones(1, n);
    else  % per-peak [n x 2]
        widthBoundsMin = wb(:,1)';
        widthBoundsMax = wb(:,2)';
    end
    hasWidthBounds = true;
else
    hasWidthBounds = false;
end

amin = zeros(size(amp));
amax = max(abs(y)) * ones(size(amp));

% Optional linear amplitude upper bounds
if isempty(baselineOpt.amplUB)
    amplUB = amax(:);
elseif isscalar(baselineOpt.amplUB)
    amplUB = baselineOpt.amplUB * ones(n,1);
else
    amplUB = baselineOpt.amplUB(:);
    if length(amplUB) ~= n
        error('baselineOpt.amplUB must be scalar or length n')
    end
end

% Olivero-Longbothum approximation coefficients for Voigt
Acoef = 0.5346; Bcoef = 0.2166;

% Build spline basis once
if baselineOpt.enable
    [B, knots] = makeSplineBasis(x, baselineOpt.knotSpacing, baselineOpt.normalize);
else
    B = zeros(np,0);
    knots = [];
end

%% parameter setup
if mode==5 % complex Voigt

    if hasWidthBounds
        widthmin = widthBoundsMin;
        widthmax = widthBoundsMax;
    else
        widthmin = max(1e-6*width, minw);
        widthmax = 2*width;
    end

    widthL = 0.4 * width;

    term = (width - Acoef*widthL).^2 - Bcoef*(widthL.^2);
    widthG = sqrt(max(0, term));

    pars(1:n) = center;
    lb(1:n) = centermin;
    ub(1:n) = centermax;

    pars(n + (1:n)) = widthL;
    lb(n + (1:n)) = widthmin;
    ub(n + (1:n)) = widthmax;

    if nargin<7 || isempty(ph_range) || length(ph_range)~=2
        pars(2*n + (1:n)) = 0;
        lb(2*n + (1:n)) = -179;
        ub(2*n + (1:n)) = 179;
    else
        pars(2*n + (1:n)) = mean(ph_range);
        lb(2*n + (1:n)) = min(ph_range);
        ub(2*n + (1:n)) = max(ph_range);
    end

    pars(3*n + (1:n)) = widthG;
    lb(3*n + (1:n)) = widthmin;
    ub(3*n + (1:n)) = widthmax;

elseif mode==7 % magnitude Voigt

    if hasWidthBounds
        widthmin = widthBoundsMin;
        widthmax = widthBoundsMax;
    else
        widthmin = max(1e-6*width, minw);
        widthmax = 2*width;
    end

    widthL = 0.4 * width;

    term = (width - Acoef*widthL).^2 - Bcoef*(widthL.^2);
    widthG = sqrt(max(0, term));

    pars(1:n) = center;
    lb(1:n) = centermin;
    ub(1:n) = centermax;

    pars(n + (1:n)) = widthL;
    lb(n + (1:n)) = widthmin;
    ub(n + (1:n)) = widthmax;

    pars(2*n + (1:n)) = widthG;
    lb(2*n + (1:n)) = widthmin;
    ub(2*n + (1:n)) = widthmax;

else
    if hasWidthBounds
        widthmin = widthBoundsMin;
        widthmax = widthBoundsMax;
    else
        widthmin = 0.5*width;
        widthmax = 2*width;
        for k = 1:n
            if (widthmin(k) < minw)
                widthmin(k) = minw;
            end
        end
    end

    pars(1:n) = center;
    lb(1:n) = centermin;
    ub(1:n) = centermax;

    pars(n + (1:n)) = width;
    lb(n + (1:n)) = widthmin;
    ub(n + (1:n)) = widthmax;

    if mode>2 && mode~=6 && mode~=7
        if nargin<7 || isempty(ph_range) || length(ph_range)~=2
            pars(2*n + (1:n)) = 0;
            lb(2*n + (1:n)) = -179;
            ub(2*n + (1:n)) = 179;
        else
            pars(2*n + (1:n)) = mean(ph_range);
            lb(2*n + (1:n)) = min(ph_range);
            ub(2*n + (1:n)) = max(ph_range);
        end
    end
end

% Optional per-peak phase bounds [n x 2] (deg), overriding the uniform ph_range
% for the phase block (modes 3/4/5).  Used e.g. to keep broad macromolecular
% components near-absorptive while leaving the narrow peaks free.
if isfield(baselineOpt,'phaseBounds') && ~isempty(baselineOpt.phaseBounds) ...
        && (mode==3 || mode==4 || mode==5)
    pb = baselineOpt.phaseBounds;
    if size(pb,1)==1, pb = repmat(pb, n, 1); end
    idx = 2*n + (1:n);
    lb(idx) = pb(:,1).';
    ub(idx) = pb(:,2).';
    pars(idx) = min(max(pars(idx), lb(idx)), ub(idx));
end

% Optional phase initialization (deg), overriding the default
% mean(ph_range).  Scalar -> applied to all peaks; vector of length n ->
% per-peak.  Bounds (lb/ub) are unchanged so phase remains free to move.
% Used by the two-stage water fit in bru2mat_svsNAD to warm-start stage 2
% at the converged stage-1 phase, so the optimizer starts in the right
% basin without having to rotate the input data.
if isfield(baselineOpt,'phaseInit') && ~isempty(baselineOpt.phaseInit) ...
        && (mode==3 || mode==4 || mode==5)
    pi0 = baselineOpt.phaseInit(:).';
    if isscalar(pi0), pi0 = pi0 * ones(1, n); end
    assert(numel(pi0) == n, 'baselineOpt.phaseInit must be scalar or length n=%d', n);
    idx = 2*n + (1:n);
    pars(idx) = min(max(pi0, lb(idx)), ub(idx));
end

% Append optional shared lineshape distortion parameters
nPeakPars = length(pars);
[lineShapePars0, lineShapeLb, lineShapeUb, lineShapeParnames] = initLineShapePars(lineShapeOpt);
if ~isempty(lineShapePars0)
    pars = [pars lineShapePars0];
    lb   = [lb   lineShapeLb];
    ub   = [ub   lineShapeUb];
end

% safety: keep every initial guess strictly within its bounds.  Broad
% components (wide widthBounds) or edge cases can otherwise seed widthL/widthG
% outside [widthmin widthmax]; this is a no-op when seeds are already valid.
pars = min(max(pars, lb), ub);

%% solve
if size(x)==size(y'), y = y'; end

oldoptions = optimset('lsqnonlin');
% Default relative tolerances â€” scaled by SSR of the data so behavior
% is independent of data amplitude.  Override with baselineOpt fields.
ssrData = sum((w .* y).^2);
relTolFun = 1e-6;  relTolX = 1e-9; % default tolerances
if isfield(baselineOpt,'TolFun'), relTolFun = baselineOpt.TolFun; end
if isfield(baselineOpt,'TolX'),   relTolX   = baselineOpt.TolX;   end
tolFun = relTolFun * ssrData;
tolX   = relTolX;
options = optimset(oldoptions, 'TolFun', tolFun,'TolX', tolX, ...
    'MaxFunEval',20000*n,'MaxIter',12000,'Display',getDisplayFlag(baselineOpt.verbose));

if (mode == 1) % Gaussian

    parnames = {'amp','pos','width'};
    ip_nonlin = lsqnonlin(@residual_gaussian,pars,lb,ub,options);
    Apeak = gaussian_basis(ip_nonlin, x);
    [ampl, beta, baseline, yfit, y_comp] = solve_linear_block(Apeak, B, y);
    pos = ip_nonlin(1:n);
    width = abs(ip_nonlin((1:n)+n));
    ip = [ampl(:).' pos width];
    for k = 1:n
        areas(k) = ampl(k) * width(k) * 0.6006 * sqrt(pi);
    end
    lbout = [amin centermin widthmin];
    ubout = [amplUB(:).' centermax widthmax];

elseif mode==2 % Lorentzian

    parnames = {'amp','pos','width'};
    ip_nonlin = lsqnonlin(@residual_lorentzian,pars,lb,ub,options);
    Apeak = lorentzian_basis(ip_nonlin, x);
    [ampl, beta, baseline, yfit, y_comp] = solve_linear_block(Apeak, B, y);
    pos = ip_nonlin(1:n);
    width = abs(ip_nonlin((1:n)+n));
    ip = [ampl(:).' pos width];
    for k = 1:n
        areas(k) = ampl(k) * width(k) * 1.5708;
    end
    lbout = [amin centermin widthmin];
    ubout = [amplUB(:).' centermax widthmax];

elseif mode==3 % complex Lorentzian

    parnames = {'amp','pos','width','ph'};
    ip_nonlin = lsqnonlin(@residual_lorentzian_complex,pars,lb,ub,options);
    Apeak = lorentzian_complex_basis(ip_nonlin, x);
    [ampl, beta, baseline, yfit, y_comp] = solve_linear_block(Apeak, B, y);
    pos = ip_nonlin(1:n);
    width = abs(ip_nonlin((1:n)+n));
    phase = ip_nonlin((1:n)+2*n);
    ip = [ampl(:).' pos width phase];
    areas = ampl(:).' .* width * pi/2;
    lbout = [amin centermin widthmin lb(2*n+(1:n))];
    ubout = [amplUB(:).' centermax widthmax ub(2*n+(1:n))];

elseif mode==4 % complex Gaussian

    parnames = {'amp','pos','width','ph'};
    ip_nonlin = lsqnonlin(@residual_gaussian_complex,pars,lb,ub,options);
    Apeak = gaussian_complex_basis(ip_nonlin, x);
    [ampl, beta, baseline, yfit, y_comp] = solve_linear_block(Apeak, B, y);
    pos = ip_nonlin(1:n);
    width = abs(ip_nonlin((1:n)+n));
    phase = ip_nonlin((1:n)+2*n);
    ip = [ampl(:).' pos width phase];
    areas = ampl(:).' .* width * 0.6006 * sqrt(pi);
    lbout = [amin centermin widthmin lb(2*n+(1:n))];
    ubout = [amplUB(:).' centermax widthmax ub(2*n+(1:n))];

elseif mode==5 % complex Voigt

    parnames = {'amp','pos','widthL','ph','widthG','width'};
    ip_nonlin = lsqnonlin(@residual_voigt_complex,pars,lb,ub,options);
    Apeak = voigt_complex_basis(ip_nonlin, x);
    [ampl, beta, baseline, yfit, y_comp] = solve_linear_block(Apeak, B, y);
    pos = ip_nonlin(1:n);
    widthL = abs(ip_nonlin((1:n)+n));
    phase = ip_nonlin((1:n)+2*n);
    widthG = abs(ip_nonlin((1:n)+3*n));
    width = Acoef * widthL + sqrt(Bcoef*widthL.^2 + widthG.^2);
    ip = [ampl(:).' pos widthL phase widthG width];
    areas = ampl(:).' .* width;
    lbout = [amin centermin widthmin lb(2*n+(1:n)) widthmin zeros(1,n)];
    ubout = [amplUB(:).' centermax widthmax ub(2*n+(1:n)) widthmax inf(1,n)];

elseif mode==6 % magnitude Lorentzian

    parnames = {'amp','pos','width'};
    ip_nonlin = lsqnonlin(@residual_lorentzian_magnitude,pars,lb,ub,options);
    Apeak = lorentzian_magnitude_basis(ip_nonlin, x);
    [ampl, beta, baseline, yfit, y_comp] = solve_linear_block(Apeak, B, y);
    pos = ip_nonlin(1:n);
    width = abs(ip_nonlin((1:n)+n));
    ip = [ampl(:).' pos width];
    areas = ampl(:).' .* width;
    lbout = [amin centermin widthmin];
    ubout = [amplUB(:).' centermax widthmax];

elseif mode==7 % magnitude Voigt

    parnames = {'amp','pos','widthL','widthG','width'};
    ip_nonlin = lsqnonlin(@residual_voigt_magnitude,pars,lb,ub,options);
    Apeak = voigt_magnitude_basis(ip_nonlin, x);
    [ampl, beta, baseline, yfit, y_comp] = solve_linear_block(Apeak, B, y);
    pos = ip_nonlin(1:n);
    widthL = abs(ip_nonlin((1:n)+n));
    widthG = abs(ip_nonlin((1:n)+2*n));
    width = Acoef * widthL + sqrt(Bcoef*widthL.^2 + widthG.^2);
    ip = [ampl(:).' pos widthL widthG width];
    areas = ampl(:).' .* width;
    lbout = [amin centermin widthmin widthmin zeros(1,n)];
    ubout = [amplUB(:).' centermax widthmax widthmax inf(1,n)];

end

if ~isempty(lineShapeParnames)
    ip = [ip ip_nonlin(nPeakPars+1:end)];
    lbout = [lbout lineShapeLb];
    ubout = [ubout lineShapeUb];
    parnames = [parnames lineShapeParnames];
end

% output
output.pars = ip;
output.lb = lbout;
output.ub = ubout;
output.areas = areas;
% Kernel-aware per-peak area: integrate the convolved component on an
% EXTENDED x-axis so peaks near the fit window edge keep their tails.
% Why extend?  trapz over the fit window alone under-counts a peak's
% area when the peak sits close to the window edge -- Voigt tails decay
% slowly and the lineshape kernel adds further shoulder mass.  Rebuilding
% the basis on a wider axis using the SAME fitted parameters + kernel
% recovers the missing tail without re-running the optimizer.
%
% The pad is `areaPadFactor * maxWidth` on each side; for Voigt/Gaussian
% the "width" used here is the full FWHM-like measure (widthL+widthG for
% Voigt) so the pad grows with the broadest fitted peak.  When the kernel
% is normalized to sum=1 (the default), in the infinite-window limit
% `areasConv` converges to the analytic peak integral -- the per-peak
% factor vs `areas` becomes a constant set only by the basis function's
% normalization, so cross-peak ratios are restored.
if isfield(lineShapeOpt,'areaPadFactor') && ~isempty(lineShapeOpt.areaPadFactor)
    padFactor = lineShapeOpt.areaPadFactor;
else
    padFactor = 10;
end
switch mode
    case 1, basisFn = @gaussian_basis;            wRefAll = abs(ip_nonlin((1:n)+n));
    case 2, basisFn = @lorentzian_basis;          wRefAll = abs(ip_nonlin((1:n)+n));
    case 3, basisFn = @lorentzian_complex_basis;  wRefAll = abs(ip_nonlin((1:n)+n));
    case 4, basisFn = @gaussian_complex_basis;    wRefAll = abs(ip_nonlin((1:n)+n));
    case 5, basisFn = @voigt_complex_basis;       wRefAll = Acoef*abs(ip_nonlin((1:n)+n)) + abs(ip_nonlin((1:n)+3*n));
    case 6, basisFn = @lorentzian_magnitude_basis;wRefAll = abs(ip_nonlin((1:n)+n));
    case 7, basisFn = @voigt_magnitude_basis;     wRefAll = Acoef*abs(ip_nonlin((1:n)+n)) + abs(ip_nonlin((1:n)+2*n));
    otherwise, basisFn = []; wRefAll = abs(width(:));
end
% Quantitation must be PHASE-INVARIANT: a peak's area is a property of the
% molecule, not of whatever phase the fit assigned it.  The complex basis
% carries a per-peak phase (modes 3/4/5 store it at ip(2n+(1:n))), and
% real(component) is scaled by ~cos(phase) -- so integrating the component
% AS FITTED makes the area shrink with residual phase (e.g. a Trp component
% landing at -27 deg loses ~11% of its area for no physical reason).  We
% instead integrate the ZERO-PHASE model peak: the same fitted position,
% widths and shared lineshape kernel (kept, because capturing kernel /
% Gaussian-Lorentzian shape is the whole point of an area rather than an
% amplitude), but with the phase set to 0.  The amplitude `ampl` is left as
% fit -- only the phase rotation is removed.  This makes areasConv the
% frequency-domain equivalent of |complex amplitude|.
ip_area = ip_nonlin;
if mode==3 || mode==4 || mode==5
    ip_area(2*n + (1:n)) = 0;    % zero the per-peak phase block
end
if ~isempty(basisFn)
    padX     = padFactor * max(wRefAll);
    dxIn     = mean(abs(diff(x)));
    xExt     = ((min(x) - padX) : dxIn : (max(x) + padX)).';
    ApeakExt = basisFn(ip_area, xExt);               % zero-phase; kernel applied inside
    y_comp_ext = ApeakExt .* ampl(:).';
    output.areasConv = trapz(xExt, real(y_comp_ext), 1);
    output.areasConvPad = padX;                      % for diagnostics
else
    % Unknown mode: fall back to in-window trapz of the post-kernel comp.
    % No basisFn to rebuild, so fall back to |complex component| to stay
    % phase-invariant (real(comp) would carry the cos(phase) scaling).
    output.areasConv = trapz(x(:), abs(y_comp), 1);
    output.areasConvPad = 0;
end
output.fit = yfit;
output.fit_peaks = yfit - baseline;
output.baseline = baseline;
output.comp = y_comp;
output.beta = beta;
output.B = B;
output.knots = knots;
output.lineShapePars = ip_nonlin(nPeakPars+1:end);
output.lineShapeKernel = buildLineShapeKernel(ip_nonlin);
output.lineShapeParnames = lineShapeParnames;
output.parnames = parnames;
output.mode = mode;
output.baselineOpt = baselineOpt;
output.lineShapeOpt = lineShapeOpt;

%% residual functions (nested â€” share parent workspace: x, y, w, n, B, amplUB, baselineOpt, lineShapeOpt, nPeakPars)

% Soft prior on peak linewidth, appended to the residual vector so lsqnonlin
% trades width error against data misfit instead of hitting a hard bound.
% A hard cap forces the width to the bound and dumps the discrepancy into the
% amplitude; a penalty lets a genuinely broad peak stay broad if the data
% insist, while stopping a peak from running away to absorb baseline.
% Compares the TOTAL Voigt FWHM (Olivero-Longbothum, the same relation used to
% seed widthG at setup) against the prior, since that is what an HSVD
% linewidth measures.
function rw = widthPriorResidual(ip)
    rw = zeros(0,1);
    if baselineOpt.lambdaWidth <= 0 || isempty(baselineOpt.priorWidth), return; end
    pw = baselineOpt.priorWidth(:);
    if numel(pw) ~= n, return; end
    switch mode
        case {5,7}   % complex / magnitude Voigt: [center, widthL, ph, widthG]
            wL = abs(ip(n   + (1:n)));
            wG = abs(ip(3*n + (1:n)));
            wt = Acoef*wL(:) + sqrt(Bcoef*wL(:).^2 + wG(:).^2);
        otherwise
            return;                      % width layout not verified for other modes
    end
    ok = isfinite(pw) & pw > 0;
    if ~any(ok), return; end
    % Normalise so lambdaWidth is dimensionless and comparable to the data
    % term: scale by the rms weighted signal, and by sqrt(npts/npeaks) so the
    % few width residuals carry the same total weight as the many data ones.
    sc = sqrt(baselineOpt.lambdaWidth * numel(y) / sum(ok)) * rmsWY;
    rw = sc * (wt(ok) - pw(ok)) ./ pw(ok);
end
function r = residual_gaussian(ip)
    Apeak = gaussian_basis(ip, x);
    [~, ~, ~, yfit] = solve_linear_block(Apeak, B, y);
    r = [w .* (yfit - y); widthPriorResidual(ip)];
end

function r = residual_lorentzian(ip)
    Apeak = lorentzian_basis(ip, x);
    [~, ~, ~, yfit] = solve_linear_block(Apeak, B, y);
    r = [w .* (yfit - y); widthPriorResidual(ip)];
end

function r = residual_lorentzian_complex(ip)
    Apeak = lorentzian_complex_basis(ip, x);
    [~, ~, ~, yfit] = solve_linear_block(Apeak, B, y);
    r = [w .* (yfit - y); widthPriorResidual(ip)];
end

function r = residual_gaussian_complex(ip)
    Apeak = gaussian_complex_basis(ip, x);
    [~, ~, ~, yfit] = solve_linear_block(Apeak, B, y);
    r = [w .* (yfit - y); widthPriorResidual(ip)];
end

function r = residual_voigt_complex(ip)
    Apeak = voigt_complex_basis(ip, x);
    [~, ~, ~, yfit] = solve_linear_block(Apeak, B, y);
    r = [w .* (yfit - y); widthPriorResidual(ip)];
end

function r = residual_lorentzian_magnitude(ip)
    Apeak = lorentzian_magnitude_basis(ip, x);
    [~, ~, ~, yfit] = solve_linear_block(Apeak, B, y);
    r = [w .* (yfit - y); widthPriorResidual(ip)];
end

function r = residual_voigt_magnitude(ip)
    Apeak = voigt_magnitude_basis(ip, x);
    [~, ~, ~, yfit] = solve_linear_block(Apeak, B, y);
    r = [w .* (yfit - y); widthPriorResidual(ip)];
end

%% linear solve for amplitudes + spline coefficients
function [ampl, beta, baseline, yfit, y_comp] = solve_linear_block(Apeak, Bmat, data)
    nloc = size(Apeak,2);
    nb = size(Bmat,2);

    % Apply weights to data rows of the linear system
    Xw = w .* [Apeak Bmat];
    rhsw = w .* data;

    if nb>0 && baselineOpt.lambda>0
        D = diff(speye(nb),2);
        Xw = [Xw; zeros(size(D,1),nloc) sqrt(baselineOpt.lambda)*D];
        rhsw = [rhsw; zeros(size(D,1),1)];
    end

    % Soft prior pulling the spline toward a supplied baseline estimate.
    % Penalises the CURVE difference ||w.*(B*beta - prior)||^2 rather than the
    % coefficient difference, so the strength does not depend on how the
    % B-spline basis happens to be scaled, and the prior needs to be
    % representable only to the accuracy the basis allows (lsqlin lands on the
    % projection of the prior onto the basis, which is exactly what we want).
    if nb>0 && baselineOpt.lambdaPrior>0 && ~isempty(baselineOpt.priorBaseline)
        pb = baselineOpt.priorBaseline(:);
        if numel(pb) == numel(w)
            sp = sqrt(baselineOpt.lambdaPrior);
            Xw   = [Xw;   sp * [zeros(numel(w),nloc), w .* Bmat]];
            rhsw = [rhsw; sp * (w .* pb)];
        end
    end

    % penalize peak amplitudes toward zero â€” prevents the solver from
    % assigning amplitude to peaks with no supporting signal
    if baselineOpt.lambdaAmpl > 0
        Xw = [Xw; sqrt(baselineOpt.lambdaAmpl)*eye(nloc) zeros(nloc,nb)];
        rhsw = [rhsw; zeros(nloc,1)];
    end

    lbin = [zeros(nloc,1); -inf(nb,1)];
    ubin = [amplUB(:); inf(nb,1)];

    opts = optimoptions('lsqlin','Display','off', ...
        'OptimalityTolerance',baselineOpt.linearTol, ...
        'ConstraintTolerance',baselineOpt.linearTol, ...
        'StepTolerance',baselineOpt.linearTol);

    x0 = [max(0, Apeak' * data); zeros(nb,1)];
    x0 = min(max(x0, lbin), ubin);

    chat = lsqlin(Xw, rhsw, [], [], [], [], lbin, ubin, x0, opts);

    ampl = chat(1:nloc);
    beta = chat(nloc+(1:nb));
    if nb>0
        baseline = Bmat * beta;
    else
        baseline = zeros(size(data));
    end

    y_comp = zeros(size(Apeak));
    for kk = 1:nloc
        y_comp(:,kk) = ampl(kk) * Apeak(:,kk);
    end
    yfit = sum(y_comp,2) + baseline;
end

%% unit-amplitude basis functions
function Apeak = gaussian_basis(ip,x)
    Apeak = zeros(length(x), n);
    pos = ip(1:n);
    width = abs(ip((1:n)+n));
    for k = 1:n
        temppars = [1 pos(k) width(k)];
        Apeak(:,k) = composite_unit(temppars, x);
    end
    Apeak = applyExtraLineShape(Apeak, ip);
end

function Apeak = lorentzian_basis(ip,x)
    Apeak = zeros(length(x), n);
    pos = ip(1:n);
    width = abs(ip((1:n)+n));
    for k = 1:n
        temppars = [1 pos(k) width(k)];
        Apeak(:,k) = compositel_unit(temppars, x);
    end
    Apeak = applyExtraLineShape(Apeak, ip);
end

function Apeak = lorentzian_complex_basis(ip,x)
    Apeak = zeros(length(x), n);
    pos = ip(1:n);
    width = abs(ip((1:n)+n));
    phase = ip((1:n)+2*n);
    for k = 1:n
        temppars = [1 pos(k) width(k) phase(k)];
        Apeak(:,k) = compositel_complex_unit(temppars, x);
    end
    Apeak = applyExtraLineShape(Apeak, ip);
end

function Apeak = gaussian_complex_basis(ip,x)
    Apeak = zeros(length(x), n);
    pos = ip(1:n);
    width = abs(ip((1:n)+n));
    phase = ip((1:n)+2*n);
    for k = 1:n
        temppars = [1 pos(k) width(k) phase(k)];
        Apeak(:,k) = compositeG_complex_unit(temppars, x);
    end
    Apeak = applyExtraLineShape(Apeak, ip);
end

function Apeak = voigt_complex_basis(ip,x)
    Apeak = zeros(length(x), n);
    pos = ip(1:n);
    widthL = abs(ip((1:n)+n));
    phase = ip((1:n)+2*n);
    widthG = abs(ip((1:n)+3*n));
    for k = 1:n
        temppars = [1 pos(k) widthL(k) phase(k) widthG(k)];
        Apeak(:,k) = compositeV_complex_unit(temppars, x);
    end
    Apeak = applyExtraLineShape(Apeak, ip);
end

function Apeak = lorentzian_magnitude_basis(ip,x)
    Apeak = zeros(length(x), n);
    pos = ip(1:n);
    width = abs(ip((1:n)+n));
    for k = 1:n
        temppars = [1 pos(k) width(k)];
        Apeak(:,k) = compositel_magnitude_unit(temppars, x);
    end
    Apeak = applyExtraLineShape(Apeak, ip);
end

function Apeak = voigt_magnitude_basis(ip,x)
    Apeak = zeros(length(x), n);
    pos = ip(1:n);
    widthL = abs(ip((1:n)+n));
    widthG = abs(ip((1:n)+2*n));
    for k = 1:n
        temppars = [1 pos(k) widthL(k) widthG(k)];
        Apeak(:,k) = compositeV_magnitude_unit(temppars, x);
    end
    Apeak = applyExtraLineShape(Apeak, ip);
end

function Aout = applyExtraLineShape(Ain, ip)
    Aout = Ain;
    kern = buildLineShapeKernel(ip);
    if length(kern) == 1
        return
    end
    for kk = 1:size(Ain,2)
        Aout(:,kk) = conv(Ain(:,kk), kern, 'same');
    end
end

function kern = buildLineShapeKernel(ip)
    if ~lineShapeOpt.enable || isempty(ip) || length(ip) <= nPeakPars
        kern = 1;
        return
    end

    nSide = lineShapeOpt.nSide;
    sidePars = ip(nPeakPars+1:end);
    sidePars = max(0, sidePars(:));

    if lineShapeOpt.asymmetric
        left = sidePars(1:nSide);
        right = sidePars(nSide+(1:nSide));
    else
        left = sidePars(1:nSide);
        right = left;
    end

    kern = [flipud(left(:)); 1; right(:)];

    if lineShapeOpt.normalize
        s = sum(kern);
        if ~isfinite(s) || s <= 0
            kern = zeros(size(kern));
            kern(nSide+1) = 1;
        else
            kern = kern / s;
        end
    end
end

end % end of main function curvefitAuto_varproBaseline

%% unit-amplitude line shape functions (local â€” no parent workspace access needed)
function yfit = composite_unit(ip,x)
    p = ip(2);
    w = abs(ip(3));
    xpw = (x-p)/w;
    yfit = exp( -((xpw/0.6006).^2) );
end

function yfit = compositel_unit(ip,x)
    p = ip(2);
    w = abs(ip(3));
    xpw = (x-p)/w;
    yfit = 1 ./ ( (4.0*(xpw .*xpw)) + 1.0 );
end

function yfit = compositel_complex_unit(ip,x)
    p = ip(2);
    w = abs(ip(3));
    ph = ip(4);
    xpw = (x-p)/w;
    L = 1 ./ (1 + 1i*2.0*xpw) * exp(-1i*pi/180*ph);
    yfit = real(L);
end

function yfit = compositeG_complex_unit(ip, x)
    p  = ip(2);
    w  = abs(ip(3));
    ph = ip(4) * pi/180;

    z = (1/0.6006) * (x - p) / w;
    G_real = exp(-z.^2);
    G_imag = (2/sqrt(pi)) * dawson(z);
    G = G_real + 1i*G_imag;
    yfit = real(G * exp(-1i*ph));
end

function yfit = compositeV_complex_unit(ip, x)
    p   = ip(2);
    wL  = abs(ip(3));
    ph  = ip(4) * pi/180;
    wG  = abs(ip(5));

    sigma = wG / (2*sqrt(2*log(2)));
    gamma = wL / 2;

    if sigma == 0 && gamma == 0
        yfit = zeros(size(x));
    elseif sigma == 0
        xpw = (x - p) / gamma;
        Vc = 1 ./ (1 + 1i*2*xpw);
        yfit = real(Vc * exp(-1i*ph));
    else
        z = ((x - p) + 1i*gamma) / (sigma*sqrt(2));
        Vc = fadf(z);
        z0 = (1i*gamma) / (sigma*sqrt(2));
        V0 = fadf(z0);
        N0 = real(V0);
        if ~isfinite(N0) || N0 == 0
            [~, idx0] = min(abs(x - p));
            N0 = real(Vc(idx0));
            if ~isfinite(N0) || N0 == 0
                N0 = 1;
            end
        end
        % conj(Vc): fadf returns w(z) = exp(-z^2) erfc(-iz), whose
        % imaginary part is POSITIVE for u = (x-p)/(sigma sqrt(2)) > 0.
        % The Lorentzian sigma=0 branch above (and compositel_complex_unit)
        % use the standard NMR dispersion convention where Im is NEGATIVE
        % for x > p.  Without the conjugation, the phase parameter sign
        % flipped discontinuously at sigma=0 and was inverted relative to
        % mode 3 for the same physical phase rotation.
        yfit = real((conj(Vc) / N0) * exp(-1i*ph));
    end
end

function yfit = compositel_magnitude_unit(ip,x)
    p = ip(2);
    w = abs(ip(3));
    xpw = (x-p)/w;
    yfit = 1 ./ sqrt(1 + 4.0*(xpw.^2));
end

function yfit = compositeV_magnitude_unit(ip,x)
    p   = ip(2);
    wL  = abs(ip(3));
    wG  = abs(ip(4));

    sigma = wG / (2*sqrt(2*log(2)));
    gamma = wL / 2;

    if sigma == 0 && gamma == 0
        yfit = zeros(size(x));
    elseif sigma == 0
        yfit = gamma ./ sqrt(gamma^2 + (x-p).^2);
    else
        z = ((x - p) + 1i*gamma) / (sigma*sqrt(2));
        Vc = fadf(z);
        z0 = (1i*gamma) / (sigma*sqrt(2));
        V0 = fadf(z0);
        N0 = abs(V0);
        if ~isfinite(N0) || N0 == 0
            [~, idx0] = min(abs(x - p));
            N0 = abs(Vc(idx0));
            if ~isfinite(N0) || N0 == 0
                N0 = 1;
            end
        end
        yfit = abs(Vc / N0);
    end
end

%% spline basis
function [Bout, knotsout] = makeSplineBasis(xin, knotSpacing, normalizeFlag)
    xin = xin(:);
    xmin = min(xin);
    xmax = max(xin);
    if knotSpacing <= 0
        error('baseline knotSpacing must be > 0')
    end

    nIntervals = max(2, round((xmax - xmin)/knotSpacing));
    breaks = linspace(xmin, xmax, nIntervals + 1);
    knotsout = augknt(breaks, 4); % cubic B-splines
    nBasis = length(knotsout) - 4;

    Bout = zeros(length(xin), nBasis);
    for jj = 1:nBasis
        coef = zeros(1, nBasis);
        coef(jj) = 1;
        sp = spmak(knotsout, coef);
        Bout(:,jj) = fnval(sp, xin).';
    end

    if normalizeFlag
        scl = max(abs(Bout), [], 1);
        scl(scl==0) = 1;
        Bout = Bout ./ scl;
    end
end

%% baseline options
function opt = parseBaselineOpt(x, opt)
    if ~isfield(opt,'enable') || isempty(opt.enable)
        opt.enable = true;
    end
    if ~isfield(opt,'style') || isempty(opt.style)
        opt.style = 'osprey';
    end
    if ~isfield(opt,'knotSpacing') || isempty(opt.knotSpacing)
        opt.knotSpacing = abs(x(end)-x(1))/10;
    end
    if ~isfield(opt,'lambda') || isempty(opt.lambda)
        switch lower(opt.style)
            case 'osprey'
                opt.lambda = 0;
            case 'lcmodel'
                opt.lambda = 1e-2;
            otherwise
                error('Unknown baseline style')
        end
    end
    if ~isfield(opt,'normalize') || isempty(opt.normalize)
        opt.normalize = true;
    end
    if ~isfield(opt,'amplUB')
        opt.amplUB = [];
    end
    if ~isfield(opt,'linearTol') || isempty(opt.linearTol)
        opt.linearTol = 1e-10;
    end
    % --- HSVD-derived soft priors (both default OFF) --------------------
    % .priorBaseline [npts x 1] a baseline estimate on the same x grid, in the
    %                same phase frame as y (e.g. the unassigned-component sum
    %                from an HSVD decomposition).  Pulled toward, not imposed.
    % .lambdaPrior   weight for it.  The penalty is ||w.*(B*beta - prior)||^2,
    %                weighted exactly like the data term, so lambdaPrior = 1
    %                means "trust the prior as much as the data", 0.1 a gentle
    %                nudge, 0 off.
    % .priorWidth    [n x 1] prior TOTAL linewidth per peak, x-units.  NaN or
    %                <=0 entries are skipped (e.g. peaks HSVD did not find).
    % .lambdaWidth   weight for it.  Normalised so lambdaWidth = 1 means a
    %                100% relative width error costs the same total
    %                sum-of-squares as a 100% relative misfit on the data.
    %                Only modes 5 and 7 (Voigt) are supported; other modes
    %                silently apply no width prior.
    if ~isfield(opt,'lambdaPrior') || isempty(opt.lambdaPrior)
        opt.lambdaPrior = 0;
    end
    if ~isfield(opt,'priorBaseline'), opt.priorBaseline = []; end
    if ~isfield(opt,'lambdaWidth') || isempty(opt.lambdaWidth)
        opt.lambdaWidth = 0;
    end
    if ~isfield(opt,'priorWidth'), opt.priorWidth = []; end
    if ~isfield(opt,'lambdaAmpl') || isempty(opt.lambdaAmpl)
        opt.lambdaAmpl = 0;
    end
    if ~isfield(opt,'verbose') || isempty(opt.verbose)
        opt.verbose = false;
    end
    if ~isfield(opt,'weights') || isempty(opt.weights)
        opt.weights = [];  % no weighting (uniform)
    end
end

%% extra lineshape options
function opt = parseLineShapeOpt(opt)
    if ~isfield(opt,'enable') || isempty(opt.enable)
        opt.enable = true;
    end
    if ~isfield(opt,'nSide') || isempty(opt.nSide)
        opt.nSide = 4;
    end
    if ~isfield(opt,'asymmetric') || isempty(opt.asymmetric)
        opt.asymmetric = true;
    end
    if ~isfield(opt,'maxSide') || isempty(opt.maxSide)
        opt.maxSide = 0.50;
    end
    if ~isfield(opt,'normalize') || isempty(opt.normalize)
        opt.normalize = true;
    end

    if opt.nSide < 0 || round(opt.nSide) ~= opt.nSide
        error('lineShapeOpt.nSide must be a nonnegative integer')
    end
    if opt.maxSide < 0
        error('lineShapeOpt.maxSide must be >= 0')
    end
end

function [pars0, lb0, ub0, names] = initLineShapePars(opt)
    if ~opt.enable || opt.nSide == 0
        pars0 = [];
        lb0 = [];
        ub0 = [];
        names = {};
        return
    end

    if opt.asymmetric
        nls = 2 * opt.nSide;
        pars0 = zeros(1, nls);
        lb0 = zeros(1, nls);
        ub0 = opt.maxSide * ones(1, nls);
        names = {'lineShapeLeft','lineShapeRight'};
    else
        nls = opt.nSide;
        pars0 = zeros(1, nls);
        lb0 = zeros(1, nls);
        ub0 = opt.maxSide * ones(1, nls);
        names = {'lineShape'};
    end

    % Optional warm-start of the side taps.  Starting the kernel at the
    % all-zero (lb) corner leaves the trust-region optimizer with a nearly
    % flat gradient there, so the kernel never engages -- the fit collapses
    % to the bare (un-convolved) peak.  Seeding the taps with a non-zero
    % pedestal profile lets the optimizer start with an engaged kernel and
    % refine it (analogous to baselineOpt.phaseInit for the phase term).
    % opt.init may be a scalar (same value for every tap), a vector of
    % length nls, or a function handle @(k) evaluated at tap offsets
    % k = 1..nSide (broadcast across both sides when asymmetric).  Values
    % are clamped to [lb, ub].  Default (absent/empty) preserves the
    % original all-zero start, so other fits are unaffected.
    if isfield(opt,'init') && ~isempty(opt.init)
        if isa(opt.init,'function_handle')
            p = opt.init((1:opt.nSide));
            if opt.asymmetric, p = [p(:).' p(:).']; end
        elseif isscalar(opt.init)
            p = opt.init * ones(1, nls);
        else
            p = opt.init(:).';
            if numel(p) == opt.nSide && opt.asymmetric, p = [p p]; end
            assert(numel(p) == nls, 'lineShapeOpt.init must be scalar, length nSide, or length %d', nls);
        end
        pars0 = min(max(p, lb0), ub0);
    end
end

function flag = getDisplayFlag(verboseFlag)
    if verboseFlag
        flag = 'iter';
    else
        flag = 'off';
    end
end
