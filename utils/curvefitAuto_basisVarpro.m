function output = curvefitAuto_basisVarpro(x, y, basisFIDs, basisInfo, timeInfo, fitOpt, baselineOpt, lineShapeOpt, parametricOpt)
% curvefitAuto_basisVarpro - LCModel/Osprey-style basis-set fitting with
% variable-projection spline baseline for 31P (or other) MR spectra.
%
% Fits a measured MR spectrum as a linear combination of simulated (or
% measured) basis FIDs. Each basis component is perturbed by an independent
% set of nonlinear parameters (frequency shift, Lorentzian broadening,
% Gaussian broadening, zero-order phase) applied in the time domain, plus
% an optional shared first-order phase term applied in the frequency
% domain. Parametric peaks (Lorentzian, specified by position/width/phase)
% can be supplied alongside the basis FIDs and are treated uniformly in the
% same per-component parameter framework. Soft amplitude-ratio constraints
% between basis components are supported via Tikhonov penalty rows.
%
% Inputs:
%   x          - full ppm axis [N x 1] aligned with fftshift(fft(fid))
%   y          - full measured spectrum [N x 1], complex preferred
%   basisFIDs  - simulated basis FIDs [Nb x nMet], complex, same dwell time
%                as the data (pass [] if only parametric peaks are used)
%   basisInfo  - struct array [nMet x 1] describing each basis:
%                  .name        - metabolite label (char)
%                  .refPpm      - nominal chemical shift (ppm, for reports)
%                  .group       - optional integer group index (for future
%                                 grouped-parameter extensions; not used)
%   timeInfo   - struct with acquisition metadata:
%                  .dwellTime   - dwell time in seconds (required)
%                  .f0          - Larmor frequency in MHz (required)
%                  .t           - optional time axis [Nb x 1]
%                  .ppmPivot    - pivot ppm for 1st-order phase (default 0)
%   fitOpt     - struct of fit options (all fields optional):
%                  .fitRange      - [min max] ppm window to fit
%                  .complex       - true/false; default: ~isreal(y)
%                  .shiftInit     - initial per-component freq shift (Hz),
%                                   scalar or [nMet x 1] (default 0)
%                  .shiftBounds   - [min max] Hz, scalar or [nMet x 2]
%                                   (default [-20 20])
%                  .lbLInit       - initial Lorentzian LB (Hz) (default 2)
%                  .lbLBounds     - [min max] (default [0 40])
%                  .lbGInit       - initial Gaussian LB (Hz)  (default 0)
%                  .lbGBounds     - [min max] (default [0 30])
%                  .phaseInit     - initial phi0 per component (deg)
%                  .phaseBounds   - [min max] deg (default [-30 30])
%                  .phase1Enable  - include global phi1 (default false)
%                  .phase1Init    - initial phi1 (deg/ppm) (default 0)
%                  .phase1Bounds  - [min max] deg/ppm (default [-45 45])
%                  .softConstraints - struct array of amplitude-ratio
%                                   soft constraints:
%                                     .idxA    - index or name of comp A
%                                     .idxB    - index or name of comp B
%                                     .ratio   - target amp(A)/amp(B)
%                                     .lambda  - penalty weight (>=0)
%                  .tieGroups     - struct array of parameter-tying groups.
%                                   Each entry shares one or more nonlinear
%                                   parameters across listed members:
%                                     .members - cell array of basis names
%                                                or indices; first listed
%                                                is the LEADER whose init
%                                                value and bounds apply.
%                                     .params  - cell of which params to
%                                                tie: any of {'shift',
%                                                'lbL','lbG','phi0'}.
%                                   Only the leader's parameter is exposed
%                                   to lsqnonlin; followers are aliased to
%                                   it. Untied parameters remain free and
%                                   independent.
%                  .maxIter, .tolFun, .tolX, .verbose
%   baselineOpt - spline baseline options (same as curvefitAuto_varproBaseline):
%                  .enable, .knotSpacing, .lambda, .normalize, .amplUB,
%                  .linearTol, .lambdaAmpl, .verbose, .weights
%   parametricOpt - optional parametric add-on peaks (treated like basis
%                   FIDs generated analytically):
%                  .enable     - true/false (default false)
%                  .names      - cell of strings
%                  .peakPpm    - peak positions (ppm)
%                  .widthHz    - initial Lorentzian FWHM (Hz)
%                  .phaseDeg   - initial phase (deg, default 0)
%                  .shiftBounds, .lbLBounds, .phaseBounds - per-peak bound
%                                overrides (rows appended in same order as
%                                .peakPpm)
%   lineShapeOpt - optional shared spectrum-domain lineshape kernel applied
%                  to peak basis spectra only (not the spline baseline),
%                  similar in spirit to LCModel/Osprey extra lineshape
%                  distortion. Allows asymmetric peak shapes that the
%                  per-component Lorentzian/Gaussian/phi0 model cannot
%                  capture. The kernel acts in pixel space.
%                  .enable      = true/false (default false)
%                  .nSide       = number of side points on each side of
%                                 the unit center tap (default 2)
%                  .asymmetric  = if true, left and right sides have
%                                 independent coefficients (default false)
%                  .maxSide     = upper bound for each side coefficient
%                                 (default 0.5)
%                  .normalize   = normalize kernel to sum to 1 (default true)
%
% Output (fields):
%   .ampl            - fitted amplitudes [nMet x 1]
%   .names           - metabolite names [nMet x 1 cell]
%   .fit             - total fitted spectrum on fit range
%   .fit_peaks       - peak-only fitted spectrum on fit range
%   .baseline        - spline baseline on fit range
%   .comp            - per-component fitted spectra on fit range [nFit x nMet]
%   .beta            - spline coefficients
%   .B               - spline basis matrix
%   .basisSpec       - modified (shift/LB/phase applied) basis spectra
%                      on fit range (before amplitude scaling)
%   .basisFIDmod     - modified basis FIDs (time domain, full length)
%   .pars            - nonlinear parameters, struct with fields:
%                        .shift, .lbL, .lbG, .phase0, .phase1
%   .lb, .ub         - nonlinear parameter bounds (struct, same fields)
%   .parnames        - cell array labeling nonlinear pars
%   .residualFit     - yfit - y on fit range
%   .softResiduals   - residuals of the soft-constraint rows
%   .fitIdx          - logical mask of points used in the fit
%   .lineShapePars   - fitted side coefficients of the shared kernel
%   .lineShapeKernel - shared kernel actually used (vector, length 2*nSide+1)
%   .fitOpt, .baselineOpt, .parametricOpt, .lineShapeOpt, .timeInfo - echoes of inputs
%
% NW - 04/2026

% ------------------------------------------------------------------- %
% input checks and defaults
% ------------------------------------------------------------------- %
x = x(:);
y = y(:);
N = length(x);
if length(y) ~= N
    error('x and y length mismatch');
end

if nargin < 3, basisFIDs = []; end
if nargin < 4, basisInfo = struct([]); end
if nargin < 5 || isempty(timeInfo)
    error('timeInfo (with .dwellTime and .f0) is required');
end
if ~isfield(timeInfo,'dwellTime') || isempty(timeInfo.dwellTime)
    error('timeInfo.dwellTime is required');
end
if ~isfield(timeInfo,'f0') || isempty(timeInfo.f0)
    error('timeInfo.f0 (MHz) is required');
end
if ~isfield(timeInfo,'ppmPivot') || isempty(timeInfo.ppmPivot)
    timeInfo.ppmPivot = 0;
end

if nargin < 6 || isempty(fitOpt),        fitOpt = struct();        end
if nargin < 7 || isempty(baselineOpt),   baselineOpt = struct();   end
if nargin < 9 || isempty(parametricOpt), parametricOpt = struct(); end
if nargin < 8 || isempty(lineShapeOpt),  lineShapeOpt = struct();  end
lineShapeOpt = parseLineShapeOpt(lineShapeOpt);

dwell = timeInfo.dwellTime;
f0    = timeInfo.f0;  % MHz

% ------------------------------------------------------------------- %
% build parametric peak FIDs (if any) and append to basis
% ------------------------------------------------------------------- %
parametricOpt = parseParametricOpt(parametricOpt);

% total time-axis length for the basis = length of basisFIDs if given,
% otherwise match data (N). If both, enforce equal length.
if ~isempty(basisFIDs)
    Nt = size(basisFIDs,1);
else
    Nt = N;
end
if isfield(timeInfo,'t') && ~isempty(timeInfo.t)
    t = timeInfo.t(:);
    if length(t) ~= Nt
        error('timeInfo.t length does not match basisFIDs');
    end
else
    t = (0:Nt-1).' * dwell;
end

nBasis = size(basisFIDs,2);
if nBasis == 0
    nameBasis = {};
else
    nameBasis = cell(nBasis,1);
    for k = 1:nBasis
        if k <= length(basisInfo) && isfield(basisInfo,'name') && ~isempty(basisInfo(k).name)
            nameBasis{k} = basisInfo(k).name;
        else
            nameBasis{k} = sprintf('basis%02d', k);
        end
    end
end

% Build parametric peak FIDs.  Each is an analytic complex Lorentzian at
% the requested ppm position with width widthHz and phase phaseDeg.
paramFIDs = [];
nameParam = {};
if parametricOpt.enable && ~isempty(parametricOpt.peakPpm)
    np = length(parametricOpt.peakPpm);
    paramFIDs = zeros(Nt, np);
    nameParam = cell(np,1);

    freqAxisFull = buildFreqAxis(N, dwell);
    for kp = 1:np
        peakHz = ppmToHzByAxis(parametricOpt.peakPpm(kp), x, freqAxisFull);
        w = parametricOpt.widthHz(kp);
        ph = parametricOpt.phaseDeg(kp) * pi/180;
        paramFIDs(:,kp) = exp(1i*2*pi*peakHz*t) .* exp(-pi*w*t) .* exp(-1i*ph);
        if kp <= length(parametricOpt.names) && ~isempty(parametricOpt.names{kp})
            nameParam{kp} = parametricOpt.names{kp};
        else
            nameParam{kp} = sprintf('param%02d', kp);
        end
    end
end

% concatenate: basis first, parametric after
if isempty(basisFIDs) && isempty(paramFIDs)
    error('At least one basis FID or one parametric peak is required');
end
allFIDs = [basisFIDs, paramFIDs];
names   = [nameBasis; nameParam];
nMet    = size(allFIDs, 2);

% ------------------------------------------------------------------- %
% fit options and bounds
% ------------------------------------------------------------------- %
nParam = size(paramFIDs,2);
fitOpt = parseFitOpt(fitOpt, nMet, nBasis, nParam, parametricOpt, y);

% restrict data to fit range
if isempty(fitOpt.fitRange)
    fitIdx = true(N,1);
else
    fitIdx = (x >= min(fitOpt.fitRange)) & (x <= max(fitOpt.fitRange));
end
xFit = x(fitIdx);
yFit = y(fitIdx);
nFit = length(xFit);
if nFit < 10
    error('fitOpt.fitRange produced <10 points; check axis units');
end

complexFit = fitOpt.complex;
if complexFit && isreal(yFit)
    warning('fitOpt.complex=true but y is purely real; fitting real part only');
    complexFit = false;
end

% ------------------------------------------------------------------- %
% baseline options
% ------------------------------------------------------------------- %
baselineOpt = parseBaselineOpt(xFit, baselineOpt);

if baselineOpt.enable
    [B, knots] = makeSplineBasis(xFit, baselineOpt.knotSpacing, baselineOpt.normalize);
else
    B = zeros(nFit, 0);
    knots = [];
end
nb = size(B,2);

% weights (applied to residual)
if ~isempty(baselineOpt.weights)
    w = baselineOpt.weights(:);
    if length(w) ~= nFit
        error('baselineOpt.weights must match number of fit points (%d)', nFit);
    end
else
    w = ones(nFit, 1);
end

% amplitude upper bounds
if isempty(baselineOpt.amplUB)
    amplUB = inf(nMet, 1);
elseif isscalar(baselineOpt.amplUB)
    amplUB = baselineOpt.amplUB * ones(nMet, 1);
else
    amplUB = baselineOpt.amplUB(:);
    if length(amplUB) ~= nMet
        error('baselineOpt.amplUB must be scalar or length nMet (%d)', nMet);
    end
end

% ------------------------------------------------------------------- %
% soft amplitude-ratio constraints
% ------------------------------------------------------------------- %
[Csoft, cRhs, cLam] = buildSoftConstraints(fitOpt.softConstraints, names, nMet);

% ------------------------------------------------------------------- %
% nonlinear parameter layout
%   ip_full layout (per-component, length 4*nMet [+1] if phi1):
%     ip_full(1:nMet)              shift (Hz)
%     ip_full(nMet+(1:nMet))       lbL   (Hz)
%     ip_full(2*nMet+(1:nMet))     lbG   (Hz)
%     ip_full(3*nMet+(1:nMet))     phi0  (deg)
%     ip_full(4*nMet+1)  (opt)     phi1  (deg/ppm)
%
%   With tieGroups, lsqnonlin operates on a REDUCED vector ip_red where
%   tied followers collapse onto their leader's slot.  Expansion is
%   ip_full = ip_red(mapToFull).
% ------------------------------------------------------------------- %
parnames = {'shift(Hz)','lbL(Hz)','lbG(Hz)','phi0(deg)'};

% per-component aliasing: each component k by default maps to its own
% slot k.  tieGroups overwrite alias(member) <- leader for the relevant
% parameter blocks, collapsing followers onto the leader.
aliasShift = (1:nMet).';
aliasLbL   = (1:nMet).';
aliasLbG   = (1:nMet).';
aliasPhi0  = (1:nMet).';

if isfield(fitOpt,'tieGroups') && ~isempty(fitOpt.tieGroups)
    tg = fitOpt.tieGroups;
    if ~isstruct(tg)
        error('fitOpt.tieGroups must be a struct array');
    end
    for gi = 1:length(tg)
        members = tg(gi).members;
        if ischar(members) || isstring(members), members = {char(members)}; end
        if isempty(members) || length(members) < 2, continue; end
        idxs = zeros(length(members),1);
        for mi = 1:length(members)
            idxs(mi) = resolveIdx(members{mi}, names);
        end
        leader = idxs(1);
        params = tg(gi).params;
        if ischar(params) || isstring(params), params = {char(params)}; end
        for pp = 1:length(params)
            switch lower(char(params{pp}))
                case 'shift'
                    aliasShift(idxs) = leader;
                case {'lbl','lb_l'}
                    aliasLbL(idxs) = leader;
                case {'lbg','lb_g'}
                    aliasLbG(idxs) = leader;
                case {'phi0','phase0','phase'}
                    aliasPhi0(idxs) = leader;
                otherwise
                    error('Unknown tieGroups param "%s"', char(params{pp}));
            end
        end
    end
end

% Reduce: for each block, find the unique leaders (in order of first
% appearance) and the position each component reads from.
[uniqShift, ~, posShift] = unique(aliasShift, 'stable');
[uniqLbL,   ~, posLbL  ] = unique(aliasLbL,   'stable');
[uniqLbG,   ~, posLbG  ] = unique(aliasLbG,   'stable');
[uniqPhi0,  ~, posPhi0 ] = unique(aliasPhi0,  'stable');

nFreeShift = length(uniqShift);
nFreeLbL   = length(uniqLbL);
nFreeLbG   = length(uniqLbG);
nFreePhi0  = length(uniqPhi0);

% Map full-vector index -> reduced-vector index.
mapToFull = [posShift(:).' ...
             nFreeShift + posLbL(:).' ...
             nFreeShift + nFreeLbL + posLbG(:).' ...
             nFreeShift + nFreeLbL + nFreeLbG + posPhi0(:).'];

pars0_red = [fitOpt.shiftInit(uniqShift).'    fitOpt.lbLInit(uniqLbL).'    fitOpt.lbGInit(uniqLbG).'    fitOpt.phaseInit(uniqPhi0).'];
lb_red    = [fitOpt.shiftBounds(uniqShift,1).' fitOpt.lbLBounds(uniqLbL,1).' fitOpt.lbGBounds(uniqLbG,1).' fitOpt.phaseBounds(uniqPhi0,1).'];
ub_red    = [fitOpt.shiftBounds(uniqShift,2).' fitOpt.lbLBounds(uniqLbL,2).' fitOpt.lbGBounds(uniqLbG,2).' fitOpt.phaseBounds(uniqPhi0,2).'];

if fitOpt.phase1Enable
    pars0_red = [pars0_red fitOpt.phase1Init];
    lb_red    = [lb_red    fitOpt.phase1Bounds(1)];
    ub_red    = [ub_red    fitOpt.phase1Bounds(2)];
    mapToFull = [mapToFull length(pars0_red)];
    parnames{end+1} = 'phi1(deg/ppm)';
end

% Append optional shared lineshape distortion parameters.  These are
% global (not per-component) and bypass the tieGroups machinery, so we
% extend mapToFull with an identity tail covering the kernel side coefs.
nFreePeak = length(pars0_red);
[lineShapePars0, lineShapeLb, lineShapeUb, lineShapeParnames] = initLineShapePars(lineShapeOpt);
nLineShape = length(lineShapePars0);
if nLineShape > 0
    pars0_red = [pars0_red lineShapePars0];
    lb_red    = [lb_red    lineShapeLb];
    ub_red    = [ub_red    lineShapeUb];
    mapToFull = [mapToFull nFreePeak + (1:nLineShape)];
    parnames  = [parnames  lineShapeParnames];
end

% Full-length echoes for output (each tied follower repeats the leader).
lb = lb_red(mapToFull);
ub = ub_red(mapToFull);

nPars = length(pars0_red);
if any(pars0_red < lb_red) || any(pars0_red > ub_red)
    warning('Some initial parameters are outside their bounds; clipping.');
    pars0_red = min(max(pars0_red, lb_red), ub_red);
end

% ------------------------------------------------------------------- %
% set up lsqnonlin
% ------------------------------------------------------------------- %
ssrData = sum(abs(w .* yFit).^2);
relTolFun = 1e-6;
relTolX   = 1e-9;
finDiffStep = 1e-8;
if isfield(fitOpt,'tolFun') && ~isempty(fitOpt.tolFun), relTolFun = fitOpt.tolFun; end
if isfield(fitOpt,'tolX')   && ~isempty(fitOpt.tolX),   relTolX   = fitOpt.tolX;   end
if isfield(fitOpt,'finDiffStep') && ~isempty(fitOpt.finDiffStep), finDiffStep = fitOpt.finDiffStep; end
tolFun = relTolFun * ssrData;
tolX   = relTolX;

maxIter = 4000;
if isfield(fitOpt,'maxIter') && ~isempty(fitOpt.maxIter), maxIter = fitOpt.maxIter; end

oldopts = optimset('lsqnonlin');
% options = optimset(oldopts, 'TolFun', tolFun, 'TolX', tolX, ...
%     'FinDiffRelStep', finDiffStep, ...
%     'MaxFunEval', 2000*nPars, 'MaxIter', maxIter, ...
%     'Display', getDisplayFlag(fitOpt.verbose));
options = optimset(oldopts, 'TolFun', tolFun, 'TolX', tolX, ...
    'MaxFunEval', 2000*nPars, 'MaxIter', maxIter, ...
    'Display', getDisplayFlag(fitOpt.verbose));

% ------------------------------------------------------------------- %
% run optimizer (operates on reduced free-parameter vector; residual
% expands to full per-component vector via mapToFull)
% ------------------------------------------------------------------- %
ip_red    = lsqnonlin(@residual, pars0_red, lb_red, ub_red, options);
ip_nonlin = ip_red(mapToFull);

% ------------------------------------------------------------------- %
% final evaluation
% ------------------------------------------------------------------- %
[Apeak, fidMod] = buildApeak(ip_nonlin);
[ampl, beta, baseline, yfit, y_comp, softRes] = solveLinearBlock(Apeak, B, yFit);

[shiftOut, lbLout, lbGout, phi0out, phi1out] = unpackPars(ip_nonlin);

% ------------------------------------------------------------------- %
% output
% ------------------------------------------------------------------- %
output.ampl           = ampl;
output.names          = names;
output.fit            = yfit;
output.fit_peaks      = yfit - baseline;
output.baseline       = baseline;
output.comp           = y_comp;
output.beta           = beta;
output.B              = B;
output.knots          = knots;
output.basisSpec      = Apeak;
output.basisFIDmod    = fidMod;
output.pars.shift     = shiftOut;
output.pars.lbL       = lbLout;
output.pars.lbG       = lbGout;
% Voigt FWHM (Hz) per component via Olivero-Longbothum approximation:
%   fwhmV ~ 0.5346*lbL + sqrt(0.2166*lbL^2 + lbG^2)
% (uses fitted addition only; ignores any baseline broadening baked into
% the basis FID at simulation time)
output.pars.fwhmV     = 0.5346*lbLout(:) + sqrt(0.2166*lbLout(:).^2 + lbGout(:).^2);
% Effective peak widths in the data: include the Lorentzian baked into
% each basis FID (basisInfo.lwHzBasis) so that lbL_total / fwhmV_total
% reflect what's actually observed, not just what the fit added on top.
lwBase = zeros(nMet,1);
if ~isempty(basisInfo) && isfield(basisInfo,'lwHzBasis')
    for kBI = 1:min(numel(basisInfo), nMet)
        if ~isempty(basisInfo(kBI).lwHzBasis)
            lwBase(kBI) = basisInfo(kBI).lwHzBasis;
        end
    end
end
output.pars.lwHzBasis  = lwBase;
output.pars.lbL_total  = lwBase + lbLout(:);
output.pars.fwhmV_total = 0.5346*output.pars.lbL_total + ...
                         sqrt(0.2166*output.pars.lbL_total.^2 + lbGout(:).^2);
output.pars.phase0    = phi0out;
output.pars.phase1    = phi1out;
output.lb.shift       = lb(1:nMet);
output.lb.lbL         = lb(nMet+(1:nMet));
output.lb.lbG         = lb(2*nMet+(1:nMet));
output.lb.phase0      = lb(3*nMet+(1:nMet));
output.ub.shift       = ub(1:nMet);
output.ub.lbL         = ub(nMet+(1:nMet));
output.ub.lbG         = ub(2*nMet+(1:nMet));
output.ub.phase0      = ub(3*nMet+(1:nMet));
if fitOpt.phase1Enable
    output.lb.phase1 = lb(4*nMet + 1);
    output.ub.phase1 = ub(4*nMet + 1);
end
% Lineshape side coefs (live at the tail of ip_full / lb / ub).
if nLineShape > 0
    output.lineShapePars   = ip_nonlin(end - nLineShape + 1 : end);
    output.lineShapeKernel = buildLineShapeKernelLocal(output.lineShapePars, lineShapeOpt);
    output.lb.lineShape    = lb(end - nLineShape + 1 : end);
    output.ub.lineShape    = ub(end - nLineShape + 1 : end);
else
    output.lineShapePars   = [];
    output.lineShapeKernel = 1;
end
output.lineShapeParnames = lineShapeParnames;
output.parnames       = parnames;
output.residualFit    = yfit - yFit;
output.softResiduals  = softRes;
output.fitIdx         = fitIdx;
output.xFit           = xFit;
output.yFit           = yFit;
output.fitOpt         = fitOpt;
output.baselineOpt    = baselineOpt;
output.parametricOpt  = parametricOpt;
output.lineShapeOpt   = lineShapeOpt;
output.timeInfo       = timeInfo;
output.nBasis         = nBasis;
output.nParam         = nParam;
output.basisInfo      = basisInfo;

% =================================================================== %
% nested functions (share parent workspace: allFIDs, t, x, xFit, yFit,
% fitIdx, nMet, B, w, amplUB, baselineOpt, fitOpt, Csoft, cRhs, cLam,
% complexFit, N, f0, timeInfo)
% =================================================================== %

    function r = residual(ipRed)
        ipFull = ipRed(mapToFull);
        Apk = buildApeak(ipFull);
        [~, ~, ~, yf] = solveLinearBlock(Apk, B, yFit);
        resid = yf - yFit;
        if complexFit
            r = [w .* real(resid); w .* imag(resid)];
        else
            r = w .* real(resid);
        end
    end

    function [Apk, fidMod] = buildApeak(ip)
        [shift, lbL, lbG, phi0, phi1, sidePars] = unpackPars(ip);

        % Time-domain modulation (broadcast)
        % fid_mod(t,k) = fid(t,k) .* exp(1i*2*pi*shift(k)*t)
        %                         .* exp(-pi*lbL(k)*t)
        %                         .* exp(-(pi*lbG(k)*t).^2 / (4*log(2)))
        %                         .* exp(1i*phi0(k)*pi/180)
        shiftRow = shift(:).';
        lbLrow   = lbL(:).';
        lbGrow   = lbG(:).';
        phi0row  = phi0(:).';

        tt = t;   % Nt x 1
        mod_shift = exp(1i * 2*pi * (tt * shiftRow));
        mod_lbL   = exp(-pi * (tt * lbLrow));
        mod_lbG   = exp(-((pi * (tt * lbGrow)).^2) / (4*log(2)));
        mod_phi0  = exp(-1i * (phi0row * pi/180));

        fidMod = allFIDs .* mod_shift .* mod_lbL .* mod_lbG .* mod_phi0;

        % FFT -> spectrum (full length)
        Nt_local = size(fidMod, 1);
        specFull = fftshift(fft(fidMod, [], 1), 1);

        % If basis length differs from data length, interpolate to x grid.
        if Nt_local == N
            specOnData = specFull;
        else
            freqBasis = buildFreqAxis(Nt_local, dwell);
            freqData  = buildFreqAxis(N, dwell);
            specOnData = zeros(N, nMet);
            for kk = 1:nMet
                specOnData(:,kk) = interp1(freqBasis, specFull(:,kk), freqData, 'linear', 0);
            end
        end

        % Apply optional global 1st-order phase in spectrum domain.
        if fitOpt.phase1Enable
            ppmRel = x - timeInfo.ppmPivot;
            phaseRamp = exp(-1i * phi1 * ppmRel * pi/180);
            specOnData = specOnData .* phaseRamp;
        end

        % Apply optional shared lineshape kernel (peak basis only — the
        % spline baseline is unaffected because it lives in B and is
        % solved separately in solveLinearBlock).  Convolution is run on
        % the full N-length spectrum so that 'same' edge effects fall
        % outside the fit window.
        if lineShapeOpt.enable && nLineShape > 0
            kern = buildLineShapeKernel(sidePars);
            if numel(kern) > 1
                for kk = 1:nMet
                    specOnData(:,kk) = conv(specOnData(:,kk), kern, 'same');
                end
            end
        end

        % Restrict to fit range
        specFit = specOnData(fitIdx, :);

        if complexFit
            Apk = specFit;
        else
            Apk = real(specFit);
        end
    end

    function [ampl, beta, baseline, yfit, y_comp, softRes] = solveLinearBlock(Apk, Bmat, data)
        nloc = size(Apk, 2);
        nbL  = size(Bmat, 2);

        if complexFit
            % Stack real and imag rows; amplitudes and spline coefs are real.
            Xtop = [real(Apk), real(Bmat)];
            Xbot = [imag(Apk), imag(Bmat)];
            Xw   = [w .* Xtop; w .* Xbot];
            rhsw = [w .* real(data); w .* imag(data)];
        else
            Xw   = w .* [real(Apk), real(Bmat)];
            rhsw = w .* real(data);
        end

        % Spline roughness penalty (second differences on beta only)
        if nbL > 0 && baselineOpt.lambda > 0
            D = diff(speye(nbL), 2);
            Xw   = [Xw; zeros(size(D,1), nloc) sqrt(baselineOpt.lambda)*D];
            rhsw = [rhsw; zeros(size(D,1), 1)];
        end

        % Amplitude Tikhonov (pulls unused amps to zero)
        if baselineOpt.lambdaAmpl > 0
            Xw   = [Xw; sqrt(baselineOpt.lambdaAmpl)*eye(nloc) zeros(nloc, nbL)];
            rhsw = [rhsw; zeros(nloc, 1)];
        end

        % Soft amplitude-ratio constraints: sqrt(lam) * C * [ampl; beta] = sqrt(lam) * rhs
        if ~isempty(Csoft)
            scaleC = sqrt(cLam(:));
            Crows  = bsxfun(@times, [Csoft zeros(size(Csoft,1), nbL)], scaleC);
            rhsC   = scaleC .* cRhs(:);
            Xw     = [Xw; Crows];
            rhsw   = [rhsw; rhsC];
        end

        lbin = [zeros(nloc,1);  -inf(nbL,1)];
        ubin = [amplUB(:);       inf(nbL,1)];

        opts = optimoptions('lsqlin', 'Display', 'off', ...
            'OptimalityTolerance', baselineOpt.linearTol, ...
            'ConstraintTolerance', baselineOpt.linearTol, ...
            'StepTolerance',       baselineOpt.linearTol);

        % initial guess: project data onto real(Apk) columns, non-negative
        if complexFit
            x0amp = max(0, real(Apk)' * real(data) + imag(Apk)' * imag(data));
        else
            x0amp = max(0, real(Apk)' * real(data));
        end
        x0 = [x0amp; zeros(nbL,1)];
        x0 = min(max(x0, lbin), ubin);

        chat = lsqlin(Xw, rhsw, [], [], [], [], lbin, ubin, x0, opts);

        ampl = chat(1:nloc);
        beta = chat(nloc+(1:nbL));

        if nbL > 0
            baseline = Bmat * beta;   % real-valued
        else
            baseline = zeros(size(data));
        end

        y_comp = zeros(size(Apk));
        for kk = 1:nloc
            y_comp(:,kk) = ampl(kk) * Apk(:,kk);
        end
        yfit = sum(y_comp, 2) + baseline;

        if nargout >= 6
            if ~isempty(Csoft)
                softRes = Csoft * ampl - cRhs(:);
            else
                softRes = [];
            end
        end
    end

    function [shift, lbL, lbG, phi0, phi1, sidePars] = unpackPars(ip)
        shift = ip(1:nMet);
        lbL   = ip(nMet + (1:nMet));
        lbG   = ip(2*nMet + (1:nMet));
        phi0  = ip(3*nMet + (1:nMet));
        if fitOpt.phase1Enable
            phi1 = ip(4*nMet + 1);
            offset = 4*nMet + 1;
        else
            phi1 = 0;
            offset = 4*nMet;
        end
        if length(ip) > offset
            sidePars = ip(offset+1:end);
        else
            sidePars = [];
        end
    end

    function kern = buildLineShapeKernel(sidePars)
        if ~lineShapeOpt.enable || isempty(sidePars)
            kern = 1;
            return
        end
        nSide = lineShapeOpt.nSide;
        sidePars = max(0, sidePars(:));
        if lineShapeOpt.asymmetric
            left  = sidePars(1:nSide);
            right = sidePars(nSide + (1:nSide));
        else
            left  = sidePars(1:nSide);
            right = left;
        end
        kern = [flipud(left); 1; right];
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

end % ---- end of main function ------------------------------------- %


% =================================================================== %
% local (non-nested) helpers
% =================================================================== %

function freq = buildFreqAxis(N, dwell)
    freq = ((-N/2):(N/2-1)).' / (N*dwell);
end

function hz = ppmToHzByAxis(ppmVal, ppmAxis, freqAxis)
    % Find Hz value corresponding to a ppm value by nearest-point lookup on
    % the data's ppm axis.  This avoids hard-coding a sign convention for
    % the ppm<->Hz map.
    [~, idx] = min(abs(ppmAxis - ppmVal));
    hz = freqAxis(idx);
end

function [C, rhs, lam] = buildSoftConstraints(scList, names, nMet)
    C = [];
    rhs = [];
    lam = [];
    if isempty(scList) || ~isstruct(scList)
        return
    end
    nC = length(scList);
    C   = zeros(nC, nMet);
    rhs = zeros(nC, 1);
    lam = zeros(nC, 1);
    hasLambda = isfield(scList, 'lambda');
    for k = 1:nC
        iA = resolveIdx(scList(k).idxA, names);
        iB = resolveIdx(scList(k).idxB, names);
        r  = scList(k).ratio;
        if hasLambda && ~isempty(scList(k).lambda)
            l = scList(k).lambda;
        else
            l = 1;
        end
        % constraint: amp(A) - ratio*amp(B) = 0
        row = zeros(1, nMet);
        row(iA) = row(iA) + 1;
        row(iB) = row(iB) - r;
        C(k,:)  = row;
        rhs(k)  = 0;
        lam(k)  = l;
    end
end

function idx = resolveIdx(v, names)
    if ischar(v) || isstring(v)
        idx = find(strcmpi(names, char(v)), 1);
        if isempty(idx)
            error('Soft-constraint name "%s" not found in basis names', char(v));
        end
    else
        idx = v;
    end
end

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
            case 'osprey',  opt.lambda = 0;
            case 'lcmodel', opt.lambda = 1e-2;
            otherwise,      error('Unknown baseline style');
        end
    end
    if ~isfield(opt,'normalize') || isempty(opt.normalize)
        opt.normalize = true;
    end
    if ~isfield(opt,'amplUB'),        opt.amplUB = [];      end
    if ~isfield(opt,'linearTol') || isempty(opt.linearTol)
        opt.linearTol = 1e-10;
    end
    if ~isfield(opt,'lambdaAmpl') || isempty(opt.lambdaAmpl)
        opt.lambdaAmpl = 0;
    end
    if ~isfield(opt,'verbose') || isempty(opt.verbose)
        opt.verbose = false;
    end
    if ~isfield(opt,'weights') || isempty(opt.weights)
        opt.weights = [];
    end
end

function opt = parseParametricOpt(opt)
    if ~isfield(opt,'enable'),     opt.enable = false; end
    if ~isfield(opt,'peakPpm'),    opt.peakPpm = [];    end
    if ~isfield(opt,'widthHz'),    opt.widthHz = [];    end
    if ~isfield(opt,'phaseDeg'),   opt.phaseDeg = [];   end
    if ~isfield(opt,'names'),      opt.names = {};      end
    if ~isempty(opt.peakPpm) && isempty(opt.widthHz)
        error('parametricOpt.widthHz required when peakPpm given');
    end
    if isempty(opt.phaseDeg) && ~isempty(opt.peakPpm)
        opt.phaseDeg = zeros(size(opt.peakPpm));
    end
    opt.peakPpm  = opt.peakPpm(:);
    opt.widthHz  = opt.widthHz(:);
    opt.phaseDeg = opt.phaseDeg(:);
    if ~isempty(opt.peakPpm) && ~opt.enable
        opt.enable = true;
    end
end

function opt = parseFitOpt(opt, nMet, nBasis, nParam, parametricOpt, y)
    if ~isfield(opt,'fitRange'),  opt.fitRange = []; end

    if ~isfield(opt,'complex') || isempty(opt.complex)
        opt.complex = ~isreal(y);
    end

    % per-component defaults (scalar expanded to nMet; [min max] expanded
    % to [nMet x 2])
    opt.shiftInit    = expandScalarVec(getOrDefault(opt,'shiftInit',0),           nMet);
    opt.shiftBounds  = expandBounds(getOrDefault(opt,'shiftBounds',[-20 20]),     nMet);
    opt.lbLInit      = expandScalarVec(getOrDefault(opt,'lbLInit',2),             nMet);
    opt.lbLBounds    = expandBounds(getOrDefault(opt,'lbLBounds',[0 40]),         nMet);
    opt.lbGInit      = expandScalarVec(getOrDefault(opt,'lbGInit',0),             nMet);
    opt.lbGBounds    = expandBounds(getOrDefault(opt,'lbGBounds',[0 30]),         nMet);
    opt.phaseInit    = expandScalarVec(getOrDefault(opt,'phaseInit',0),           nMet);
    opt.phaseBounds  = expandBounds(getOrDefault(opt,'phaseBounds',[-30 30]),     nMet);

    if ~isfield(opt,'phase1Enable') || isempty(opt.phase1Enable)
        opt.phase1Enable = false;
    end
    if ~isfield(opt,'phase1Init') || isempty(opt.phase1Init)
        opt.phase1Init = 0;
    end
    if ~isfield(opt,'phase1Bounds') || isempty(opt.phase1Bounds)
        opt.phase1Bounds = [-45 45];
    end

    % per-peak overrides for parametric add-ons (last nParam rows)
    if nParam > 0
        if isfield(parametricOpt,'shiftBounds') && ~isempty(parametricOpt.shiftBounds)
            opt.shiftBounds(nBasis+(1:nParam),:) = expandBounds(parametricOpt.shiftBounds, nParam);
        end
        if isfield(parametricOpt,'lbLBounds') && ~isempty(parametricOpt.lbLBounds)
            opt.lbLBounds(nBasis+(1:nParam),:)   = expandBounds(parametricOpt.lbLBounds,   nParam);
        end
        if isfield(parametricOpt,'phaseBounds') && ~isempty(parametricOpt.phaseBounds)
            opt.phaseBounds(nBasis+(1:nParam),:) = expandBounds(parametricOpt.phaseBounds, nParam);
        end
    end

    if ~isfield(opt,'softConstraints'), opt.softConstraints = []; end
    if ~isfield(opt,'verbose') || isempty(opt.verbose), opt.verbose = false; end
end

function v = getOrDefault(s, f, d)
    if isfield(s, f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = d;
    end
end

function out = expandScalarVec(v, n)
    v = v(:);
    if isscalar(v)
        out = v * ones(n,1);
    elseif length(v) == n
        out = v;
    else
        error('Vector length %d does not match nMet=%d', length(v), n);
    end
end

function out = expandBounds(b, n)
    % Accept [1 x 2] scalar bounds or [n x 2] per-component bounds.
    if isvector(b) && numel(b) == 2
        out = repmat(b(:).', n, 1);
    elseif size(b,1) == n && size(b,2) == 2
        out = b;
    else
        error('Bounds must be [1 x 2] or [%d x 2]', n);
    end
end

function [Bout, knotsout] = makeSplineBasis(xin, knotSpacing, normalizeFlag)
    xin = xin(:);
    xmin = min(xin);
    xmax = max(xin);
    if knotSpacing <= 0
        error('baseline knotSpacing must be > 0');
    end
    nIntervals = max(2, round((xmax - xmin)/knotSpacing));
    breaks = linspace(xmin, xmax, nIntervals + 1);
    knotsout = augknt(breaks, 4);
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

function flag = getDisplayFlag(verboseFlag)
    if verboseFlag
        flag = 'iter';
    else
        flag = 'off';
    end
end

function opt = parseLineShapeOpt(opt)
    if ~isfield(opt,'enable') || isempty(opt.enable)
        opt.enable = false;
    end
    if ~isfield(opt,'nSide') || isempty(opt.nSide)
        opt.nSide = 2;
    end
    if ~isfield(opt,'asymmetric') || isempty(opt.asymmetric)
        opt.asymmetric = false;
    end
    if ~isfield(opt,'maxSide') || isempty(opt.maxSide)
        opt.maxSide = 0.50;
    end
    if ~isfield(opt,'normalize') || isempty(opt.normalize)
        opt.normalize = true;
    end
    if opt.nSide < 0 || round(opt.nSide) ~= opt.nSide
        error('lineShapeOpt.nSide must be a nonnegative integer');
    end
    if opt.maxSide < 0
        error('lineShapeOpt.maxSide must be >= 0');
    end
end

function [pars0, lb0, ub0, names] = initLineShapePars(opt)
    if ~opt.enable || opt.nSide == 0
        pars0 = []; lb0 = []; ub0 = []; names = {};
        return
    end
    if opt.asymmetric
        nls = 2 * opt.nSide;
        pars0 = zeros(1, nls);
        lb0   = zeros(1, nls);
        ub0   = opt.maxSide * ones(1, nls);
        names = {'lineShapeLeft','lineShapeRight'};
    else
        nls = opt.nSide;
        pars0 = zeros(1, nls);
        lb0   = zeros(1, nls);
        ub0   = opt.maxSide * ones(1, nls);
        names = {'lineShape'};
    end
end

function kern = buildLineShapeKernelLocal(sidePars, opt)
    % Mirror of the nested buildLineShapeKernel for use after lsqnonlin
    % returns (when nested-function workspace is no longer accessible).
    if ~opt.enable || isempty(sidePars)
        kern = 1; return
    end
    nSide = opt.nSide;
    sidePars = max(0, sidePars(:));
    if opt.asymmetric
        left  = sidePars(1:nSide);
        right = sidePars(nSide + (1:nSide));
    else
        left  = sidePars(1:nSide);
        right = left;
    end
    kern = [flipud(left); 1; right];
    if opt.normalize
        s = sum(kern);
        if ~isfinite(s) || s <= 0
            kern = zeros(size(kern)); kern(nSide+1) = 1;
        else
            kern = kern / s;
        end
    end
end
