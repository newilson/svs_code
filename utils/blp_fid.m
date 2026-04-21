function [fid_ext, nPredicted] = blp_fid(fid, nPredict, modelOrder)
% BLP_FID  Backward linear prediction to reconstruct missing initial FID points.
%
%   [fid_ext, nPredicted] = blp_fid(fid, nPredict, modelOrder)
%
%   Reconstructs the first few points of a free induction decay (FID) that
%   were lost due to a readout delay (e.g., finite excitation pulse duration
%   in FID-based MRS). This eliminates first-order phase errors and the
%   associated rolling baseline artifact without the need for post-hoc
%   phase correction.
%
%   The algorithm models the FID as a sum of damped complex exponentials
%   and uses truncated-SVD linear prediction to extrapolate backward in
%   time.
%
%   INPUTS
%       fid        - Complex column or row vector, the acquired FID
%                    (1024 points typical for 31P MRS).
%       nPredict   - Number of points to predict backward (e.g., 3-5).
%                    Default: 3.
%       modelOrder - Number of singular values (poles) retained in the
%                    SVD model. Controls the complexity of the signal model.
%                    Rule of thumb: 2x the number of expected resonances.
%                    Set to 0 or [] for automatic selection via singular
%                    value knee detection.
%                    Default: 0 (auto).
%
%   OUTPUTS
%       fid_ext    - Extended FID with predicted points prepended.
%                    Length = length(fid).
%       nPredicted - Actual number of points prepended (same as nPredict).
%
%   USAGE EXAMPLE
%       % Load or define your FID (1024 complex points)
%       fid_raw = load_my_fid();  % your data loading function
%
%       % Predict 3 missing initial points with auto model order
%       [fid_ext, nPred] = blp_fid(fid_raw, 3);
%
%       % Or specify model order manually
%       [fid_ext, nPred] = blp_fid(fid_raw, 3, 25);
%
%       % Now Fourier transform -- should need only zero-order phase correction
%       spec = fftshift(fft(fid_ext));
%
%       % If feeding to LCModel, you can either:
%       %   (a) supply the extended FID and adjust the header accordingly, or
%       %   (b) truncate back to 1024 points (the first nPred points are the
%       %       reconstructed ones, so the original data starts at index nPred+1)
%       fid_for_lcmodel = fid_ext(1:1024);  % option (b): keep original length
%
%   NOTES
%       - The model order is the most important tuning parameter. If the
%         baseline is still wavy, try increasing it. If the predicted points
%         look noisy or unstable, decrease it.
%       - For 31P MRS at 7T (PCr, 3xATP, Pi, PME, PDE, etc.), a model
%         order of 20-40 is typically appropriate.
%       - Auto model order selection works by building the LP matrix at a
%         generous candidate order, computing the SVD, and finding the
%         "knee" in the singular value curve where signal-dominated values
%         transition to the noise floor. It uses a maximum curvature method
%         on the log-scaled singular values.
%       - This function does NOT apply any apodization or zero-filling.
%         Apply those after backward LP if desired.
%
%   REFERENCE
%       Barache D, Antoine JP, Dereppe JM. "The continuous wavelet transform,
%       an alternative to Fourier transform..." J Magn Reson. 1997.
%       (General LP-SVD methodology in NMR context)
%
%   See also: svd, fft, fftshift

% -------------------------------------------------------------------------
%  Input handling
% -------------------------------------------------------------------------
if nargin < 2 || isempty(nPredict),   nPredict   = 3; end
if nargin < 3 || isempty(modelOrder), modelOrder = 0; end  % 0 = auto

autoSelect = (modelOrder == 0);

% Force column vector
fid = fid(:);
N   = length(fid);

% Sanity checks
assert(~isreal(fid), 'FID should be complex-valued.');
if ~autoSelect
    assert(modelOrder < N/2, ...
        'Model order (%d) must be much less than FID length (%d).', modelOrder, N);
end
assert(nPredict >= 1 && nPredict <= 20, ...
    'nPredict=%d seems unreasonable. Expected 1-20.', nPredict);

% -------------------------------------------------------------------------
%  Build the backward LP system
% -------------------------------------------------------------------------
%  Backward LP relation:  x[n] = b_1*x[n+1] + b_2*x[n+2] + ... + b_P*x[n+P]
%
%  We set up an overdetermined system  A * b = d  where:
%    d(n) = fid(n)           for n = 1, 2, ..., N-P
%    A(n,:) = [fid(n+1), fid(n+2), ..., fid(n+P)]
%
%  Then solve via truncated SVD.
%
%  For auto model order selection, we build the matrix at a generous
%  candidate order (Pmax), compute the full SVD, find the knee, then
%  truncate.

if autoSelect
    % Candidate order: N/4 is generous but keeps the system well
    % overdetermined (matrix has N - Pmax rows vs Pmax columns).
    Pmax = min(floor(N/4), 256);
else
    Pmax = modelOrder;
end

nRows = N - Pmax;

% Target vector (left-hand side of the LP equations)
d = fid(1:nRows);

% Prediction matrix
A = zeros(nRows, Pmax);
for k = 1:Pmax
    A(:, k) = fid((1:nRows) + k);
end

% -------------------------------------------------------------------------
%  SVD and model order selection
% -------------------------------------------------------------------------
[U, S, V] = svd(A, 'econ');
svals = diag(S);

if autoSelect
    % -----------------------------------------------------------------
    %  Automatic knee detection on log-scaled singular values
    % -----------------------------------------------------------------
    %  Method: find the point of maximum curvature on the curve
    %  (index, log(sval)). This corresponds to the sharpest bend where
    %  signal singular values give way to the noise floor.
    
    logS = log10(svals);
    nSV  = length(logS);
    
    % Normalise x-axis to [0, 1] so curvature is scale-invariant
    x = (0:nSV-1)' / (nSV - 1);
    y = logS;
    
    % First and second derivatives via finite differences
    % (use 3-point stencil, pad edges)
    dy  = gradient(y, x);
    d2y = gradient(dy, x);
    
    % Curvature: kappa = |y''| / (1 + y'^2)^(3/2)
    kappa = abs(d2y) ./ (1 + dy.^2).^1.5;
    
    % Search for the knee only in a sensible range: between index 5 and
    % 80% of Pmax. Too early can catch initial curvature; too late is
    % noise.
    searchMin = max(5, 1);
    searchMax = round(0.8 * nSV);
    [~, kneeIdx] = max(kappa(searchMin:searchMax));
    kneeIdx = kneeIdx + searchMin - 1;
    
    % Apply a minimum floor so we don't pick something unreasonably small
    P = max(kneeIdx, 10);
    
    fprintf('Auto model order: knee detected at %d (of %d candidate SVs).\n', ...
        P, Pmax);
else
    P = Pmax;  % user-specified, already equals modelOrder
end

% Truncate SVD to the chosen model order
nKeep = P;

Uk = U(:, 1:nKeep);
Sk = S(1:nKeep, 1:nKeep);
Vk = V(:, 1:nKeep);

% LP coefficients
b = Vk * (Sk \ (Uk' * d));

% -------------------------------------------------------------------------
%  Backward extrapolation
% -------------------------------------------------------------------------
fid_ext = zeros(N + nPredict, 1);
fid_ext(nPredict+1 : end) = fid;   % original data in the latter portion

for m = nPredict:-1:1
    % Index m in fid_ext corresponds to time point -(nPredict - m + 1)
    % relative to the original fid(1).
    %
    % The backward LP formula:
    %   fid_ext(m) = b(1)*fid_ext(m+1) + b(2)*fid_ext(m+2) + ... + b(P)*fid_ext(m+P)
    
    idx = m + (1:Pmax);
    fid_ext(m) = b.' * fid_ext(idx);
end

fid_ext = fid_ext(1:N); % keep same number of points

nPredicted = nPredict;

% -------------------------------------------------------------------------
%  Diagnostic output
% -------------------------------------------------------------------------
fprintf('Backward LP: predicted %d points, model order %d%s, FID %d -> %d pts.\n', ...
    nPredict, P, ternary(autoSelect, ' (auto)', ''), N, N + nPredict);

end

% -------------------------------------------------------------------------
%  Helper: inline ternary (for older MATLAB without string expressions)
% -------------------------------------------------------------------------
function out = ternary(cond, valTrue, valFalse)
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end