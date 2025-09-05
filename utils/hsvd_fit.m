function [freq_Hz, damp_Hz, amp, y_comp, y_fit, lambda] = hsvd_fit(t, y, K, lambda, rangeHz, do_svds)
% HSVDFIT  HSVD decomposition of an NMR/FID signal with ridge (Tikhonov) LS
%
% Inputs
%   t        : Nx1 time vector (s)
%   y        : Nx1 complex signal vector (FID)
%   K        : number of damped sinusoids to retain over the entire range
%   lambda   : (optional) Tikhonov parameter (>=0). If omitted or empty,
%              a small heuristic based on svd(U1) is used.
%   rangeHz  : prune all frequencies outside this range
%   do_svds  : use the svds algorithm else use svd( ,'econ')
%
% Outputs
%   freq_Hz    : K×1 frequencies (Hz)
%   damp_Hz    : K×1 damping factors (1/s, positive = decaying)
%   amp        : K×1 complex amplitudes
%   phase_rad  : K×1 initial phases (rad)
%   y_fit      : Nx1 complex FID reconstructed from the K sinusoids
%   lambda     : scalar lambda actually used
%
% Example
%   [f,d,a,p,fit,lam] = hsvd_fit(t, y, 30);
%   plot(t, real(y), t, real(fit)); legend('data','HSVD fit');

% --- 0.  Basic checks ----------------------------------------------------
t = t(:); y = y(:);
N = length(y);
if length(t) ~= N, error('t and y must have the same length.'); end
dt = t(2) - t(1);

if nargin < 5
    rangeHz = []; % full spectral range
end

% --- 1.  Build the Hankel matrix H ---------------------------------------
L  = floor(N/2);
H  = hankel(y(1:L), y(L:end));   % L×(N-L+1)

% --- 1B. decide on algorithm -----------------------------------------
if nargin < 6 || isempty(do_svds)
    do_svds = K<0.15*L && N>2000;
end

% --- 2. SVD --------------------------------------------------
if do_svds
    [Uk,Sk,~] = svds(H,K,'largest','Tolerance',1e-6,'MaxIterations',2000);
else
    [U,S,~] = svd(H,'econ');         
    Uk = U(:,1:K);
end

% --- 3.  Shifted subspaces and ridge LS for A ----------------------------
U1 = Uk(1:end-1,:);              % (L-1)×K
U2 = Uk(2:end,:);                % (L-1)×K

% Choose lambda if not provided
if nargin < 4 || isempty(lambda)
    s = svd(U1,'econ');
    sK = s(min(K, numel(s)));
    lambda = (0.05 * sK)^2;      % heuristic (try 0.01–0.2 of sK)^2
end

% Augmented system (ridge) using backslash (avoids normal equations)
A = [U1; sqrt(lambda)*eye(K)] \ [U2; zeros(K,K)];   % K×K

% --- 4.  Poles -----------------------------------------------------------
z = eig(A);
z = z( abs(z) <= 1 ); % remove exponentially increasing poles
K = numel(z); % update K
lnz     = log(z);          % ln z = (−d + j2πf)·dt
damp_Hz = -real(lnz) / dt;   % keep original units convention
freq_Hz =  imag(lnz) / (2*pi*dt);

% --- 4B. Prune frequency range
if ~isempty(rangeHz)
    inds = find( freq_Hz > min(rangeHz) & freq_Hz < max(rangeHz) );
    if isempty(inds)
        warning('In HSVD_FIT.m : There are no peaks in the prescribed frequency range')
    end
    lnz = lnz(inds);
    freq_Hz = freq_Hz(inds);
    damp_Hz = damp_Hz(inds);
    K = numel(lnz);
end

% --- 5.  Amplitudes ------------------------------------------------------
n = 0:N-1;
E = exp( lnz * n );          % K×N, since lnz is K×1 and n is 1×N

amp = (E.' \ y);                 % K×1
phase_rad = angle(amp);
y_fit = E.' * amp;

% indidvidual components
y_comp = zeros(length(n),K);
for ii=1:K
    y_comp(:,ii) = E(ii,:) * amp(ii);
end

% --- 6.  Sort (optional) -------------------------------------------------
[~,idx] = sort(freq_Hz);
freq_Hz   = freq_Hz(idx);
damp_Hz   = damp_Hz(idx);
amp       = amp(idx);
phase_rad = phase_rad(idx);
y_comp = y_comp(:,idx);
end
