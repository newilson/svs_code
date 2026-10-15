function T = fitSatRecT1(tau, S, fixK)
% FITSATRECT1  Fit a saturation-recovery curve  S(tau) = S0*(1 - k*exp(-tau/T1)).
%   tau   : recovery delays (ms)
%   S     : signal areas (same length)
%   fixK  : true -> k fixed at 1 (perfect saturation), default false
% Returns T with fields S0, T1 (ms), k, their 95% CI half-widths, R2, rmse,
% Sfit (model at tau), tau, S.

if nargin < 3 || isempty(fixK), fixK = false; end

tau = tau(:); S = S(:);
ok = isfinite(tau) & isfinite(S);
tau = tau(ok); S = S(ok);

sc = max(abs(S));
if sc == 0, sc = 1; end
s = S / sc;

T1_0 = median(tau);
opts = optimoptions('lsqcurvefit','Display','off', ...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12,'MaxFunctionEvaluations',5000);

if fixK
    model = @(b,x) b(1) .* (1 - exp(-x ./ b(2)));
    b0 = [max(s), T1_0];  lb = [0, 1e-3];  ub = [Inf, Inf];
else
    model = @(b,x) b(1) .* (1 - b(3)*exp(-x ./ b(2)));
    b0 = [max(s), T1_0, 1];  lb = [0, 1e-3, 0];  ub = [Inf, Inf, 3];
end

[b,resnorm,~,~,~,~,J] = lsqcurvefit(model, b0, tau, s, lb, ub, opts);

np   = numel(b);
dof  = max(numel(s) - np, 1);
mse  = resnorm / dof;
J    = full(J);
ci   = tinv95(dof) * sqrt(diag(mse * inv(J.'*J))); %#ok<MINV>

T = struct();
T.S0    = b(1)*sc;  T.S0_ci = ci(1)*sc;
T.T1    = b(2);     T.T1_ci = ci(2);
if fixK
    T.k = 1;        T.k_ci = 0;
else
    T.k = b(3);     T.k_ci = ci(3);
end
T.R2   = 1 - resnorm / max(sum((s - mean(s)).^2), eps);
T.rmse = sqrt(mse)*sc;
T.tau  = tau;
T.S    = S;
T.Sfit = model(b, tau)*sc;
T.fixK = fixK;
end

function t = tinv95(dof)
tbl = [1 12.71; 2 4.303; 3 3.182; 4 2.776; 5 2.571; 6 2.447; 7 2.365; ...
       8 2.306; 9 2.262; 10 2.228; 15 2.131; 20 2.086; 30 2.042; 60 2.000];
if dof >= 60
    t = 1.96;
elseif dof <= 1
    t = 12.71;
else
    t = interp1(tbl(:,1), tbl(:,2), dof, 'linear');
end
end
