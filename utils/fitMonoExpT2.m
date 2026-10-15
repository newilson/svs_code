function T = fitMonoExpT2(TE, S)
% FITMONOEXPT2  Fit a mono-exponential decay  S(TE) = S0 * exp(-TE/T2).
%   TE : echo times (ms)
%   S  : signal areas (same length)
% Returns T with fields S0, T2 (ms), S0_ci, T2_ci (95% half-widths), R2,
% rmse, Sfit (model at TE), TE, S.

TE = TE(:); S = S(:);
ok = isfinite(TE) & isfinite(S);
TE = TE(ok); S = S(ok);

sc = max(abs(S));
if sc == 0, sc = 1; end
s = S / sc;

pos = s > 0;
if nnz(pos) >= 2
    p = polyfit(TE(pos), log(s(pos)), 1);
    T2_0 = -1/p(1);
    S0_0 = exp(p(2));
    if ~(T2_0 > 0 && isfinite(T2_0)), T2_0 = max(TE); end
else
    S0_0 = max(s); T2_0 = max(TE);
end

model = @(b,te) b(1) .* exp(-te ./ b(2));
opts = optimoptions('lsqcurvefit','Display','off', ...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12,'MaxFunctionEvaluations',5000);
[b,resnorm,~,~,~,~,J] = lsqcurvefit(model, [S0_0, T2_0], TE, s, [0, 1e-3], [Inf, Inf], opts);

dof  = max(numel(s) - 2, 1);
mse  = resnorm / dof;
J    = full(J);
covb = mse * inv(J.'*J); %#ok<MINV>
ci   = tinv95(dof) * sqrt(diag(covb));

T = struct();
T.S0    = b(1)*sc;  T.S0_ci = ci(1)*sc;
T.T2    = b(2);     T.T2_ci = ci(2);
T.R2    = 1 - resnorm / max(sum((s - mean(s)).^2), eps);
T.rmse  = sqrt(mse)*sc;
T.TE    = TE;
T.S     = S;
T.Sfit  = model(b, TE)*sc;
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
