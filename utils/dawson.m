function D = dawson(z)
%DAWSON  Dawson's integral, D(z) = exp(-z^2) * integral_0^z exp(t^2) dt.
%
%   D = dawson(z) evaluates elementwise for real or complex z; D has the
%   size of z.
%
%   Computed from the Faddeeva function w(z) (fadf.m) via
%       w(z) = exp(-z^2) + (2i/sqrt(pi)) * D(z)
%   For real z this reduces to D = sqrt(pi)/2 * imag(w(z)), which is used
%   directly because it avoids cancelling the two exp(-z^2) real parts.
%   Complex arguments are evaluated in the upper half plane (fadf conjugates
%   lower-half-plane input without undoing it) and reflected back using
%   D(conj(z)) = conj(D(z)).
%
%   Shadows the Symbolic Math Toolbox DAWSON, which does not accept double
%   input.

sz = size(z);
zv = z(:);

if isreal(zv)
    D = sqrt(pi)/2 * imag(fadf(zv));
else
    lo = imag(zv) < 0;
    zv(lo) = conj(zv(lo));
    D = 0.5i*sqrt(pi) * (exp(-zv.^2) - fadf(zv));
    D(lo) = conj(D(lo));
end

D = reshape(D, sz);
end
