function L = locateRefExtLine31P(x, y, nominalPpm)
% locateRefExtLine31P - external-phantom line position from the spectrum.
%
%   L = locateRefExtLine31P(x, y, nominalPpm)
%
%   x           ppm axis (PCr = 0), ascending
%   y           complex spectrum on x, pre-phased and in the absorptive frame
%   nominalPpm  position used when no line is found
%
% The real part minus its 0.8 ppm running mean is searched for its maximum in
% 7.3-10 ppm (parabolic refinement between bins).  SNR = that maximum over the
% SD of the same high-passed spectrum (real and imaginary) in 14-19 ppm.
%
%   L.ppm       located position if found, otherwise nominalPpm
%   L.ppmFound  located position (reported even when not found)
%   L.snr       SNR of the located maximum
%   L.found     L.snr >= 5
%   L.nominal   nominalPpm
win      = [7.3 10];
noiseWin = [14 19];
hwPpm    = 0.8;
snrMin   = 5;

x = x(:); y = y(:);
dx = abs(x(2) - x(1));
hp = y - movmean(y, round(hwPpm/dx));
nz = x > noiseWin(1) & x < noiseWin(2);
sd = std([real(hp(nz)); imag(hp(nz))]);
m  = real(hp);
w  = find(x > win(1) & x < win(2));
[pk, k] = max(m(w));
k = w(k);
d = 0;
a = m(k-1); b = m(k); c = m(k+1);
if a - 2*b + c < 0
    d = 0.5*(a - c)/(a - 2*b + c);
end
L = struct('ppm', nominalPpm, 'ppmFound', x(k) + d*dx, 'snr', pk/sd, ...
           'found', pk/sd >= snrMin, 'nominal', nominalPpm);
if L.found
    L.ppm = L.ppmFound;
end
end
