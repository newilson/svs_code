function [pc0, pc1, p] = peakLinearPhase(ctr, phDeg, w, firstOrder, ppm, center_freq)
% Weighted phase through a set of peaks.  Zeroth order is the w-weighted
% circular mean.  First order puts each phase on the branch nearest that mean,
% then fits a w-weighted line in ppm.  p is the line [slope intercept] in ppm.
% No peaks -> no correction.
    ctr = ctr(:);  phDeg = phDeg(:);  w = w(:);
    pc1 = 0;
    if isempty(ctr)
        pc0 = 0;  p = [0 0];
        return
    end
    pc0 = angle(sum(w .* exp(1i*pi/180*phDeg))) * 180/pi;
    p   = [0 pc0];
    if ~firstOrder || numel(ctr) < 2 || (max(ctr) - min(ctr)) < 10*eps
        return
    end
    phU = pc0 + mod(phDeg - pc0 + 180, 360) - 180;
    V   = [ctr, ones(numel(ctr),1)];
    p   = ((V' * (w .* V)) \ (V' * (w .* phU)))';
    pc0 = polyval(p, center_freq);
    % shiftSpectrumPhase convention: pc1 is the total 1st-order phase across
    % half the bandwidth
    pc1 = p(1) * abs(ppm(end) - ppm(1)) / 2;
