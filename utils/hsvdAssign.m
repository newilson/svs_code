function [assigned, found, sig, phDeg, ctr, wdt, nComp, selMat] = hsvdAssign( ...
        ppmC, lwHz, phC, amp, f0, peakRange, lwMin, lwMax, phTol)
% Assign HSVD components to peaks.  A component qualifies for peak kk when its
% center lies in that peak's ppm window AND its linewidth and phase are inside
% the per-peak caps.  Everything left unassigned is baseline.
%
% Note what "unassigned" means: the test is against the PEAK TABLE, not against
% any region window.  A component sitting inside a region but under no peak is
% baseline, and so is a component centered far outside -- which is exactly what
% lets a region's baseline carry the tails of resonances outside it.
%
% ctr/wdt/phDeg are only meaningful where found(kk) is true; the caller keeps
% its own fallback seeds for the rest.
    np       = size(peakRange, 1);
    assigned = false(numel(ppmC), 1);
    found    = false(np, 1);
    sig      = zeros(np, 1);
    phDeg    = zeros(np, 1);
    ctr      = mean(peakRange, 2);
    wdt      = nan(np, 1);
    nComp    = zeros(np, 1);
    selMat   = false(numel(ppmC), np);
    for kk = 1:np
        sel = ppmC >= min(peakRange(kk,:)) & ppmC <= max(peakRange(kk,:)) & ...
              lwHz >= lwMin(kk)            & lwHz <= lwMax(kk)            & ...
              phC  >= -phTol(kk)           & phC  <= phTol(kk);
        selMat(:,kk) = sel;
        nComp(kk)    = sum(sel);
        if ~any(sel), continue; end
        assigned  = assigned | sel;
        found(kk) = true;
        A = sum(amp(sel));
        w = abs(amp(sel));
        sig(kk)   = abs(A);
        phDeg(kk) = angle(A) * 180/pi;
        ctr(kk)   = sum(w .* ppmC(sel)) / sum(w);
        wdt(kk)   = sum(w .* lwHz(sel)) / sum(w) / f0;
    end
