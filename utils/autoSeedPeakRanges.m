function newRange = autoSeedPeakRanges(rspec, ppm, peaks, opt)
% AUTOSEEDPEAKRANGES  Data-driven tightening of peak center ranges.
%
%   newRange = autoSeedPeakRanges(rspec, ppm, peaks, opt)
%
% For each peak listed in opt.peaks, locate the strongest local maximum of the
% background-subtracted real spectrum within a search window and replace that
% peak's row in peaks.range with [center-halfWidth, center+halfWidth].  Peaks
% not listed (or with no positive peak in their window) keep their original
% range.  Used by the blood 1H pipeline so per-dataset frequency drift can't
% leave a metabolite component stranded in a valley (which the >=0 amplitude
% solve in curvefitAuto_varproBaseline then nulls).
%
% Inputs:
%   rspec - real spectrum (vector, any ppm ordering)
%   ppm   - ppm axis (same length as rspec)
%   peaks - struct with .name (cell) and .range ([n x 2])
%   opt   - struct:
%             .bgSmoothPPM - broad-background moving-average width (ppm), def 0.30
%             .peaks       - struct array with fields:
%                              .name      - must match an entry in peaks.name
%                              .search    - [lo hi] ppm search window
%                              .halfWidth - half-width (ppm) of the new range
%
% Output:
%   newRange - peaks.range with matched rows replaced.

newRange = peaks.range;
rspec = rspec(:); ppm = ppm(:);
[ppmS, si] = sort(ppm);
rS = rspec(si);
dppm = median(abs(diff(ppmS)));
if dppm <= 0 || ~isfinite(dppm), return; end

if ~isfield(opt,'bgSmoothPPM') || isempty(opt.bgSmoothPPM), opt.bgSmoothPPM = 0.30; end
wbg = max(3, round(opt.bgSmoothPPM/dppm));
bg  = movmean(rS, wbg);
d   = rS - bg;                                  % peaks above broad background
d   = movmean(d, max(1, round(0.02/dppm)));     % light denoise for picking

specList = opt.peaks;
for i = 1:numel(specList)
    nm = specList(i).name;
    k  = find(strcmpi(peaks.name, nm), 1);
    if isempty(k), continue; end
    sw = specList(i).search;
    hw = specList(i).halfWidth;
    m  = ppmS >= min(sw) & ppmS <= max(sw);
    if ~any(m), continue; end
    dd = d; dd(~m) = -inf;
    [pk, ix] = max(dd);
    if ~isfinite(pk) || pk <= 0, continue; end   % no positive peak -> keep original
    c = ppmS(ix);
    newRange(k,:) = [c-hw, c+hw];
end
end
