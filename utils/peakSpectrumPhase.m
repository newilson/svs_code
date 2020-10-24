function [spec_ph,ph] = peakSpectrumPhase(spec,ppm,t,lb,peaks,ranges)
%
% [spec_ph,ph] = peakSpectrumPhase(spec,ppm,t,lb,peaks,ranges)
%
% spec is [npts x whatever]
% ppm is [npts x 1]
% t is [npts x 1]
% peaks is [npeaks x 1]
% ranges is [npeaks x 2]
%
% See also peakSpectrumShift.m

si = size(spec);
npts = si(1);
npeaks = length(peaks);
if ~isequal(npts,length(ppm),length(t))
    error('inconsistent dimension')
end
if ~isequal(npeaks,size(ranges,1))
    error('inconsistent dimensions')
end
if ~isequal(size(ranges,2),2)
    error('inconsistent range')
end

spec = reshape(spec,npts,[]);
fid = ifft(ifftshift(spec,1),[],1);
fid = expFilter(t,lb,fid);
speclb = fftshift(fft(fid,[],1),1);
specmag = abs(speclb);
ph = zeros(size(spec));

for ii=1:size(specmag,2)
    peakphs = zeros(npeaks,1);
    peakinds = zeros(npeaks,1);
    for jj=1:npeaks
        ind = find((ppm)>(peaks(jj)-ranges(jj,1)) & (ppm)<peaks(jj)+ranges(jj,2));
        [~,maxind] = max(col(specmag(ind,ii)));
        peakinds(jj) = maxind + ind(1) - 1;
        peakphs(jj) = angle(speclb(peakinds(jj)));
    end
    p = polyfit(peakinds,peakphs,1);
    ph(:,ii) = polyval(p,1:npts); 
end
spec_ph = shiftSpectrumPhase(spec,ph);

        

