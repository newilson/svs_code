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

if isempty(peaks)
    spec_ph = spec;
    ph = [0 0];
    return;
end

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
if npeaks>1
    ph = zeros(size(spec));
else
    ph = zeros(1,size(specmag,2));
end

for ii=1:size(specmag,2)
    peakphs = zeros(npeaks,1);
    peakinds = zeros(npeaks,1);
    for jj=1:npeaks
        ind = find((ppm)>(peaks(jj)-ranges(jj,1)) & (ppm)<peaks(jj)+ranges(jj,2));
        % original - [~,maxind] = max(col(specmag(ind,ii)));
        [~,maxind] = findpeaks(specmag(ind,ii),'SortStr','descend','NPeaks',1);
        peakinds(jj) = maxind + ind(1) - 1;
        peakphs(jj) = angle(speclb(peakinds(jj)));
    end
    peakphs = unwrap(peakphs);
    if npeaks==1
        ph(ii) = peakphs;
    else
        p = polyfit(peakinds,peakphs,1);
%         p(:,1) = 0; % is 0 order phase correction necessary here?
        disp(p)
        ph(:,ii) = polyval(p,1:npts);
    end
end
spec_ph = shiftSpectrumPhase(spec,ph);

        

