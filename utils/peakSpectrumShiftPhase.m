function [spec_sh,hzshift,phfull] = peakSpectrumShiftPhase(spec,hz,f0,lb,peaks,ranges)
%
% 
% spec_sh is [npts x whatever]
% ppm is [npts x 1]
% peaks is center frequencies for singlets and is [npeaks x 1]
% ranges is the lb and ub for singlets and is [npeaks x 2]
%

if ~isequal(length(hz),size(spec,1))
    error('inconsistent dimensions')
end
if ~isequal(length(peaks),size(ranges,1))
    error('inconsistent dimensions')
end
if ~isequal(size(ranges,2),2)
    error('inconsistent range')
end

si = size(spec);
npts = si(1);
npeaks = length(peaks);
bw = (hz(2)-hz(1))*npts;
t = (0:npts-1)*1/bw;
center = 4.72;
ppm0 = hz/f0;

peaks = peaks - center;

spec = reshape(spec,npts,[]);
fid = ifft(ifftshift(spec,1),[],1);
fid = expFilter(t,lb,fid);
specbroad = fftshift(fft(fid,[],1),1);
spec_sh = zeros(size(spec));

hzshift = zeros(size(spec_sh,2),1);
phfull = zeros(size(spec_sh));
for ii=1:size(specbroad,2)
    zpsh = zeros(npeaks,1);
    ph = zeros(npeaks,1);
    for jj=1:npeaks
        ind = find((ppm0)>(peaks(jj)-ranges(jj,1)) & (ppm0)<peaks(jj)+ranges(jj,2));
        [~,zptemp] = max(col(abs(specbroad(ind,ii))));
        zpsh(jj) = ppm0(ind(1)-1+zptemp) - peaks(jj);
        ph(jj) = angle(specbroad(ind(1)-1+zptemp,ii));
    end
    zp = mean(zpsh);
    hzshift(ii) = zp*f0;
    spec_sh(:,ii) = shiftSpectrumFrequency(spec(:,ii),hzshift(ii),t);
    
    p = polyfit(zpsh(:)+peaks(:),ph(:),1);
    phfull(:,ii) = polyval(p,ppm0);
    
    spec_sh(:,ii) = spec_sh(:,ii) .* exp(1i*-phfull(:,ii));
end
