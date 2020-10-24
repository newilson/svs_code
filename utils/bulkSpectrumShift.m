function [spec_sh, hzshift] = bulkSpectrumShift(spec,hz)
%
% spec_sh = bulkSpectrumShift(spec,hz)
%
% spec is [npts x whatever]
% hz is [npts x 1]
%
% water should be at 0 hz

if ~isequal(length(hz),size(spec,1))
    error('inconsistent dimensions')
end

si = size(spec);
npts = si(1);
bw = (hz(2)-hz(1))*npts;
t = (0:npts-1)*1/bw;

spec = reshape(spec,npts,[]);
specmag = abs(spec);
spec_sh = zeros(size(spec));

hzshift = zeros(size(specmag,2),1);
for ii=1:size(specmag,2)
%     [~,hzshift(ii)] = findpeaks(specmag(:,ii),hz,'MinPeakHeight',0.98*max(specmag(:,ii)),'Npeaks',1);
    [~,ind] = max(specmag(:,ii),[],1);
    hzshift(ii) = hz(ind);
    spec_sh(:,ii) = shiftSpectrumFrequency(spec(:,ii),hzshift(ii),t);
end

