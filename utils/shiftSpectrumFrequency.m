function spec_sh = shiftSpectrumFrequency(spec,hzshift,t)
%
% spec is [npts x whatever]
% hzshift is shift in Hz = shift in ppm * f0

if isequal(hzshift,zeros(size(hzshift)))
    spec_sh = spec;
    return;
end

if ~isequal(length(t),size(spec,1))
    error('inconsistent dimensions')
end

fid = ifft(ifftshift(spec,1),[],1);
si = size(fid);
fid = reshape(fid,si(1),[]);
if length(hzshift>1)
    if ~isequal(length(hzshift),size(fid,2))
        error('inconsistent dimensions')
    else
        for ii=1:size(fid,2)
            fid(:,ii) = fid(:,ii) .* exp(-1i*hzshift(ii)*2*pi*t(:));
        end
    end
else
    for ii=1:size(fid,2)
        fid(:,ii) = fid(:,ii) .* exp(-1i*hzshift*2*pi*t(:));
    end
end

spec_sh = fftshift(fft(fid,[],1),1);
