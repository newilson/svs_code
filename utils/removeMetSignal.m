function fid = removeMetSignal(fid,freqShiftHz,time,WSopts)
%
% removes specified metabolite signal by frequency demodulation followed by water removal 
% fid must be [npts x whatever]
% freqShiftHz is the metabolite frequency in Hz relative to water
% time is the time values for the fid
% WSopts are options to remove residual water. See removeResidualWater.m


if ~isequal(length(time),size(fid,1))
	error('mismatched size between fid and time axis')
end

if nargin<4
	WSopts = [];
end
	
% Shift metabolite to 0 frequency
spec = fftshift(fft(fid,[],1),1);
spec = shiftSpectrumFrequency(spec,freqShiftHz,time);
fid = ifft(ifftshift(spec,1),[],1);

% Remove metabolite signal
fid = removeResidualWater(fid,WSopts);

% Shift back to original frequency
spec = fftshift(fft(fid,[],1),1);
spec = shiftSpectrumFrequency(spec, -1*freqShiftHz, time);
fid = ifft(ifftshift(spec,1),[],1);
