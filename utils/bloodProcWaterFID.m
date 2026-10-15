function spec = bloodProcWaterFID(fid, t, lb, f_water_Hz)
% BLOODPROCWATERFID  Process one blood water FID the way bru2mat_svsNAD
% processes the water reference: Klose ECC against itself, water restored to
% its nominal offset, exponential line broadening, FFT.
%
%   fid        - complex FID (column)
%   t          - time axis (s), same length as fid
%   lb         - exponential line broadening (Hz)
%   f_water_Hz - water offset from the carrier (4.7*f0 - O1)

[~, fid] = eddyCurrentCorrection(fid, fid, 0, 1);
spec = fftshift(fft(fid, [], 1), 1);
spec = shiftSpectrumFrequency(spec, -f_water_Hz, t);
fid  = ifft(ifftshift(spec, 1), [], 1);
fid  = expFilter(t, lb, fid);
spec = fftshift(fft(fid, [], 1), 1);
end
