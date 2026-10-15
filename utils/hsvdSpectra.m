function Sc = hsvdSpectra(fHz, dHz, amp, N, dt)
%HSVDSPECTRA  Full-length spectrum of every HSVD component, one column each.
%
%   Sc = hsvdSpectra(fHz, dHz, amp, N, dt)
%
%   fHz  component frequencies (Hz, relative to the axis center)
%   dHz  component damping rates (rad/s; Lorentzian FWHM in Hz is dHz/pi)
%   amp  complex component amplitudes (FID units)
%   N    points on the spectrum axis
%   dt   dwell time (s)
%
% Evaluating all N points (not just a window) is what allows a component
% centered outside a region to contribute its tail inside that region.
%
% Used by the coarse HSVD in dat2mat_svsNAD and by compare_hsvd_varpro /
% plot_hsvd_varpro_fits to rebuild an HSVD model spectrum from the saved
% component list.
    tf = (0:N-1)' * dt;
    Yc = exp(tf * (-dHz(:).' + 1i*2*pi*fHz(:).')) .* amp(:).';
    Sc = fftshift(fft(Yc, [], 1), 1);
end
