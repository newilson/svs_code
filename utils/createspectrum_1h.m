function fid = createspectrum_1h(cs_ppm, ampl, phase, lb, noiselevel, directory, filename, acqparam)
% createspectrum_1h  Generate a simulated FID and write LCModel .RAW file
%
%   For 1H basis spectra. Adapted from createspectrum.m (PGH, Aug 2001)
%   used in the psimulation 31P framework.
%
%   Inputs:
%     cs_ppm     - chemical shifts in ppm (1 x N)
%     ampl       - amplitudes (1 x N)
%     phase      - phases in radians (1 x N)
%     lb         - Lorentzian line broadening in Hz
%     noiselevel - noise level (0 for basis spectra)
%     directory  - output directory (with trailing filesep)
%     filename   - metabolite name (no extension)
%     acqparam   - struct with .nbpoints, .sw_hz, .hzpppm

plotflag = 1;

if length(cs_ppm) ~= length(ampl)
    error('Lengths of chemical shift and amplitude arrays must be equal');
end

nblines = length(cs_ppm);
fid = 0;

% Time axis
t = (0:acqparam.nbpoints-1) / acqparam.sw_hz;

% Generate FID as sum of damped complex exponentials
% Convention matches createspectrum.m: fid = re - 1i*im (LCModel convention)
% Frequencies are relative to receiveroffset_ppm
offset = acqparam.receiveroffset_ppm;
for ii = 1:nblines
    freq_hz = (cs_ppm(ii) - offset) * acqparam.hzpppm;
    re = ampl(ii) * cos(2*pi*t*freq_hz + phase(ii)) .* exp(-t*pi*lb);
    im = ampl(ii) * sin(2*pi*t*freq_hz + phase(ii)) .* exp(-t*pi*lb);
    fid = fid + re - 1i*im;
end

% Add noise
if noiselevel > 0
    noise = noiselevel * (randn(1, acqparam.nbpoints) + 1i * randn(1, acqparam.nbpoints));
    fid = fid + noise;
end

% Write .RAW file for LCModel
disp(['Creating ' directory filename '.RAW']);

fileid = fopen([directory filename '.RAW'], 'w');
fprintf(fileid, ' $NMID ID=''%s'', FMTDAT=''(2E14.5)''\n', filename);
fprintf(fileid, ' TRAMP=1., VOLUME=1. $END\n');
fprintf(fileid, '%14.5E%14.5E\n', [real(fid); imag(fid)]);
fclose(fileid);

% Display
if plotflag
    figure;

    subplot(2,1,1);
    plot(t*1e3, real(fid), 'b'); hold on;
    plot(t*1e3, imag(fid), 'r');
    xlabel('Time (ms)'); title(filename);
    legend('Real','Imag');

    subplot(2,1,2);
    ft = fftshift(fft(fid));
    fmax = acqparam.sw_hz / 2;
    freq = linspace(fmax, -fmax, acqparam.nbpoints);
    ppm_axis = freq / acqparam.hzpppm + offset;
    plot(ppm_axis, real(ft));
    set(gca, 'Xdir', 'reverse');
    xlabel('Chemical shift (ppm)');
    xlim([6 11]);
end
