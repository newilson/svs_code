function [ws, nws, ppm] = dat2mat_svs(ws_dat,nws_dat)
%
% [ws, nws, ppm] = dat2mat_svs(ws_dat,nws_dat)
%
% Steps:
% eddy current/first order correction
% filtering
% coil combination
% block averaging
% residual water removal
% Fourier Transform
% H2O frequency correction
% spectrum-based phase correction/frequency shifting
% resolve averages
% baseline correction

plt = true;
addpath('mapVBVD')

if nargin<2
    [file,path] = uigetfile('*.dat','Choose NWS dat file');
    nws_dat = fullfile(path,file);
    [file,path] = uigetfile('*.dat','Choose WS dat file');
    ws_dat = fullfile(path,file);
end

% non water suppressed
nws_obj = mapVBVD(nws_dat);
if iscell(nws_obj)
    nws_obj = nws_obj{end};
end
nws_obj.image.flagDoAverage = true;
bw = 1e9/nws_obj.hdr.Config.DwellTimeSig;
f0 = 1e-6*nws_obj.hdr.Config.Frequency;
nch = nws_obj.image.NCha;
npts = nws_obj.image.NCol;

nws = nws_obj.image{''};

% water suppressed
ws_obj = mapVBVD(ws_dat);
if iscell(ws_obj)
    ws_obj = ws_obj{end};
end
if ~isequal(bw,1e9/ws_obj.hdr.Config.DwellTimeSig) || ~isequal(f0,1e-6*ws_obj.hdr.Config.Frequency) || ~isequal(nch,ws_obj.image.NCha) || ~isequal(npts,ws_obj.image.NCol)
    error('inconsistent protocols')
end
nave = ws_obj.image.NAve;

ws = ws_obj.image{''};

% time axis
t = (0:npts-1)*1/bw;

% freq axis
hz = (-1/2:1/npts:1/2-1/npts)*bw;
ppm = 4.72 + hz/f0;

disp('eddy current correction')
eccopt = 0;
[ws, nws] = eddyCurrentCorrection(nws,ws,eccopt);

disp('filtering')
lb = 7;
[ws,filt] = expFilter(t,lb,ws);
nws = expFilter(t,lb,nws);

disp('coil combination')
ccopt.minsig_frac = 0.2;
[weights, ws, nws] = coilCombinationNoPC(nws,ccopt.minsig_frac,ws);
if plt
    figure, plot(weights,'*')
end

disp('block averaging')
block_size = 32;
ws = blockAverage(ws,block_size);

disp('residual water removal')
if plt
    figure, plot(ppm,real(fftshift(fft(sum(ws,2),[],1),1)),'b')
    hold on
end
wsopts.type = 'hsvd';
wsopts.hsvd.bounds = [-0.025 0.025]; % normalized bounds for water
wsopts.hsvd.nsin = 25; % number of decaying sinusoids
wsopts.wavelet.zf = 0; % zero filling
wsopts.wavelet.scale = 7;
wsopts.wavelet.type = 'Daubechies';
wsopts.wavelet.par = 10;
ws = removeResidualWater(ws,wsopts);
if plt
    plot(ppm,real(fftshift(fft(sum(ws,2),[],1),1)),'r')
end

disp('Fourier transformation')
ws = fftshift(fft(ws,[],1),1);
nws = fftshift(fft(nws,[],1),1);

disp('water frequency correction')
[nws, hzshift] = bulkSpectrumShift(nws,hz);
ws = shiftSpectrumFrequency(ws,hzshift,t);

if plt
    figure, plot(ppm,real(sum(ws,2)),'b');
    hold on
end
disp('individual spectra freq/ph correction')
% peaks = [3.03 3.21 3.91 3.42]; % cr-cho-cr-tau
peaks = [3.03 3.21 3.91 2.02]; % cr-cho-cr-naa
ranges = 0.1*ones(length(peaks),2);
addlb = 4; % additional line broadening to locate peaks
[ws,hzshift] = peakSpectrumShift(ws,hz,f0,addlb,peaks,ranges);
if plt
    plot(ppm,real(sum(ws,2)),'k')
end
ws = peakSpectrumPhase(ws,ppm,t,addlb,peaks,ranges/2);
if plt
    plot(ppm,real(sum(ws,2)),'r')
    figure, plot(hzshift,'-*')
end

disp('resolve block averages')
ws = mean(ws,2);

disp('baseline correction')
baseopts.method = 'pchip';
baseopts.stepsize = 0.5;
baseopts.windowsize = 0.5;
if plt
    figure, plot(ppm,real(ws),'b')
    hold on
end
% ws = baselineCorrect(ws,ppm,baseopts); % real part only
if plt
    plot(ppm,real(ws),'r')
end

    