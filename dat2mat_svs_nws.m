function [ws, nws, ppm] = dat2mat_svs_nws(nws_dat,pars)
%
% [ws, nws, ppm] = dat2mat_svs_nws(nws_dat)
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
addpath('utils')

if nargin<1 || isempty(nws_dat)
    [file,path] = uigetfile('*.dat','Choose NWS dat file');
    nws_dat = fullfile(path,file);
end

if nargin<2 || isempty(pars)
    pars.plt = false;
end
if ~isfield(pars,'removeOS')
    pars.removeOS = false;
end
if ~isfield(pars,'lb')
    pars.lb = 1;
end
if ~isfield(pars,'ccopt')
    pars.ccopt.minsig_frac = 0.05;
end
if ~isfield(pars,'block_size')
    pars.block_size = []; % no block averaging
end

% non water suppressed
nws_obj = mapVBVD(nws_dat);
if iscell(nws_obj)
    nws_obj = nws_obj{end};
end
nws_obj.image.flagRemoveOS = false; % oversampling
nws_obj.image.flagDoAverage = false;
nws_obj.image.flagAverageSets = false; % cmrr sequence uses 'sets' variable for avgs
bw = 1e9/nws_obj.hdr.Config.DwellTimeSig;
f0 = 1e-6*nws_obj.hdr.Config.Frequency;
nch = nws_obj.image.NCha;
npts = nws_obj.image.NCol;
vectorSize = 2*nws_obj.hdr.MeasYaps.sSpecPara.lVectorSize; % 2 for oversampling
if ~isequal(npts,vectorSize)
    warning(['Actual points is: ' num2str(npts) newline 'Prescribed points is: ' num2str(vectorSize)])
end
seqname = nws_obj.hdr.MeasYaps.tSequenceFileName;

nave = nws_obj.image.NAve * nws_obj.image.NSet; % cmrr uses 'sets' for 'avgs'

nws = nws_obj.image{''};

% timing correction - needed for eja_slaser_svs
if contains(seqname,'eja_svs_slaser')
% nwstemp = sum(nws,2);
% nwstemp = movmean(nwstemp,32);
% [~,delay_ind] = max(abs(nwstemp)); % 192
%     delay_ind = npts - vectorSize; % 256
    delay_ind = nws_obj.noise.iceParam(5); % 230
else
    delay_ind = 0;
end
nws = nws(delay_ind+1:end,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
nws = nws(1:vectorSize,:,:,:,:,:,:,:,:,:,:,:,:,:,:);

% remove OS
if pars.removeOS
    keepOS = [1:vectorSize/4, 1+vectorSize*3/4:vectorSize];
    nws = ifft(nws,[],1);
    nws = fft(nws(keepOS,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:),[],1);
    bw = bw/2;
    npts = vectorSize/2;
else
    npts = vectorSize;
end

% time axis
t = (0:npts-1)*1/bw;

% freq axis
hz = (-1/2:1/npts:1/2-1/npts)*bw;
ppm = 4.72 + hz/f0;

disp('zero order phase correction')
pars.eccopt = -1;
[~, nws] = eddyCurrentCorrection(nws,nws,pars.eccopt);

disp('filtering')
[nws,filt] = expFilter(t,pars.lb,nws);

disp('coil combination')
[weights, ~, nws] = coilCombinationNoPC(nws,pars.ccopt.minsig_frac,nws);
if pars.plt
    figure, plot(weights,'*'), title('coil weights')
end

disp('block averaging')
nws = blockAverage(nws,pars.block_size);

disp('water signal separation')

wsopts.type = 'none';
if pars.plt && ~strcmp(wsopts.type,'none')
    figure, plot(ppm,real(fftshift(fft(sum(ws(:,:),2),[],1),1)),'b'), title('water removal')
    hold on
end
wsopts.hsvd.bounds = [-0.015 0.015]; % normalized bounds for water
wsopts.hsvd.nsin = 25; % number of decaying sinusoids
wsopts.wavelet.zf = 0; % zero filling
wsopts.wavelet.scale = 7;
wsopts.wavelet.type = 'Daubechies';
wsopts.wavelet.par = 10;
wsopts.filt.N = 30;
ws = removeResidualWater(ws,wsopts);
if plt && ~strcmp(wsopts.type,'none')
    plot(ppm,real(fftshift(fft(sum(ws(:,:),2),[],1),1)),'r')
end

disp('Fourier transformation')
ws = fftshift(fft(ws,[],1),1);
nws = fftshift(fft(nws,[],1),1);

[nws, hzshift] = bulkSpectrumShift(nws,hz);
hzshift = 0; 
disp(['water frequency correction: ' num2str(hzshift) ' hz'])
ws = shiftSpectrumFrequency(ws,hzshift,t);

if plt && ~isequal(hzshift,0*hzshift)
    figure, plot(ppm,real(sum(ws(:,:),2)),'b'); title('frequency/phase correction')
    hold on
end
disp('individual spectra freq/ph correction')
% peaks = [3.03 3.21 3.91 3.42]; % cr-cho-cr-tau
% peaks = [3.03 3.21 3.91 2.02]; % cr-cho-cr-naa
peaks = [];
ranges = 0.1*ones(length(peaks),2);
addlb = 3; % additional line broadening to locate peaks
[ws,hzshift] = peakSpectrumShift(ws,hz,f0,addlb,peaks,ranges);
for ii=1:length(hzshift)
    disp(['      freq shift: ' num2str(hzshift(ii)) ' hz'])
end
if plt && ~isequal(hzshift,0*hzshift)
    plot(ppm,real(sum(ws(:,:),2)),'k')
end
ws = peakSpectrumPhase(ws,ppm,t,addlb,peaks,ranges/2);
if plt && ~isequal(hzshift,0*hzshift)
    plot(ppm,real(sum(ws(:,:),2)),'r'), legend({'orig','freq corr','phase corr'})
%     figure, plot(hzshift,'-*')
end

% Additional denoising
den = 'none';
switch den
    case 'svd'        
        % svd denoising
        if size(ws,2)>1
            f = NWsvdTS(ws,ppm,'SVD denoising: save and close when done');
            waitfor(f);
            ws = out; clear out
        end        
    case 'wav'
        % wavelet filtering
        if exist('wmaxlev')
            f = NWwavden(ws,ppm,'Wavelet denoising: save and close when done');
            waitfor(f);
            ws = out; clear out
        end
end

disp('resolve block averages')
ws = squeeze(mean(ws,2));

disp('baseline correction')
baseopts.algo = 'none';
baseopts.method = 'pchip';
baseopts.stepsize = 0.5;
baseopts.windowsize = 0.5;
if plt && ~strcmp(baseopts.algo,'none')
    figure, plot(ppm,real(ws),'b'), title('baseline correction')
    hold on
end
ws = baselineCorrect(ws,ppm,baseopts); % real part only
if plt && ~strcmp(baseopts.algo,'none')
    plot(ppm,real(ws),'r')
end
    