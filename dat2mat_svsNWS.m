function [ws, nws, ppm] = dat2mat_svsNWS(nws_dat,pars)
%
% [ws, nws, ppm] = dat2mat_svs(nws_dat)
%
% Steps:
% filtering
% residual water removal
% eddy current/first order correction
% coil combination
% block averaging
% Fourier Transform
% H2O frequency correction
% resolve averages
% baseline correction

plt = true;
addpath('mapVBVD')
addpath('utils')

if nargin<1 || isempty(nws_dat)
    [file,path] = uigetfile('*.dat','Choose dat file');
    nws_dat = fullfile(path,file);
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
keepOS = [1:vectorSize/4, 1+vectorSize*3/4:vectorSize];
nws = ifft(nws,[],1);
nws = fft(nws(keepOS,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:),[],1);
npts = vectorSize/2;
bw = bw/2;

nave = nws_obj.image.NAve * nws_obj.image.NSet; % cmrr uses 'sets' for 'avgs'

% time axis
t = (0:npts-1)*1/bw;

% freq axis
hz = (-1/2:1/npts:1/2-1/npts)*bw;
ppm = 4.72 + hz/f0;

disp('filtering')
lb = 3;
nws = expFilter(t,lb,nws);


disp('eddy current correction')
eccopt = -1;
nws = eddyCurrentCorrection(nws,nws,eccopt);


disp('coil combination')
mean_nws = mean(reshape(nws,npts,nch,[]),3);
ccopt.minsig_frac = 0.15;
[weights, nws] = coilCombinationNoPC(mean_nws,ccopt.minsig_frac,nws);
if plt
    figure, plot(weights,'*'), title('coil weights')
end

% block averaging before residual water removal to save time - move after
% if concerns about motion
disp('block averaging')
block_size = 2;
nws = blockAverage(nws,block_size);


disp('residual water removal')
if plt
    figure, plot(ppm,real(fftshift(fft(sum(nws,2),[],1),1)),'b'), title('water removal')
    hold on
end
wsopts.type = 'hsvd';
wsopts.hsvd.bounds = [-0.015 0.015]; % normalized bounds for water
wsopts.hsvd.nsin = 25; % number of decaying sinusoids
wsopts.wavelet.zf = 0; % zero filling
wsopts.wavelet.scale = 7;
wsopts.wavelet.type = 'Daubechies';
wsopts.wavelet.par = 10;
wsopts.filt.N = 30;
ws = removeResidualWater(nws,wsopts);
nws = nws - ws; % only water
if plt
    plot(ppm,real(fftshift(fft(sum(ws,2),[],1),1)),'r')
end


disp('Fourier transformation')
ws = fftshift(fft(ws,[],1),1);
nws = fftshift(fft(nws,[],1),1);


[nws, hzshift] = bulkSpectrumShift(nws,hz);
disp(['water frequency correction (hz):'])
disp(hzshift(:))
ws = shiftSpectrumFrequency(ws,hzshift,t);



% Additional denoising
den = 'svd';
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


% Baseline Correction
disp('baseline correction')
baseopts.algo = 'none';
baseopts.method = 'pchip';
baseopts.stepsize = 0.5;
baseopts.windowsize = 0.5;
if plt
    figure, plot(ppm,real(ws),'b'), title('baseline correction')
    hold on
end
ws = baselineCorrect(ws,ppm,baseopts); % real part only
if plt
    plot(ppm,real(ws),'r')
end