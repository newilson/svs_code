function [ws, nws, ppm] = dat2mat_svs(ws_dat,nws_dat,pars)
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
addpath('utils')

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
nws_obj.image.flagRemoveOS = false; % oversampling
nws_obj.image.flagDoAverage = true;
nws_obj.image.flagAverageSets = true; % cmrr sequence uses 'sets' variable for avgs
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
if ndims(nws)>2
    nws = nws(:,:,1);
%     nws = mean(nws,3);
end

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
nws = nws(delay_ind+1:end,:);
nws = nws(1:vectorSize,:);

% remove OS
removeOS = false;
if removeOS
    keepOS = [1:vectorSize/4, 1+vectorSize*3/4:vectorSize];
    nws = ifft(nws,[],1);
    nws = fft(nws(keepOS,:),[],1);
    bw = bw/2;
    npts = vectorSize/2;
else
    npts = vectorSize;
end

% water suppressed
ws_obj = mapVBVD(ws_dat);
if iscell(ws_obj)
    ws_obj = ws_obj{end};
end
ws_obj.image.flagRemoveOS = false; % oversampling
ws_obj.image.flagAverageSets = false; % cmrr sequence uses 'sets' variable for avgs

if ~isequal(f0,1e-6*ws_obj.hdr.Config.Frequency) || ~isequal(nch,ws_obj.image.NCha) || ~isequal(nws_obj.image.NCol,ws_obj.image.NCol)
    error('inconsistent protocols')
end
nave = ws_obj.image.NAve * ws_obj.image.NSet; % cmrr uses 'sets' for 'avgs'

ws = ws_obj.image{''};
ws = ws(delay_ind+1:end,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
ws = ws(1:vectorSize,:,:,:,:,:,:,:,:,:,:,:,:,:,:);

% remove OS
if removeOS
    ws = ifft(ws,[],1);
    ws = fft(ws(keepOS,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:),[],1);
end

% time axis
t = (0:npts-1)*1/bw;

% freq axis
hz = (-1/2:1/npts:1/2-1/npts)*bw;
ppm = 4.72 + hz/f0;

disp('eddy current correction')
eccopt = -1;
[ws, nws] = eddyCurrentCorrection(nws,ws,eccopt);

disp('exponential filtering')
lb = 3;
[ws,filt] = expFilter(t,lb,ws);
nws = expFilter(t,lb,nws);

disp('mlsvd denoising')
if isfield(pars,'domlsvd')
    if isfield(pars.domlsvd,'status')
        if pars.domlsvd.status==true % domlsvd
            addpath('..\tensorlab_2016-03-28');
            if isfield(pars.mlsvd,'singularValues') && length(pars.mlsvd.singularValues)==ndims(ws)
                [Ufull,Sfull,svfull] = mlsvd(ws);
                for ii=1:length(Ufull)
                    Ut{ii} = Ufull{ii}(1:pars.mlsvd.singularValues(ii));
                end
                idx = arrayfun(@(k) 1:size_core(k), 1:numel(size_core), 'UniformOutput', false);
                St  = Sfull(idx{:});     % equivalent to Sfull(1:n1, 1:n2, ..., 1:nN)
                ws = lmlragen(Ut,St);
            else % automatic rank estimation
                [~,~,~,~,Ut,St] = mlsvd_wRankEstNW(ws);
                ws = lmlragen(Ut,St);
            end
            rmpath('..\tensorlab_2016-03-28');
        end
    end
end


disp('coil combination')
ccopt.minsig_frac = 0.05;
[weights, ws, nws] = coilCombinationNoPC(nws,ccopt.minsig_frac,ws);
if plt
    figure, plot(weights,'*'), title('coil weights')
end

disp('block averaging')
block_size = [];
ws = blockAverage(ws,block_size);

disp('residual water removal')

wsopts.type = 'none';
if plt && ~strcmp(wsopts.type,'none')
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
    