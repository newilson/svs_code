%% demoWaterRemoval
clear,clc
addpath('utils')

%%
% filename = fullfile(pwd,'testData\20200914_DYER_MR_ES_CC017_19.rda');
% [output.hdr, output.complex_fid] = Read_rda_file(filename);
% npts = output.hdr.VectorSize;
% bw = 1e6/output.hdr.DwellTime;
% f0 = output.hdr.MRFrequency;

filename = fullfile(pwd,'testData\S008_NW.MR.RESEARCH_ARMBRUSTER.0003.0001.2021.11.01.12.31.24.125164.1240580336.IMA');
[output.complex_fid, output.hdr, output.short_hdr] = readSiemensDicomSpectrum(filename);
npts = output.short_hdr.npts;
bw = output.short_hdr.bw;
f0 = output.short_hdr.f0;

fid = output.complex_fid;
spec = fftshift(fft(fid,[],1),1);

% time axis
t = (0:npts-1)*1/bw;

% freq axes
hz = (-1/2:1/npts:1/2-1/npts)*bw;
ppm = 4.72 - hz/f0;


% remove residual water
doWS{1} = 'hsvd';
doWS{2} = 'filtA11';
doWS{3} = 'filtA10';
doWS{4} = 'filtA21';

fidWS = zeros(length(fid),length(doWS));
specWS = 0*fidWS;
tE = zeros(length(doWS),1);
for ii=1:length(doWS)
    wsopts.plt = false;
    if strcmp(doWS{ii},'hsvd')
        wsopts.type = 'hsvd';
        wsopts.hsvd.bounds = [-0.025 0.025]; % normalized bounds for water
        wsopts.hsvd.nsin = 25; % number of decaying sinusoids
    elseif strcmp(doWS{ii},'wav')
        wsopts.type = 'wav';
        wsopts.wavelet.zf = 0; % zero filling
        wsopts.wavelet.scale = 7;
        wsopts.wavelet.type = 'Daubechies';
        wsopts.wavelet.par = 10;
    elseif startsWith(doWS{ii},'filtA')
        wsopts.type = 'filt';
        wsopts.filt.type = 'average';
        wsopts.filt.N = str2double(doWS{ii}(6:end));
        % disp([doWS{ii} num2str(wsopts.filt.N)])
    elseif startsWith(doWS{ii},'filtG')
        wsopts.type = 'filt';
        wsopts.filt.type = 'gaussian';
        wsopts.filt.N = str2double(doWS{ii}(6:end));
        wsopts.filt.sigma = wsopts.filt.N/2;
    elseif strcmp(doWS{ii},'filtG')
        wsopts.type = 'filt';
        wsopts.filt.type = 'gaussian';
        wsopts.filt.N = 10;
        wsopts.filt.sigma = 5;
    end

    tS = tic;
    fidWS(:,ii) = removeResidualWater(fid,wsopts);
    tE(ii) = toc(tS)
    specWS(:,ii) = fftshift(fft(fidWS(:,ii),[],1),1);

end

% NWplayplot(real(specWS));
% NWplayplot(real(fftshift(fft(fid,[],1),1)),[],'water removal',real(fftshift(fft(fidws,[],1),1)),[]);

figure, hold on
for ii=1:length(doWS)
    plot(ppm,real(specWS(:,ii)))
end
legend(doWS)
plot(ppm,real(spec),'-r')

rmpath('utils')

