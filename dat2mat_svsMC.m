% dat2mat_svsMC(dat_file);
clear,clc

% EDIT THESE AS NEEDED

% fullfilepath = [];
fullfilepath = 'C:\Users\z0043etb\Desktop\WILSON_NEIL_20_01_24-16_36_06-STD-1_3_12_2_1107_5_2_0_79106\RESEARCH_PBAGGA_20200124_163752_091000\meas_MID00216_FID03634_metcycst_2stepPC.dat'; % NEIL brain
% fullfilepath = 'C:\Users\z0043etb\Desktop\meas_MID00272_FID01824_metcycst_upfieldTest2step.dat'; % PHANTOM

addpath('mapVBVD')

if isempty(fullfilepath)
    [file,path] = uigetfile('*.dat','Choose dat file');
    fullfilepath = fullfile(path,file);
    disp(fullfilepath);
end

twix_obj = mapVBVD(fullfilepath);
if iscell(twix_obj)
    twix_obj = twix_obj{end};
end

twix_obj.image.flagRemoveOS = false;

if isfield(twix_obj.hdr.Config,'DwellTime')
    bw = 1e9/twix_obj.hdr.Config.DwellTime;
elseif isfield(twix_obj.hdr.Config,'DwellTimeSig')
    bw = 1e9/twix_obj.hdr.Config.DwellTimeSig;
else
    error('unknown dwell time')
end
if twix_obj.image.flagRemoveOS
    bw = bw/2;
end
f0 = 1e-6*twix_obj.hdr.Config.Frequency;
nch = twix_obj.image.NCha;
nave = twix_obj.image.NAve; % MC dimension

fid = twix_obj.image{''}; % pts-ch-avgs-reps

si = size(fid);
npts = si(1);

% time axis
t = (0:npts-1)*1/bw;

% freq axis
hz = (-1/2:1/npts:1/2-1/npts)*bw;
ppm = hz/f0 + 4.7;

% preliminary coil signal scaling - for DEBUGGING ONLY
coilsig1 = sum(fid(1,:,1:2:end,1),3).*conj(sum(fid(1,:,1:2:end,1),3));
coilsig2 = sum(fid(1,:,2:2:end,1),3).*conj(sum(fid(1,:,2:2:end,1),3));
coilsig = 0.5 * (coilsig1 + coilsig2);
% figure, plot(1:length(coilsig),coilsig1,'o',1:length(coilsig),coilsig2,'*',1:length(coilsig),coilsig,'s')
coilweight = coilsig / sum(coilsig);
coilweight = sqrt( coilweight / max(coilweight) );
coilweight( coilweight < 0.05 ) = 0; % ignore coils with less than 5% of max signal
maxcoil = find(coilweight==1);
[sorted_coilweights,sorted_coilinds] = sort(coilweight,'descend');

% some debug plotting
NWplayplot(real(fftshift(fft(squeeze(fid(:,maxcoil,:)),[],1),1)));
figure, plot(ppm,real(fftshift(fft(squeeze(sum(fid(:,sorted_coilinds(5),:),3)),[],1),1)));

% noise per coil per average
noisepts = 400;
noise = std(fid(end-noisepts:end,:,:,:),1); % last 400 complex points in fid

% apodize here if desired
lb = 0;
fid = expFilter(t,lb,fid);

% pointwise SNR
pointSNR = abs(fid)./noise;
inds = pointSNR<100; % empiric SNR threshold

% linear fit fid phase
sigpts = 400;
fid0 = fid(1:sigpts,:,:,:);
ang = unwrap(angle(fid0),pi,1);
p = NWpolyfitim(1,t(1:sigpts),permute(ang,[2 3 4 1]));

% debug fitting
ang_fit = NWpolyvalim(p,t(1:sigpts));
ang_fit = permute(ang_fit,[ndims(ang_fit) 1:ndims(ang_fit)-1]);
figure
for ii=1:length(sorted_coilinds)
    plot(1:si(1),unwrap(angle(fid(:,sorted_coilinds(ii),1)),pi,1),'b',1:sigpts,ang_fit(:,sorted_coilinds(ii),1,:),'r');
    pause;
end
t_p = repmat(t(:),[1 size(p,1) size(p,2)]);
t_p = permute(t_p,[2:ndims(t_p) 1]);
phcorr = exp(-1i*(p(:,:,1).*t_p + p(:,:,2)));
phcorr = permute(phcorr,[ndims(phcorr) 1:ndims(phcorr)-1]);
fidcorr = fid .* phcorr;
fidodd = -1*fidcorr(:,:,1:2:end,:);
fideven = fidcorr(:,:,2:2:end,:);

% figure, plot(1:si(1),angle(fidcorr(:,maxcoil,1)))
figure,
for ii=1:length(sorted_coilinds)
    plot(1:si(1),real(fftshift(fft(fid(:,sorted_coilinds(ii),1),[],1),1)),'b',1:si(1),real(fftshift(fft(fidcorr(:,sorted_coilinds(ii),1),[],1),1)),'r')
    pause;
end
% figure, plot(ppm,real(fftshift(fft(squeeze(sum(fidcorr(:,sorted_coilinds(5),:),3)),[],1),1)));

% Scale even/odd
scpts = 32;
mag1 = sum(abs(fidodd(1:scpts,:,:,:)),1);
mag2 = sum(abs(fideven(1:scpts,:,:,:)),1);
ratio = mag1 ./ mag2;
fideven = fideven .* ratio;

% Weight and Combine coils
coilthresh = 0.05;
coilsig = 0.5 * ( sum(mag1.^2,3) + sum(mag2.^2,3) );
coilweight = coilsig / sum(coilsig(:));
coilweight = sqrt( coilweight / max(coilweight) );
coilweight( coilweight < coilthresh ) = 0; % ignore coils with less than threshold of max signal

fidodd = sum(fidodd .* coilweight,2);
fideven = sum(fideven .* coilweight,2);

% Resolve averages
fidodd = sum(fidodd,3);
fideven = sum(fideven,3);
metfid = fidodd + fideven;
watfid = fideven - fidodd;

% Fourier transform
watspec = fftshift(fft(watfid,[],1),1);
metspec = fftshift(fft(metfid,[],1),1);

% final plotting
% figure, subplot(1,2,1), plot(ppm,real(watspec)), axis tight,
% subplot(1,2,2), plot(ppm,real(metspec)), axis tight
figure, plot(ppm,real(fftshift(fft(squeeze(sum(fid(:,sorted_coilinds(1),:),3)),[],1),1))); 
% xlim([1 4.5])
figure, plot(ppm,real(metspec)), 
% xlim([1 4.5])

% water_extent = 0.5;
% range = 0 * ppm;
% range(ppm<water_extent & ppm>-water_extent) = 1;
% lb = 0.80; ub = 1.2;
% switch coilweight % scale even scans
%     case 'individual'
%         F = zeros(1,npairs);
%         for ii = 1:npairs
%             curodd = odd(:,ii); cureven = even(:,ii);
%             objfun = @(F)abs(curodd - F*cureven).*range;
%             F(ii) = fminbnd(objfun,lb,ub);
%             even(:,ii) = even(:,ii) * F(ii);
%         end
%     case 'all'
%         oddsum = sum(odd,2);
%         evensum = sum(even,2);
%         objfun = @(F)abs(oddsum - F*evensum).*range;
%         F = fminbnd(objfun,lb,ub);
%         even = even * F;
% end



