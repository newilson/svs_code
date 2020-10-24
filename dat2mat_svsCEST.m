function [sig, hdr, offset] = dat2mat_svsCEST(fullfilepath)

addpath('mapVBVD')

if nargin<1 || isempty(fullfilepath)
    [file,path] = uigetfile('*.dat','Choose dat file');
    fullfilepath = fullfile(path,file);
    disp(fullfilepath);
end

[hdr, fid] = mapraw_cestsvs(fullfilepath);

% get signal
firstpt = fid(10,:,:,:,:,:);
TDsig = squeeze(abs(firstpt));
TDsig = reshape(TDsig,size(TDsig,1),[]);
coil_weight = mean(TDsig,2);
coil_weight = coil_weight/sum(coil_weight);

spec = fftshift(fft(fid,[],1),1);

% combine coils - weighted sum of squares
spectemp = 0*squeeze(spec(:,1,:,:,:,:));
for ii=1:size(spec,2)
    spectemp(:,:,:,:,:) = spectemp(:,:,:,:,:) + squeeze(sqrt(coil_weight(ii).*spec(:,ii,:,:,:,:).*conj(spec(:,ii,:,:,:,:))));
end
spec = spectemp; clear spectemp

% freq axis
si = size(spec);
hz = (-1/2:1/si(1):1/2-1/si(1))*hdr.MRS_BW;

spec = reshape(spec,si(1),[]);
sig = zeros(size(spec,2),1);
offset = 0*sig;
for ii=1:size(spec,2)
%     findpeaks(real(spec(:,ii)),hz,'sortstr','descend','npeaks',1,'Annotate','extents');
%     pause
   [amp,loc,width] = findpeaks(real(spec(:,ii)),hz,'sortstr','descend','npeaks',1);
   sig(ii) = amp*width;
   offset(ii) = loc; % check this, maybe should be -loc
end
if length(si)>2
    sig = reshape(sig,si(2:end));
    offset = reshape(offset,si(2:end));
end

