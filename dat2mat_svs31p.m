% dat2mat_svs31p(dat_file);
clear,clc

% EDIT THESE AS NEEDED
mag_mode = false; % use absolute value spectrum only?
npeaks = 6;

fullfilepath = [];

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
nch = twix_obj.image.NCha; % should be 1
if ~isequal(nch,1)
    warning(['Number of coils is ' num2str(nch) '. Should be 1.'])
end

nave = twix_obj.image.NAve;

fid = twix_obj.image{''};

% zero fill
fid = cat(1,fid,0*fid);

si = size(fid);
npts = si(1);
fid = reshape(fid,npts,[]); % makes fid 2D - check this depending on combination of other variables

% time axis
t = (0:npts-1)*1/bw;

% freq axis
hz = (-1/2:1/npts:1/2-1/npts)*bw;
ppm = hz/f0;

% apodization
lb = str2double(inputdlg('Enter line broadening (Hz): ','Apodization'));
if isempty(lb) || isnan(lb), lb = 0; end
[fid,filt] = expFilter(t,lb,fid);

% moving average
avg_block = round(str2double(inputdlg('Enter number of scans to average: ', 'Moving average')));
if isempty(avg_block) || isnan(avg_block), avg_block = 1; end
if avg_block>1 && avg_block<=size(fid,2)
    fidtemp = zeros(npts,size(fid,2)-avg_block+1);
    for ii=1:size(fid,2)-avg_block+1
        fidtemp(:,ii) = mean(fid(:,ii:ii+avg_block-1),2);
    end
    fid = fidtemp;
end

% FT
spec = fftshift(fft(fid,[],1),1);
if mag_mode
    spec = abs(spec);
end

den = 'svd';
switch den
    case 'svd'        
        % hsvd denoising
        if size(fid,2)>1
            f = NWsvdTS(spec,ppm,'SVD denoising: save and close when done');
            waitfor(f);
            spec = out; clear out
        end        
    case 'wav'
        % wavelet filtering
        if exist('wmaxlev')
            f = NWwavden(spec,ppm,'Wavelet denoising: save and close when done');
            waitfor(f);
            spec = out; clear out
        end
end

if ~mag_mode
    % phase correction
    f = NWman_phase(spec,ppm,'Manual phase correction: save and close when done');
    waitfor(f);
    spec = out; clear out
end

% baseline correction
f = NWsemiman_base(spec,ppm,'Semi-manual baseline correction: save and close when done');
waitfor(f);
spec = out; clear out

% semi automatic fitting
if mag_mode
    fitmode = 0;
else
    fitmode = 2; % 0: magnitude, 1: gaussians, 2: lorentzians
end
count = 0;
keepFitting = true;
while keepFitting
    if count>0
        for ii=1:size(spec,2)
            if ii==1
                [tempfit,n,names,tempheight,temppos,tempwidth,tempampl] = curvefitIMCL(ppm,real(spec(:,ii)),0,fitmode);
                specfit = zeros(size(spec)); specfit(:,1) = tempfit;
                height = zeros(n,size(spec,2)); height(:,1) = tempheight;
                pos = zeros(n,size(spec,2)); pos(:,1) = temppos;
                width = zeros(n,size(spec,2)); width(:,1) = tempwidth;
                ampl = zeros(n,size(spec,2)); ampl(:,1) = tempampl;
            else
                [specfit(:,ii),n,names,height(:,ii),pos(:,ii),width(:,ii),ampl(:,ii)] = curvefitIMCL(ppm,real(spec(:,ii)),0,fitmode);
            end
        end
    else
        ampl0 = zeros(npeaks,size(spec,2));
        width0 = 0*ampl0;
        pos0 = 0*ampl0;
        count = count+1;
        specfit = 0*spec;
        ampl = 0*ampl0;
        pos = 0*pos0;
        width = 0*ampl;
        height = 0*ampl;
        for ii=1:size(spec,2)
            noise = std(spec(1:50,ii));
            %             findpeaks(real(spec(:,ii)),ppm,'WidthReference','halfheight','MinPeakHeight',2*noise,'NPeaks',npeaks,'SortStr','descend','MinPeakDistance',0.6,'annotate','extent');
            [pks,locs,wid,p] = findpeaks(real(spec(:,ii)),ppm,'WidthReference','halfheight','MinPeakHeight',2*noise,'NPeaks',npeaks,'SortStr','descend','MinPeakDistance',0.6);
            width0(:,ii) = wid(:)/2; % hwhm
            ampl0(:,ii) = pks(:) .* wid(:)/2;
            pos0(:,ii) = locs(:);
        end
        for ii=1:size(spec,2)
            [specfit(:,ii),ampl(:,ii),pos(:,ii)] = mycurvefit(ppm(:),spec(:,ii),ampl0(:,ii),pos0(:,ii),width0(:,ii),fitmode);
        end
    end
    
    f1 = NWplayplot(real(spec),ppm,'fit and spectrum',specfit,ppm);
    f2 = NWplayplot(real(spec)-specfit,ppm,'fit and residual',specfit,ppm);
    keepFitting = false;
    answer = questdlg('Try again?','Fit','Yes','No','No');
    if strcmp(answer,'Yes')
        keepFitting = true;
    else
        keepFitting = false;
    end
    close(f1);
    close(f2);
end

% reshaping
spec = reshape(spec,si);
specfit = reshape(specfit,si);

% write to file
if exist('writematrix')
    M = [pos(:) ampl(:)];
    writematrix(M,fullfile(path,'autofit.csv'));
else
    T = table(pos(:),ampl(:),'VariableNames',{'ppm','signal'});
    writetable(T,fullfile(path,'autofit.csv'));
end