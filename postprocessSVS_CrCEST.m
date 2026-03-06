function output = postprocessSVS_CrCEST(pathname,pars)
%
% CrCEST experiment using prepSVS sequence
% We look at water signal change with MTR asymmetry analysis
%
% NW - 08/2025

% check for parallel pool
doParFor = false;
try
    checkParTool = canUseParallelPool;
catch
    checkParTool = false;
end

if checkParTool
    if ~isempty(gcp('nocreate')) || Z>1
        doParFor = true;
    end
end

% in case current directory is not this one
thisFile = mfilename('fullpath');
thisPath = fileparts(thisFile);
addpath([thisPath filesep 'utils']);

if nargin<2 || ~isfield(pars,'fitmode')
    pars.fitmode = 1; % 0: 1st time pt, 1: max spectrum, 2: HSVD
end

if nargin<2 || ~isfield(pars,'peak')
    pars.peak = 0; % ppm - usually looking at water at 0 ppm
end

if nargin<2 || ~isfield(pars,'dummyShots')
    pars.dummyShots = 0;
end

if nargin<2 || ~isfield(pars,'denoiseSeries')
    pars.denoiseSeries = false;
end

if nargin<1 || isempty(pathname)
    pathname = uigetdir('Choose Spectroscopy Dicom Directory');
end
output.pathname = pathname;

% read files
files = dir(pathname);
fid = [];
count = 0;
for ii=1:length(files)
    fname = fullfile(files(ii).folder, files(ii).name);
    [~,~,ext] = fileparts(fname);
    if (~files(ii).isdir && (strcmp(ext,'.IMA')||strcmp(ext,'.dcm')) )
        [tempFid, hdr, shortHdr] = readSiemensDicomSpectrum(fname);
        fid = cat(2,fid,tempFid);        
        count = count+1;
        if count>1
            if ~isequal(npts,shortHdr.npts) || ~isequal(bw,shortHdr.bw) || ~isequal(f0,shortHdr.f0) || ~isequal(B1,shortHdr.WIPlong(9)) || ~isequal(satDur,shortHdr.WIPlong(6)) ...
                    || ~isequal(startPPM,shortHdr.WIPdbl(3)) || ~isequal(stopPPM,shortHdr.WIPdbl(4)) || ~isequal(stepPPM, shortHdr.WIPdbl(5))
                error('dicoms do not appear to be from the same scan')
            end
        end
        npts = shortHdr.npts;
    	bw = shortHdr.bw;
    	f0 = shortHdr.f0;
        satDur = shortHdr.WIPlong(6);
        B1 = shortHdr.WIPlong(9);
        startPPM = shortHdr.WIPdbl(3);
        stopPPM = shortHdr.WIPdbl(4);
        stepPPM = shortHdr.WIPdbl(5);
        ppmlist = startPPM:stepPPM:stopPPM;
    end
end

output.fid = fid;
output.satDur = satDur;
output.B1 = B1;
output.ppmlist = ppmlist;
output.bw = bw;
output.f0 = f0;

% time axis
t = (0:npts-1)*1/bw;

% freq axes
hz = (-1/2:1/npts:1/2-1/npts)*bw;
ppm = hz/f0; % centered on 0 ppm

% SVD over scans - frequency domain
if isnumeric(pars.denoiseSeries)
    spec = fftshift(fft(fid,[],1),1);
    [U,S,V] = svds(spec,pars.denoiseSeries,'largest','MaxIterations',250);
    spec = U*S*V';
    fid = ifft(ifftshift(spec,1),[],1);
elseif pars.denoiseSeries
    spec = fftshift(fft(fid,[],1),1);
    f = NWsvdTS(spec,ppm,'SVD denoising: save and close when done');
    waitfor(f);
    spec = out; clear out
    fid = ifft(ifftshift(spec,1),[],1);
end

% get signal
switch pars.fitmode
    case 0 % 1st time point
        sig = fid(1,:);

    case 1 % max spectrum
        spec = fftshift(fft(fid,[],1),1);
        inds = abs(ppm-pars.peak) < 0.2; % peak should be within this range
        sig = max(abs(spec(inds,:)),[],1);

    case 1.5 % findpeaks
        spec = abs(fftshift(fft(fid,[],1),1));
        inds = abs(ppm) < 1; % restricted range from 3.7-5.7 ppm
        sig = zeros(1,size(spec,2)); width = 0*sig;
        for jj=1:size(spec,2)
            tempspec = spec(:,jj);
            [~,loc,width,prom] = findpeaks(tempspec(inds),SortStr="descend",NPeaks=1);
            sig(jj) = width * prom;
            width(jj) = width;
        end

    case 2 % hsvd
        rangeHz = f0 * [-0.15 0.15]; % +/- 0.15 ppm range
        do_svds = false;
        sig = zeros(1,size(fid,2)); width = 0*sig;
        for jj=1:size(fid,2)
            tempfid = fid(:,jj);
            [freqHz, damp, amp, y_comp, y_fit,lambda] = hsvd_fit(t,tempfid,30,[],rangeHz,do_svds);
            sig(jj) = abs(y_fit(1));
            width(jj) = abs(1/damp);
        end

    case 2.1 % hsvd - 1 peak only
        rangeHz = f0 * [-0.15 0.15]; % +/- 0.15 ppm range
        do_svds = true;
        sig = zeros(1,size(fid,2)); width = 0*sig;
        for jj=1:size(fid,2)
            tempfid = fid(:,jj);
            [freqHz, damp, amp, y_comp, y_fit,lambda] = hsvd_fit(t,tempfid,1,[],rangeHz,do_svds);
            sig(jj) = abs(y_fit(1));
            width(jj) = abs(1/damp);
        end
        
    case 3 % 1 peak Lorentzian - frequency domain 
        warning('not yet implemented')

end
output.sig = sig;

% +/- only
if length(ppmlist)==1
    if ppmlist>0
        pos = sig(1:2:end);
        neg = sig(2:2:end);
    else
        pos = sig(2:2:end);
        neg = sig(1:2:end);
    end
else
    warning('not yet implemented')
    return;
end

mtr = 100 * (neg - pos) ./ neg;
if pars.dummyShots>length(mtr)
    warning('too many dummy shots')
elseif pars.dummyShots>0
    mtr = mtr(1+pars.dummyShots:end);
end

output.pos = pos;
output.neg = neg;
output.CrCEST = mtr;
output.pars = pars;
output.ppm = ppm;
output.time = t;
if exist('width','var')
    output.width = width;
end

rmpath([thisPath filesep 'utils']);

        
