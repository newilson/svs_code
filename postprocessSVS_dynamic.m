function output = postprocessSVS_dynamic(pathname,pars)
%
% Dynamic spectroscopy acquisition processing
%
% Inputs:
%   pathname - path to directory containing DICOM spectroscopy files
%   pars     - structure with processing parameters:
%       .fitmode      - fitting method (default 2, see options below)
%       .peaklistPPM  - expected peak positions in ppm (default 0)
%       .dummyShots   - number of dummy acquisitions to skip (default 0)
%       .denoiseSeries- SVD denoising: false, true (interactive), or rank (numeric)
%       .lb           - line broadening in Hz (optional)
%       .nComp        - number of HSVD components (for fitmode 10+, default 30)
%
%   For fitmode 1-7 (curvefitAuto):
%       .ppmRange     - [min max] fitting range in ppm (default: full range)
%       .amp0         - initial amplitudes (auto-detect if empty)
%       .width0       - initial linewidths in ppm (auto-detect if empty)
%       .ph_range     - [min max] phase bounds in degrees (default [-180 180])
%       .minw         - minimum linewidth (default 0)
%
% Fitmode options:
%   0: 1st time point, 0.1: max spectrum, 0.2: findpeaks
%   1: Gaussian, 2: Lorentzian, 3: complex Lorentzian
%   4: complex Gaussian, 5: complex Voigt
%   6: magnitude Lorentzian, 7: magnitude Voigt
%   10: HSVD (full range), 10.1: HSVD (1 peak), 10.2: HSVD (multi-peak)
%
% Output:
%   output.sig        - extracted signal time series
%   output.fid        - raw FID data
%   output.ppm        - ppm axis
%   output.time       - time axis
%   output.pars       - parameters used
%   output.pathname   - input pathname
%   output.bw         - bandwidth
%   output.f0         - Larmor frequency
%   (for HSVD modes): .freqHz, .damp, .amp, .fit
%   (for curvefitAuto modes): .fit, .fitResults, .ppm_fit
%
% NW - 01/2026

% check for parallel pool
% doParFor = false;
% try
%     checkParTool = canUseParallelPool;
% catch
%     checkParTool = false;
% end
% 
% if checkParTool
%     if ~isempty(gcp('nocreate')) || Z>1
%         doParFor = true;
%     end
% end

DEBUG = true;

% in case current directory is not this one
thisFile = mfilename('fullpath');
thisPath = fileparts(thisFile);
addpath([thisPath filesep 'utils']);

if nargin<2 || ~isfield(pars,'fitmode')
    pars.fitmode = 2;
end

if nargin<2 || ~isfield(pars,'peaklistPPM')
    pars.peaklistPPM = 0; % 0 ppm
end

if nargin<2 || ~isfield(pars,'dummyShots')
    pars.dummyShots = 0;
end

if nargin<2 || ~isfield(pars,'denoiseSeries')
    pars.denoiseSeries = false;
end

if pars.fitmode>=10 && pars.fitmode<11 % HSVD fitting
    if nargin<2 || ~isfield(pars,'nComp')
        pars.nComp = 30;
    end
end

% curvefitAuto parameters (fitmode 1-7)
if pars.fitmode>=1 && pars.fitmode<10
    if nargin<2 || ~isfield(pars,'ppmRange')
        pars.ppmRange = []; % default fitting range in ppm
    end
    if nargin<2 || ~isfield(pars,'amp0')
        pars.amp0 = []; % auto-detect if empty
    end
    if nargin<2 || ~isfield(pars,'width0')
        pars.width0 = []; % auto-detect if empty
    end
    if nargin<2 || ~isfield(pars,'ph_range')
        pars.ph_range = [-180 180]; % phase bounds for complex modes
    end
    if nargin<2 || ~isfield(pars,'minw')
        pars.minw = 0; % minimum linewidth
    end
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
            if ~isequal(npts,shortHdr.npts) || ~isequal(bw,shortHdr.bw) || ~isequal(f0,shortHdr.f0)
                error('dicoms do not appear to be from the same scan')
            end
        end
        npts = shortHdr.npts;
    	bw = shortHdr.bw;
        nucleus = shortHdr.Nucleus;
    	f0 = shortHdr.f0;
    end
end

output.fid = fid;
output.bw = bw;
output.f0 = f0;

% time axis
t = (0:npts-1)*1/bw;

% freq axes
hz = (-1/2:1/npts:1/2-1/npts)*bw;
if strcmpi(nucleus,'1H')
    if contains(shortHdr.SeqName,'svssel_latest') && ~contains(shortHdr.SeqName,'svssel_latest0') % older versions up to svssel_latest021224 do not include this
        center_freq = shortHdr.WIPdbl(4);
    else
        center_freq = 4.72;
    end
    ppm = center_freq - hz/f0;
else
    ppm = hz/f0; % centered on 0 ppm
end

% filtering
if isfield(pars,'lb')
    if DEBUG && pars.lb>0
        figure, plot(t,abs(fid)/abs(fid(1)),'b'), hold on
    end
    [fid, filt] = expFilter(t,pars.lb,fid);
    if DEBUG && pars.lb>0
        plot(t,filt,'r')
    end
end

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

if DEBUG
    f = figure;
end

% get signal
switch pars.fitmode
    case 0 % 1st time point
        sig = fid(1,:);

    case 0.1 % max spectrum
        spec = fftshift(fft(fid,[],1),1);
        inds = abs(ppm-pars.peaklistPPM) < 0.2; % peak should be within this range
        sig = max(abs(spec(inds,:)),[],1);

    case 0.2 % findpeaks
        spec = abs(fftshift(fft(fid,[],1),1));
        inds = abs(ppm) < 1; % restricted range from 3.7-5.7 ppm
        sig = zeros(1,size(spec,2));
        for jj=1:size(spec,2)
            tempspec = spec(:,jj);
            [~,loc,width,prom] = findpeaks(tempspec(inds),SortStr="descend",NPeaks=1);
            sig(jj) = width * prom;
        end

    case 10 % hsvd
        rangeHz = []; % full range
        do_svds = false;
        sig = zeros(1,size(fid,2));
        freqHz = 0*sig; damp = 0*sig; amp = 0*sig; fit = 0*fid;
        for jj=1:size(fid,2)
            tempfid = fid(:,jj);
            [temp_freqHz, temp_damp, temp_amp, y_comp, temp_y_fit,lambda] = hsvd_fit(t,tempfid,pars.nComp,[],rangeHz,do_svds);
            sig(jj) = abs(temp_y_fit(1));
            if jj==1
                freqHz = zeros(length(temp_freqHz),size(fid,2));
                damp = 0*freqHz; amp = 0*freqHz; 
            end
            freqHz(:,jj) = temp_freqHz;
            damp(:,jj) = temp_damp;
            amp(:,jj) = temp_amp;
            fit(:,jj) = temp_y_fit;
            if DEBUG
                disp(['Using ' num2str(length(temp_freqHz)) ' components'])
                debug_plot_time(f,ppm,tempfid,temp_y_fit);
                pause(1)
            end
        end

    case 10.1 % hsvd - 1 peak only centered at 0 - works for PCr or H20
        if strcmpi(nucleus,'31P')
            rangeHz = f0 * [-1 1]; % +/- 1 ppm - suitable for PCr
        else
            rangeHz = f0 * [-0.15 0.15]; % +/- 0.15 ppm range - suitable for H20
        end
        do_svds = true;
        sig = zeros(1,size(fid,2));
        freqHz = 0*sig; damp = 0*sig; amp = 0*sig; fit = 0*fid;
        for jj=1:size(fid,2)
            tempfid = fid(:,jj);
            [temp_freqHz, temp_damp, temp_amp, y_comp, temp_y_fit,lambda] = hsvd_fit(t,tempfid,1,[],rangeHz,do_svds);
            sig(jj) = abs(temp_y_fit(1));
            freqHz(jj) = temp_freqHz;
            damp(jj) = temp_damp;
            amp(jj) = temp_amp;
            fit(:,jj) = temp_y_fit;
            if DEBUG
                disp(['Using ' num2str(length(temp_freqHz)) ' component'])
                debug_plot_time(f,ppm,tempfid,temp_y_fit);
                pause(1)
            end
        end
        
    case 10.2 % hsvd - 1 peak only centered at 0 - works for PCr or H20. This one allows for multiple components.
        if strcmpi(nucleus,'31P')
            rangeHz = f0 * [-0.5 0.5]; % +/- 0.5 ppm - suitable for PCr
        else
            rangeHz = f0 * [-0.15 0.15]; % +/- 0.15 ppm range - suitable for H20
        end
        do_svds = false;
        sig = zeros(1,size(fid,2));
        freqHz = 0*sig; damp = 0*sig; amp = 0*sig; fit = 0*fid;
        for jj=1:size(fid,2)
            tempfid = fid(:,jj);
            [temp_freqHz, temp_damp, temp_amp, y_comp, temp_y_fit,lambda] = hsvd_fit(t,tempfid,pars.nComp,[],rangeHz,do_svds);
            sig(jj) = abs(temp_y_fit(1));
            if jj==1
                freqHz = zeros(length(temp_freqHz),size(fid,2));
                damp = 0*freqHz; amp = 0*freqHz; 
            end
            freqHz(:,jj) = temp_freqHz;
            damp(:,jj) = temp_damp;
            amp(:,jj) = temp_amp;
            fit(:,jj) = temp_y_fit;
            if DEBUG
                title(['Frame ' num2str(jj) ' of ' num2str(size(fid,2))]);
                disp(['Using ' num2str(length(temp_freqHz)) ' components'])
                debug_plot_time(f,ppm,tempfid,temp_y_fit);
                pause(1)
            end
        end

    otherwise  % frequency domain fitting (fitmode 1-7, see curvefitAuto.m)
        % fitmode 1: Gaussian, 2: Lorentzian, 3: complex Lorentzian,
        % 4: complex Gaussian, 5: complex Voigt, 6: mag Lorentzian, 7: mag Voigt

        if pars.fitmode < 1 || pars.fitmode > 7
            error('Invalid fitmode. Use 0-0.2 for simple, 1-7 for curvefitAuto, 10-10.2 for HSVD');
        end

        % compute spectrum
        spec = fftshift(fft(fid,[],1),1);

        % restrict to ppm range
        if ~isempty(pars.ppmRange)
            ppmInds = ppm >= min(pars.ppmRange) & ppm <= max(pars.ppmRange);
            ppm_fit = ppm(ppmInds);
        else
            ppmInds = 1:length(ppm);
            ppm_fit = ppm;
        end

        % determine if we need real or magnitude spectrum based on mode
        % Note: "complex" modes (3,4,5) fit real data with a phase parameter
        if pars.fitmode == 6 || pars.fitmode == 7
            useSpec = abs(spec(ppmInds,:));  % magnitude modes
        else
            useSpec = real(spec(ppmInds,:)); % all other modes fit real part
        end

        % get initial guesses from first spectrum if not provided
        nPeaks = length(pars.peaklistPPM);
        if isempty(pars.amp0) || isempty(pars.width0)
            % find peaks in spectrum
            [~, detectedCenter, detectedWidth, detectedAmp] = findpeaks(abs(spec(ppmInds,1)), ppm_fit);

            if isempty(detectedAmp)
                % no peaks found, use defaults at expected positions
                tempAmp = repmat(max(abs(spec(ppmInds,1))), 1, nPeaks);
                tempWidth = repmat(0.05, 1, nPeaks); % default width in ppm
            else
                % match detected peaks to expected positions in peaklistPPM
                tempAmp = zeros(1, nPeaks);
                tempWidth = zeros(1, nPeaks);
                for kk = 1:nPeaks
                    [~, closestIdx] = min(abs(detectedCenter - pars.peaklistPPM(kk)));
                    tempAmp(kk) = detectedAmp(closestIdx);
                    tempWidth(kk) = detectedWidth(closestIdx);
                end
            end
            if isempty(pars.amp0), pars.amp0 = tempAmp; end
            if isempty(pars.width0), pars.width0 = tempWidth; end
        end

        % initialize tracking variables
        amp_current = pars.amp0;
        center_current = pars.peaklistPPM(:)';
        width_current = pars.width0;
        sig = zeros(1, size(fid,2));
        fit = zeros(length(ppm_fit), size(fid,2));
        fitResults = cell(1, size(fid,2));

        for jj = 1:size(fid,2)
            tempSpec = useSpec(:,jj);

            % call curvefitAuto
            fitOut = curvefitAuto(ppm_fit(:), tempSpec(:), pars.fitmode, ...
                amp_current, center_current, width_current, pars.ph_range, pars.minw);

            fitResults{jj} = fitOut;
            fit(:,jj) = fitOut.fit;

            % extract signal - use sum of areas or first peak amplitude
            if nPeaks == 1
                sig(jj) = fitOut.areas(1);
            else
                sig(jj) = sum(fitOut.areas);
            end

            % update tracking variables for next iteration
            amp_current = fitOut.pars(1:nPeaks);
            center_current = fitOut.pars(nPeaks + (1:nPeaks));
            width_current = fitOut.pars(2*nPeaks + (1:nPeaks));

            if DEBUG
                debug_plot_freq(f, ppm_fit, tempSpec, fitOut.fit);
                title(['Frame ' num2str(jj) ' of ' num2str(size(fid,2))]);
                pause(0.5)
            end
        end

        output.fit = fit;
        output.fitResults = fitResults;
        output.ppm_fit = ppm_fit;

end
output.sig = sig;

output.pars = pars;
output.ppm = ppm;
output.time = t;

rmpath([thisPath filesep 'utils']);

end

function debug_plot_time(figId,x,FTy1,FTy2)
% debug plot for time domain fitting

y1 = fftshift(fft(FTy1));
y2 = fftshift(fft(FTy2));

figure(figId)
subplot(2,1,1), plot(x,abs(y1),'k',x,abs(y2),'r',x,abs(y1)-abs(y2),'g');
subplot(2,1,2), plot(x,real(y1),'k',x,real(y2),'r',x,real(y1)-real(y2),'g');
legend('data','fit','residual');
xlabel('ppm');

end

function debug_plot_freq(figId,x,y1,y2)
% debug plot for frequency domain fitting

figure(figId)
plot(x,y1,'k',x,y2,'r',x,y1-y2,'g');
legend('data','fit','residual');
xlabel('ppm');

end
