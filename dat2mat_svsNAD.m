function [output,pars,fullHdr] = dat2mat_svsNAD(pars,ws_dat,nws_dat)
%
% [output,pars,fullHdr] = dat2mat_svsNAD(pars,ws_dat,nws_dat)
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



% in case current directory is not this one
thisFile = mfilename('fullpath');
thisPath = fileparts(thisFile);
addpath([thisPath filesep 'utils']);
addpath([thisPath filesep 'mapVBVD']);

% defaults
if ~isfield(pars,'plt') || isempty(pars.plt), pars.plt = false; end
if ~isfield(pars,'removeOS') || isempty(pars.removeOS), pars.removeOS = false; end
if ~isfield(pars,'lb') || isempty(pars.lb), pars.lb = 0; end
if ~isfield(pars,'den') || isempty(pars.den), pars.den = false; end
if ~isfield(pars,'eccopt') || isempty(pars.eccopt), pars.eccopt = -1; end % uses eddy current correction to frequency correct signal from large voxels
if ~isfield(pars,'doWatSupPost') || isempty(pars.doWatSupPost), pars.doWatSupPost = false; end
if ~isfield(pars,'doWatSupPre') || isempty(pars.doWatSupPre), pars.doWatSupPre = false; end
if ~isfield(pars,'peaks'), pars.peaks = []; end
if (~isfield(pars,'ccopts') || ~isfield(pars.ccopts,'minsig_frac')) || isempty(pars.ccopt.minsig_frac), pars.ccopt.minsig_frac = 0.05; end 
if ~isfield(pars,'block_size'), pars.block_size = []; end
if ~isfield(pars,'base'), pars.base = nan; end
if ~isfield(pars,'PC'), pars.PC = 0; end
if ~isfield(pars,'dofit'), pars.dofit = false; end
if pars.dofit
    if ~isfield(pars.fit,'mode'), pars.fit.mode = 3; end % complex Lorentzian
end

if pars.doWatSupPost
    pars.doWatSupPre = false;
end

if pars.doWatSupPre
    pars.doWatSupPost = false;
end

if nargin<2
    [file,path] = uigetfile('*.dat','Choose NWS dat file');
    nws_dat = fullfile(path,file);
    [file,path] = uigetfile('*.dat','Choose WS dat file');
    ws_dat = fullfile(path,file);
elseif nargin<3
    nws_dat = [];
end


%% non water suppressed
if isempty(nws_dat)
    pars.eccopt = -1; % can't do eddy current correction without water reference data
    nws = [];

elseif ~isempty(nws_dat)

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

nch = nws_obj.image.NCha;
nave = nws_obj.image.NAve;
npts = nws_obj.image.NCol;
reps = nws_obj.image.NRep;
avgDim = find(strcmp(nws_obj.image.sqzDims,'Ave'));
repDim = find(strcmp(nws_obj.image.sqzDims,'Rep'));

wiplong = nws_obj.hdr.MeasYaps.sWipMemBlock.alFree;
empInd = cellfun('isempty',wiplong);
wiplong(empInd) = {0};
wiplong = cell2mat(wiplong);
wipdbl = nws_obj.hdr.MeasYaps.sWipMemBlock.adFree;
empInd = cellfun('isempty',wipdbl);
wipdbl(empInd) = {0};
wipdbl = cell2mat(wipdbl);

% read in data
nws = nws_obj.image{''};

% averages
if ~isempty(avgDim)
    nws = squeeze(mean(nws,avgDim));
end

if ndims(nws)>2 % take first point if there are extra, unexpected dimensions besides read-channels
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
elseif contains(seqname,'svssel_latest') % check this for different versions of svssel_latest
    delay_ind = 2 * wiplong(8);
else
    delay_ind = 0;
end
nws = nws(delay_ind+1:end,:);
nws = nws(1:vectorSize,:);

% remove OS
if pars.removeOS
    keepOS = [1:vectorSize/4, 1+vectorSize*3/4:vectorSize];
    nws = ifft(nws,[],1);
    nws = fft(nws(keepOS,:),[],1);
    bw = bw/2;
    npts = vectorSize/2;
else
    npts = vectorSize;
end


end

%% water suppressed
ws_obj = mapVBVD(ws_dat);
if iscell(ws_obj)
    ws_obj = ws_obj{end};
end
ws_obj.image.flagRemoveOS = false; % oversampling
ws_obj.image.flagAverageSets = false; % cmrr sequence uses 'sets' variable for avgs

fullHdr = ws_obj.hdr;

if exist('f0','var')
    if ~isequal(f0,1e-6*ws_obj.hdr.Config.Frequency) || ~isequal(nch,ws_obj.image.NCha) || ~isequal(nws_obj.image.NCol,ws_obj.image.NCol)
        error('inconsistent protocols')
    end
end
nave = ws_obj.image.NAve * ws_obj.image.NSet; % cmrr uses 'sets' for 'avgs'

seqname = ws_obj.hdr.MeasYaps.tSequenceFileName;
examMemoryUID = ws_obj.hdr.Config.ExamMemoryUID;
fparts = split(examMemoryUID,'_');
for ii=1:length(fparts)
    if length(fparts{ii})==8 & strcmp(fparts{ii}(1:3),'202')
        scandate = datetime(fparts{ii},'InputFormat','yyyyMMdd');
    end
end
bw = 1e9/ws_obj.hdr.Config.DwellTimeSig;
f0 = 1e-6*ws_obj.hdr.Config.Frequency;
nch = ws_obj.image.NCha;
nave = ws_obj.image.NAve;
npts = ws_obj.image.NCol;
reps = ws_obj.image.NRep;
avgDim = find(strcmp(ws_obj.image.sqzDims,'Ave'));
repDim = find(strcmp(ws_obj.image.sqzDims,'Rep'));
vectorSize = 2*ws_obj.hdr.MeasYaps.sSpecPara.lVectorSize; % 2 for oversampling

wiplong = ws_obj.hdr.MeasYaps.sWipMemBlock.alFree;
empInd = cellfun('isempty',wiplong);
wiplong(empInd) = {0};
wiplong = cell2mat(wiplong);
wipdbl = ws_obj.hdr.MeasYaps.sWipMemBlock.adFree;
empInd = cellfun('isempty',wipdbl);
wipdbl(empInd) = {0};
wipdbl = cell2mat(wipdbl);

% timing correction - needed for eja_slaser_svs
if contains(seqname,'eja_svs_slaser')
    delay_ind = ws_obj.noise.iceParam(5); % 230
elseif contains(seqname,'svssel_latest') % check this for different versions of svssel_latest
    delay_ind = 2 * wiplong(8);
else
    delay_ind = 0;
end


ws = ws_obj.image{''};
ws = ws(delay_ind+1:end,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
ws = ws(1:vectorSize,:,:,:,:,:,:,:,:,:,:,:,:,:,:);

% remove OS
if pars.removeOS
    ws = ifft(ws,[],1);
    ws = fft(ws(keepOS,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:),[],1);
end

% time axis
t = (0:vectorSize-1)*1/bw;

% freq axis
hz = (-1/2:1/vectorSize:1/2-1/vectorSize)*bw;
if contains(seqname,'svssel_latest') && ~contains(seqname,'svssel_latest0')  % older versions up to svssel_latest021224 do not include this
    if scandate < datetime('20240212','InputFormat','yyyyMMdd')
        center_freq = 4.72;
    else
        center_freq = wipdbl(4);
    end
else
    center_freq = 4.72;
end
ppm = center_freq + hz/f0;

disp('eddy current correction')

if isempty(nws)
    pars.eccopt = -1;
    ref_temp = squeeze(mean(ws,avgDim));
    [ws, ~, phcorr] = eddyCurrentCorrection(ref_temp,ws,pars.eccopt);
    % phcorr = squeeze(phcorr(1,:,1));
else
    [ws, nws, phcorr] = eddyCurrentCorrection(nws,ws,pars.eccopt);
    % phcorr = squeeze(phcorr(1,:,1));
end

disp('filtering')
[ws,filt] = expFilter(t,pars.lb,ws);
if ~isempty(nws)
    nws = expFilter(t,pars.lb,nws);
end

if pars.doWatSupPre
    wsopts = pars.wsopts;

    disp('water signal removal')
    if ~strcmp(wsopts.type,'filt')
        disp('    changing to average-based water signal removal')
        wsopts.type = 'filt';
        wsopts.filt.type = 'average';
        wsopts.filt.N = 11;
        wsopts.plt = true;
    end
    
    if center_freq==4.72
        ws = removeResidualWater(double(ws),wsopts);
    else
        freqShiftHz = f0*(center_freq - 4.72);
        ws = removeMetSignal(double(ws),freqShiftHz,t,wsopts);
    end
end

% should i squeeze here first?
if strcmp(pars.den,'mlsvdPre')
    disp('MLSVD of the fid')
    addpath('..\tensorlab_2016-03-28');
    if isfield(pars.den,'MLSingularValues')
        [Ufull,Sfull,Vfull] = mlsvd(ws);
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

disp('block averaging')
% block_size = min(64,size(ws,2));
if isempty(pars.block_size)
    pars.block_size = size(ws,3);
end
nd = ndims(ws);
if nd<4
    ws = permute(ws,[1 3 2]);
else
    ws = permute(ws,[1 3 2 4:nd]);
end
ws = blockAverage(ws,pars.block_size); % average dimension MUST be 2nd in ws here
if nd<4
    ws = permute(ws,[1 3 2]);
else
    ws = permute(ws,[1 3 2 4:nd]);
end


disp('coil combination')
if isempty(nws)
    [weights,ws] = coilCombinationNoPC(ws,pars.ccopt.minsig_frac,ws);
else
    [weights, ws, nws] = coilCombinationNoPC(nws,pars.ccopt.minsig_frac,ws);
end
if pars.plt
    % figure, plot(weights,'*'), title('coil weights')
end


wsopts.type = 'none'; % not needed for NAD spectroscopy
if ~strcmp(wsopts.type,'none')
    disp('residual water removal')
end
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
if pars.plt && ~strcmp(wsopts.type,'none')
    plot(ppm,real(fftshift(fft(sum(ws(:,:),2),[],1),1)),'r')
end

if strcmp(pars.den,'hsvd')
    disp('HSVD of the fid')
    f = NWsvdden(ws,t,bw);
    waitfor(f);
    if exist('out','var')
        ws = out; clear out
        output.hdr.flag.hsvdden = true;
    end
end

% Correct for zero order phase discrepancy between nws and ws - sequence
% imperfections? Should ideally be redundant.
phcorr = exp(-1i*squeeze(angle(ws(1,:,:,:))));
ws = ws .* repmat(phcorr,[size(ws,1) 1 1 1 ]);

% resolve remaining averaging
ws = squeeze(mean(ws,2));

% output
output.met.fid = ws;
output.wat.fid = nws;

disp('Fourier transformation')
ws = fftshift(fft(ws,[],1),1);
nws = fftshift(fft(nws,[],1),1);

%%%%%%%%%%% Spectrum %%%%%%%%%%%%%%%%%%%%%

% already taken care of by ECC
% [nws, hzshift] = bulkSpectrumShift(nws,hz);
% hzshift = 0; 
% disp(['water frequency correction: ' num2str(hzshift) ' hz'])
% ws = shiftSpectrumFrequency(ws,hzshift,t);

% if plt && ~isequal(hzshift,0*hzshift)
%     figure, plot(ppm,real(sum(ws(:,:),2)),'b'); title('frequency/phase correction')
%     hold on
% end


% phase correction
if ~isfield(pars,'PC') || isempty(pars.PC)
    f = NWman_phase(ws,ppm,'Manual phase correction: save and close when done');
    waitfor(f);
    if exist('out','var')
        ws = out; clear out
        output.short_hdr.flag.pc = true;
        pars.PC = parsPC;
    end
elseif isnumeric(pars.PC)
    if pars.PC~=0
        disp('individual spectra freq/ph correction')
        if ~isfield(pars,'PCpivot')
            pars.PCpivot = 4.72;
        end
        ws = shiftSpectrumPhase(ws,pars.PC,ppm,pars.PCpivot);
    end
end

% removed for now
if false % ~isempty(pars.peaks)
    peaks = pars.peaks;

    ranges = 0.1*ones(length(peaks),2);
    addlb = 2; % additional line broadening to locate peaks
    hzshift = 0;
    % [ws,hzshift] = peakSpectrumShift(ws,hz,f0,addlb,peaks,ranges);
    % for ii=1:length(hzshift)
    %     disp(['      freq shift: ' num2str(hzshift(ii)) ' hz'])
    % end
    % if plt && ~isequal(hzshift,0*hzshift)
    %     plot(ppm,real(sum(ws(:,:),2)),'k')
    % end
    ws = peakSpectrumPhase(ws,ppm,t,addlb,peaks,ranges/2);
    if pars.plt && ~isequal(hzshift,0*hzshift)
        plot(ppm,real(sum(ws(:,:),2)),'r'), legend({'orig','freq corr','phase corr'})
    %     figure, plot(hzshift,'-*')
    end
end

% Additional denoising
switch pars.den
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


if isfield(pars,'base') && strcmpi(pars.base,'full')
    disp('baseline correction')
    f = NWsemiman_base(ws,ppm,'Semi manual baseline correction: save and close when done');
    waitfor(f);
    if exist('out','var')
        output.met.baseline = ws - out.spec;
        ws = out.spec;
        pars.flag.baseCorr = true;
        pars.base = out.pars;
        clear out
    end
end

% baseopts.algo = 'none';
% baseopts.method = 'pchip';
% baseopts.stepsize = 0.5;
% baseopts.windowsize = 0.5;
% if plt && ~strcmp(baseopts.algo,'none')
%     figure, plot(ppm,real(ws),'b'), title('baseline correction')
%     hold on
% end
% ws = baselineCorrect(ws,ppm,baseopts); % real part only
% if plt && ~strcmp(baseopts.algo,'none')
%     plot(ppm,real(ws),'r')
% end

% fitting
pars.waterScaling = false;
if pars.dofit
    if strcmpi(pars.dofit,'lcmodel')

        warning('no yet implemented. skipping...')

    elseif strcmpi(pars.dofit,'auto')

        % scale to water first
        pars.waterScaling = true;
        nws = nws / max(abs(nws(:)));
        ws = ws / max(abs(nws(:)));

        if isfield(pars.fit,'ppm_range')
            inds = find( ppm > min(pars.fit.ppm_range) & ppm < max(pars.fit.ppm_range));
            x = ppm(inds);
            y = double(real(ws(inds)));
        else
            x = ppm;
            y = double(real(ws));
        end

        minw = 10/f0;
        [yfit,n,names,ampl,pos,width,integral,ip,pars] = curvefitAuto(x,y,minw,pars.fit.mode);


    elseif strcmpi(pars.dofit,'man')

        % scale to water first
        pars.waterScaling = true;
        nws = nws / max(abs(nws(:)));
        ws = ws / max(abs(nws(:))); 

        if isfield(pars.fit,'ppm_range')
            inds = find( ppm > min(pars.fit.ppm_range) & ppm < max(pars.fit.ppm_range));
            x = ppm(inds);
            y = double(real(ws(inds)));
        else
            x = ppm;
            y = double(real(ws));
        end

        if isfield(pars,'base') && strcmpi(pars.base,'fit') % baseline correction on fitted range only
            f = NWsemiman_base(y,x,'Semi manual baseline correction: save and close when done');
            waitfor(f);
            if exist('out','var')
                output.met.fit.baseline = y - out.spec;
                y = out.spec;
                pars.flag.baseCorr = true;
                pars.base = out.pars;
                clear out
            end
        end

        minw = 10/f0; % 10 Hz minimum linewidth
        [yfit,names,areas,ip,ip0,lb,ub] = curvefitMan(x,double(real(y)),minw,pars.fit.mode,pars.fit.peaks);

        if pars.plt
            figure, plot(x,real(y),'k',x,yfit,'r',x,real(y)-yfit,'g'), legend({'data','fit','residual'}), set(gca,'xdir','reverse')
        end

        % fit output
        output.met.fit.spec = y;
        output.met.fit.spec_fit = yfit;
        output.met.fit.pars = ip;
        output.met.fit.names = names;
        output.met.fit.areas = areas;
        output.met.fit.ppm = x;

        % water fit - automatic
        if ~isempty(nws)
            % inds = find( ppm > 3.7 & ppm < 5.7 );
            % x = ppm(inds);
            % y = double(real(nws(inds)));
            % [yfit,~,areas,ip] = curvefitMan(x,y,minw,pars.fit.mode,{'wat'});
            x = ppm;
            y = double(real(nws));
            [yfit,areas,ip] = curvefitWat(x,y,minw,pars.fit.mode);
            output.wat.fit.spec = y;
            output.wat.fit.ppm = x;
            output.wat.fit.spec_fit = yfit;
            output.wat.fit.pars = ip;
            output.wat.fit.areas = areas;
        end

    else
        warning('unknown fitting method. skipping...')
    end
end



% output
output.met.spec = ws;
output.time = t;
output.hz = hz;
output.ppm = ppm;
output.wat.spec = nws;


rmpath([thisPath filesep 'mapVBVD']);
rmpath([thisPath filesep 'utils']);
