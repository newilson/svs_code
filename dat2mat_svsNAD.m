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

% Warning only -- unknown fields are never fatal, and the schema below WILL fall behind as options are added.
checkParsFields(pars);

% defaults
if ~isfield(pars,'plt') || isempty(pars.plt), pars.plt = false; end % debugging/fitting plots
if ~isfield(pars,'pltSpec') || isempty(pars.pltSpec), pars.pltSpec = true; end % plot final spectrum before fitting
if ~isfield(pars,'removeOS') || isempty(pars.removeOS), pars.removeOS = false; end
if ~isfield(pars,'lb') || isempty(pars.lb), pars.lb = 0; end
if ~isfield(pars,'lb_met') || isempty(pars.lb_met), pars.lb_met = pars.lb; end
if ~isfield(pars,'lb_wat') || isempty(pars.lb_wat), pars.lb_wat = pars.lb; end
if ~isfield(pars,'den') || isempty(pars.den), pars.den = false; end
if ~isfield(pars,'eccopt') || isempty(pars.eccopt), pars.eccopt = -1; end % Brown method
if ~isfield(pars,'WatSupPost') || isempty(pars.WatSupPost) || isempty(pars.WatSupPost.opts.type), pars.WatSupPost.opts.type = 'none'; end % water removal AFTER coil combination and averaging
if ~isfield(pars,'WatSupPre') || isempty(pars.WatSupPre) || isempty(pars.WatSupPre.opts.type), pars.WatSupPre.opts.type = 'none'; end % water removal BEFORE coil combination and averaging
if ~isfield(pars,'peaks'), pars.peaks = []; end
if (~isfield(pars,'ccopt') || ~isfield(pars.ccopt,'minsig_frac')) || isempty(pars.ccopt.minsig_frac), pars.ccopt.minsig_frac = 0.05; end 
if ~isfield(pars,'block_align') || ~isfield(pars.block_align,'size'),   pars.block_align.size = 64; end
if ~isfield(pars,'block_align') || ~isfield(pars.block_align,'method'), pars.block_align.method = 'fd'; end % 'td' (fidfreqshift) or 'fd' (fidfreqshift_fd)
if ~isfield(pars.block_align,'Nfit'),     pars.block_align.Nfit = 256; end       % TD: # FID points fit
if ~isfield(pars,'base'), pars.base = nan; end
if ~isfield(pars,'hsvd') || ~isfield(pars.hsvd,'flag'), pars.hsvd.flag = false; end
if ~isfield(pars.hsvd,'sv'), pars.hsvd.sv = []; end % HSVD singular values
if ~isfield(pars.hsvd,'cad'), pars.hsvd.cad = 3; end % HSVD Cadzow iterations
if isfield(pars,'pc'), pars.PC = pars.pc; end
if ~isfield(pars,'PC'), pars.PC = 0; end
if ~isfield(pars,'dofit'), pars.dofit = false; end
if pars.dofit
    if ~isfield(pars.fit,'mode'), pars.fit.mode = 3; end % complex Lorentzian
end

% Can have both?
% if pars.WatSupPost.flag
%     pars.WatSupPre.flag = false;
% end
% 
% if pars.WatSupPre.flag
%     pars.WatSupPost.flag = false;
% end

if nargin<2
    [file,path] = uigetfile('*.dat','Choose NWS dat file');
    nws_dat = fullfile(path,file);
    [file,path] = uigetfile('*.dat','Choose WS dat file');
    ws_dat = fullfile(path,file);
elseif nargin<3
    nws_dat = [];
end

[~,fname,~] = fileparts(ws_dat); % filename 
fname = strrep(fname,'_', ' '); % replace underscores with spaces


%% non water suppressed
if isempty(nws_dat)
    pars.eccopt = -1; % can't do eddy current correction without water reference data
    nws = [];

elseif ~isempty(nws_dat)

nws_obj = mapVBVD(nws_dat);
if iscell(nws_obj)
    nws_obj = nws_obj{end};
end
output.wat.refvolt = nws_obj.hdr.Dicom.flTransRefAmpl;
nws_obj.image.flagRemoveOS = false; % oversampling
nws_obj.image.flagDoAverage = false;
nws_obj.image.flagAverageSets = false; % cmrr sequence uses 'sets' variable for avgs
seqname = nws_obj.hdr.MeasYaps.tSequenceFileName;
examMemoryUID = nws_obj.hdr.Config.ExamMemoryUID;
fparts = split(examMemoryUID,'_');
for ii=1:length(fparts)
    if length(fparts{ii})==8 & strcmp(fparts{ii}(1:3),'202')
        scandate = datetime(fparts{ii},'InputFormat','yyyyMMdd');
    end
end
if isfield(nws_obj.hdr.Config,'DwellTimeSig')
    bw = 1e9/nws_obj.hdr.Config.DwellTimeSig;
elseif isfield(nws_obj.hdr.Config,'DwellTime')
    bw = 1e9/nws_obj.hdr.Config.DwellTime;
end
f0 = 1e-6*nws_obj.hdr.Config.Frequency;
nch = nws_obj.image.NCha;
npts = nws_obj.image.NCol;
vectorSize = 2*nws_obj.hdr.MeasYaps.sSpecPara.lVectorSize; % 2 for oversampling
% if ~isequal(npts,vectorSize)
%     warning(['Actual points is: ' num2str(npts) newline 'Prescribed points is: ' num2str(vectorSize)])
% end

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
if length(nws_obj.image.sqzDims)>1
    nws = nws_obj.image{''};
else
    nws = nws_obj.image();
end

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
elseif contains(seqname,'svssel') % check this for different versions of svssel_latest
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

% time axis
t = (0:vectorSize-1)*1/bw;

% freq axis
hz = (-1/2:1/vectorSize:1/2-1/vectorSize)*bw;
if contains(seqname,'svssel') && ~contains(seqname,'svssel_latest0')  % older versions up to svssel_latest021224 do not include this
    if length(wipdbl)<10 % scandate < datetime('20240212','InputFormat','yyyyMMdd') % this date needs to be confirmed
        center_freq = 4.72;
    else
        center_freq = wipdbl(4);
    end
else
    center_freq = 4.72;
end
ppm = center_freq + hz/f0;

output.wat.ppm = ppm;
output.wat.hz = hz;
output.wat.seqname = seqname;
output.wat.TR_ms = nws_obj.hdr.MeasYaps.alTR{1}/1000; 
output.wat.TE_ms = nws_obj.hdr.MeasYaps.alTE{1}/1000;
output.wat.nAvg = nws_obj.image.NAve * nws_obj.image.NSet;
output.wat.scandate = scandate;


end

%% water suppressed
ws_obj = mapVBVD(ws_dat);
if iscell(ws_obj)
    ws_obj = ws_obj{end};
end
output.met.refvolt = ws_obj.hdr.Dicom.flTransRefAmpl;
if ~isequal(output.met.refvolt, output.wat.refvolt)
    warning('Reference voltage mismatch between WATREF (%.1f V) and MET (%.1f V) scans', output.wat.refvolt, output.met.refvolt);
end
ws_obj.image.flagRemoveOS = false; % oversampling
ws_obj.image.flagAverageSets = false; % cmrr sequence uses 'sets' variable for avgs

fullHdr = ws_obj.hdr;

if exist('f0','var')
    if ~isequal(f0,1e-6*ws_obj.hdr.Config.Frequency) || ~isequal(nch,ws_obj.image.NCha) || ~isequal(nws_obj.image.NCol,ws_obj.image.NCol)
        error('inconsistent protocols')
    end
end
% nave = ws_obj.image.NAve * ws_obj.image.NSet; % cmrr uses 'sets' for 'avgs'

seqname = ws_obj.hdr.MeasYaps.tSequenceFileName;
examMemoryUID = ws_obj.hdr.Config.ExamMemoryUID;
fparts = split(examMemoryUID,'_');
for ii=1:length(fparts)
    if length(fparts{ii})==8 & strcmp(fparts{ii}(1:3),'202')
        scandate = datetime(fparts{ii},'InputFormat','yyyyMMdd');
    end
end
if isfield(ws_obj.hdr.Config,'DwellTimeSig')
    bw = 1e9/ws_obj.hdr.Config.DwellTimeSig;
elseif isfield(ws_obj.hdr.Config,'DwellTime')
    bw = 1e9/ws_obj.hdr.Config.DwellTime;
end
f0 = 1e-6*ws_obj.hdr.Config.Frequency;
nch = ws_obj.image.NCha;
nave = ws_obj.image.NAve;
npts = ws_obj.image.NCol;
reps = ws_obj.image.NRep;
avgDim = find(strcmp(ws_obj.image.sqzDims,'Ave'));
repDim = find(strcmp(ws_obj.image.sqzDims,'Rep'));
vectorSize = 2*ws_obj.hdr.MeasYaps.sSpecPara.lVectorSize; % 2 for oversampling
% if ~isequal(npts,vectorSize)
%     warning(['Actual points is: ' num2str(npts) newline 'Prescribed points is: ' num2str(vectorSize)])
% end

wiplong = ws_obj.hdr.MeasYaps.sWipMemBlock.alFree;
empInd = cellfun('isempty',wiplong);
wiplong(empInd) = {0};
wiplong = cell2mat(wiplong);
wipdbl = ws_obj.hdr.MeasYaps.sWipMemBlock.adFree;
empInd = cellfun('isempty',wipdbl);
wipdbl(empInd) = {0};
wipdbl = cell2mat(wipdbl);

% timing correction - needed for eja_slaser_svs and svssel_latest
if contains(seqname,'eja_svs_slaser')
    delay_ind = ws_obj.noise.iceParam(5); % 230
elseif contains(seqname,'svssel') % check this for different versions of svssel_latest
    delay_ind = 2 * wiplong(8);
else
    delay_ind = 0;
end


ws = ws_obj.image{''};
ws = ws(delay_ind+1:end,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
ws = ws(1:vectorSize,:,:,:,:,:,:,:,:,:,:,:,:,:,:);

% remove OS
if pars.removeOS
    keepOS = [1:vectorSize/4, 1+vectorSize*3/4:vectorSize];
    ws = ifft(ws,[],1);
    ws = fft(ws(keepOS,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:),[],1);
    bw = bw/2;
    npts = vectorSize/2;
else
    npts = vectorSize;
end

% time axis
t = (0:vectorSize-1)*1/bw;

% freq axis
hz = (-1/2:1/vectorSize:1/2-1/vectorSize)*bw;
if contains(seqname,'svssel') && ~contains(seqname,'svssel_latest0')  % older versions up to svssel_latest021224 do not include this
    if length(wipdbl)<10 % scandate < datetime('20240212','InputFormat','yyyyMMdd') % this date needs to be confirmed
        center_freq = 4.72;
        exc_freq = wipdbl(4);
        exc_bw = wipdbl(5);
    else
        center_freq = wipdbl(4);
        exc_freq = wipdbl(5);
        exc_bw = wipdbl(6);
    end
    if ~isfield(pars.block_align,'ppmRange'), pars.block_align.ppmRange = [8.5 10]; end    % FD: [ppm_lo ppm_hi] match window (tuned for 9.7ppm NAD)
else
    center_freq = 4.72;
    exc_freq = 4.72;
    exc_bw = NaN;
    if ~isfield(pars.block_align,'ppmRange'), pars.block_align.ppmRange = [8.5 10]; end    % FD: [ppm_lo ppm_hi] match window
end
ppm = center_freq + hz/f0;

output.met.exc_freq = exc_freq;
output.met.exc_bw = exc_bw;
output.met.ppm = ppm;
output.met.hz = hz;
output.met.seqname = seqname;
output.met.TR_ms = ws_obj.hdr.MeasYaps.alTR{1}/1000; 
output.met.TE_ms = ws_obj.hdr.MeasYaps.alTE{1}/1000; 
output.met.nAvg = ws_obj.image.NAve * ws_obj.image.NSet;
output.met.scandate = scandate;

disp('eddy current correction')

if isempty(nws)
    pars.eccopt = -1;
    ref_temp = squeeze(mean(ws,avgDim));
    [ws, ~, phcorr] = eddyCurrentCorrection(ref_temp,ws,pars.eccopt,nch);
    % phcorr = squeeze(phcorr(1,:,1));
else
    [ws, nws, phcorr] = eddyCurrentCorrection(nws,ws,pars.eccopt,nch);
    % phcorr = squeeze(phcorr(1,:,1));
end

disp('filtering')
ws = expFilter(t,pars.lb_met,ws);
if ~isempty(nws)
    nws = expFilter(t,pars.lb_wat,nws);
end

if ~strcmpi(pars.WatSupPre.opts.type,'none')
    wsopts = pars.WatSupPre.opts;

    disp('water signal removal (pre)')
    if ~strcmp(wsopts.type,'filt') && nch > 1 && nave > 1
        disp('    changing to average-based water signal removal for computational speed ')
        wsopts.type = 'filt';
        wsopts.filt.type = 'average';
        wsopts.filt.N = 11;
        wsopts.plt = true;
    end
    
    if center_freq==4.72
        ws = removeResidualWater(double(ws),wsopts);
    else
        freqShiftHz = f0*(4.72 - center_freq);
        ws = removeMetSignal(double(ws),freqShiftHz,t,wsopts);
    end
end

% should i squeeze here first?
if strcmp(pars.den,'mlsvdPre')
    disp('MLSVD of the fid')
    addpath('..\tensorlab_2016-03-28');
    if isfield(pars.den,'MLSingularValues')
        [Ufull,Sfull,Vfull] = mlsvd(ws);
        sv = min(pars.den.MLSingularValues(:).', size(Sfull));
        Ut = cell(1,numel(Ufull));
        for ii=1:length(Ufull)
            Ut{ii} = Ufull{ii}(:,1:sv(ii));
        end
        idx = arrayfun(@(k) 1:sv(k), 1:numel(sv), 'UniformOutput', false);
        St  = Sfull(idx{:});     % equivalent to Sfull(1:n1, 1:n2, ..., 1:nN)
        ws = lmlragen(Ut,St);
        clear Ut St Ufull Sfull
    else % automatic rank estimation
        [~,~,~,size_core,Ut,St] = mlsvd_wRankEstNW(ws);
        ws = lmlragen(Ut,St);
        pars.den.MLSingularValues = size_core;
        clear Ut St
    end
    rmpath('..\tensorlab_2016-03-28');
end

if nch>1
    disp('coil combination')
    if isempty(nws)
        ref_temp = squeeze(mean(ws,avgDim));
        [weights,ws] = coilCombinationNoPC(ref_temp,pars.ccopt.minsig_frac,ws);
    else
        [weights, ws, nws] = coilCombinationNoPC(nws,pars.ccopt.minsig_frac,ws);
    end
    output.met.cc_weights = weights;
    if pars.plt
        % figure, plot(weights,'*'), title('coil weights')
    end
    if ~isempty(avgDim)
        avgDim = avgDim - 1; % coils are always inside of averages and reps
    end
    if ~isempty(repDim)
        repDim = repDim -1; % coils are always inside of averages and reps
    end
end

if nave>1
    disp('block averaging')
    if isempty(pars.block_align.size)
        pars.block_align.size = size(ws,avgDim); % defaults to 1 block consisting of all the averages
    end
    nd = ndims(ws);
    if avgDim>2 % make avgDim 2nd
        if avgDim==nd
            ws = permute(ws,[1 avgDim 2:avgDim-1]);
        else
            ws = permute(ws,[1 avgDim 2:avgDim-1 avgDim+1:nd]);
        end
    end
    ws = blockAverage(ws,pars.block_align.size); % average dimension MUST be 2nd in ws here
    if avgDim>2 % put back in order
        if avgDim==nd
            ws = permute(ws,[1 avgDim 2:avgDim-1]);
        else
            ws = permute(ws,[1 avgDim 2:avgDim-1 avgDim+1:nd]);
        end
    end
end


if ~strcmpi(pars.WatSupPost.opts.type,'none')
    wsopts = pars.WatSupPost.opts;

    disp('water signal removal (post)')   
    if pars.plt
        temp = real(fftshift(fft(ws(:,:),[],1),1));
    end

    if center_freq==4.72
        ws = removeResidualWater(double(ws),wsopts);
    else
        freqShiftHz = f0*(4.72 - center_freq);
        ws = removeMetSignal(double(ws),freqShiftHz,t,wsopts);
    end

    if pars.plt
        temp2 = real(fftshift(fft(ws(:,:),[],1),1));
        if isvector(ws)
            figure('Name','1H residual water removal'), plot(ppm, temp,'-k',ppm,temp2,'-r'), title('residual water removal'), set(gca,'xdir','reverse')
        else
            NWplayplot(temp2,ppm,'residual water removal',temp, ppm);
        end
    end
end

if pars.block_align.size < nave
    disp('block alignment and resolve averages')
    if isfield(pars.block_align,'saveDebug') && pars.block_align.saveDebug
        output.met.blocks_pre = ws;  % FIDs before alignment, per block
    end
    switch lower(pars.block_align.method)
        case 'fd'
            ws = fidfreqshift_fd(t, ws, ppm, pars.block_align.ppmRange);
        otherwise % 'td'
            ws = fidfreqshift(t, ws, pars.block_align.Nfit, 0);
    end
    if isfield(pars.block_align,'saveDebug') && pars.block_align.saveDebug
        output.met.blocks_post = ws; % FIDs after alignment, per block
    end
    ws = squeeze(mean(ws,2));
end


% Correct for zero order phase discrepancy between nws and ws - sequence
% imperfections? Should ideally be redundant but is not. Must be after averaging.
phcorr = exp(-1i*squeeze(angle(ws(1,:,:,:))));
ws = ws .* repmat(phcorr,[size(ws,1) 1 1 1 ]);


if pars.hsvd.flag
    disp('HSVD denoising...')
    if ~isempty(pars.hsvd.sv)
        for ii=1:size(ws,2)
            ws(:,ii) = NWhsvd(ws(:,ii),pars.hsvd.sv,pars.hsvd.cad);
        end
    else
        f = NWhsvdden(ws,t,bw);
        waitfor(f);
        ws = out; clear out
    end
end

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
    % f = NWman_phase(ws,ppm,'Manual phase correction: save and close when done');
    f = NWman_phase(ws,ppm,['Phase correction: ' fname]);
    waitfor(f);
    if exist('out','var')
        ws = out; clear out
        output.short_hdr.flag.pc = true;
        pars.PC = parsPC;
        output.met.PC = parsPC;
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

if pars.pltSpec
    figure('Name','1H processed spectrum'), plot(ppm,real(ws),'-k'), title(fname), set(gca,'xdir','reverse')
end


% output - pre scaling for fitting
output.met.spec = ws;
output.wat.spec = nws;

% fitting
if pars.dofit
    if strcmpi(pars.dofit,'lcmodel')

        % creates .RAW files for LCModel fitting
        rawfile = regexprep(ws_dat,'datfiles','rawfiles');
        rawfile = regexprep(rawfile,'.dat','.raw');
        fid = ifft(ifftshift(ws,1),[],1);
        create_lcmodelRAW(rawfile,fid,'NUNFIL',length(fid),'DELTAT',1/bw,'HZPPPM',f0);

        if ~isempty(nws)
            % creates .RAW files for LCModel water referencing
            rawfile = regexprep(nws_dat,'datfiles','rawfiles');
            rawfile = regexprep(rawfile,'.dat','.raw');
            fid = ifft(ifftshift(nws,1),[],1);
            create_lcmodelRAW(rawfile,fid,'NUNFIL',length(fid),'DELTAT',1/bw,'HZPPPM',f0);
        end

    elseif strcmpi(pars.dofit,'varpro')

        % VARPRO fitting 
        % pars.peaks should be a struct with:
        %   .name  - cell array of peak names, e.g. {'NAD_H2','NAD_H4'}
        %   .range - [n x 2] ppm ranges for each peak, e.g. [9.0 9.5; 8.6 9.0]
        % pars.fit fields:
        %   .mode      - lineshape mode (default 5 = Voigt)
        %   .ppm_range - [min max] overall fitting range in ppm
        %   .ph_range  - [min max] phase bounds in degrees (default [-45 45])
        %   .baselineOpt - struct passed to curvefitAuto_varproBaseline
        %   .lineShapeOpt - struct passed to curvefitAuto_varproBaseline (optional)
        %   .coarseFirstOrder - fit the 1st-order term of the coarse phase
        %       (default false = zeroth order only; both coarse estimators)
        %   .hsvdClean - struct for HSVD nuisance signal removal before fitting
        %       .enable - true/false
        %       .K      - number of HSVD sinusoids (default 60)
        %       .range  - [min max] ppm range to preserve (default: ppm_range)

        % default expanded peak table (de Graaf et al., MRM 2017)
        % NAD+ H2/H6/H4, NAA amide, and 8 unspecified singlets from
        % adenosine, ATP, macromolecules in the downfield region
        
        if isempty(pars.peaks) || ~isfield(pars.peaks,'range') || ~isfield(pars.peaks,'name')
            pars.peaks.name  = {'NADH2','NADH6','NADH4','NAA_NH', ...
                                'S1','S2','S3','S4','S5','S6','S7','S8'};
            pars.peaks.range = [9.30 9.36; 9.11 9.17; 8.80 8.86; 7.82 7.88; ...
                                7.24 7.30; 7.48 7.54; 7.90 7.96; 8.01 8.07; ...
                                8.13 8.19; 8.24 8.30; 8.36 8.42; 8.47 8.53];
            disp('Using default de Graaf peak table (12 peaks)')
        end

        % defaults
        if ~isfield(pars.fit,'ph_range'), pars.fit.ph_range = [-45 45]; end
        if ~isfield(pars.fit,'baselineOpt'), pars.fit.baselineOpt = struct(); end
        if ~isfield(pars.fit,'ppm_range_coarse')
            if isfield(pars.fit,'ppm_range')
                pars.fit.ppm_range_coarse = pars.fit.ppm_range; % default: same as fine fit range
            end
        end
        % Linewidth bounds may be given in Hz via baselineOpt.widthBoundsHz,
        % which is converted to ppm here using the run's own f0.  Hz is the
        % field-independent way to express a linewidth, and it lets the driver
        % set per-peak caps without knowing f0 (which only exists after the
        % header is read).  Same shape rules as widthBounds: a [min max] row
        % applied to all peaks, or an [npeaks x 2] matrix.  An explicit
        % widthBounds (ppm) still wins if both are supplied.
        if isfield(pars.fit.baselineOpt,'widthBoundsHz') && ...
                ~isempty(pars.fit.baselineOpt.widthBoundsHz) && ...
                ~isfield(pars.fit.baselineOpt,'widthBounds')
            pars.fit.baselineOpt.widthBounds = pars.fit.baselineOpt.widthBoundsHz / f0;
        end
        % default linewidth bounds: 9-36 Hz (de Graaf), applied uniformly to all peaks
        if ~isfield(pars.fit.baselineOpt,'widthBounds')
            pars.fit.baselineOpt.widthBounds = [10/f0 60/f0]; % ppm
        end

        % --- Broad component handling (consistent with the 31P pipeline) -----
        % Peaks whose name starts with 'Broad' are wide macromolecular/baseline
        % components. They get wide linewidth bounds and a broad initial width, and are excluded
        % from the coarse phase fit (they are nuisance signal, not phase
        % references).  Everything else keeps the default narrow bounds.
        npksAll = numel(pars.peaks.name);
        wb0 = pars.fit.baselineOpt.widthBounds;
        if size(wb0,1) == 1
            peakWB = repmat(wb0, npksAll, 1);   % expand global -> per-peak
        else
            peakWB = wb0;
        end
        % per-peak starting Gaussian fraction (Voigt modes), pars.peaks.name order
        peakEta = [];
        if isfield(pars.fit.baselineOpt,'etaInit') && ~isempty(pars.fit.baselineOpt.etaInit)
            peakEta = pars.fit.baselineOpt.etaInit(:).';
            if isscalar(peakEta), peakEta = peakEta * ones(1,npksAll); end
            assert(numel(peakEta) == npksAll, ...
                'baselineOpt.etaInit must be scalar or length %d (pars.peaks.name)', npksAll);
        end

        isBroad = startsWith(pars.peaks.name(:), 'Broad');
        if isfield(pars.fit,'broadWidthBounds') && ~isempty(pars.fit.broadWidthBounds)
            broadWB = pars.fit.broadWidthBounds;
        else
            broadWB = [65/f0 180/f0];   % ~0.22-0.61 ppm: broader than the narrow peaks. Optimized for 7T.
        end
        if isfield(pars.fit,'broadWidthInit') && ~isempty(pars.fit.broadWidthInit)
            broadW0 = pars.fit.broadWidthInit;
        else
            broadW0 = 150/f0;           % ~0.5 ppm init (cf. 31P 100 Hz)
        end

        % broad components are kept near-absorptive
        if isfield(pars.fit,'broadPhaseBounds') && ~isempty(pars.fit.broadPhaseBounds)
            broadPB = pars.fit.broadPhaseBounds;
        else
            broadPB = [-15 15];
        end
        if any(isBroad)
            peakWB(isBroad,:) = repmat(broadWB, sum(isBroad), 1);
            fprintf('Broad component(s): %s  (LW %.3g-%.3g ppm, init %.3g)\n', ...
                strjoin(pars.peaks.name(isBroad), ', '), broadWB(1), broadWB(2), broadW0);
        end
        % focus weighting for Step 2: emphasize a sub-range of the fit
        if ~isfield(pars.fit,'focusRange')
            if isfield(pars.fit,'ppm_range')
                pars.fit.focusRange = pars.fit.ppm_range; % default: original narrow range
            else
                pars.fit.focusRange = [];
            end
        end
        if ~isfield(pars.fit,'focusWeight'), pars.fit.focusWeight = 1; end


        % HSVD nuisance signal removal of everything EXCEPT the NAD+ region 
        if isfield(pars.fit,'hsvdClean') && isfield(pars.fit.hsvdClean,'enable') && pars.fit.hsvdClean.enable
            if ~isfield(pars.fit.hsvdClean,'K'), pars.fit.hsvdClean.K = 60; end
            if ~isfield(pars.fit.hsvdClean,'range')
                pars.fit.hsvdClean.range = [8.75 12];
            end

            disp('HSVD nuisance signal removal...')
            fid_temp = ifft(ifftshift(ws, 1), [], 1);
            excludeRangeHz = (pars.fit.hsvdClean.range - center_freq) * f0;
            [f_nuis,~,~,~, y_nuisance] = hsvd_fit(t(:), double(fid_temp), ...
                pars.fit.hsvdClean.K, [], [], [], excludeRangeHz);
            fid_temp = fid_temp - y_nuisance;
            ws_clean = fftshift(fft(fid_temp, [], 1), 1);

            output.met.fit.hsvdClean.nuisance_freqHz = f_nuis;
            fprintf('  %d HSVD components removed (outside [%.1f - %.1f] ppm)\n', ...
                numel(f_nuis), min(pars.fit.hsvdClean.range), max(pars.fit.hsvdClean.range))

            if pars.plt
                figure('Name','1H HSVD nuisance removal')
                plot(ppm, real(ws), 'k', ppm, real(ws_clean), 'r')
                set(gca, 'xdir', 'reverse')
                legend({'original','HSVD cleaned'})
                title('HSVD nuisance removal')
            end

            ws = ws_clean;
        end

        % ============================================================
        % Region setup (multi-region capable)
        % ============================================================
        % pars.fit.regions (optional) is a struct array; each element:
        %   .fitRange - [min max] ppm window (required)
        %   .peaks    - cell of names from pars.peaks.name active in this
        %               region (default: peaks whose center is in fitRange)
        %   .name     - optional label
        %   .mode/.ph_range/.baselineOpt/.lineShapeOpt/.focusRange/.focusWeight
        %             - optional per-region overrides of the shared settings
        % When unset, a single implicit region spans the coarse range. 
        % pars.fit.coarseMode controls how coarse fit + phase correction are
        % derived for multiple regions:
        %   'perRegion' (default) - independent coarse->phase->fine per region
        %   'union'               - one coarse fit + phase over the union of
        %                           all region windows, then per-region fine
        useRegions = isfield(pars.fit,'regions') && ~isempty(pars.fit.regions);

        minw = 10/f0;

        if isfield(pars.fit,'ppm_range_coarse')
            baseRange = pars.fit.ppm_range_coarse;
        elseif isfield(pars.fit,'ppm_range')
            baseRange = pars.fit.ppm_range;
        else
            baseRange = [min(ppm) max(ppm)];
        end

        if useRegions
            regionList = pars.fit.regions(:);
        else
            regionList = struct('fitRange', {[min(baseRange) max(baseRange)]}, ...
                                'peaks', {pars.peaks.name}, 'name', {'all'});
        end
        nRegions = numel(regionList);

        if isfield(pars.fit,'coarseMode') && ~isempty(pars.fit.coarseMode)
            coarseMode = lower(pars.fit.coarseMode);
        else
            coarseMode = 'perregion';
        end
        if ~ismember(coarseMode, {'perregion','union'})
            error('pars.fit.coarseMode must be ''perRegion'' or ''union''')
        end
        % 'union' re-partitions ONE HSVD component set per region
        % (unionCoarseForRegion).  The VARPRO coarse path produces no
        % components to re-partition, so the combination is unsupported --
        % fail loudly rather than silently degrade.
        if strcmp(coarseMode,'union') && isfield(pars.fit,'hsvdCoarse') && ...
                isfield(pars.fit.hsvdCoarse,'enable') && ...
                ~isempty(pars.fit.hsvdCoarse.enable) && ~pars.fit.hsvdCoarse.enable
            error(['coarseMode ''union'' requires hsvdCoarse.enable = true; ' ...
                   'the VARPRO coarse path has no component set to re-partition. ' ...
                   'Use coarseMode ''perRegion'' with hsvdCoarse.enable = false.'])
        end

        coarseFO = false;
        if isfield(pars.fit,'coarseFirstOrder') && ~isempty(pars.fit.coarseFirstOrder)
            coarseFO = logical(pars.fit.coarseFirstOrder);
        end

        % Per-peak enable for the HSVD width prior (lambdaWidth).  Logical
        % vector in pars.peaks.name order, or [] for all peaks.
        widthPriorPeaks = [];
        if isfield(pars.fit.baselineOpt,'widthPriorPeaks') && ...
                ~isempty(pars.fit.baselineOpt.widthPriorPeaks)
            widthPriorPeaks = logical(pars.fit.baselineOpt.widthPriorPeaks(:))';
            assert(numel(widthPriorPeaks) == numel(pars.peaks.name), ...
                'baselineOpt.widthPriorPeaks must be logical, length %d (pars.peaks.name)', ...
                numel(pars.peaks.name));
        end

        % resolve each region's peak indices into pars.peaks.name
        regSel = cell(nRegions,1);
        peakCenters = mean(pars.peaks.range, 2);
        for r = 1:nRegions
            reg = regionList(r);
            if isfield(reg,'peaks') && ~isempty(reg.peaks)
                nm = reg.peaks;
                if ischar(nm), nm = {nm}; end
                sel = zeros(1, numel(nm));
                for q = 1:numel(nm)
                    ii = find(strcmpi(pars.peaks.name, nm{q}), 1);
                    if isempty(ii)
                        error('region %d peak "%s" not found in pars.peaks.name', r, nm{q})
                    end
                    sel(q) = ii;
                end
            else
                sel = find(peakCenters(:).' > min(reg.fitRange) & peakCenters(:).' < max(reg.fitRange));
                if isempty(sel)
                    error('region %d (%g-%g ppm) contains no peaks; set .peaks or widen fitRange', ...
                        r, min(reg.fitRange), max(reg.fitRange))
                end
            end
            regSel{r} = sel;
        end

        % ============================================================
        % Step 1: Coarse VARPRO fit(s) + phase correction
        % ============================================================
        npks = numel(pars.peaks.name);
        center0_all = zeros(npks,1);
        width0_all  = zeros(npks,1);
        ph0_all     = zeros(npks,1);
        pcInfo      = repmat(struct('pc0',0,'pc1',0), nRegions, 1);
        coarseStore = cell(nRegions,1);
        wsP_store   = cell(nRegions,1);

        % HSVD coarse-seeding options (see nadCoarsePhase).  Per-peak vectors
        % are given in pars.peaks.name order and subset per region below.
        hsvdC = struct();
        if isfield(pars.fit,'hsvdCoarse') && ~isempty(pars.fit.hsvdCoarse)
            hsvdC = pars.fit.hsvdCoarse;
        end

        % Peaks the coarse HSVD actually found a qualifying component for.
        foundAll = false(1, npks);

        if strcmp(coarseMode,'union')
            allSel = unique(cat(2, regSel{:}));
            allSelC = allSel(~isBroad(allSel));      % broad peaks excluded from phase fit
            if isempty(allSelC), allSelC = allSel; end
            allRanges = cat(1, regionList.fitRange);
            unionRange = [min(allRanges(:)) max(allRanges(:))];
            disp('Step 1: Coarse VARPRO (union)...')
            boC = pars.fit.baselineOpt;  boC.widthBounds = peakWB(allSelC,:);
            if ~isempty(peakEta), boC.etaInit = peakEta(allSelC); end
            wiC = []; if isfield(pars.peaks,'widthInit') && ~isempty(pars.peaks.widthInit), wiC = pars.peaks.widthInit(allSelC); end
            [wsP, c0, w0, p0, pc0, pc1, cf] = nadCoarsePhase(ws, ppm, f0, center_freq, ...
                unionRange, pars.peaks.range(allSelC,:), pars.fit.mode, ...
                boC, minw, pars.plt, 'union', wiC, subsetHsvdOpt(hsvdC, allSelC, npks), ...
                coarseFO);
            center0_all(allSelC) = c0;  width0_all(allSelC) = w0;  ph0_all(allSelC) = p0;
            if isstruct(cf) && isfield(cf,'hsvd') && isfield(cf.hsvd,'found')
                foundAll(allSelC) = logical(cf.hsvd.found(:)');
            end
            allSelB = allSel(isBroad(allSel));        % seed broad peaks for the fine fit
            if ~isempty(allSelB)
                center0_all(allSelB) = mean(pars.peaks.range(allSelB,:), 2);
                width0_all(allSelB)  = broadW0;
                ph0_all(allSelB)     = 0;
            end
            % One coarse fit and one phase for every region, but the baseline
            % partition is redone per region: each region's own peaks are the
            % ones removed, everything else (including components centred
            % outside the region) stays in that region's baseline.  Also puts
            % baselinePhased on the region's own grid, which is what the
            % lambdaPrior hand-off below checks the length of.
            for r = 1:nRegions
                pcInfo(r).pc0 = pc0;  pcInfo(r).pc1 = pc1;
                wsP_store{r}  = wsP;
                selR  = regSel{r};
                selRC = selR(~isBroad(selR));
                if isempty(selRC), selRC = selR; end
                coarseStore{r} = unionCoarseForRegion(cf, ppm, f0, center_freq, ...
                    regionList(r).fitRange, pars.peaks.range(selRC,:), ...
                    subsetHsvdOpt(hsvdC, selRC, npks), pc0, pc1);
            end
            fprintf('Phase correction (union): pc0 = %.1f deg, pc1 = %.1f deg\n', pc0, pc1)
        else
            disp('Step 1: Coarse VARPRO (per region)...')
            for r = 1:nRegions
                sel  = regSel{r};
                selC = sel(~isBroad(sel));        % broad peaks excluded from phase fit
                if isempty(selC), selC = sel; end
                lbl = regionLabel(regionList, r);
                boC = pars.fit.baselineOpt;  boC.widthBounds = peakWB(selC,:);
                if ~isempty(peakEta), boC.etaInit = peakEta(selC); end
                wiC = []; if isfield(pars.peaks,'widthInit') && ~isempty(pars.peaks.widthInit), wiC = pars.peaks.widthInit(selC); end
                [wsP_r, c0, w0, p0, pc0, pc1, cf] = nadCoarsePhase(ws, ppm, f0, center_freq, ...
                    regionList(r).fitRange, pars.peaks.range(selC,:), pars.fit.mode, ...
                    boC, minw, pars.plt, lbl, wiC, subsetHsvdOpt(hsvdC, selC, npks), ...
                    coarseFO);
                center0_all(selC) = c0;  width0_all(selC) = w0;  ph0_all(selC) = p0;
                if isstruct(cf) && isfield(cf,'hsvd') && isfield(cf.hsvd,'found')
                    foundAll(selC) = logical(cf.hsvd.found(:)');
                end
                pcInfo(r).pc0 = pc0;  pcInfo(r).pc1 = pc1;
                wsP_store{r}  = wsP_r;  coarseStore{r} = cf;
                selB = sel(isBroad(sel));         % seed broad peaks for the fine fit
                if ~isempty(selB)
                    center0_all(selB) = mean(pars.peaks.range(selB,:), 2);
                    width0_all(selB)  = broadW0;
                    ph0_all(selB)     = 0;
                end
                fprintf('  region %d (%s): pc0 = %.1f deg, pc1 = %.1f deg\n', r, lbl, pc0, pc1)
            end
        end

        disp('Coarse VARPRO parameters:')
        for kk = 1:npks
            fprintf('  %s: pos=%.3f ppm, width=%.4f ppm, phase=%.1f deg\n', ...
                pars.peaks.name{kk}, center0_all(kk), width0_all(kk), ph0_all(kk))
        end

        % Optional per-peak width seed override for the FINE fit.  Use when a
        % peak parks in a narrow local minimum during coarse (e.g. broad Trp
        % swamped by the NAD-dominant LS).  Vector aligned with
        % pars.peaks.name; entries <=0 or NaN keep the coarse-fitted width.
        if isfield(pars.peaks,'widthInit') && ~isempty(pars.peaks.widthInit)
            wi = pars.peaks.widthInit(:);
            assert(numel(wi) == npks, ...
                'pars.peaks.widthInit must be length numel(pars.peaks.name) = %d', npks);
            override = wi > 0 & ~isnan(wi);
            if any(override)
                fprintf('Per-peak widthInit override (fine-fit seed):\n');
                for kk = find(override(:))'
                    fprintf('  %s: %.4f -> %.4f ppm\n', ...
                        pars.peaks.name{kk}, width0_all(kk), wi(kk));
                end
                width0_all(override) = wi(override);
            end
        end

        % ============================================================
        % Step 2: Fine VARPRO fit per region (on phased spectrum)
        % ============================================================
        disp('Step 2: Fine VARPRO fit...')
        nppm         = numel(ppm);
        areasAll     = zeros(1, npks);   % base-Voigt area (ampl .* width), kernel-blind
        areasConvAll = zeros(1, npks);   % kernel-aware (trapz of convolved component)
        compAll      = nan(nppm, npks);
        fitFull  = nan(nppm, 1);
        baseFull = nan(nppm, 1);
        specFull = nan(nppm, 1);
        modesUsed   = zeros(nRegions,1);
        regFitStore = cell(nRegions,1);

        for r = 1:nRegions
            reg = regionList(r);
            sel = regSel{r};
            inds_r = find(ppm > min(reg.fitRange) & ppm < max(reg.fitRange));
            x_r = ppm(inds_r);
            y_r = double(real(wsP_store{r}(inds_r)));

            modeR = pars.fit.mode;
            if isfield(reg,'mode') && ~isempty(reg.mode), modeR = reg.mode; end
            phR = pars.fit.ph_range;
            if isfield(reg,'ph_range') && ~isempty(reg.ph_range), phR = reg.ph_range; end
            baseR = pars.fit.baselineOpt;
            if isfield(reg,'baselineOpt') && isstruct(reg.baselineOpt)
                fn = fieldnames(reg.baselineOpt);
                for ii = 1:numel(fn), baseR.(fn{ii}) = reg.baselineOpt.(fn{ii}); end
            end
            % per-peak linewidth bounds (broad components wide), unless the
            % region explicitly overrode widthBounds itself
            if ~(isfield(reg,'baselineOpt') && isfield(reg.baselineOpt,'widthBounds'))
                baseR.widthBounds = peakWB(sel,:);
            end
            if ~isempty(peakEta)
                baseR.etaInit = peakEta(sel);
            end
            % per-peak phase bounds: broad components near-absorptive, narrow
            % peaks keep the region phase range (phR)
            pbR = repmat(phR(:).', numel(sel), 1);
            ib  = isBroad(sel);
            if any(ib), pbR(ib,:) = repmat(broadPB, sum(ib), 1); end
            baseR.phaseBounds = pbR;

            % ---- HSVD-derived soft priors ------------------------------
            % The coarse HSVD step already produces the unassigned-component sum
            % (a direct measurement of the baseline under these peaks) and a
            % per-peak linewidth.  Feed both in as PRIORS, not constraints --
            % baselineOpt.lambdaPrior / .lambdaWidth control how hard the fit
            % is pulled toward them, and 0 reproduces the old behavior.
            cfR = coarseStore{r};
            baseR.priorBaseline = [];
            if isstruct(cfR) && isfield(cfR,'baselinePhased') && ...
                    numel(cfR.baselinePhased) == numel(inds_r)
                baseR.priorBaseline = cfR.baselinePhased;
            elseif isfield(baseR,'lambdaPrior') && ~isempty(baseR.lambdaPrior) ...
                    && baseR.lambdaPrior > 0
                % Do not let a missing or mis-sized prior quietly turn itself
                % off -- the fit would still run, just with a free baseline,
                % and nothing downstream would show it.
                if isstruct(cfR) && isfield(cfR,'baselinePhased')
                    why = sprintf('coarse baseline is %d points but the region grid is %d', ...
                                  numel(cfR.baselinePhased), numel(inds_r));
                else
                    why = ['the coarse estimator supplied no baseline ' ...
                           '(hsvdCoarse.enable = false has no component decomposition)'];
                end
                warning('dat2mat_svsNAD:priorGridMismatch', ...
                    'region %d (%s): %s; lambdaPrior=%g will have NO effect here.', ...
                    r, regionLabel(regionList, r), why, baseR.lambdaPrior);
            end
            % Width prior: the HSVD linewidth seed, but only for peaks HSVD actually found
            pwR = width0_all(sel);
            pwR(~foundAll(sel)) = NaN;      % NaN entries are skipped downstream
            if ~isempty(widthPriorPeaks)
                pwR(~widthPriorPeaks(sel)) = NaN;
            end
            baseR.priorWidth = pwR;

            lsR = [];
            if isfield(pars.fit,'lineShapeOpt'), lsR = pars.fit.lineShapeOpt; end
            if isfield(reg,'lineShapeOpt') && ~isempty(reg.lineShapeOpt), lsR = reg.lineShapeOpt; end
            focusR = pars.fit.focusRange;
            focusW = pars.fit.focusWeight;
            if isfield(reg,'focusRange'),  focusR = reg.focusRange;  end
            if isfield(reg,'focusWeight'), focusW = reg.focusWeight; end

            outR = nadFineFit(x_r, y_r, center0_all(sel), width0_all(sel), ...
                pars.peaks.range(sel,:), baseR, focusR, focusW, modeR, phR, minw, lsR);

            regFitStore{r} = outR;
            modesUsed(r)   = modeR;
            areasAll(sel)  = outR.areas;
            if isfield(outR, 'areasConv')
                areasConvAll(sel) = outR.areasConv;
            else
                areasConvAll(sel) = outR.areas;
            end
            compAll(inds_r, sel) = outR.comp;
            fitFull(inds_r)  = outR.fit;
            baseFull(inds_r) = outR.baseline;
            specFull(inds_r) = y_r;

            if pars.plt
                figure('Name', sprintf('1H fine fit - region %d %s', r, regionLabel(regionList,r)))
                plot(x_r, y_r, 'k', x_r, outR.fit, 'r', x_r, y_r - outR.fit, 'g')
                hold on, plot(x_r, outR.baseline, 'b--')
                legend({'data','fit','residual','baseline'})
                set(gca,'xdir','reverse')
                title(sprintf('Step 2: Fine VARPRO fit (region %d: %s)', r, regionLabel(regionList,r)))

                % per-component diagnostic plots (overlay + per-peak gallery),
                rLbl    = regionLabel(regionList, r);
                pkNames = pars.peaks.name(sel);
                ttl     = sprintf('1H fit - region %d: %s', r, rLbl);
                fOv = plot_varpro_1h(x_r, y_r, outR, pkNames, ttl);
                set(fOv, 'Name', sprintf('1H varpro fit - region %d %s', r, rLbl));
                % fGl = plot_componentSpectra_1h(x_r, y_r, outR, pkNames, ttl);
                % set(fGl, 'Name', sprintf('1H component spectra - region %d %s', r, rLbl));
            end
        end

        % union window axis for combined storage
        allRanges = cat(1, regionList.fitRange);
        unionInds = find(ppm > min(allRanges(:)) & ppm < max(allRanges(:)));
        x = ppm(unionInds);

        % Each region is fit independently. Downstream code expects one vector of pars.
        if all(modesUsed == modesUsed(1))
            switch modesUsed(1)
                case {1,2,6}, nBlk = 3;
                case {3,4},   nBlk = 4;
                case 7,       nBlk = 5;
                case 5,       nBlk = 6;
                otherwise,    nBlk = 3;
            end
            parsFlat = zeros(1, nBlk*npks);
            lbFlat   = zeros(1, nBlk*npks);
            ubFlat   = zeros(1, nBlk*npks);
            for r = 1:nRegions
                sel = regSel{r}; n_r = numel(sel);
                for b = 1:nBlk
                    blk = (b-1)*n_r + (1:n_r);
                    parsFlat((b-1)*npks + sel) = regFitStore{r}.pars(blk);
                    lbFlat((b-1)*npks + sel)   = regFitStore{r}.lb(blk);
                    ubFlat((b-1)*npks + sel)   = regFitStore{r}.ub(blk);
                end
            end
            parnamesFlat = regFitStore{1}.parnames(1:nBlk);
        else
            parsFlat = []; lbFlat = []; ubFlat = []; parnamesFlat = {};
        end

        % store output (peak-ordered to match pars.peaks.name)
        if useRegions
            output.met.fit.spec_full = ws;
        else
            output.met.fit.spec_full = wsP_store{1};
        end
        output.met.fit.spec      = specFull(unionInds);
        output.met.fit.spec_fit  = fitFull(unionInds);
        output.met.fit.baseline  = baseFull(unionInds);
        output.met.fit.comp      = compAll(unionInds, :);
        output.met.fit.pars      = parsFlat;
        output.met.fit.parnames  = parnamesFlat;
        output.met.fit.areas     = areasAll;
        % Kernel-aware areas
        output.met.fit.areasConv = areasConvAll;
        output.met.fit.names     = pars.peaks.name;
        output.met.fit.ppm       = x;
        output.met.fit.lb        = lbFlat;
        output.met.fit.ub        = ubFlat;
        output.met.fit.init.center0 = center0_all;
        output.met.fit.init.width0  = width0_all;
        output.met.fit.init.ph0     = ph0_all;
        output.met.fit.coarse.fit      = coarseStore{1}.fit;
        output.met.fit.coarse.baseline = coarseStore{1}.baseline;
        output.met.fit.phaseCorr.pc0   = pcInfo(1).pc0;
        output.met.fit.phaseCorr.pc1   = pcInfo(1).pc1;
        output.met.fit.phaseCorr.pivot = center_freq;
        output.met.fit.coarseMode      = coarseMode;
        for r = 1:nRegions
            output.met.fit.region(r).name       = regionLabel(regionList, r);
            output.met.fit.region(r).fitRange   = regionList(r).fitRange;
            output.met.fit.region(r).peakSelect = regSel{r};
            output.met.fit.region(r).peakNames  = pars.peaks.name(regSel{r});
            output.met.fit.region(r).fit        = regFitStore{r};
            output.met.fit.region(r).coarse     = coarseStore{r};
            output.met.fit.region(r).pc0        = pcInfo(r).pc0;
            output.met.fit.region(r).pc1        = pcInfo(r).pc1;
        end

        % --- Water area estimation ---
        if ~isempty(nws) && isfield(pars,'watfit')
            ppm = output.wat.ppm;
            watfit = pars.watfit;
            if ~isfield(watfit,'method'), watfit.method = 'integrate'; end

            % extract water range
            if ~strcmpi(watfit.method,'fid')
                wat_inds = find( ppm > min(watfit.ppm_range) & ppm < max(watfit.ppm_range));
                xw = ppm(wat_inds);
                yw = double(real(nws(wat_inds)));
            end

            if strcmpi(watfit.method, 'fid')
                % --- FID first-point method ---
                % fid(1) = (1/N) * sum(spec) by MATLAB FFT convention.
                % Spectral integral in ppm:
                %   trapz(ppm, spec) ≈ sum(spec) * dppm = fid(1) * N * dppm
                disp('Water area: FID first-point method')
                nws_fid = output.wat.fid;
                N = size(nws_fid, 1);
                dppm = abs(ppm(2) - ppm(1));
                waterArea = real(nws_fid(1))/2 * N * dppm; % halve first time point to reduce DC offset

                output.wat.fit.method = 'fid';
                output.wat.fit.areaTotal = waterArea;
                fprintf('  Water area (FID): %.4g\n', waterArea)

            elseif strcmpi(watfit.method, 'integrate')
                % --- Numerical integration with linear baseline ---
                disp('Water area: numerical integration')
                bl = linspace(yw(1), yw(end), length(yw))';
                waterArea = trapz(xw, yw - bl);

                output.wat.fit.method = 'integrate';
                output.wat.fit.areaTotal = waterArea;
                output.wat.fit.spec = yw;
                output.wat.fit.baseline = bl;
                output.wat.fit.ppm = xw;
                fprintf('  Water area (integrate): %.4g\n', waterArea)

                if pars.plt
                    figure('Name','1H water integration')
                    plot(xw, yw, 'k', xw, bl, 'b--')
                    hold on
                    area(xw, yw - bl, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
                    legend({'data','baseline','integrated area'})
                    set(gca,'xdir','reverse')
                    title('Water integration')
                end

            elseif strcmpi(watfit.method, 'fit')
                % --- Multi-component VARPRO fit (original approach) ---
                disp('Water area: VARPRO fit')
                if ~isfield(watfit,'ph_range'), watfit.ph_range = pars.fit.ph_range; end
                if ~isfield(watfit,'baselineOpt'), watfit.baselineOpt = pars.fit.baselineOpt; end
                if ~isfield(watfit,'mode'), watfit.mode = pars.fit.mode; end
                watfit.baselineOpt.weights = []; % make sure no weighting for water fit

                nWatComp = 3;
                if isfield(watfit,'nComp'), nWatComp = watfit.nComp; end
                wat_peaks_name  = arrayfun(@(k) sprintf('Water%d',k), 1:nWatComp, 'uni', false);
                wat_peaks_range = repmat(watfit.ppm_range, nWatComp, 1);

                w0 = max(0.05, minw);
                if nWatComp == 1
                    amp0_w    = 1;
                    center0_w = 4.72;
                    width0_w  = w0;
                else
                    offsets   = linspace(-0.05, 0.05, nWatComp);
                    amp0_w    = ones(1, nWatComp);
                    center0_w = 4.72 + offsets;
                    width0_w  = w0 * ones(1, nWatComp);
                end

                watBaseOpt = watfit.baselineOpt;
                % The metabolite fit may carry PER-PEAK linewidth bounds. Those rows are meaningless here.
                if isfield(watBaseOpt,'widthBoundsHz')
                    watBaseOpt = rmfield(watBaseOpt,'widthBoundsHz');
                end
                % same for the per-peak Voigt seed and width-prior mask
                for fq = {'etaInit','widthPriorPeaks'}
                    if isfield(watBaseOpt, fq{1})
                        watBaseOpt = rmfield(watBaseOpt, fq{1});
                    end
                end
                if isfield(watfit,'etaInit') && ~isempty(watfit.etaInit)
                    watBaseOpt.etaInit = watfit.etaInit;
                end
                if isfield(watBaseOpt,'widthBounds') && ...
                        size(watBaseOpt.widthBounds,1) > 1 && ...
                        size(watBaseOpt.widthBounds,1) ~= nWatComp
                    wbW = watBaseOpt.widthBounds;
                    watBaseOpt.widthBounds = [min(wbW(:,1)) max(wbW(:,2))];
                end
                watBaseOpt.TolFun = 1e-8;
                watBaseOpt.TolX   = 1e-8;
                watBaseOpt.centerBounds = wat_peaks_range;
                watBaseOpt.knotSpacing = max(watBaseOpt.knotSpacing, diff(watfit.ppm_range));
                localInds_w = xw > min(wat_peaks_range(:)) & xw < max(wat_peaks_range(:));
                if any(localInds_w)
                    watBaseOpt.amplUB = max(abs(real(yw(localInds_w))));
                else
                    watBaseOpt.amplUB = max(abs(real(yw)));
                end

                if isfield(watfit,'lineShapeOpt')
                    watFitOutput = curvefitAuto_varproBaseline(xw, yw, watfit.mode, amp0_w, center0_w, width0_w, watfit.ph_range, minw, watBaseOpt, watfit.lineShapeOpt);
                else
                    watFitOutput = curvefitAuto_varproBaseline(xw, yw, watfit.mode, amp0_w, center0_w, width0_w, watfit.ph_range, minw, watBaseOpt);
                end

                if pars.plt
                    figure('Name','1H water VARPRO fit')
                    plot(xw, real(yw), 'k', xw, watFitOutput.fit, 'r', xw, real(yw) - watFitOutput.fit, 'g')
                    hold on
                    plot(xw, watFitOutput.baseline, 'b--')
                    legend({'data','fit','residual','baseline'})
                    set(gca,'xdir','reverse')
                    title('VARPRO water fit')
                end

                output.wat.fit.method = 'fit';
                output.wat.fit.spec = yw;
                output.wat.fit.spec_fit = watFitOutput.fit;
                output.wat.fit.baseline = watFitOutput.baseline;
                output.wat.fit.comp = watFitOutput.comp;
                output.wat.fit.pars = watFitOutput.pars;
                output.wat.fit.parnames = watFitOutput.parnames;
                output.wat.fit.areas = watFitOutput.areas;
                if isfield(watFitOutput, 'areasConv')
                    output.wat.fit.areasConv = watFitOutput.areasConv;
                else
                    output.wat.fit.areasConv = watFitOutput.areas;
                end
                % areaTotal = canonical kernel-aware water area
                output.wat.fit.areaTotal    = sum(output.wat.fit.areasConv);
                output.wat.fit.areaTotalWin = trapz(xw(:), watFitOutput.fit(:) - watFitOutput.baseline(:));
                output.wat.fit.areaBase     = sum(watFitOutput.areas);
                output.wat.fit.names = wat_peaks_name;
                output.wat.fit.ppm = xw;
                output.wat.fit.lb = watFitOutput.lb;
                output.wat.fit.ub = watFitOutput.ub;
                fprintf('  Water area (fit, ext-axis): %.4g  (in-window: %.4g, base-only: %.4g)\n', ...
                    output.wat.fit.areaTotal, output.wat.fit.areaTotalWin, output.wat.fit.areaBase)

            else
                error('Unknown watfit.method: %s. Use ''integrate'', ''fid'', or ''fit''.', watfit.method)
            end
        end

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

        % water fit - automatic
        if ~isempty(nws)
            x = output.wat.ppm;
            y = double(real(nws));
            [yfit,areas,ip] = curvefitWat(x,y,minw,pars.fit.mode);
            output.wat.fit.spec = y;
            output.wat.fit.ppm = x;
            output.wat.fit.spec_fit = yfit;
            output.wat.fit.pars = ip;
            output.wat.fit.areas = areas;
        end


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
            figure('Name','1H metabolite fit (manual)'), plot(x,real(y),'k',x,yfit,'r',x,real(y)-yfit,'g'), legend({'data','fit','residual'}), set(gca,'xdir','reverse')
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
            x = output.wat.ppm;
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

fclose('all');

rmpath([thisPath filesep 'mapVBVD']);
rmpath([thisPath filesep 'utils']);


% =================================================================== %
% local helpers for the multi-region VARPRO pipeline
% =================================================================== %

function lbl = regionLabel(regionList, r)
% Region display label: its .name if set, else 'region<r>'.
    lbl = sprintf('region%d', r);
    if isfield(regionList, 'name') && ~isempty(regionList(r).name)
        lbl = regionList(r).name;
    end

function o = subsetHsvdOpt(o, sel, npks)
% Subset per-peak HSVD option vectors to the peaks of one region.  Fields
% given as length-npks vectors (pars.peaks.name order) are indexed by sel;
% scalars pass through untouched so a single cap can cover every peak.
    if isempty(o), o = struct(); return; end
    fn = {'lwMinHz','lwMaxHz','phaseTolDeg'};
    for ii = 1:numel(fn)
        if isfield(o, fn{ii}) && numel(o.(fn{ii})) == npks && npks > 1
            o.(fn{ii}) = o.(fn{ii})(sel);
        end
    end

function [assigned, found, sig, phDeg, ctr, wdt] = hsvdAssign( ...
        ppmC, lwHz, phC, amp, f0, peakRange, lwMin, lwMax, phTol)
% Assign HSVD components to peaks.  A component qualifies for peak kk when its
% centre lies in that peak's ppm window AND its linewidth and phase are inside
% the per-peak caps.  Everything left unassigned is baseline.
%
% Note what "unassigned" means: the test is against the PEAK TABLE, not against
% any region window.  A component sitting inside a region but under no peak is
% baseline, and so is a component centred far outside -- which is exactly what
% lets a region's baseline carry the tails of resonances outside it.
%
% ctr/wdt/phDeg are only meaningful where found(kk) is true; the caller keeps
% its own fallback seeds for the rest.
    np       = size(peakRange, 1);
    assigned = false(numel(ppmC), 1);
    found    = false(np, 1);
    sig      = zeros(np, 1);
    phDeg    = zeros(np, 1);
    ctr      = mean(peakRange, 2);
    wdt      = nan(np, 1);
    for kk = 1:np
        sel = ppmC >= min(peakRange(kk,:)) & ppmC <= max(peakRange(kk,:)) & ...
              lwHz >= lwMin(kk)            & lwHz <= lwMax(kk)            & ...
              phC  >= -phTol(kk)           & phC  <= phTol(kk);
        if ~any(sel), continue; end
        assigned  = assigned | sel;
        found(kk) = true;
        A = sum(amp(sel));
        w = abs(amp(sel));
        sig(kk)   = abs(A);
        phDeg(kk) = angle(A) * 180/pi;
        ctr(kk)   = sum(w .* ppmC(sel)) / sum(w);
        wdt(kk)   = sum(w .* lwHz(sel)) / sum(w) / f0;
    end

function Sc = hsvdSpectra(fHz, dHz, amp, N, dt)
% Full-length spectrum of every HSVD component, one column each.  Evaluating
% all N points (not just a window) is what allows a component centred outside
% a region to contribute its tail inside that region.
    tf = (0:N-1)' * dt;
    Yc = exp(tf * (-dHz(:).' + 1i*2*pi*fHz(:).')) .* amp(:).';
    Sc = fftshift(fft(Yc, [], 1), 1);

function cf = hsvdPartition(Sc, assigned, inds, ppm, center_freq, pc0, pc1)
% Split a component set into peaks (.fit) and baseline (.baseline) over the
% points inds, unphased and phased by [pc0 pc1].  The phased baseline is the
% one to hand the fine fit as a prior: it is in the same phase frame as the
% data the fine fit will see.
    cf.fit      = real(sum(Sc(inds,  assigned), 2));
    cf.baseline = real(sum(Sc(inds, ~assigned), 2));
    phvec = exp(-1i*pi/180*(pc0 + 2*pc1/abs(ppm(end)-ppm(1)) * (ppm(:)-center_freq)));
    cf.baselinePhased = real(sum(Sc(inds, ~assigned) .* phvec(inds), 2));
    cf.fitPhased      = real(sum(Sc(inds,  assigned) .* phvec(inds), 2));

function cfR = unionCoarseForRegion(cf, ppm, f0, center_freq, fitRange, ...
                                    peakRange, hsvdOpt, pc0, pc1)
% Re-partition the ONE global HSVD component set for a single region.
%
% coarseMode 'union' runs a single coarse fit so that both regions share one
% phase.  The baseline, though, must stay region-specific: the components to
% subtract are the ones THIS region models, so for the Trp region the NAD
% components stay in the baseline where they belong.  Assigning against the
% union peak table instead would strip the NAD shoulder out of the Trp
% baseline while region 2's model has no NAD peaks to put it back.
%
% No HSVD is re-run -- fHz/dHz/amp come straight off the stored struct.
    h  = cf.hsvd;
    np = size(peakRange, 1);
    lwMin = optVec(hsvdOpt,'lwMinHz',      10,  np);
    lwMax = optVec(hsvdOpt,'lwMaxHz',      90,  np);
    phTol = optVec(hsvdOpt,'phaseTolDeg', 180,  np);

    assignedR = hsvdAssign(h.ppm, h.lwHz, h.phaseDeg, h.amp, f0, ...
                           peakRange, lwMin, lwMax, phTol);

    inds = find(ppm > min(fitRange) & ppm < max(fitRange));
    Sc   = hsvdSpectra(h.fHz, h.dHz, h.amp, h.N, h.dt);
    cfR  = hsvdPartition(Sc, assignedR, inds, ppm, center_freq, pc0, pc1);

    cfR.pars        = [];
    cfR.phResidDeg  = cf.phResidDeg;      % phase is global, so this is too
    cfR.hsvd        = h;
    cfR.hsvd.assigned = assignedR;        % but the partition is this region's

function v = optVec(o, fld, dflt, n)
% Per-peak option vector, accepting a scalar, a length-n vector, or nothing.
    v = dflt;
    if isstruct(o) && isfield(o, fld) && ~isempty(o.(fld)), v = o.(fld); end
    v = v(:);
    if isscalar(v), v = v * ones(n,1); end
    assert(numel(v) == n, '%s must be scalar or length %d', fld, n);

function [wsP, center0, width0, ph0, pc0, pc1, coarseFit] = nadCoarsePhase( ...
        ws, ppm, f0, center_freq, fitRange, peakRange, mode, baselineOpt, minw, plt, label, widthInit, hsvdOpt, firstOrder)
% Coarse HSVD decomposition + phase correction derived from this region's peaks.
% Returns the phase-corrected full spectrum (wsP), coarse seeds for the
% region's peaks (center0/width0/ph0), the applied phase (pc0/pc1), and a
% coarse fit struct (.fit = assigned components, .baseline = unassigned).
%
% firstOrder (optional, default false) fits the 1st-order phase term; honoured
% by both coarse estimators.
%
% hsvdOpt fields (all optional):
%   .enable      true.  false selects the previous coarse estimator, a VARPRO
%                Voigt fit of the region (local function coarseVarpro).  That
%                path supplies no component decomposition, so it cannot feed
%                baselineOpt.lambdaPrior and cannot be combined with
%                coarseMode 'union'.
%   .K           60    HSVD order              (SpecTickle opt.hsvd_order)
%   .npts        1024  FID points decomposed   (SpecTickle opt.hsvd_np)
%   .lwMinHz     10    scalar or per-peak      (peaktable linewidth_min)
%   .lwMaxHz     90    scalar or per-peak      (peaktable linewidth_max)
%   .phaseTolDeg 180   scalar or per-peak      (peaktable phase_tol)
%
% Optional widthInit (length np vector): per-peak width fallback, used only
% for peaks where HSVD finds no qualifying component.  Entries <=0 or NaN
% keep the default 30/f0.
    inds = find(ppm > min(fitRange) & ppm < max(fitRange));
    x  = ppm(inds);
    y  = double(real(ws(inds)));
    np = size(peakRange, 1);

    if nargin < 13 || isempty(hsvdOpt), hsvdOpt = struct(); end
    if ~isfield(hsvdOpt,'enable')      || isempty(hsvdOpt.enable),      hsvdOpt.enable = true;   end
    if ~isfield(hsvdOpt,'K')           || isempty(hsvdOpt.K),           hsvdOpt.K = 60;          end
    if ~isfield(hsvdOpt,'npts')        || isempty(hsvdOpt.npts),        hsvdOpt.npts = 1024;     end
    if ~isfield(hsvdOpt,'lwMinHz')     || isempty(hsvdOpt.lwMinHz),     hsvdOpt.lwMinHz = 10;    end
    if ~isfield(hsvdOpt,'lwMaxHz')     || isempty(hsvdOpt.lwMaxHz),     hsvdOpt.lwMaxHz = 90;    end
    if ~isfield(hsvdOpt,'phaseTolDeg') || isempty(hsvdOpt.phaseTolDeg), hsvdOpt.phaseTolDeg = 180; end
    lwMin = hsvdOpt.lwMinHz(:);     if isscalar(lwMin), lwMin = lwMin*ones(np,1); end
    lwMax = hsvdOpt.lwMaxHz(:);     if isscalar(lwMax), lwMax = lwMax*ones(np,1); end
    phTol = hsvdOpt.phaseTolDeg(:); if isscalar(phTol), phTol = phTol*ones(np,1); end

    % fallback seeds (also used for any peak HSVD does not find)
    center0 = mean(peakRange, 2);
    width0  = (30/f0) * ones(np, 1);
    ph0     = zeros(np, 1);
    if nargin >= 12 && ~isempty(widthInit)
        wi = widthInit(:);
        assert(numel(wi) == np, 'widthInit length must match peakRange rows (%d)', np);
        m = wi > 0 & ~isnan(wi);
        width0(m) = wi(m);
    end

    % hsvdOpt.enable=false selects the PREVIOUS coarse estimator: a VARPRO
    % Voigt fit over the region (see coarseVarpro).  Until 2026-08 this flag
    % was parsed but never tested, and the VARPRO code it named was commented
    % out, so setting it false silently changed nothing.
    if nargin < 14 || isempty(firstOrder), firstOrder = false; end

    if ~hsvdOpt.enable
        [wsP, center0, width0, ph0, pc0, pc1, coarseFit] = coarseVarpro( ...
            ws, x, y, ppm, center_freq, peakRange, mode, baselineOpt, minw, ...
            plt, label, center0, width0, firstOrder);
        return
    end

    % ---- HSVD on the FID underlying this spectrum -----------------------
    N   = numel(ppm);
    bw  = N * f0 * abs(ppm(2) - ppm(1));   % Hz; ppm step = bw/(N*f0)
    dt  = 1/bw;
    fid = ifft(ifftshift(ws(:)));
    nu  = min(N, hsvdOpt.npts);
    tv  = (0:nu-1)' * dt;

    [fHz, dHz, amp] = hsvd_fit(tv, fid(1:nu), hsvdOpt.K);

    ppmC = center_freq + fHz(:)/f0;      % component centre, ppm
    lwHz = dHz(:)/pi;                    % Lorentzian FWHM, Hz
    phC  = angle(amp(:)) * 180/pi;       % component phase, deg

    % ---- assign components to peaks ---------------
    [assigned, found, sig, phA, ctrA, wdtA] = hsvdAssign(ppmC, lwHz, phC, amp, ...
                                                 f0, peakRange, lwMin, lwMax, phTol);
    ph0(found)     = phA(found);
    center0(found) = ctrA(found);
    width0(found)  = wdtA(found);
    if ~any(found)
        warning('nadCoarsePhase:noHSVDpeaks', ...
            '%s: HSVD found no qualifying components; phase left uncorrected.', label);
    end

    % ---- linear phase (pc0 + pc1) through this region's peaks -----------
    % Weighted by signal^2 so the slope is set by the well-determined peaks.
    use = find(found);
    if numel(use) < 2 || (max(center0(use)) - min(center0(use))) < 10*eps
        pc0 = 0;
        if ~isempty(use), pc0 = ph0(use(1)); end
        pc1 = 0;
        p   = [0 pc0];
    else
        [cU, ord] = sort(center0(use));
        phU = unwrap(ph0(use(ord)) * pi/180) * 180/pi;   % avoid +/-180 wrap
        wU  = sig(use(ord)).^2;                          % signal^2 weights
        if ~firstOrder
            % signal^2-weighted mean phase over the peaks
            pc0 = sum(wU(:) .* phU(:)) / sum(wU);
            pc1 = 0;
            p   = [0 pc0];
        else
            V   = [cU(:), ones(numel(cU),1)];
            p   = ((V' * (wU .* V)) \ (V' * (wU .* phU(:))))';
            pc0 = polyval(p, center_freq);           % 0th order at pivot
            % shiftSpectrumPhase convention: pc1 is total 1st-order phase
            % across half the bandwidth
            pc1 = p(1) * abs(ppm(end) - ppm(1)) / 2;
        end
    end
    wsP = shiftSpectrumPhase(ws, [pc0 pc1], ppm, center_freq);

    % ---- synthesize a coarseFit struct from the components --------------
    % .fit = peak-assigned components, .baseline = everything else
    Sc        = hsvdSpectra(fHz, dHz, amp, N, dt);
    coarseFit = hsvdPartition(Sc, assigned, inds, ppm, center_freq, pc0, pc1);
    coarseFit.pars     = [];
    % fHz/dHz/dt/N are kept so a caller holding this struct can re-partition
    % the SAME components for another window without re-running the HSVD
    % (see unionCoarseForRegion).
    coarseFit.hsvd     = struct('ppm',ppmC, 'lwHz',lwHz, 'phaseDeg',phC, ...
                                'amp',amp(:), 'fHz',fHz(:), 'dHz',dHz(:), ...
                                'dt',dt, 'N',N, 'assigned',assigned, ...
                                'found',found, 'signal',sig, ...
                                'polyCoef',p, 'K',numel(ppmC));
    % residual per-peak phase AFTER correction -- should sit near 0
    coarseFit.phResidDeg = ph0 - polyval(p, center0);

    fprintf('  %s: HSVD %d comps, %d/%d peaks seeded; resid phase %s deg\n', ...
        label, numel(ppmC), sum(found), np, ...
        strtrim(sprintf('%.1f ', coarseFit.phResidDeg(found))));

    if plt
        figure('Name', sprintf('1H coarse HSVD - %s', label))
        plot(x, y, 'k', x, coarseFit.fit + coarseFit.baseline, 'r', x, y - coarseFit.fit - coarseFit.baseline, 'g')
        hold on, plot(x, coarseFit.baseline, 'b--')
        legend({'data','assigned comps','residual','unassigned (baseline)'})
        set(gca,'xdir','reverse')
        title(sprintf('Coarse HSVD: %s', label))
    end

function [wsP, center0, width0, ph0, pc0, pc1, coarseFit] = coarseVarpro( ...
        ws, x, y, ppm, center_freq, peakRange, mode, baselineOpt, minw, ...
        plt, label, center0, width0, firstOrder)
% Coarse estimator used before 041082f (2026-07-23): a loose VARPRO fit of the
% region with the SAME lineshape model as the fine fit, from which the seeds
% and the linear phase are read off.  Selected by hsvdCoarse.enable = false.
%
% Differs from the HSVD path in what it can and cannot supply: there is no
% component decomposition, so coarseFit carries no .hsvd and no
% .baselinePhased.  Callers must therefore treat both as optional --
% baselineOpt.lambdaPrior has no prior to work with here, and coarseMode
% 'union' cannot re-partition per region (dat2mat_svsNAD errors on that
% combination rather than silently dropping the prior).
%
% center0/width0 arrive pre-seeded with the caller's fallbacks so the two
% paths start from the same place.
    np = size(peakRange, 1);

    bo = baselineOpt;
    bo.centerBounds = peakRange;
    bo.TolFun = 1e-3;
    bo.TolX   = 1e-3;
    aub = zeros(np, 1);
    for kk = 1:np
        rng = peakRange(kk,:);
        li  = x > min(rng) & x < max(rng);
        if any(li)
            aub(kk) = max(abs(real(y(li))));
        else
            aub(kk) = max(abs(real(y)));
        end
    end
    bo.amplUB = aub;

    coarseFit = curvefitAuto_varproBaseline(x, y, mode, ones(np,1), center0, width0, ...
        [-179 179], minw, bo);

    center0 = coarseFit.pars(np + (1:np))';
    if mode == 5 || mode == 7
        ph0    = coarseFit.pars(3*np + (1:np))';
        width0 = coarseFit.pars(5*np + (1:np))';   % total Voigt width
    elseif mode == 3 || mode == 4
        ph0    = coarseFit.pars(3*np + (1:np))';
        width0 = coarseFit.pars(2*np + (1:np))';
    else
        ph0    = zeros(np, 1);
        width0 = coarseFit.pars(2*np + (1:np))';
    end

    % linear phase (pc0 + pc1) through this region's peaks
    if nargin < 14 || isempty(firstOrder), firstOrder = false; end
    if np == 1
        pc0 = ph0;
        pc1 = 0;
    elseif ~firstOrder
        pc0 = mean(ph0);                           % zeroth order only (default)
        pc1 = 0;
    else
        p   = polyfit(center0(:), ph0(:), 1);
        pc0 = polyval(p, center_freq);             % 0th order at pivot
        % shiftSpectrumPhase convention: pc1 is total 1st-order phase
        % across half the bandwidth
        pc1 = p(1) * abs(ppm(end) - ppm(1)) / 2;
    end
    % NEGATED, unlike the HSVD path.  ph0 here is the phase parameter of the
    % Voigt model itself, and that model applies conj(Vc) (see
    % compositeV_complex_unit): its phase runs opposite to shiftSpectrumPhase.
    % Applying +pc0 would ADD the measured phase instead of removing it,
    % roughly doubling the error -- which is exactly what happened between
    % 1a4640b and 041082f, when this estimator was live alongside the conj.
    % The HSVD path is unaffected: it reads phase from component amplitudes,
    % which already run in shiftSpectrumPhase's direction.
    wsP = shiftSpectrumPhase(ws, [-pc0 -pc1], ppm, center_freq);

    fprintf('  %s: coarse VARPRO; pc0 = %.1f deg, pc1 = %.1f deg (applied negated)\n', ...
            label, pc0, pc1);

    if plt
        figure('Name', sprintf('1H coarse fit - %s', label))
        plot(x, y, 'k', x, coarseFit.fit, 'r', x, y - coarseFit.fit, 'g')
        hold on, plot(x, coarseFit.baseline, 'b--')
        legend({'data','coarse fit','residual','baseline'})
        set(gca,'xdir','reverse')
        title(sprintf('Coarse VARPRO fit: %s', label))
    end

function out = nadFineFit(x, y, center0, width0, peakRange, baselineOpt, ...
        focusRange, focusWeight, mode, ph_range, minw, lineShapeOpt)
% Fine VARPRO fit for one region: per-peak amplitude upper bounds and center
% bounds from peakRange, optional focus weighting, and optional shared extra
% lineshape kernel.  Phase bounds straddle 0 (spectrum already corrected).
    np = size(peakRange, 1);
    bo = baselineOpt;
    bo.centerBounds = peakRange;
    aub = zeros(np, 1);
    for kk = 1:np
        rng = peakRange(kk,:);
        li = x > min(rng) & x < max(rng);
        if any(li)
            aub(kk) = max(abs(real(y(li))));
        else
            aub(kk) = max(abs(real(y)));
        end
    end
    bo.amplUB = aub;

    wts = ones(size(x(:)));
    if ~isempty(focusRange) && focusWeight > 1
        fi = x(:) > min(focusRange) & x(:) < max(focusRange);
        wts(fi) = focusWeight;
    end
    bo.weights = wts;

    amp0 = ones(np, 1);
    if ~isempty(lineShapeOpt)
        out = curvefitAuto_varproBaseline(x, y, mode, amp0, center0, width0, ...
            ph_range, minw, bo, lineShapeOpt);
    else
        out = curvefitAuto_varproBaseline(x, y, mode, amp0, center0, width0, ...
            ph_range, minw, bo);
    end

% =====================================================================
function checkParsFields(pars)
% Warn about pars fields that nothing in this pipeline reads.
%
% Catches two failure modes that are otherwise silent:
%   1. Version skew -- a caller copied from a machine running newer code sets
%      an option this file has never heard of.  The fit still runs, just with
%      a different algorithm than the caller asked for.
%   2. Typos and legacy names (pars.fitmode vs pars.fit.mode, pars.block_size
%      vs pars.block_align.size), which have the same silent-no-op effect.
%
% WARNING ONLY, by design.  The schema below is a hand-maintained mirror of
% what this file and the fitters it calls consume, so it will drift as options
% are added. If you add a pars option, add it here too.
%
% Sub-structs handed wholesale to other routines (WatSupPre.opts,
% WatSupPost.opts) are deliberately NOT descended into: their option sets
% belong to those routines, not here.

baselineFields = {'enable','knotSpacing','lambda','lambdaAmpl','lambdaPrior', ...
                  'lambdaWidth','TolFun','TolX','linearTol','normalize','verbose', ...
                  'weights','amplUB','widthBounds','widthBoundsHz','centerBounds', ...
                  'phaseBounds','phaseInit','priorBaseline','priorWidth','etaInit', ...
                  'widthPriorPeaks'};

KNOWN = { ...
  '', {'plt','pltSpec','removeOS','lb','lb_met','lb_wat','den','eccopt', ...
       'WatSupPre','WatSupPost','peaks','ccopt','block_align','base','hsvd', ...
       'PC','PCpivot','dofit','fit','flag','watfit','waterScaling'}; ...
  'fit', {'mode','ph_range','ppm_range','ppm_range_coarse','baselineOpt','regions', ...
          'coarseMode','coarseFirstOrder','hsvdCoarse','hsvdClean','lineShapeOpt','focusRange', ...
          'focusWeight','broadWidthBounds','broadWidthInit','broadPhaseBounds','peaks'}; ...
  'fit.baselineOpt',         baselineFields; ...
  'fit.regions',             {'fitRange','peaks','name','baselineOpt','peakSelect'}; ...
  'fit.regions.baselineOpt', baselineFields; ...
  'fit.hsvdCoarse',          {'enable','K','npts','lwMinHz','lwMaxHz','phaseTolDeg'}; ...
  'fit.hsvdClean',           {'enable','K','range','nuisance_freqHz'}; ...
  'fit.lineShapeOpt',        {'enable','nSide','asymmetric','maxSide','init', ...
                              'normalize','areaPadFactor'}; ...
  'peaks',                   {'name','range','widthInit'}; ...
  'watfit',                  {'method','mode','nComp','ppm_range','ph_range', ...
                              'baselineOpt','lineShapeOpt','etaInit'}; ...
  'watfit.baselineOpt',      baselineFields; ...
  'block_align',             {'size','method','Nfit','ppmRange','saveDebug'}; ...
  'ccopt',                   {'minsig_frac'}; ...
  'hsvd',                    {'flag','sv','cad'}; ...
  'den',                     {'MLSingularValues'}; ...
  'flag',                    {'baseCorr'}; ...
};

% Struct-array elements (pars.fit.regions) all carry the same field set, so a
% stray field is seen once per element -- report it once.
bad = unique(scanParsLevel(pars, '', KNOWN, {}), 'stable');

if ~isempty(bad)
    warning('dat2mat_svsNAD:unknownParsField', ...
        ['%d pars field(s) are not read by this code and have NO effect:\n' ...
         '    %s\n' ...
         'If you expected these to do something, this copy of dat2mat_svsNAD.m ' ...
         'is probably older than the caller.'], ...
        numel(bad), strjoin(bad, sprintf('\n    ')));
end

% =====================================================================
function bad = scanParsLevel(s, prefix, KNOWN, bad)
% Recurse one level of pars, collecting field paths absent from the schema.
% Levels with no schema entry are not checked at all (their options belong to
% some other routine), so a missing entry can only ever under-report.

    if ~isstruct(s) || isempty(s), return; end

    idx = find(strcmp(KNOWN(:,1), prefix), 1);
    if isempty(idx), return; end
    allowed = KNOWN{idx,2};

    fn = fieldnames(s);
    for k = 1:numel(fn)
        if isempty(prefix)
            here = fn{k};
        else
            here = [prefix '.' fn{k}];
        end

        if ~any(strcmp(fn{k}, allowed))
            bad{end+1} = ['pars.' here]; %#ok<AGROW>
            continue
        end

        % Struct arrays (pars.fit.regions) get every element checked.
        v = s(1).(fn{k});
        if isstruct(v)
            for e = 1:numel(v)
                bad = scanParsLevel(v(e), here, KNOWN, bad);
            end
        end
    end
