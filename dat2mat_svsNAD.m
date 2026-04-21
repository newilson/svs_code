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
if ~isfield(pars,'plt') || isempty(pars.plt), pars.plt = false; end % debugging plots
if ~isfield(pars,'pltSpec') || isempty(pars.pltSpec), pars.pltSpec = true; end % plot final spectrum
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

% time axis
t = (0:vectorSize-1)*1/bw;

% freq axis
hz = (-1/2:1/vectorSize:1/2-1/vectorSize)*bw;
if contains(seqname,'svssel_latest') && ~contains(seqname,'svssel_latest0')  % older versions up to svssel_latest021224 do not include this
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
    warning('Mismatch in reference voltage between WAT REF and MET scans');
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
if contains(seqname,'svssel_latest') && ~contains(seqname,'svssel_latest0')  % older versions up to svssel_latest021224 do not include this
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
            figure, plot(ppm, temp,'-k',ppm,temp2,'-r'), title('residual water removal'), set(gca,'xdir','reverse')
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
    figure, plot(ppm,real(ws)), title(fname), set(gca,'xdir','reverse')
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
        % default linewidth bounds: 9-36 Hz (de Graaf), applied uniformly to all peaks
        if ~isfield(pars.fit.baselineOpt,'widthBounds')
            pars.fit.baselineOpt.widthBounds = [10/f0 60/f0]; % ppm
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
                figure
                plot(ppm, real(ws), 'k', ppm, real(ws_clean), 'r')
                set(gca, 'xdir', 'reverse')
                legend({'original','HSVD cleaned'})
                title('HSVD nuisance removal')
            end

            ws = ws_clean;
        end

        % extract fitting range (used for both steps)
        if isfield(pars.fit,'ppm_range_coarse')
            inds = find( ppm > min(pars.fit.ppm_range_coarse) & ppm < max(pars.fit.ppm_range_coarse));
            x = ppm(inds);
            y = double(real(ws(inds)));
        elseif isfield(pars.fit,'ppm_range')
            inds = find( ppm > min(pars.fit.ppm_range) & ppm < max(pars.fit.ppm_range));
            x = ppm(inds);
            y = double(real(ws(inds)));
        else
            x = ppm;
            y = double(real(ws));
        end

        minw = 10/f0;
        npeaks = size(pars.peaks.range, 1);

        % ============================================================
        % Step 1: Coarse VARPRO fit (initialization + phase estimation)
        % ============================================================
        % Naive seeds: center of each peak range, default width, phase 0.
        % Wide phase bounds and relaxed tolerances — just needs to land
        % in the right basin to estimate phase for correction.
        disp('Step 1: Coarse VARPRO initialization...')
        center0 = mean(pars.peaks.range, 2);
        width0 = (30/f0) * ones(npeaks, 1);

        coarseBaseOpt = pars.fit.baselineOpt;
        coarseBaseOpt.centerBounds = pars.peaks.range;
        coarseBaseOpt.TolFun = 1e-3;
        coarseBaseOpt.TolX   = 1e-3;

        % per-peak amplitude upper bounds
        coarseAmplUB = zeros(npeaks, 1);
        for kk = 1:npeaks
            rng = pars.peaks.range(kk,:);
            localInds = x > min(rng) & x < max(rng);
            if any(localInds)
                coarseAmplUB(kk) = max(abs(real(y(localInds))));
            else
                coarseAmplUB(kk) = max(abs(real(y)));
            end
        end
        coarseBaseOpt.amplUB = coarseAmplUB;

        % wide phase bounds, no lineshape kernel, no focus weighting
        coarseFit = curvefitAuto_varproBaseline(x, y, pars.fit.mode, ...
            ones(npeaks,1), center0, width0, [-179 179], minw, ...
            coarseBaseOpt);

        % extract coarse parameters
        center0 = coarseFit.pars(npeaks + (1:npeaks))';
        if pars.fit.mode == 5 || pars.fit.mode == 7
            ph0 = coarseFit.pars(3*npeaks + (1:npeaks))';
            width0 = coarseFit.pars(5*npeaks + (1:npeaks))'; % total Voigt width
        elseif pars.fit.mode == 3 || pars.fit.mode == 4
            ph0 = coarseFit.pars(3*npeaks + (1:npeaks))';
            width0 = coarseFit.pars(2*npeaks + (1:npeaks))';
        else
            ph0 = zeros(npeaks, 1);
            width0 = coarseFit.pars(2*npeaks + (1:npeaks))';
        end

        disp('Coarse VARPRO parameters:')
        for kk = 1:npeaks
            fprintf('  %s: pos=%.3f ppm, width=%.4f ppm, phase=%.1f deg\n', ...
                pars.peaks.name{kk}, center0(kk), width0(kk), ph0(kk))
        end

        if pars.plt
            figure
            plot(x, real(y), 'k', x, coarseFit.fit, 'r', x, real(y) - coarseFit.fit, 'g')
            hold on, plot(x, coarseFit.baseline, 'b--')
            legend({'data','coarse fit','residual','baseline'})
            set(gca,'xdir','reverse')
            title('Step 1: Coarse VARPRO fit')
        end

        % ============================================================
        % Phase correction from coarse fit
        % ============================================================
        % Fit linear phase (pc0 + pc1) through per-peak phases vs ppm,
        % then apply to full complex spectrum before the fine fit.
        peakPositions_ppm = center0(:);
        peakPhases_deg = ph0(:);
        if npeaks == 1
            pc0 = peakPhases_deg;
            pc1 = 0;
        else
            p = polyfit(peakPositions_ppm, peakPhases_deg, 1);
            pc1_per_ppm = p(1); % deg/ppm
            pc0 = polyval(p, center_freq); % 0th order at pivot
            % shiftSpectrumPhase convention: pc1 is total 1st-order
            % phase across half the bandwidth
            pc1 = pc1_per_ppm * abs(ppm(end) - ppm(1)) / 2;
        end
        fprintf('Phase correction: pc0 = %.1f deg, pc1 = %.1f deg\n', pc0, pc1)

        ws = shiftSpectrumPhase(ws, [pc0 pc1], ppm, center_freq);
        output.met.fit.phaseCorr.pc0 = pc0;
        output.met.fit.phaseCorr.pc1 = pc1;
        output.met.fit.phaseCorr.pivot = center_freq;

        % re-extract fitting range from phase-corrected spectrum
        y = double(real(ws(inds)));

        % ============================================================
        % Step 2: Fine VARPRO fit on phase-corrected spectrum
        % ============================================================
        disp('Step 2: Fine VARPRO fit...')
        pars.fit.baselineOpt.centerBounds = pars.peaks.range;
        amp0 = ones(size(center0));

        % per-peak amplitude upper bounds from phase-corrected data
        amplUB = zeros(npeaks, 1);
        for kk = 1:npeaks
            rng = pars.peaks.range(kk,:);
            localInds = x > min(rng) & x < max(rng);
            if any(localInds)
                amplUB(kk) = max(abs(real(y(localInds))));
            else
                amplUB(kk) = max(abs(real(y)));
            end
        end
        pars.fit.baselineOpt.amplUB = amplUB;

        % build focus weight vector for Step 2
        fitWeights = ones(size(x(:)));
        if ~isempty(pars.fit.focusRange) && pars.fit.focusWeight > 1
            focusInds = x > min(pars.fit.focusRange) & x < max(pars.fit.focusRange);
            fitWeights(focusInds) = pars.fit.focusWeight;
            fprintf('Focus weighting: %.1fx on [%.1f - %.1f] ppm (%d/%d points)\n', ...
                pars.fit.focusWeight, min(pars.fit.focusRange), max(pars.fit.focusRange), ...
                sum(focusInds), length(x))
        end
        pars.fit.baselineOpt.weights = fitWeights;

        % phase bounds straddle 0 — spectrum is now corrected
        if isfield(pars.fit,'lineShapeOpt')
            fitOutput = curvefitAuto_varproBaseline(x, y, pars.fit.mode, ...
                amp0, center0, width0, pars.fit.ph_range, minw, ...
                pars.fit.baselineOpt, pars.fit.lineShapeOpt);
        else
            fitOutput = curvefitAuto_varproBaseline(x, y, pars.fit.mode, ...
                amp0, center0, width0, pars.fit.ph_range, minw, ...
                pars.fit.baselineOpt);
        end

        if pars.plt
            figure
            plot(x, real(y), 'k', x, fitOutput.fit, 'r', x, real(y) - fitOutput.fit, 'g')
            hold on
            plot(x, fitOutput.baseline, 'b--')
            legend({'data','fit','residual','baseline'})
            set(gca,'xdir','reverse')
            title('Step 2: Fine VARPRO fit')
        end

        % store output
        output.met.fit.spec_full = ws;
        output.met.fit.spec = y;
        output.met.fit.spec_fit = fitOutput.fit;
        output.met.fit.baseline = fitOutput.baseline;
        output.met.fit.comp = fitOutput.comp;
        output.met.fit.pars = fitOutput.pars;
        output.met.fit.parnames = fitOutput.parnames;
        output.met.fit.areas = fitOutput.areas;
        output.met.fit.names = pars.peaks.name;
        output.met.fit.ppm = x;
        output.met.fit.lb = fitOutput.lb;
        output.met.fit.ub = fitOutput.ub;
        output.met.fit.init.center0 = center0;
        output.met.fit.init.width0 = width0;
        output.met.fit.init.ph0 = peakPhases_deg;
        output.met.fit.coarse.fit = coarseFit.fit;
        output.met.fit.coarse.baseline = coarseFit.baseline;
        output.met.fit.weights = fitWeights;

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
                    figure
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
                    figure
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
                output.wat.fit.areaTotal = sum(watFitOutput.areas);
                output.wat.fit.names = wat_peaks_name;
                output.wat.fit.ppm = xw;
                output.wat.fit.lb = watFitOutput.lb;
                output.wat.fit.ub = watFitOutput.ub;
                fprintf('  Water area (fit): %.4g\n', output.wat.fit.areaTotal)

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
