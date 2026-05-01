function [output,pars,fullHdr] = dat2mat_svs31p(pars,dat_file)
%
% [output,pars,fullHdr] = dat2mat_svs31p(pars,dat_file)
%
% Steps:
% block averaging
% filtering
% hsvd denoising
% Fourier Transform
% phase correction
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
if ~isfield(pars,'block_size'), pars.block_size = []; end
if ~isfield(pars,'zf') || isempty(pars.zf), pars.zf = false; end
if ~isfield(pars,'peaks'), pars.peaks = [0 5]; end % PCr, Pi
if ~isfield(pars,'peakRanges'), pars.peakRanges = []; end
if ~isfield(pars,'peakAddLB'), pars.peakAddLB = 3; end
if ~isfield(pars,'PC'), pars.PC = 0; end
if ~isfield(pars,'base'), pars.base = nan; end
if ~isfield(pars,'hsvd') || ~isfield(pars.hsvd,'flag'), pars.hsvd.flag = false; end
if ~isfield(pars.hsvd,'sv'), pars.hsvd.sv = []; end % HSVD singular values
if ~isfield(pars.hsvd,'cad'), pars.hsvd.cad = 3; end % HSVD Cadzow iterations
if ~isfield(pars,'dofit'), pars.dofit = false; end
if ~isfield(pars,'blp'), pars.blp.flag = false; end % backwards linear prediction for TE delay
if ~isfield(pars.blp,'modelOrder'), pars.blp.modelOrder = 0; end % auto select
if ~isfield(pars,'block_align') || ~isfield(pars.block_align,'size'),   pars.block_align.size = 64; end
if ~isfield(pars,'block_align') || ~isfield(pars.block_align,'method'), pars.block_align.method = 'fd'; end % 'td' (fidfreqshift) or 'fd' (fidfreqshift_fd)
if ~isfield(pars.block_align,'Nfit'),     pars.block_align.Nfit = 256; end       % TD: # FID points fit
if pars.dofit
    if ~isfield(pars.fit,'mode'), pars.fit.mode = 2; end % complex Lorentzian
end

if nargin<2
    [file,path] = uigetfile('*.dat','Choose dat file');
    dat_file = fullfile(path,file);
end


%% read data
twix_obj = mapVBVD(dat_file);
if iscell(twix_obj)
    twix_obj = twix_obj{end};
end

twix_obj.image.flagRemoveOS = false;
twix_obj.image.flagAverageSets = false;

fullHdr = twix_obj.hdr;

if isfield(twix_obj.hdr.Config,'DwellTime')
    bw = 1e9/twix_obj.hdr.Config.DwellTime;
elseif isfield(twix_obj.hdr.Config,'DwellTimeSig')
    bw = 1e9/twix_obj.hdr.Config.DwellTimeSig;
else
    error('unknown dwell time')
end
f0 = 1e-6*twix_obj.hdr.Config.Frequency;
nch = twix_obj.image.NCha; % should be 1
if ~isequal(nch,1)
    warning(['Number of coils is ' num2str(nch) '. Should be 1.'])
end
TEus = fullHdr.MeasYaps.alTE{1};

nave = twix_obj.image.NAve;
npts = twix_obj.image.NCol;
reps = twix_obj.image.NRep;
avgDim = find(strcmp(twix_obj.image.sqzDims,'Ave'));
repDim = find(strcmp(twix_obj.image.sqzDims,'Rep'));

seqname = twix_obj.hdr.MeasYaps.tSequenceFileName;

wiplong = twix_obj.hdr.MeasYaps.sWipMemBlock.alFree;
empInd = cellfun('isempty',wiplong);
wiplong(empInd) = {0};
wiplong = cell2mat(wiplong);
wipdbl = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree;
empInd = cellfun('isempty',wipdbl);
wipdbl(empInd) = {0};
wipdbl = cell2mat(wipdbl);

fid = twix_obj.image{''};

si = size(fid);

% remove OS
if pars.removeOS
    bw = bw/2;
    npts = npts/2;
end

% zero fill
if pars.zf
    fid = cat(1,fid,0*fid);
    npts = 2*npts;
end

fid = reshape(fid,npts,[]); % makes fid 2D

% time axis
t = (0:npts-1)*1/bw;

% freq axis
hz = (-1/2:1/npts:1/2-1/npts)*bw;
ppm = hz/f0;

adcshift = fullHdr.Phoenix.sSpecPara.dDeltaFrequency;
ppm = ppm + adcshift;

if nave>1
    disp('block averaging')
    if isempty(pars.block_align.size)
        pars.block_align.size = size(fid,avgDim); % defaults to 1 block consisting of all the averages
    end
    nd = ndims(fid);
    if avgDim>2 % make avgDim 2nd
        if avgDim==nd
            fid = permute(fid,[1 avgDim 2:avgDim-1]);
        else
            fid = permute(fid,[1 avgDim 2:avgDim-1 avgDim+1:nd]);
        end
    end
    fid = blockAverage(fid,pars.block_align.size); % average dimension MUST be 2nd in ws here
    if avgDim>2 % put back in order
        if avgDim==nd
            fid = permute(fid,[1 avgDim 2:avgDim-1]);
        else
            fid = permute(fid,[1 avgDim 2:avgDim-1 avgDim+1:nd]);
        end
    end
end


disp('filtering')
fid = expFilter(t,pars.lb,fid);

if pars.hsvd.flag
    if ~isempty(pars.hsvd.sv)
        for ii=1:size(fid,2)
            fid(:,ii) = NWhsvd(fid(:,ii),pars.hsvd.sv,pars.hsvd.cad);
        end
    else
        f = NWhsvdden(fid,t,bw);
        waitfor(f);
        fid = out; clear out
    end
end

if pars.block_align.size < nave
    disp('block alignment and resolve averages')
    if isfield(pars.block_align,'saveDebug') && pars.block_align.saveDebug
        output.met.blocks_pre = fid;  % FIDs before alignment, per block
    end
    switch lower(pars.block_align.method)
        case 'fd'
            fid = fidfreqshift_fd(t, fid, ppm, pars.block_align.ppmRange);
        otherwise % 'td'
            fid = fidfreqshift(t, fid, pars.block_align.Nfit, 0);
    end
    if isfield(pars.block_align,'saveDebug') && pars.block_align.saveDebug
        output.met.blocks_post = fid; % FIDs after alignment, per block
    end
    fid = squeeze(mean(fid,2));
end


% output fid pre backwards linear prediction
output.rawfid = fid;

if pars.blp.flag
    % Correct for TE delay
    TEs = TEus * 1e-6;
    nPredict = round(TEs/(1/bw)); % points to backwards predict
    fid = blp_fid(fid,nPredict,pars.blp.modelOrder);
end


% Fourier Transform
spec = fftshift(fft(fid,[],1),1);

% peak-based frequency/phase correction
if ~isempty(pars.peaks)
    peaks = pars.peaks;
    if isempty(pars.peakRanges)
        pars.peakRanges = 0.5*ones(length(peaks),2);
    end
    ranges = pars.peakRanges;
    addlb = pars.peakAddLB;

    disp('peak-based frequency correction')
    % line-broadened magnitude spectrum for peak detection
    fidlb = expFilter(t,addlb,fid);
    specmag = abs(fftshift(fft(fidlb,[],1),1));

    hzshift = zeros(size(spec,2),1);
    for ii=1:size(specmag,2)
        zpsh = zeros(length(peaks),1);
        for jj=1:length(peaks)
            ind = find(ppm > (peaks(jj)-ranges(jj,1)) & ppm < (peaks(jj)+ranges(jj,2)));
            [~,zpind] = max(specmag(ind,ii));
            zpsh(jj) = ppm(ind(1)-1+zpind) - peaks(jj);
        end
        hzshift(ii) = mean(zpsh)*f0;
        spec(:,ii) = shiftSpectrumFrequency(spec(:,ii),hzshift(ii),t);
    end
    for ii=1:length(hzshift)
        disp(['      freq shift: ' num2str(hzshift(ii)) ' hz'])
    end

    disp('peak-based phase correction')
    % recompute line-broadened spectrum after frequency correction
    speclb = fftshift(fft(fidlb,[],1),1);
    specmag = abs(speclb);

    npeaks = length(peaks);
    if npeaks>1
        ph = zeros(size(spec));
    else
        ph = zeros(1,size(specmag,2));
    end
    for ii=1:size(specmag,2)
        peakphs = zeros(npeaks,1);
        peakinds = zeros(npeaks,1);
        for jj=1:npeaks
            ind = find(ppm > (peaks(jj)-ranges(jj,1)/2) & ppm < (peaks(jj)+ranges(jj,2)/2));
            [~,maxind] = findpeaks(specmag(ind,ii),'SortStr','descend','NPeaks',1);
            if isempty(maxind)
                [~,maxind] = max(specmag(ind,ii));
            end
            peakinds(jj) = maxind + ind(1) - 1;
            peakphs(jj) = angle(speclb(peakinds(jj),ii));
        end
        peakphs = unwrap(peakphs);
        if npeaks==1
            ph(ii) = peakphs;
        else
            p = polyfit(peakinds,peakphs,1);
            ph(:,ii) = polyval(p,(1:npts)');
        end
    end
    spec = shiftSpectrumPhase(spec,ph);

    if pars.plt
        figure, plot(ppm,real(sum(spec(:,:),2))), title('after peak-based correction'), set(gca,'xdir','reverse')
    end
end

% phase correction
if ~isfield(pars,'PC') || isempty(pars.PC)
    pivot0 = 0; 
    f = NWman_phase(spec,ppm,'Manual phase correction: save and close when done',pivot0);
    waitfor(f);
    if exist('out','var')
        spec = out; clear out
        output.flag.pc = true;
        pars.PC = parsPC;
    end
elseif isstruct(pars.PC)
    if ~isfield(pars.PC,'pivot')
        pars.PC.pivot = 0;
    end
    spec = shiftSpectrumPhase(spec,[pars.PC.pc0 pars.PC.pc1],ppm,pars.PC.pivot);
end


% switch pars.den
%     case 'svd'
%         % svd denoising
%         if size(spec,2)>1
%             f = NWsvdTS(spec,ppm,'SVD denoising: save and close when done');
%             waitfor(f);
%             spec = out; clear out
%         end
%     case 'wav'
%         % wavelet filtering
%         if exist('wmaxlev')
%             f = NWwavden(spec,ppm,'Wavelet denoising: save and close when done');
%             waitfor(f);
%             spec = out; clear out
%         end
% 
% end


% baseline correction
if isfield(pars,'base') && strcmpi(pars.base,'full')
    disp('baseline correction')
    f = NWsemiman_base(spec,ppm,'Semi manual baseline correction: save and close when done');
    waitfor(f);
    if exist('out','var')
        output.baseline = spec - out.spec;
        spec = out.spec;
        pars.flag.baseCorr = true;
        pars.base = out.pars;
        clear out
    end
end



% fitting
if pars.dofit

    if ischar(pars.fit.mode) && strcmpi(pars.fit.mode,'lcmodel') 

        % creates .RAW files for LCModel fitting
        rawfile = regexprep(dat_file,'datfile','rawfile');
        rawfile = regexprep(rawfile,'.dat','.raw');
        fid = ifft(ifftshift(spec,1),[],1);
        create_lcmodelRAW(rawfile,fid,'NUNFIL',length(fid),'DELTAT',1/bw,'HZPPPM',f0);

    elseif ischar(pars.fit.mode) && strcmpi(pars.fit.mode,'basisvarpro')

        % ----------------------------------------------------------
        % Basis-set varpro fit (LCModel/Osprey-style) in frequency
        % domain via curvefitAuto_basisVarpro.
        % ----------------------------------------------------------

        % default metabolite subset = NAD/UDP region
        if ~isfield(pars.fit,'metabs') || isempty(pars.fit.metabs)
            pars.fit.metabs = {'aATP','NADH','NADplus', ...
                               'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'};
        end
        if ~isfield(pars.fit,'ppm_range') || isempty(pars.fit.ppm_range)
            pars.fit.ppm_range = [-10.5 -6.5];
        end
        if ~isfield(pars.fit,'phase1Enable') || isempty(pars.fit.phase1Enable)
            pars.fit.phase1Enable = true;
        end
        if ~isfield(pars.fit,'softConstraints')
            pars.fit.softConstraints = [];
        end
        if ~isfield(pars.fit,'tieGroups')
            pars.fit.tieGroups = [];
        end

        % build basis FIDs at data length / dwell
        [basisFIDs, basisInfo] = make31P_basisFIDs(t(:), f0, ...
            'metabSubset', pars.fit.metabs, ...
            'adcshift',    adcshift);

        % fit input: prefer first-spectrum column for now
        x = ppm(:);
        y = double(spec(:,1));

        % --- timing / fit options
        % Anchor phi1's linear-phase ramp at whichever ppm we just used for
        % zero-order phase correction (so phi0 per component only has to
        % absorb small residuals, not the global phase).
        if isfield(pars.fit,'ppmPivot') && ~isempty(pars.fit.ppmPivot)
            ppmPivot = pars.fit.ppmPivot;
        elseif ~isempty(pars.peaks)
            ppmPivot = pars.peaks(1);
        else
            ppmPivot = 0;
        end

        timeInfo = struct( ...
            'dwellTime', 1/bw, ...
            'f0',        f0, ...
            't',         t(:), ...
            'ppmPivot',  ppmPivot);

        fitOpt = struct();
        fitOpt.fitRange       = pars.fit.ppm_range;
        fitOpt.complex        = true;
        fitOpt.phase1Enable   = pars.fit.phase1Enable;
        fitOpt.phaseBounds    = [-30 30];     % narrow; global phase comes from peak-based PC
        fitOpt.shiftBounds    = [-25 25];     % Hz, ~0.2 ppm at 7T
        fitOpt.lbLInit        = 20;           % typical 7T 31P brain LW
        fitOpt.lbLBounds      = [0 60];
        fitOpt.lbGInit        = 5;
        fitOpt.lbGBounds      = [0 30];
        fitOpt.softConstraints = pars.fit.softConstraints;
        fitOpt.tieGroups       = pars.fit.tieGroups;
        if isfield(pars.fit,'verbose'), fitOpt.verbose = pars.fit.verbose; end

        % allow user override of any fitOpt field via pars.fit.fitOpt
        if isfield(pars.fit,'fitOpt') && isstruct(pars.fit.fitOpt)
            fn = fieldnames(pars.fit.fitOpt);
            for ii_fn = 1:numel(fn)
                fitOpt.(fn{ii_fn}) = pars.fit.fitOpt.(fn{ii_fn});
            end
        end

        % --- baseline options (gentle spline; main job is broad envelope)
        baselineOpt = struct();
        baselineOpt.enable      = true;
        baselineOpt.style       = 'lcmodel';
        baselineOpt.knotSpacing = 1.5;        % ppm
        baselineOpt.lambda      = 1e-2;
        if isfield(pars.fit,'baselineOpt') && isstruct(pars.fit.baselineOpt)
            fn = fieldnames(pars.fit.baselineOpt);
            for ii_fn = 1:numel(fn)
                baselineOpt.(fn{ii_fn}) = pars.fit.baselineOpt.(fn{ii_fn});
            end
        end

        % --- run the fit
        outFit = curvefitAuto_basisVarpro(x, y, basisFIDs, basisInfo, ...
            timeInfo, fitOpt, baselineOpt);

        % --- report
        % Displayed amplitudes are scaled by n (number of 31P nuclei) so
        % that printed values are proportional to molecular concentration
        % and comparable across species with different n.  Raw amplitudes
        % remain in outFit.ampl(k).
        nVec = ones(length(outFit.names),1);
        if isfield(outFit,'basisInfo') && ~isempty(outFit.basisInfo)
            for k = 1:min(numel(outFit.basisInfo), length(outFit.names))
                if isfield(outFit.basisInfo(k),'n') && ~isempty(outFit.basisInfo(k).n) ...
                        && ~isnan(outFit.basisInfo(k).n) && outFit.basisInfo(k).n > 0
                    nVec(k) = outFit.basisInfo(k).n;
                end
            end
        end
        amplScaled = outFit.ampl ./ nVec;

        fprintf('\n*** basisVarpro 31P fit (ampl_sc shown; raw in outFit.ampl) ***\n');
        fprintf('    lbL_tot = basis lwHz (baked in) + fitted lbL; fwhmV_tot is the effective Voigt FWHM.\n');
        for k = 1:length(outFit.names)
            fprintf(['  %-12s amplSc= %-10.4g shift= %+5.2f Hz  ' ...
                     'lbL_tot= %6.1f Hz (=%.0f+%.1f)  lbG= %5.1f Hz  fwhmV_tot= %6.1f Hz  phi0= %+6.1f deg\n'], ...
                outFit.names{k}, amplScaled(k), outFit.pars.shift(k), ...
                outFit.pars.lbL_total(k), outFit.pars.lwHzBasis(k), outFit.pars.lbL(k), ...
                outFit.pars.lbG(k), outFit.pars.fwhmV_total(k), ...
                outFit.pars.phase0(k));
        end
        if outFit.fitOpt.phase1Enable
            fprintf('  global phi1 = %+5.2f deg/ppm\n', outFit.pars.phase1);
        end

        % --- totals (n-scaled, i.e., concentration-proportional)
        nadNames = {'NADH','NADplus'};
        udpNames = {'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'};
        totalNAD = 0;  nadFound = {};
        totalUDP = 0;  udpFound = {};
        for k = 1:length(outFit.names)
            if any(strcmp(outFit.names{k}, nadNames))
                totalNAD = totalNAD + amplScaled(k);
                nadFound{end+1} = outFit.names{k}; %#ok<AGROW>
            elseif any(strcmp(outFit.names{k}, udpNames))
                totalUDP = totalUDP + amplScaled(k);
                udpFound{end+1} = outFit.names{k}; %#ok<AGROW>
            end
        end
        fprintf('  -----\n');
        if isempty(nadFound)
            fprintf('  totalNAD     = (no NAD components in basis)\n');
        else
            fprintf('  totalNAD     = %.4g  (%s)\n', totalNAD, strjoin(nadFound,'+'));
        end
        if isempty(udpFound)
            fprintf('  totalUDP     = (no UDP components in basis)\n');
        else
            fprintf('  totalUDP     = %.4g  (%s)\n', totalUDP, strjoin(udpFound,'+'));
        end

        outFit.amplScaled            = amplScaled;
        outFit.nVec                  = nVec;
        outFit.totals.NAD            = totalNAD;
        outFit.totals.UDP            = totalUDP;
        outFit.totals.NAD_components = nadFound;
        outFit.totals.UDP_components = udpFound;

        if pars.plt
            [~,fname,~] = fileparts(dat_file);
            plot_basisvarpro_31p(outFit, fname);
            plot_basisSpectra_31p(outFit, fname);
        end

        output.basisVarproResults = outFit;

    elseif ischar(pars.fit.mode) && strcmpi(pars.fit.mode,'amares')
        
        % OXSA->AMARES
        addpath(genpath('Z:\code\OXSA'))

        % Pass spectrum and let AMARES reconstruct the FID via specInvFft
        % (which correctly handles the 2x t=0 point convention).
        % Do NOT pass signals — let AMARES use the spectra path.
        amaresStruct.samples = size(spec,1);
        amaresStruct.dwellTime = 1 / bw;
        amaresStruct.imagingFrequency = f0;
        amaresStruct.timeAxis = t(:);
        amaresStruct.ppmAxis = ppm(:);
        amaresStruct.signals = {[]};              % empty so AMARES uses spectra path
        amaresStruct.spectra = {spec(:,1)};

        % DEBUG
        % debug_pk;

        % Fit
        % beginTime = time from pulse centroid to first ADC sample (s)
        beginTime = TEus * 1e-6;

        % --- Stage 1: Fit major peaks with basic PK ---
        PK1 = AMARES.priorKnowledge.PK_7T_31P_Basic_NW();
        for iPK = 1:length(PK1.initialValues)
            PK1.initialValues(iPK).chemShift = PK1.initialValues(iPK).chemShift - adcshift;
            if ~isempty(PK1.bounds(iPK).chemShift)
                PK1.bounds(iPK).chemShift = PK1.bounds(iPK).chemShift - adcshift;
            end
        end
        Results1 = AMARES.amares(amaresStruct,1,1,beginTime,0.0,PK1,0,'fixOffset',0);

        % --- Stage 2: Fit full model, seeding amplitudes from Stage 1 ---
        PK2 = AMARES.priorKnowledge.PK_7T_31P_Ren_NW();
        for iPK = 1:length(PK2.initialValues)
            PK2.initialValues(iPK).chemShift = PK2.initialValues(iPK).chemShift - adcshift;
            if ~isempty(PK2.bounds(iPK).chemShift)
                PK2.bounds(iPK).chemShift = PK2.bounds(iPK).chemShift - adcshift;
            end
        end

        % Seed matching peaks from Stage 1 results
        basicNames = {'PCr','ATPb','ATPa','ATPg','Pi','PE','PC','GPE','GPC'};
        for iPK = 1:length(PK2.initialValues)
            pn = PK2.initialValues(iPK).peakName;
            if iscell(pn), pn = pn{1}(1:end-1); else, pn = pn; end
            matchIdx = find(strcmp(basicNames, pn));
            if ~isempty(matchIdx)
                PK2.initialValues(iPK).amplitude = Results1.Amplitudes(matchIdx);
                PK2.initialValues(iPK).linewidth = Results1.Linewidths(matchIdx);
                PK2.initialValues(iPK).phase = mod(Results1.Phases(matchIdx), 360);
                PK2.initialValues(iPK).chemShift = Results1.ChemicalShifts(matchIdx);
            end
        end

        % For NAD/UDP, seed amplitudes from literature ratios relative
        % to PCr, and seed phase from Stage 1 average
        PCrAmp = Results1.Amplitudes(1);  % PCr is compound 1 in basic PK
        avgPhase = mod(mean(Results1.Phases), 360);
        for iPK = 1:length(PK2.initialValues)
            pn = PK2.initialValues(iPK).peakName;
            if iscell(pn), pn = pn{1}(1:end-1); else, pn = pn; end
            if ~any(strcmp(basicNames, pn))
                % NAD+ ~10-15% of PCr, NADH ~1-2% of PCr
                % UDP(G) total ~10% of PCr (Ren 2021 Table 1):
                %   UDPGal ~16%, UDPGlc ~30%, UDPGalNAc ~24%, UDPGlcNAc ~30%
                if strcmp(pn,'NADH')
                    PK2.initialValues(iPK).amplitude = PCrAmp * 0.015;
                elseif strncmp(pn,'NADp',4)
                    PK2.initialValues(iPK).amplitude = PCrAmp * 0.08;
                elseif strncmp(pn,'UDPS',4)
                    PK2.initialValues(iPK).amplitude = PCrAmp * 0.05;  % combined UDP
                elseif strncmp(pn,'UDPGal_',7)
                    PK2.initialValues(iPK).amplitude = PCrAmp * 0.016; % ~16% of total UDP
                elseif strncmp(pn,'UDPGlc_',7)
                    PK2.initialValues(iPK).amplitude = PCrAmp * 0.030; % ~30% of total UDP
                elseif strncmp(pn,'UDPGaNAc',8)
                    PK2.initialValues(iPK).amplitude = PCrAmp * 0.024; % ~24% of total UDP
                elseif strncmp(pn,'UDPGcNAc',8)
                    PK2.initialValues(iPK).amplitude = PCrAmp * 0.030; % ~30% of total UDP
                elseif strcmp(pn,'Broad10_5') || strcmp(pn,'Broad8_5') || strcmp(pn,'BroadMP')
                    PK2.initialValues(iPK).amplitude = PCrAmp * 0.03;
                else
                    PK2.initialValues(iPK).amplitude = PCrAmp * 0.015;
                end
                PK2.initialValues(iPK).phase = avgPhase;
            end
        end

        Results = AMARES.amares(amaresStruct,1,1,beginTime,0.0,PK2,0, ...
            'fixOffset',0,'skipLinearInit',true,'MaxIter',500,'MaxFunEvals',5000);

        if pars.plt
            plot_amares_full(Results1);
            plot_amares_full(Results);
            plot_amares_nadudp(Results);
        end

        output.AmaresResults = Results;

        rmpath(genpath('Z:\code\OXSA'))
       
    else

        if isfield(pars.fit,'ppm_range')
            inds = find( ppm > min(pars.fit.ppm_range) & ppm < max(pars.fit.ppm_range));
            x = ppm(inds);
            y = double(real(spec(inds)));
        else
            x = ppm;
            y = double(real(spec));
        end

        if strcmpi(pars.dofit,'auto')

            minw = 10/f0;
            [yfit,n,names,ampl,pos,width,integral,ip,pars] = curvefitAuto(x,y,minw,pars.fit.mode);

        elseif strcmpi(pars.dofit,'man')

            minw = 10/f0; % 10 Hz minimum linewidth
            [yfit,names,areas,ip,ip0,lb,ub] = curvefitMan(x,double(real(y)),minw,pars.fit.mode,pars.fit.peaks);

            if pars.plt
                figure, plot(x,real(y),'k',x,yfit,'r',x,real(y)-yfit,'g'), legend({'data','fit','residual'}), set(gca,'xdir','reverse')
            end

            % fit output
            output.fit.spec = y;
            output.fit.spec_fit = yfit;
            output.fit.pars = ip;
            output.fit.names = names;
            output.fit.areas = areas;
            output.fit.ppm = x;

        else
            warning('unknown fitting method. skipping...')
        end

    end
end


% output
output.fid = fid;
output.spec = spec;
output.time = t;
output.hz = hz;
output.ppm = ppm;


rmpath([thisPath filesep 'mapVBVD']);
rmpath([thisPath filesep 'utils']);
