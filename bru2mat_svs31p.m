function [output,pars,fullHdr] = bru2mat_svs31p(pars,bru_dir)
%
% [output,pars,fullHdr] = bru2mat_svs31p(pars,bru_dir)
%
% Bruker high-field 31P SVS processing. See dat2mat_svs31p.m.
%
% Steps:
% read Bruker data
% filtering
% hsvd denoising
% Fourier Transform
% phase correction
% baseline correction
% fitting



% in case current directory is not this one
thisFile = mfilename('fullpath');
thisPath = fileparts(thisFile);
addpath([thisPath filesep 'utils']);
addpath([thisPath filesep 'bruker']);


% defaults
if ~isfield(pars,'plt') || isempty(pars.plt), pars.plt = false; end
if ~isfield(pars,'removeOS') || isempty(pars.removeOS), pars.removeOS = false; end
if ~isfield(pars,'lb') || isempty(pars.lb), pars.lb = 0; end
if ~isfield(pars,'den') || isempty(pars.den), pars.den = false; end
if ~isfield(pars,'zf') || isempty(pars.zf), pars.zf = false; end
if ~isfield(pars,'peaks'), pars.peaks = [0 5]; end % PCr, Pi
if ~isfield(pars,'peakRanges'), pars.peakRanges = []; end
if ~isfield(pars,'peakAddLB'), pars.peakAddLB = 3; end
if ~isfield(pars,'PC'), pars.PC = 0; end
if ~isfield(pars,'initPhDeg'), pars.initPhDeg = 0; end % initial phase guess in degrees
if ~isfield(pars,'initPPMShift'), pars.initPPMShift = 0; end % initial shift guess in PPM
if ~isfield(pars,'base'), pars.base = nan; end
if ~isfield(pars,'hsvd') || ~isfield(pars.hsvd,'flag'), pars.hsvd.flag = false; end
if ~isfield(pars.hsvd,'sv'), pars.hsvd.sv = []; end % HSVD singular values
if ~isfield(pars.hsvd,'cad'), pars.hsvd.cad = 3; end % HSVD Cadzow iterations
if ~isfield(pars,'dofit'), pars.dofit = false; end
if pars.dofit
    if ~isfield(pars.fit,'mode'), pars.fit.mode = 2; end % complex Lorentzian
end

if nargin<2
    bru_dir = uigetdir('','Choose Bruker experiment directory');
end


%% read data

[bruParams, fid] = bruker_read(bru_dir);
fid = bruker_remove_digital_filter(bruParams, fid);
fid = fid.'; 
fullHdr = bruParams;
bw     = bruParams.acqus.SW_h;                    % Hz
f0     = bruParams.acqus.BF1;                     % MHz (31P at high field)
npts   = size(fid,1);
nave   = bruParams.acqus.NS;                       % number of scans (already summed in fid)
TRus   = bruParams.acqus.D(2) * 1e6;              % TR delay (s -> us)

adcshift = bruParams.acqus.O1;          
seqname = bruParams.acqus.PULPROG;
scanname = bruParams.acqus.EXP;

% Notes:
%   - Bruker `fid` is typically already sum-of-NS averages; if the
%     experiment was acquired as a 'ser' file we may need to mean over
%     the indirect dimension.
%   - DSPFVS/GRPDLY-based group delay must be removed before processing
%     (bruker_remove_digital_filter handles this).
% -------------------------------------------------------------------------

% if loaded as 2D (npts x nave), collapse to 1D average — no block averaging
if ~isempty(fid) && size(fid,2) > 1
    disp('averaging across scans')
    fid = mean(fid,2);
end

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

% conjugate for Bruker
% fid = conj(fid); 

% time axis
t = (0:npts-1)*1/bw;

% freq axis
hz = (-1/2:1/npts:1/2-1/npts)*bw;
% hz = (1/2:-1/npts:-1/2+1/npts)*bw;
ppm = hz/f0;
ppm = ppm + adcshift;

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


% output fid pre backwards linear prediction
output.rawfid = fid;

% Fourier Transform
spec = fftshift(fft(fid,[],1),1);

% initial shift and phase correction
if pars.initPPMShift~=0
    hzshift = pars.initPPMShift * f0;
    for ii=1:size(spec,2)
        spec(:,ii) = shiftSpectrumFrequency(spec(:,ii),hzshift,t);
    end
end
if pars.initPhDeg~=0
    spec = shiftSpectrumPhase(spec,pi/180 * pars.initPhDeg);
end

fid = ifft(ifftshift(spec,1),[],1);

% peak-based frequency/phase correction
if ~isempty(pars.peaks)
    peaks = pars.peaks;
    if isempty(pars.peakRanges)
        pars.peakRanges = 0.5*ones(length(peaks),2);
    end
    ranges = pars.peakRanges;
    addlb = pars.peakAddLB;

    disp('peak-based frequency correction')
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

    % --- TODO: high-field 31P fitting adjustments -----------------------
    % Starting point: dat2mat_svs31p uses basisVarpro and AMARES with
    % linewidths/bounds tuned for 7T. At higher field (e.g. 9.4/11.7 T):
    %   - rescale Hz-based bounds (lbL, lbG, shift) since chemical shift
    %     dispersion in Hz grows with f0
    %   - basis FIDs need to be regenerated at the new f0
    %   - AMARES PK_*_NW prior knowledge files may need a high-field variant
    % For now: leave the fit dispatch as a stub.
    % --------------------------------------------------------------------

    if ischar(pars.fit.mode) && strcmpi(pars.fit.mode,'lcmodel')

        rawfile = fullfile(bru_dir,'lcmodel.raw');
        fid_for_raw = ifft(ifftshift(spec,1),[],1);
        create_lcmodelRAW(rawfile,fid_for_raw,'NUNFIL',length(fid_for_raw),'DELTAT',1/bw,'HZPPPM',f0);

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
            plot_basisvarpro_31p(outFit, scanname);
            plot_basisSpectra_31p(outFit, scanname);
        end

        output.basisVarproResults = outFit;

    elseif ischar(pars.fit.mode) && strcmpi(pars.fit.mode,'amares')

        % --- TODO: high-field 31P AMARES prior knowledge -----------------
        warning('AMARES fit not yet wired for Bruker high-field data.')

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

            minw = 10/f0;
            [yfit,names,areas,ip,ip0,lb,ub] = curvefitMan(x,double(real(y)),minw,pars.fit.mode,pars.fit.peaks);

            if pars.plt
                figure, plot(x,real(y),'k',x,yfit,'r',x,real(y)-yfit,'g'), legend({'data','fit','residual'}), set(gca,'xdir','reverse')
            end

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


rmpath([thisPath filesep 'utils']);
rmpath([thisPath filesep 'bruker']);
