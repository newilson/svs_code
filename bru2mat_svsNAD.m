function [output,pars,fullHdr] = bru2mat_svsNAD(pars,ws_dir,nws_dir)
%
% [output,pars,fullHdr] = bru2mat_svsNAD(pars,ws_dir,nws_dir)
%
% Bruker high-field 1H NAD downfield SVS processing.  Mirrors
% dat2mat_svsNAD but removes coil combination and block averaging
% (single-channel high-field data is typically a single FID already summed
% over NS).  The fitting branch (varpro/auto/man/lcmodel) is the same
% multi-region VARPRO + watfit pipeline as the Siemens function.
%
% Inputs
%   pars     - struct of processing/fitting options (see dat2mat_svsNAD.m
%              for full field documentation). 
%   ws_dir   - Bruker experiment directory holding the water-suppressed
%              acquisition (contains fid/ser + acqus)
%   nws_dir  - Bruker experiment directory for the water reference
%              (optional; pass [] to skip)
%
% Outputs
%   output   - struct with .met / .wat sub-structs (spec, fid, fit, ...)
%   pars     - input pars with any auto-filled fields written back
%   fullHdr  - WS acqus parameter struct (bruker_read output)
%
% Steps:
%   read Bruker WS / NWS data (single coil)
%   eddy current correction
%   filtering
%   residual water removal (pre)
%   (no coil combination; no block averaging; no block alignment)
%   residual water removal (post)
%   Fourier Transform
%   phase correction
%   additional denoising
%   baseline correction
%   fitting (varpro / auto / man / lcmodel)



% in case current directory is not this one
thisFile = mfilename('fullpath');
thisPath = fileparts(thisFile);
addpath([thisPath filesep 'utils']);
addpath([thisPath filesep 'bruker']);

% defaults
if ~isfield(pars,'plt') || isempty(pars.plt), pars.plt = false; end
if ~isfield(pars,'pltSpec') || isempty(pars.pltSpec), pars.pltSpec = true; end
if ~isfield(pars,'removeOS') || isempty(pars.removeOS), pars.removeOS = false; end
if ~isfield(pars,'lb') || isempty(pars.lb), pars.lb = 0; end
if ~isfield(pars,'lb_met') || isempty(pars.lb_met), pars.lb_met = pars.lb; end
if ~isfield(pars,'lb_wat') || isempty(pars.lb_wat), pars.lb_wat = pars.lb; end
if ~isfield(pars,'den') || isempty(pars.den), pars.den = false; end
if ~isfield(pars,'eccopt') || isempty(pars.eccopt), pars.eccopt = -1; end % Brown method
if ~isfield(pars,'WatSupPost') || isempty(pars.WatSupPost) || isempty(pars.WatSupPost.opts.type), pars.WatSupPost.opts.type = 'none'; end
if ~isfield(pars,'WatSupPre')  || isempty(pars.WatSupPre)  || isempty(pars.WatSupPre.opts.type),  pars.WatSupPre.opts.type  = 'none'; end
if ~isfield(pars,'peaks'), pars.peaks = []; end
if ~isfield(pars,'base'), pars.base = nan; end
if ~isfield(pars,'hsvd') || ~isfield(pars.hsvd,'flag'), pars.hsvd.flag = false; end
if ~isfield(pars.hsvd,'sv'), pars.hsvd.sv = []; end
if ~isfield(pars.hsvd,'cad'), pars.hsvd.cad = 3; end
if ~isfield(pars,'PC'), pars.PC = 0; end
if ~isfield(pars,'initPhDeg')    || isempty(pars.initPhDeg),    pars.initPhDeg    = 0; end % initial zero-order phase guess (degrees)
if ~isfield(pars,'initPPMShift') || isempty(pars.initPPMShift), pars.initPPMShift = 0; end % initial chemical-shift offset (ppm)
if ~isfield(pars,'dofit'), pars.dofit = false; end
if pars.dofit
    if ~isfield(pars.fit,'mode'), pars.fit.mode = 3; end % complex Lorentzian
end

if nargin<2
    ws_dir  = uigetdir('','Choose WS Bruker experiment directory');
    nws_dir = uigetdir('','Choose NWS Bruker experiment directory');
elseif nargin<3
    nws_dir = [];
end

[~,fname,~] = fileparts(ws_dir);
fname = strrep(fname,'_',' ');

% single coil for Bruker high-field probes
nch = 1;

%% non water suppressed
if isempty(nws_dir)
    pars.eccopt = -1; % can't do eddy current correction without water reference data
    nws = [];

else

    [nwsParams, nws] = bruker_read(nws_dir);
    nws = bruker_remove_digital_filter(nwsParams, nws);
    nws = nws.';                                % col-major: (npts x nave)

    bw       = nwsParams.acqus.SW_h;            % Hz
    f0       = nwsParams.acqus.BF1;             % MHz, 1H Larmor
    npts     = size(nws,1);
    nave     = nwsParams.acqus.NS;              % NWS averages (sum-of-NS in fid)
    seqname  = nwsParams.acqus.PULPROG;
    TR_ms    = nwsParams.acqus.D(2) * 1e3;
    TE_ms    = 2 * nwsParams.acqus.D(17) * 1e3 + 2 * nwsParams.acqus.P(17) / 1e3 + nwsParams.acqus.P(13) / 1e3;                             % no TE field in 1D acqus
    rcv_gain = nwsParams.acqus.RG;
    center_freq_Hz = nwsParams.acqus.O1;
    center_freq = center_freq_Hz/f0;

    try
        scandate = datetime(nwsParams.acqus.DATE, 'ConvertFrom','posixtime');
    catch
        scandate = NaT;
    end

    % Normalize to Receiver Gain and averages
    nws = nws / rcv_gain / nave;

    % collapse 2D ser to 1D average (no block averaging)
    if ~isempty(nws) && size(nws,2) > 1
        disp('averaging NWS across scans')
        nws = mean(nws, 2);
    end

    % Bruker is typically NOT 2x oversampled like Siemens — leave as no-op
    if pars.removeOS
        warning('removeOS requested but Bruker data is usually not 2x oversampled — skipping.')
    end

    vectorSize = size(nws,1);

    t   = (0:vectorSize-1)/bw;
    hz  = (-1/2:1/vectorSize:1/2-1/vectorSize)*bw;
    ppm = (center_freq_Hz + hz)/f0;

    output.wat.ppm      = ppm;
    output.wat.hz       = hz;
    output.wat.center_freq_Hz = center_freq_Hz;
    output.wat.center_freq = center_freq; % ppm
    output.wat.seqname  = seqname;
    output.wat.TR_ms    = TR_ms;
    output.wat.TE_ms    = TE_ms;
    output.wat.nAvg     = nave;
    output.wat.scandate = scandate;
end

%% water suppressed

[wsParams, ws] = bruker_read(ws_dir);
ws = bruker_remove_digital_filter(wsParams, ws);
ws = ws.';

fullHdr = wsParams;

% consistency check vs NWS protocol
if exist('f0','var')
    if ~isequal(f0, wsParams.acqus.BF1) || ~isequal(bw, wsParams.acqus.SW_h) ...
            || ~isequal(npts, size(ws,1))
        warning('inconsistent protocols between NWS and WS Bruker datasets')
    end
end

bw       = wsParams.acqus.SW_h;
f0       = wsParams.acqus.BF1;
npts     = size(ws,1);
nave     = wsParams.acqus.NS;
seqname  = wsParams.acqus.PULPROG;
TR_ms    = wsParams.acqus.D(2) * 1e3;
TE_ms    = 2 * wsParams.acqus.D(17) * 1e3 + 2 * wsParams.acqus.P(17) / 1e3 + wsParams.acqus.P(13) / 1e3;  % no TE field in 1D acqus
rcv_gain = wsParams.acqus.RG;
center_freq_Hz = wsParams.acqus.O1;

try
    scandate = datetime(wsParams.acqus.DATE, 'ConvertFrom','posixtime');
catch
    scandate = NaT;
end

% Normalize to Receiver Gain and Averages
ws = ws / rcv_gain / nave;

% collapse averages — no block averaging
if ~isempty(ws) && size(ws,2) > 1
    disp('averaging WS across scans')
    ws = mean(ws, 2);
end

if pars.removeOS
    warning('removeOS requested but Bruker data is usually not 2x oversampled — skipping.')
end

vectorSize = size(ws,1);

t   = (0:vectorSize-1)/bw;
hz  = (-1/2:1/vectorSize:1/2-1/vectorSize)*bw;
ppm = (center_freq_Hz + hz)/f0;
center_freq = center_freq_Hz/f0;
% disp(['center freq = ' num2str(center_freq) ' ppm'])

% % Siemens-style WIP excitation fields don't exist on Bruker; leave as defaults
% exc_freq = center_freq;
% exc_bw   = NaN;
if ~isfield(pars,'block_align'), pars.block_align = struct(); end
if ~isfield(pars.block_align,'ppmRange'), pars.block_align.ppmRange = [8.5 10]; end

% output.met.exc_freq = exc_freq;
% output.met.exc_bw   = exc_bw;
output.met.ppm      = ppm;
output.met.hz       = hz;
output.met.center_freq_Hz = center_freq_Hz;
output.met.center_freq = center_freq; % ppm
output.met.seqname  = seqname;
output.met.TR_ms    = TR_ms;
output.met.TE_ms    = TE_ms;
output.met.nAvg     = nave;
output.met.scandate = scandate;

disp('eddy current correction')
if isempty(nws)
    pars.eccopt = -1;
    [ws, ~, ~] = eddyCurrentCorrection(ws, ws, pars.eccopt, nch);
else
    % Klose ECC (eccopt==0) divides out the entire reference phase including
    % the linear ramp 2*pi*f_water*t, so water is forced to 0 Hz on output.
    % Measure water's pre-ECC frequency from nws and restore it afterward so
    % the ppm axis (built from center_freq_Hz) keeps its meaning.
    if pars.eccopt == 0
        % Water resonates at 4.7 ppm; its offset from the transmitter is fixed
        % by ppm = (center_freq_Hz + hz)/f0, so f_water_Hz = 4.7*f0 - center_freq_Hz.
        watPPM = 4.7;
        f_water_Hz = watPPM*f0 - center_freq_Hz;
        fprintf('  Klose ECC: water expected at %.2f Hz (%.2f ppm); restoring after correction\n', ...
            f_water_Hz, watPPM);
    end
    [ws, nws, ~] = eddyCurrentCorrection(nws, ws, pars.eccopt, nch);
    if pars.eccopt == 0
        % shiftSpectrumFrequency multiplies FID by exp(-1i*hzshift*2*pi*t),
        % which shifts the spectrum by -hzshift Hz. To move water from 0 Hz
        % back to +f_water_Hz, pass hzshift = -f_water_Hz.
        ws_spec  = fftshift(fft(ws,  [], 1), 1);
        nws_spec = fftshift(fft(nws, [], 1), 1);
        ws_spec  = shiftSpectrumFrequency(ws_spec,  -f_water_Hz, t);
        nws_spec = shiftSpectrumFrequency(nws_spec, -f_water_Hz, t);
        ws  = ifft(ifftshift(ws_spec,  1), [], 1);
        nws = ifft(ifftshift(nws_spec, 1), [], 1);
        % output.met.ecc_f_water_Hz_restored = f_water_Hz;
        % output.wat.ecc_f_water_Hz_restored = f_water_Hz;
    end
end

disp('filtering')
ws = expFilter(t, pars.lb_met, ws);
if ~isempty(nws)
    nws = expFilter(t, pars.lb_wat, nws);
end

if ~strcmpi(pars.WatSupPre.opts.type,'none')
    wsopts = pars.WatSupPre.opts;
    disp('water signal removal (pre)')
    if center_freq==4.72
        ws = removeResidualWater(double(ws),wsopts);
    else
        freqShiftHz = f0*(4.72 - center_freq);
        ws = removeMetSignal(double(ws),freqShiftHz,t,wsopts);
    end
end

% --- no coil combination (single coil) ---
% --- no block averaging ------------------

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

% --- no block alignment (Bruker single-shot fid already summed) ---

% Zero-order phase discrepancy between nws and ws
phcorr = exp(-1i*squeeze(angle(ws(1,:,:,:))));
ws = ws .* repmat(phcorr,[size(ws,1) 1 1 1]);


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
ws  = fftshift(fft(ws,[],1),1);
nws = fftshift(fft(nws,[],1),1);

% initial shift / phase guess (mirrors bru2mat_svs31p).  Applied before the
% coarse VARPRO phase fit so the optimizer starts inside the right basin
% when the carrier or zero-order phase is far from default.
if pars.initPPMShift ~= 0
    hzshift = pars.initPPMShift * f0;
    for ii = 1:size(ws,2)
        ws(:,ii) = shiftSpectrumFrequency(ws(:,ii), hzshift, t);
    end
end
if pars.initPhDeg ~= 0
    ws = shiftSpectrumPhase(ws, pi/180 * pars.initPhDeg);
end


% phase correction
if ~isfield(pars,'PC') || isempty(pars.PC)
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
        if ~isfield(pars,'PCpivot')
            pars.PCpivot = 4.72;
        end
        ws = shiftSpectrumPhase(ws,pars.PC,ppm,pars.PCpivot);
    end
end

% Additional denoising
switch pars.den
    case 'svd'
        if size(ws,2)>1
            f = NWsvdTS(ws,ppm,'SVD denoising: save and close when done');
            waitfor(f);
            ws = out; clear out
        end
    case 'wav'
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

if pars.pltSpec
    figure('Name','1H processed spectrum'), plot(output.met.ppm,real(ws)), title(fname), set(gca,'xdir','reverse')
end


% output - pre scaling for fitting
output.met.spec = ws;
output.wat.spec = nws;

% =====================================================================
% Fitting  (mirrors dat2mat_svsNAD.m — full multi-region VARPRO + watfit
%           branch, plus auto/man/lcmodel fallbacks.  Bruker tweaks:
%           LCModel .raw path is <expDir>/lcmodel.raw, not regexp'd off a
%           datfiles/rawfiles convention.)
% =====================================================================
if pars.dofit
    ppm = output.met.ppm;

    if strcmpi(pars.dofit,'lcmodel')

        % creates .RAW files for LCModel fitting
        rawfile = fullfile(ws_dir,'lcmodel.raw');
        fid = ifft(ifftshift(ws,1),[],1);
        create_lcmodelRAW(rawfile,fid,'NUNFIL',length(fid),'DELTAT',1/bw,'HZPPPM',f0);

        if ~isempty(nws)
            rawfile = fullfile(nws_dir,'lcmodel.raw');
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

        % Optional data-driven re-seeding of peak center ranges (gated; used by
        % the blood 1H pipeline).  Tightens pars.peaks.range around peaks
        % detected in the processed spectrum so per-dataset frequency drift
        % can't strand a component in a valley (which the >=0 amplitude solve
        % then nulls).  Off by default => brain/calf pipelines unaffected.
        if isfield(pars.peaks,'autoSeed') && isstruct(pars.peaks.autoSeed) && ...
                isfield(pars.peaks.autoSeed,'enable') && pars.peaks.autoSeed.enable
            rng0 = pars.peaks.range;
            pars.peaks.range = autoSeedPeakRanges(real(ws), ppm, pars.peaks, pars.peaks.autoSeed);
            for kk = 1:numel(pars.peaks.name)
                if ~isequal(rng0(kk,:), pars.peaks.range(kk,:))
                    fprintf('autoSeed: %-8s range %.3f-%.3f -> %.3f-%.3f ppm\n', ...
                        pars.peaks.name{kk}, rng0(kk,1), rng0(kk,2), ...
                        pars.peaks.range(kk,1), pars.peaks.range(kk,2));
                end
            end
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

        % --- Broad component handling (consistent with the 31P pipeline) -----
        % Peaks whose name starts with 'Broad' are wide macromolecular/baseline
        % components (cf. 'BroadMP','Broad10_5' in processU01_Brain_31P).  They
        % get wide linewidth bounds and a broad initial width, and are excluded
        % from the coarse phase fit (they are nuisance signal, not phase
        % references).  Everything else keeps the default narrow bounds.
        npksAll = numel(pars.peaks.name);
        wb0 = pars.fit.baselineOpt.widthBounds;
        if size(wb0,1) == 1
            peakWB = repmat(wb0, npksAll, 1);   % expand global -> per-peak
        else
            peakWB = wb0;
        end
        isBroad = startsWith(pars.peaks.name(:), 'Broad');
        if isfield(pars.fit,'broadWidthBounds') && ~isempty(pars.fit.broadWidthBounds)
            broadWB = pars.fit.broadWidthBounds;
        else
            broadWB = [65/f0 180/f0];   % ~0.22-0.61 ppm: broader than the narrow
                                        % peaks but a local hump, not a window-wide
                                        % blob that competes with the spline
        end
        if isfield(pars.fit,'broadWidthInit') && ~isempty(pars.fit.broadWidthInit)
            broadW0 = pars.fit.broadWidthInit;
        else
            broadW0 = 150/f0;           % ~0.5 ppm init (cf. 31P 100 Hz)
        end
        % broad components are kept near-absorptive: a smooth macromolecular
        % baseline should not carry dispersive (phased) character, otherwise it
        % spreads broadly negative and the spline arches above the data to
        % cancel it (a spline<->broad degeneracy).
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
        % When unset, a single implicit region spans the coarse range with
        % every peak, reproducing the original single-region behaviour.
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
        pcInfo      = repmat(struct('pc0',0,'pc1',0,'pivot',0), nRegions, 1);
        coarseStore = cell(nRegions,1);
        wsP_store   = cell(nRegions,1);

        if strcmp(coarseMode,'union')
            allSel = unique(cat(2, regSel{:}));
            allSelC = allSel(~isBroad(allSel));      % broad peaks excluded from phase fit
            if isempty(allSelC), allSelC = allSel; end
            allRanges = cat(1, regionList.fitRange);
            unionRange = [min(allRanges(:)) max(allRanges(:))];
            disp('Step 1: Coarse VARPRO (union)...')
            boC = pars.fit.baselineOpt;  boC.widthBounds = peakWB(allSelC,:);
            wiC = []; if isfield(pars.peaks,'widthInit') && ~isempty(pars.peaks.widthInit), wiC = pars.peaks.widthInit(allSelC); end
            [wsP, c0, w0, p0, pc0, pc1, cf, pv] = nadCoarsePhase(ws, ppm, f0, center_freq, ...
                unionRange, pars.peaks.range(allSelC,:), pars.fit.mode, ...
                boC, minw, pars.plt, 'union', wiC);
            center0_all(allSelC) = c0;  width0_all(allSelC) = w0;  ph0_all(allSelC) = p0;
            allSelB = allSel(isBroad(allSel));        % seed broad peaks for the fine fit
            if ~isempty(allSelB)
                center0_all(allSelB) = mean(pars.peaks.range(allSelB,:), 2);
                width0_all(allSelB)  = broadW0;
                ph0_all(allSelB)     = 0;
            end
            for r = 1:nRegions
                pcInfo(r).pc0 = pc0;  pcInfo(r).pc1 = pc1;  pcInfo(r).pivot = pv;
                wsP_store{r}  = wsP;  coarseStore{r} = cf;
            end
            fprintf('Phase correction (union): pc0 = %.1f deg, pc1 = %.1f deg (pivot %.2f ppm)\n', pc0, pc1, pv)
        else
            disp('Step 1: Coarse VARPRO (per region)...')
            for r = 1:nRegions
                sel  = regSel{r};
                selC = sel(~isBroad(sel));        % broad peaks excluded from phase fit
                if isempty(selC), selC = sel; end
                lbl = regionLabel(regionList, r);
                boC = pars.fit.baselineOpt;  boC.widthBounds = peakWB(selC,:);
                wiC = []; if isfield(pars.peaks,'widthInit') && ~isempty(pars.peaks.widthInit), wiC = pars.peaks.widthInit(selC); end
                [wsP_r, c0, w0, p0, pc0, pc1, cf, pv] = nadCoarsePhase(ws, ppm, f0, center_freq, ...
                    regionList(r).fitRange, pars.peaks.range(selC,:), pars.fit.mode, ...
                    boC, minw, pars.plt, lbl, wiC);
                center0_all(selC) = c0;  width0_all(selC) = w0;  ph0_all(selC) = p0;
                pcInfo(r).pc0 = pc0;  pcInfo(r).pc1 = pc1;  pcInfo(r).pivot = pv;
                wsP_store{r}  = wsP_r;  coarseStore{r} = cf;
                selB = sel(isBroad(sel));         % seed broad peaks for the fine fit
                if ~isempty(selB)
                    center0_all(selB) = mean(pars.peaks.range(selB,:), 2);
                    width0_all(selB)  = broadW0;
                    ph0_all(selB)     = 0;
                end
                fprintf('  region %d (%s): pc0 = %.1f deg, pc1 = %.1f deg (pivot %.2f ppm)\n', r, lbl, pc0, pc1, pv)
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
            % per-peak phase bounds: broad components near-absorptive, narrow
            % peaks keep the region phase range (phR)
            pbR = repmat(phR(:).', numel(sel), 1);
            ib  = isBroad(sel);
            if any(ib), pbR(ib,:) = repmat(broadPB, sum(ib), 1); end
            baseR.phaseBounds = pbR;
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
                % mirroring the 31P basisVarpro plots
                rLbl    = regionLabel(regionList, r);
                pkNames = pars.peaks.name(sel);
                ttl     = sprintf('1H fit - region %d: %s', r, rLbl);
                fOv = plot_varpro_1h(x_r, y_r, outR, pkNames, ttl);
                set(fOv, 'Name', sprintf('1H varpro fit - region %d %s', r, rLbl));
                fGl = plot_componentSpectra_1h(x_r, y_r, outR, pkNames, ttl);
                set(fGl, 'Name', sprintf('1H component spectra - region %d %s', r, rLbl));
            end
        end

        % union window axis for combined storage
        allRanges = cat(1, regionList.fitRange);
        unionInds = find(ppm > min(allRanges(:)) & ppm < max(allRanges(:)));
        x = ppm(unionInds);

        % reconstruct a flat peak-ordered pars/lb/ub when all regions share
        % one lineshape mode (peak blocks only; per-region lineshape kernels
        % and full detail live in output.met.fit.region(r).fit)
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
        % Kernel-aware areas: integrates the convolved per-peak component
        % (see curvefitAuto_varproBaseline).  When pars.fit.lineShapeOpt is
        % enabled, `areas` (base-Voigt only) under-counts by the kernel sum;
        % `areasConv` is the right metric to feed absoluteQuant_svs against
        % a kernel-aware water area.
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
        output.met.fit.phaseCorr.pivot = pcInfo(1).pivot;
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
            output.met.fit.region(r).pivot      = pcInfo(r).pivot;
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

                % Initial peak width: default 0.1 ppm, matched to the
                % typical visible water FWHM at 500 MHz Bruker.  Bigger
                % seeds (e.g. 0.2) push the single-component Voigt to fit
                % the broad pedestal at the cost of the narrow core;
                % smaller seeds (e.g. 0.05) leave the optimizer stuck on a
                % too-narrow core that misses most of the area.
                % For multi-component fits (nComp > 1) the initial widths
                % are spread LOGARITHMICALLY from w0/4 to w0*4 so the
                % optimizer starts with one narrow core + one broad
                % pedestal rather than near-degenerate components all at
                % w0.  Override with a scalar (single w0 for all) or
                % vector (per-component widths) via pars.watfit.widthInit.
                if isfield(watfit,'widthInit') && ~isempty(watfit.widthInit)
                    if isscalar(watfit.widthInit)
                        w0       = watfit.widthInit;
                        width0_w = w0 * ones(1, nWatComp);
                    else
                        width0_w = watfit.widthInit(:).';
                        assert(numel(width0_w) == nWatComp, ...
                            'pars.watfit.widthInit must be scalar or length nComp=%d', nWatComp);
                        w0 = mean(width0_w);
                    end
                else
                    w0 = max(0.1, minw);
                    if nWatComp == 1
                        width0_w = w0;
                    else
                        width0_w = w0 * logspace(-log10(4), log10(4), nWatComp);
                    end
                end
                % center0_w: default splits multi-component seeds by +/-0.05
                % ppm (one peak per offset).  For a co-located narrow+broad
                % water model (both components ON the water line, separated
                % only by linewidth) override with a vector via
                % pars.watfit.centerInit so the optimizer does NOT seed two
                % distinct peaks.
                if isfield(watfit,'centerInit') && ~isempty(watfit.centerInit)
                    center0_w = watfit.centerInit(:).';
                    assert(numel(center0_w) == nWatComp, ...
                        'pars.watfit.centerInit must be length nComp=%d', nWatComp);
                    amp0_w    = ones(1, nWatComp);
                elseif nWatComp == 1
                    amp0_w    = 1;
                    center0_w = 4.72;
                else
                    offsets   = linspace(-0.05, 0.05, nWatComp);
                    amp0_w    = ones(1, nWatComp);
                    center0_w = 4.72 + offsets;
                end

                watBaseOpt = watfit.baselineOpt;
                watBaseOpt.TolFun = 1e-8;
                watBaseOpt.TolX   = 1e-8;
                % centerBounds: default is the full ppm_range for every
                % component (loose -- lets components wander into separate
                % peaks).  For a co-located narrow+broad model, supply tight
                % per-component bounds [nComp x 2] via pars.watfit.centerBounds
                % so both stay pinned on the water line.
                if isfield(watfit,'centerBounds') && ~isempty(watfit.centerBounds)
                    watBaseOpt.centerBounds = watfit.centerBounds;
                else
                    watBaseOpt.centerBounds = wat_peaks_range;
                end
                watBaseOpt.knotSpacing = max(watBaseOpt.knotSpacing, diff(watfit.ppm_range));
                localInds_w = xw > min(wat_peaks_range(:)) & xw < max(wat_peaks_range(:));
                if any(localInds_w)
                    watBaseOpt.amplUB = max(abs(real(yw(localInds_w))));
                else
                    watBaseOpt.amplUB = max(abs(real(yw)));
                end

                % Optional two-stage phase warm-start: when the lineshape
                % kernel is enabled, the kernel-side Jacobian columns can
                % keep the trust-region step in a basin where the phase
                % parameter never moves (kernel takes a bite of the
                % residual, residual changes shape, phase gradient drops
                % below the step-acceptance threshold).  Mirrors the
                % brain/calf coarse->fine pattern: stage 1 fits with the
                % kernel disabled to converge phase cleanly; phase is then
                % applied to the complex water spectrum slice, and stage 2
                % refits the pre-phased data with the kernel enabled and a
                % tight ph_range so phase has nowhere left to drift.
                % Opt-in only; complex-peak modes (3/4/5) only.
                doTwoStage = isfield(watfit,'twoStagePhase') && watfit.twoStagePhase ...
                            && isfield(watfit,'lineShapeOpt') && isfield(watfit.lineShapeOpt,'enable') ...
                            && watfit.lineShapeOpt.enable ...
                            && any(watfit.mode == [3 4 5]);

                if doTwoStage
                    % Stage 1: kernel disabled, phase free over watfit.ph_range.
                    % Loose tolerances (1e-3) match the metabolite coarse
                    % stage in nadCoarsePhase -- stage 1 just needs to land
                    % in the right phase basin, not converge tightly.
                    lsOff = watfit.lineShapeOpt;  lsOff.enable = false;
                    watBaseOpt_s1 = watBaseOpt;
                    watBaseOpt_s1.TolFun = 1e-3;
                    watBaseOpt_s1.TolX   = 1e-3;
                    stage1 = curvefitAuto_varproBaseline(xw, yw, watfit.mode, amp0_w, center0_w, width0_w, watfit.ph_range, minw, watBaseOpt_s1, lsOff);
                    % phase indices in output.pars: per peak, [amp, pos, width, ph]
                    nWat = nWatComp;
                    phIdx = 4*(1:nWat);    % the 'ph' column for each peak
                    ph_stage1 = stage1.pars(phIdx);
                    fprintf('  watfit two-stage: stage-1 ph = %s deg (warm-starting stage 2)\n', ...
                        mat2str(ph_stage1(:).', 4));

                    % Stage 2: kernel enabled, full watfit.ph_range, but
                    % phase initialized at the stage-1 converged value via
                    % baselineOpt.phaseInit so the optimizer starts in the
                    % right basin.  Data is NOT modified -- the original yw
                    % (and the original complex spectrum) are unchanged, so
                    % the plot, output.wat.fit.spec, and any downstream
                    % consumer see the raw water spectrum.
                    watBaseOpt_s2 = watBaseOpt;
                    watBaseOpt_s2.phaseInit = ph_stage1;
                    watFitOutput = curvefitAuto_varproBaseline(xw, yw, watfit.mode, amp0_w, center0_w, width0_w, watfit.ph_range, minw, watBaseOpt_s2, watfit.lineShapeOpt);

                    watFitOutput.phaseCorr.ph_stage1 = ph_stage1;
                elseif isfield(watfit,'lineShapeOpt')
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
                % areaTotal = canonical kernel-aware water area: sum of the
                % per-component areasConv (extended-axis trapz of the
                % convolved component, see curvefitAuto_varproBaseline.m,
                % padded by lineShapeOpt.areaPadFactor * max(width)).  This
                % recovers the slow tails -- critical for Lorentzian water
                % (mode 3) whose 1/x^2 wings leak well past the fit window.
                % areaTotalWin = in-window trapz of (fit - baseline) kept
                % for back-compat / sanity check (will under-count Lorentzian
                % by a few percent unless lineShapeOpt.areaPadFactor is
                % large).  areaBase = sum(ampl .* width) for legacy callers.
                output.wat.fit.areaTotal    = sum(output.wat.fit.areasConv);
                output.wat.fit.areaTotalWin = trapz(xw(:), watFitOutput.fit(:) - watFitOutput.baseline(:));
                output.wat.fit.areaBase     = sum(watFitOutput.areas);   % for reference
                output.wat.fit.names = wat_peaks_name;
                output.wat.fit.ppm = xw;
                output.wat.fit.lb = watFitOutput.lb;
                output.wat.fit.ub = watFitOutput.ub;
                if isfield(watFitOutput,'phaseCorr')
                    output.wat.fit.phaseCorr = watFitOutput.phaseCorr;
                end
                fprintf('  Water area (fit, ext-axis): %.4g  (in-window: %.4g, base-only: %.4g)\n', ...
                    output.wat.fit.areaTotal, output.wat.fit.areaTotalWin, output.wat.fit.areaBase)

            else
                error('Unknown watfit.method: %s. Use ''integrate'', ''fid'', or ''fit''.', watfit.method)
            end
        end

    elseif strcmpi(pars.dofit,'auto')

        % scale to water first
        pars.waterScaling = true;
        if ~isempty(nws)
            nws = nws / max(abs(nws(:)));
            ws  = ws  / max(abs(nws(:)));
        end

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
        if ~isempty(nws)
            nws = nws / max(abs(nws(:)));
            ws  = ws  / max(abs(nws(:)));
        end

        if isfield(pars.fit,'ppm_range')
            inds = find( ppm > min(pars.fit.ppm_range) & ppm < max(pars.fit.ppm_range));
            x = ppm(inds);
            y = double(real(ws(inds)));
        else
            x = ppm;
            y = double(real(ws));
        end

        if isfield(pars,'base') && strcmpi(pars.base,'fit')
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

        minw = 10/f0;
        [yfit,names,areas,ip,ip0,lb,ub] = curvefitMan(x,double(real(y)),minw,pars.fit.mode,pars.fit.peaks);

        if pars.plt
            figure('Name','1H metabolite fit (manual)'), plot(x,real(y),'k',x,yfit,'r',x,real(y)-yfit,'g'), legend({'data','fit','residual'}), set(gca,'xdir','reverse')
        end

        output.met.fit.spec = y;
        output.met.fit.spec_fit = yfit;
        output.met.fit.pars = ip;
        output.met.fit.names = names;
        output.met.fit.areas = areas;
        output.met.fit.ppm = x;

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

rmpath([thisPath filesep 'bruker']);
rmpath([thisPath filesep 'utils']);


% =================================================================== %
% local helpers for the multi-region VARPRO pipeline
% (copied from dat2mat_svsNAD.m; keep in sync)
% =================================================================== %

function lbl = regionLabel(regionList, r)
% Region display label: its .name if set, else 'region<r>'.
    lbl = sprintf('region%d', r);
    if isfield(regionList, 'name') && ~isempty(regionList(r).name)
        lbl = regionList(r).name;
    end

function [wsP, center0, width0, ph0, pc0, pc1, coarseFit, pivotPPM] = nadCoarsePhase( ...
        ws, ppm, f0, center_freq, fitRange, peakRange, mode, baselineOpt, minw, plt, label, widthInit) %#ok<INUSL>
% Coarse VARPRO fit over a window + phase correction derived from its peaks.
% Returns the phase-corrected full spectrum (wsP), coarse seeds for the
% region's peaks (center0/width0/ph0), the applied phase (pc0/pc1) about
% pivot pivotPPM = mean(fitRange), and the coarse fit struct.  Naive seeds:
% center of each peak range, default width, phase 0; wide phase bounds and
% relaxed tolerances just to land in the right basin for phase estimation.
% (center_freq is accepted for backwards-compatible caller signatures but is
% no longer used as the pivot.)
%
% Optional widthInit (length np vector): per-peak override of the default
% width0 = 30/f0 seed.  Entries that are <=0 or NaN keep the default.  Use
% when a broad peak (e.g. Trp on blood) would otherwise park narrow during
% coarse and corrupt the polyfit-derived phase trend.
    inds = find(ppm > min(fitRange) & ppm < max(fitRange));
    x = ppm(inds);
    y = double(real(ws(inds)));
    np = size(peakRange, 1);
    center0 = mean(peakRange, 2);
    width0  = (30/f0) * ones(np, 1);
    if nargin >= 12 && ~isempty(widthInit)
        wi = widthInit(:);
        assert(numel(wi) == np, 'widthInit length must match peakRange rows (%d)', np);
        m = wi > 0 & ~isnan(wi);
        width0(m) = wi(m);
    end

    bo = baselineOpt;
    bo.centerBounds = peakRange;
    bo.TolFun = 1e-3;
    bo.TolX   = 1e-3;
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

    % Phase pivot = centroid of the fit window, NOT the carrier.  Pivoting
    % a polyfit-derived linear phase at the carrier (~4.72 ppm) amplifies
    % any small pc1 error fivefold by the time it reaches NAD at ~9 ppm.
    % Pivoting inside the region makes pc1 the residual across the region,
    % not a slew across the whole spectrum.
    pivotPPM = mean(fitRange);
    if np == 1
        pc0 = ph0;
        pc1 = 0;
    else
        p   = polyfit(center0(:), ph0(:), 1);
        pc0 = polyval(p, pivotPPM);                % 0th order at pivot
        % shiftSpectrumPhase convention: pc1 is total 1st-order phase
        % across half the bandwidth
        pc1 = p(1) * abs(ppm(end) - ppm(1)) / 2;
    end
    wsP = shiftSpectrumPhase(ws, [pc0 pc1], ppm, pivotPPM);

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
