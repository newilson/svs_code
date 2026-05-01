function [output,pars,fullHdr] = bru2mat_svsNAD(pars,ws_dir,nws_dir)
%
% [output,pars,fullHdr] = bru2mat_svsNAD(pars,ws_dir,nws_dir)
%
% Bruker high-field 1H NAD downfield SVS processing. Mirrors
% dat2mat_svsNAD but removes coil combination and block averaging
% (single-channel high-field data) and adjusts the fitting.
%
% Steps:
% read Bruker WS / NWS data
% eddy current/first order correction
% filtering
% residual water removal (pre)
% (no coil combination; no block averaging)
% residual water removal (post)
% Fourier Transform
% phase correction
% baseline correction
% fitting



% in case current directory is not this one
thisFile = mfilename('fullpath');
thisPath = fileparts(thisFile);
addpath([thisPath filesep 'utils']);

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


%% non water suppressed
if isempty(nws_dir)
    pars.eccopt = -1; % can't do eddy current correction without water reference data
    nws = [];

else

    % --- TODO: Bruker NWS read (helper functions in utils/bruker_*) -----
    % Pseudo-code:
    %   [nwsParams, nws] = bruker_read(nws_dir);
    %   nws = bruker_remove_digital_filter(nws, nwsParams);
    %   bw    = nwsParams.acqus.SW_h;
    %   f0    = nwsParams.acqus.BF1;       % MHz, 1H Larmor
    %   npts  = size(nws,1);
    %   nave  = nwsParams.acqus.NS;        % NWS averages (sum-of-NS in fid)
    %   seqname = nwsParams.acqus.PULPROG;
    %   center_freq = nwsParams.acqus.O1 / nwsParams.acqus.BF1 + 4.72;  % verify
    %   refvolt = [];                      % no Siemens-style refvolt for Bruker
    %   scandate = NaT;                    % parse from acqus.DATE if needed
    % Notes:
    %   - Bruker single-coil 1H — nch == 1, no coil combination needed.
    %   - If 'ser' present, mean across the indirect dim (no block averaging).
    % --------------------------------------------------------------------

    bw          = NaN;   % placeholder
    f0          = NaN;   % placeholder
    npts        = NaN;   % placeholder
    nave        = NaN;   % placeholder
    nch         = 1;     % single coil
    seqname     = '';    % placeholder
    center_freq = 4.72;  % water reference; verify from O1/BF1
    nws         = [];    % placeholder
    scandate    = NaT;
    TR_ms       = NaN;
    TE_ms       = NaN;

    % collapse to 1D average if loaded as (npts x nave)
    if ~isempty(nws) && size(nws,2) > 1
        nws = mean(nws,2);
    end

    % remove OS — Bruker typically does NOT have Siemens-style 2x OS; leave as no-op
    if pars.removeOS
        warning('removeOS requested but Bruker data is usually not 2x oversampled — skipping.')
    end

    vectorSize = npts;

    t = (0:vectorSize-1)*1/bw;
    hz = (-1/2:1/vectorSize:1/2-1/vectorSize)*bw;
    ppm = center_freq + hz/f0;

    output.wat.ppm = ppm;
    output.wat.hz = hz;
    output.wat.seqname = seqname;
    output.wat.TR_ms = TR_ms;
    output.wat.TE_ms = TE_ms;
    output.wat.nAvg = nave;
    output.wat.scandate = scandate;
end

%% water suppressed

% --- TODO: Bruker WS read (helper functions in utils/bruker_*) ----------
% Pseudo-code:
%   [wsParams, ws] = bruker_read(ws_dir);
%   ws = bruker_remove_digital_filter(ws, wsParams);
%   fullHdr = wsParams;
%   sanity-check that bw/f0/npts match the NWS dataset.
% ------------------------------------------------------------------------

fullHdr     = struct();   % placeholder
ws          = [];         % placeholder
% on the WS side, also pull from wsParams:
%   bw, f0, npts, nave, seqname, scandate, TR_ms, TE_ms, center_freq
%   exc_freq, exc_bw if available

if exist('f0','var') && ~isnan(f0)
    % consistency checks would go here (compare to ws params)
end

% collapse averages — no block averaging
if ~isempty(ws) && size(ws,2) > 1
    ws = mean(ws,2);
end

if pars.removeOS
    warning('removeOS requested but Bruker data is usually not 2x oversampled — skipping.')
end

vectorSize = size(ws,1);
if vectorSize == 0, vectorSize = npts; end

t = (0:vectorSize-1)*1/bw;
hz = (-1/2:1/vectorSize:1/2-1/vectorSize)*bw;
ppm = center_freq + hz/f0;

% placeholders for fields the Siemens version pulls from WIP memory
exc_freq = center_freq;
exc_bw   = NaN;
if ~isfield(pars,'block_align'), pars.block_align = struct(); end
if ~isfield(pars.block_align,'ppmRange'), pars.block_align.ppmRange = [8.5 10]; end

output.met.exc_freq = exc_freq;
output.met.exc_bw   = exc_bw;
output.met.ppm      = ppm;
output.met.hz       = hz;
output.met.seqname  = seqname;
output.met.TR_ms    = NaN;
output.met.TE_ms    = NaN;
output.met.nAvg     = nave;
output.met.scandate = scandate;

disp('eddy current correction')
if isempty(nws)
    pars.eccopt = -1;
    [ws, ~, ~] = eddyCurrentCorrection(ws,ws,pars.eccopt,nch);
else
    [ws, nws, ~] = eddyCurrentCorrection(nws,ws,pars.eccopt,nch);
end

disp('filtering')
ws = expFilter(t,pars.lb_met,ws);
if ~isempty(nws)
    nws = expFilter(t,pars.lb_wat,nws);
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
            figure, plot(ppm, temp,'-k',ppm,temp2,'-r'), title('residual water removal'), set(gca,'xdir','reverse')
        else
            NWplayplot(temp2,ppm,'residual water removal',temp, ppm);
        end
    end
end

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
ws = fftshift(fft(ws,[],1),1);
nws = fftshift(fft(nws,[],1),1);


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
    figure, plot(ppm,real(ws)), title(fname), set(gca,'xdir','reverse')
end


% output - pre scaling for fitting
output.met.spec = ws;
output.wat.spec = nws;

% fitting
if pars.dofit

    % --- TODO: high-field NAD fitting adjustments -----------------------
    % Starting point: dat2mat_svsNAD uses VARPRO + de Graaf 12-peak table
    % for downfield NAD region. For Bruker high-field data:
    %   - widthBounds in ppm scale should narrow (Hz LW does not grow as
    %     fast as ppm dispersion at higher field), so default
    %     [10/f0 60/f0] in ppm may need re-tuning
    %   - center_freq / pivot for shiftSpectrumPhase still 4.72 ppm (1H)
    %   - hsvdClean range may need updating for high-field shift dispersion
    % For now: leave the fit dispatch as a stub.
    % --------------------------------------------------------------------

    if strcmpi(pars.dofit,'lcmodel')

        rawfile = fullfile(ws_dir,'lcmodel.raw');
        fid_for_raw = ifft(ifftshift(ws,1),[],1);
        create_lcmodelRAW(rawfile,fid_for_raw,'NUNFIL',length(fid_for_raw),'DELTAT',1/bw,'HZPPPM',f0);

        if ~isempty(nws)
            rawfile = fullfile(nws_dir,'lcmodel.raw');
            fid_for_raw = ifft(ifftshift(nws,1),[],1);
            create_lcmodelRAW(rawfile,fid_for_raw,'NUNFIL',length(fid_for_raw),'DELTAT',1/bw,'HZPPPM',f0);
        end

    elseif strcmpi(pars.dofit,'varpro')

        % --- TODO: re-tune widthBounds, focusRange, hsvdClean range ----
        warning('VARPRO fit not yet wired for Bruker high-field NAD data.')

    elseif strcmpi(pars.dofit,'auto')

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
            figure, plot(x,real(y),'k',x,yfit,'r',x,real(y)-yfit,'g'), legend({'data','fit','residual'}), set(gca,'xdir','reverse')
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

rmpath([thisPath filesep 'utils']);
