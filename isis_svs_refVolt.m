  function [refVolt,output] = isis_svs_refVolt(dat_files,pars)
%
% [refVolt,output] = isis_svs_refVolt(dat_files,pars)
%
% Transmitter reference voltage calibration from ISIS_SVS flip angle
% series.  dat_files is a cell array of .dat file names, each one an
% acquisition at a different prescribed excitation flip angle.  With the
% ISIS HSn inversion pulses overdriven, the only B1+ dependence left is
% assumed in the excitation, so
%
%   S(FA) = a * sin(relB1 * FA),   relB1 = achieved / prescribed flip
%   refVolt = TxVolt / relB1
%
% Peak signals come from HSVD.  Poles (frequency + damping) are estimated
% once from the summed series and the complex amplitudes are then re-solved
% for each flip angle, so every point in the curve is a linear fit against
% a common lineshape.  Signal per peak is the area under the magnitude
% spectrum of the poles inside that peak's ppm window, summed coherently and
% integrated.  Peaks listed in pars.mainPoleOnly (31P default PCr) count
% only their strongest pole; other poles in their window stay in the fit,
% unassigned.
%
% The HSVD poles and amplitudes are fit on FID points skipPts+1 to
% skipPts+hsvd.maxPts only, dropping the leading points corrupted by the
% excitation transient and the noise-only tail.  The basis is still
% referenced to t = 0, so the amplitudes are back-extrapolated.
%
% Every peak is fit separately; refVolt is combined across peaks per
% pars.combine.
%
% Steps:
% read and concatenate the flip angle series (mapVBVD per file)
% protocol consistency and RF pulse clipping checks
% nucleus dependent defaults
% zero order phase of the reference peak (frequency alignment only if pars.align.freq)
% HSVD pole estimation on the summed spectrum
% per flip angle amplitude solve
% per peak sine fit
%
% ***** DEFAULT values depend on the observed nucleus (31P or 1H only for now) ********
%
% 08/2026 - NW

% in case current directory is not this one
thisFile = mfilename('fullpath');
thisPath = fileparts(thisFile);
addpath([thisPath filesep 'utils']);
addpath([thisPath filesep 'mapVBVD']);

if nargin<1 || isempty(dat_files)
    % [file,path] = uigetfile('*.dat','Choose the flip angle series dat files','MultiSelect','on');
    [file,path] = uigetfile('/mnt/rtexport/RTexport_Current/*.dcm','Choose the flip angle series dicom files','MultiSelect','on');
    if isequal(file,0), return; end
    if ~iscell(file), file = {file}; end
    dat_files = fullfile(path,file);
end
if ~iscell(dat_files), dat_files = cellstr(dat_files); end
if nargin<2, pars = struct(); end

% nucleus independent defaults
if ~isfield(pars,'plt') || isempty(pars.plt), pars.plt = true; end
if ~isfield(pars,'excludeFA'), pars.excludeFA = []; end
if ~isfield(pars,'combine') || isempty(pars.combine), pars.combine = 'median'; end
if ~isfield(pars,'lb') || isempty(pars.lb), pars.lb = 0; end
if ~isfield(pars,'skipPts') || isempty(pars.skipPts), pars.skipPts = 6; end   % leading FID points left out of the HSVD
if ~isfield(pars,'nucleus'), pars.nucleus = ''; end
if ~isfield(pars,'pulseModel'), pars.pulseModel = ''; end
if ~isfield(pars,'pulseDur'), pars.pulseDur = []; end          % us, [] = read WIP
if ~isfield(pars,'excCtr'), pars.excCtr = []; end              % ppm, [] = read WIP
if ~isfield(pars,'isisBW'), pars.isisBW = []; end              % Hz, [] = read WIP
if ~isfield(pars,'wip') || isempty(pars.wip)
    pars.wip = struct('durIdx',6,'excCtrIdx',11,'isisBWIdx',3);
end
if ~isfield(pars,'hsvd'), pars.hsvd = struct(); end
if ~isfield(pars.hsvd,'maxPts') || isempty(pars.hsvd.maxPts), pars.hsvd.maxPts = 512; end   % FID points used after skipPts
if ~isfield(pars,'align'), pars.align = struct(); end
if ~isfield(pars.align,'freq') || isempty(pars.align.freq), pars.align.freq = false; end   % shift each acquisition onto the reference peak

%% read the series
hdr = struct([]);
for ii=1:numel(dat_files)
    fname = dat_files{ii};
    assert(isfile(fname),'file not found: %s',fname);

    [~,shortname,ext] = fileparts(fname);
    if any(strcmpi(ext,{'.ima','.dcm'}))
        cur = readSpecDicom(fname,pars);
    else
        cur = readSpecTwix(fname,pars);
    end

    jj = length(hdr)+1;
    hdr(jj).name    = [shortname ext];
    hdr(jj).file    = fname;
    hdr(jj).source  = cur.source;
    hdr(jj).fa      = cur.fa;
    hdr(jj).TxVolt  = cur.TxVolt;
    hdr(jj).nucleus = cur.nucleus;
    hdr(jj).rfVolt  = cur.rfVolt;
    hdr(jj).npts    = cur.npts;
    hdr(jj).nCha    = cur.nCha;
    hdr(jj).nAve    = cur.nAve;
    hdr(jj).bw      = cur.bw;
    hdr(jj).f0      = cur.f0;
    hdr(jj).adcshift= cur.adcshift;
    hdr(jj).pulseDur= cur.pulseDur;
    hdr(jj).excCtr  = cur.excCtr;
    hdr(jj).isisBW  = cur.isisBW;
    hdr(jj).fid     = cur.fid;
end
if ~isempty(pars.excludeFA)
    drop = ismember([hdr.fa],pars.excludeFA);
    if any(drop)
        fprintf('excluding flip angle(s) %s by request\n',mat2str([hdr(drop).fa]));
        hdr = hdr(~drop);
    end
end

[~,ord] = sort([hdr.fa]);
hdr = hdr(ord);
nFlip = numel(hdr);
assert(nFlip>=3,'Need at least 3 flip angles, found %d',nFlip);

% protocol consistency
for ii=2:nFlip
    if ~isequal(hdr(ii).npts,hdr(1).npts) || ~isequal(hdr(ii).bw,hdr(1).bw) || ...
       ~isequal(hdr(ii).f0,hdr(1).f0)     || ~isequal(hdr(ii).TxVolt,hdr(1).TxVolt) || ...
       ~isequal(hdr(ii).nCha,hdr(1).nCha) || ~strcmp(hdr(ii).nucleus,hdr(1).nucleus) || ...
       ~isequaln(hdr(ii).pulseDur,hdr(1).pulseDur) || ~isequaln(hdr(ii).excCtr,hdr(1).excCtr)
        error('%s does not match the rest of the series',hdr(ii).name)
    end
end

fa      = [hdr.fa]';
TxVolt  = hdr(1).TxVolt;
f0      = hdr(1).f0;
adcshift= hdr(1).adcshift;
npts    = hdr(1).npts;
nCha    = hdr(1).nCha;
bw      = hdr(1).bw;
t       = (0:npts-1)'/bw;
ppm     = (-1/2:1/npts:1/2-1/npts)'*bw/f0 + adcshift;

if isempty(pars.pulseDur), pars.pulseDur = hdr(1).pulseDur; end
if isempty(pars.excCtr),   pars.excCtr   = hdr(1).excCtr;   end
if isempty(pars.isisBW),   pars.isisBW   = hdr(1).isisBW;   end
if isempty(pars.pulseModel)
    if ~isempty(pars.pulseDur) && ~isempty(pars.excCtr)
        pars.pulseModel = 'rect';
    else
        pars.pulseModel = 'sine';
    end
end
if strcmpi(pars.pulseModel,'rect')
    assert(~isempty(pars.pulseDur) && ~isempty(pars.excCtr), ...
        ['pulseModel ''rect'' needs the excitation pulse duration and centre - ' ...
         'set pars.pulseDur (us) and pars.excCtr (ppm), or pars.wip indices']);
end

if isempty(pars.nucleus), pars.nucleus = hdr(1).nucleus; end
is1H = strcmpi(pars.nucleus,'1H');

fprintf('\n%s: %d flip angles: %s deg   (prescribed ref voltage %g V)\n', ...
    pars.nucleus, nFlip, mat2str(fa'), TxVolt);

%% nucleus dependent defaults
if is1H
    % water reference only, ppm window matching pars.watfit.ppm_range in dat2mat_svsNAD
    if ~isfield(pars,'peaks') || isempty(pars.peaks)
        % name, ppm, half width of the assignment window (ppm)
        pars.peaks = { 'Water', 4.7, 0.75 };
    end
    if ~isfield(pars,'mainPoleOnly'), pars.mainPoleOnly = {}; end
    if ~isfield(pars.hsvd,'K')     || isempty(pars.hsvd.K),     pars.hsvd.K = 7; end
    if ~isfield(pars.hsvd,'maxLW') || isempty(pars.hsvd.maxLW), pars.hsvd.maxLW = 100; end
    if ~isfield(pars.hsvd,'range') || isempty(pars.hsvd.range), pars.hsvd.range = [3.2 6.2]; end
    if ~isfield(pars.align,'ppm')  || isempty(pars.align.ppm),  pars.align.ppm = 4.7; end
    if ~isfield(pars.align,'win')  || isempty(pars.align.win),  pars.align.win = 0.5; end
    if ~isfield(pars.align,'lb')   || isempty(pars.align.lb),   pars.align.lb = 2; end
else
    if ~isfield(pars,'peaks') || isempty(pars.peaks)
        % name, ppm (PCr = 0), half width of the assignment window (ppm)
        % GPC window also covers GPE; PCr is narrow as the frequency is set on
        % it; bATP is left out as it lies beyond the first null of the rect
        % excitation
        pars.peaks = { 'PME',   6.7,  0.5
                       'Pi',    4.85, 0.3
                       'GPC',   3.0,  0.6
                       'PCr',   0.0,  0.08
                       'gATP', -2.5,  0.6
                       'aATP', -7.6,  0.8 };
    end
    % PCr is fit by a main pole plus a small one ~19 Hz away; only the main one counts
    if ~isfield(pars,'mainPoleOnly'), pars.mainPoleOnly = {'PCr'}; end
    if ~isfield(pars.hsvd,'K')     || isempty(pars.hsvd.K),     pars.hsvd.K = 22; end
    if ~isfield(pars.hsvd,'maxLW') || isempty(pars.hsvd.maxLW), pars.hsvd.maxLW = 120; end
    if ~isfield(pars.hsvd,'range') || isempty(pars.hsvd.range), pars.hsvd.range = [-25 15]; end
    if ~isfield(pars.align,'ppm')  || isempty(pars.align.ppm),  pars.align.ppm = 0; end
    if ~isfield(pars.align,'win')  || isempty(pars.align.win),  pars.align.win = 1.5; end
    if ~isfield(pars.align,'lb')   || isempty(pars.align.lb),   pars.align.lb = 20; end
end
assert(pars.align.ppm>min(ppm) && pars.align.ppm<max(ppm), ...
    ['reference peak at %g ppm is outside the acquired range [%.1f %.1f] ppm - ' ...
     'set pars.align.ppm and pars.peaks for this acquisition'], ...
    pars.align.ppm, min(ppm), max(ppm));

% check for pulse clipping
rfVolt = [hdr.rfVolt]';
if all(isfinite(rfVolt)) && any(rfVolt~=0)
    [p,S] = polyfit(fa,rfVolt,1);
    output.rfVolt = rfVolt;
    if isfield(S,'rsquared')
        output.rfVoltRsq = S.rsquared;
        if S.rsquared < 1 - 1e-3
            warning('excitation pulse voltage is not linear in flip angle (Rsq = %.4f) - possible clipping',S.rsquared);
            if pars.plt
                figure('Name','RF pulse voltage')
                plot(0:max(fa),polyval(p,0:max(fa)),'-b',fa,rfVolt,'o')
                xlabel('Flip Angle'), ylabel('Voltage (V)')
                legend(sprintf('Rsq = %3.4f',S.rsquared),'Location','NorthWest')
            end
        end
    end
end

%% preprocess
% coil weights are taken from the whole series so every flip angle is
% combined identically
if nCha>1
    p1 = zeros(nFlip,nCha);
    for ii=1:nFlip, p1(ii,:) = hdr(ii).fid(1,:); end
    [~,iref] = max(sum(abs(p1),2));
    w = conj(p1(iref,:));
    w = w/norm(w);
    for ii=1:nFlip, hdr(ii).fid = hdr(ii).fid * w.'; end
    fprintf('combined %d coils with common weights from the %g deg acquisition\n',nCha,fa(iref));
end

spec    = zeros(npts,nFlip);
hzshift = zeros(nFlip,1);
ph0     = zeros(nFlip,1);
alignInd = find(ppm > pars.align.ppm-pars.align.win & ppm < pars.align.ppm+pars.align.win);
% the reference peak phase is read at its maximum, so on a 16x zero filled
% axis; without frequency alignment the peak does not sit on a bin
nzf   = 16*npts;
ppmzf = (-1/2:1/nzf:1/2-1/nzf)'*bw/f0 + adcshift;
alignIndZf = find(ppmzf > pars.align.ppm-pars.align.win & ppmzf < pars.align.ppm+pars.align.win);
for ii=1:nFlip
    fid   = expFilter(t,pars.lb,hdr(ii).fid);
    fidlb = expFilter(t,pars.align.lb,hdr(ii).fid);

    % onto the common ppm axis if this acquisition used a different ADC centre
    dshift = (adcshift - hdr(ii).adcshift)*f0;

    % frequency alignment of the reference peak
    speclb = shiftSpectrumFrequency(fftshift(fft(fidlb,[],1),1),dshift,t);
    hzshift(ii) = dshift;
    if pars.align.freq
        [~,im] = max(abs(speclb(alignInd)));
        hzshift(ii) = dshift + (ppm(alignInd(im)) - pars.align.ppm)*f0;
        speclb = shiftSpectrumFrequency(speclb,hzshift(ii)-dshift,t);
    end

    % zero order phase of the reference peak
    speczf = fftshift(fft(ifft(ifftshift(speclb,1),[],1),nzf,1),1);
    [~,im] = max(abs(speczf(alignIndZf)));
    ph0(ii) = angle(speczf(alignIndZf(im)));

    spec(:,ii) = shiftSpectrumPhase( ...
        shiftSpectrumFrequency(fftshift(fft(fid,[],1),1),hzshift(ii),t), ph0(ii));
end
if pars.align.freq
    fprintf('reference peak (%g ppm) shifts: %s Hz\n',pars.align.ppm,mat2str(round(hzshift',2)));
else
    fprintf('frequency alignment off\n');
end

%% HSVD
% poles from the summed series, amplitudes re-solved per flip angle, both on
% the fitInd window of the FID
assert(pars.skipPts>=0 && pars.skipPts<npts/2, ...
    'pars.skipPts must be between 0 and %d',npts/2-1);
fitInd  = (pars.skipPts+1 : min(npts,pars.skipPts+pars.hsvd.maxPts))';
rangeHz = (pars.hsvd.range - adcshift) * f0;
fidAll  = ifft(ifftshift(spec,1),[],1);
fidSum  = sum(fidAll,2);
[fHz,dHz] = hsvd_fit(t(1:numel(fitInd)),fidSum(fitInd),pars.hsvd.K,[],rangeHz,false);

lw   = dHz/pi;
keep = lw>0 & lw<=pars.hsvd.maxLW;
assert(any(keep),'no HSVD poles survived the linewidth cutoff');
poleppm = fHz(keep)/f0 + adcshift;
lw      = lw(keep);
lnz     = (-dHz(keep) + 1i*2*pi*fHz(keep)) * (t(2)-t(1));
E       = exp( lnz * (0:npts-1) ).';

fprintf('HSVD: %d poles kept of %d (K=%d, max linewidth %g Hz), FID points %d to %d\n', ...
    numel(lw), numel(dHz), pars.hsvd.K, pars.hsvd.maxLW, fitInd(1), fitInd(end));

names = pars.peaks(:,1);
nPeaks    = numel(names);
poleSel = false(numel(lw),nPeaks);
for jj=1:nPeaks
    poleSel(:,jj) = abs(poleppm - pars.peaks{jj,2}) <= pars.peaks{jj,3};
end

% strongest pole only for the pars.mainPoleOnly peaks, judged on the summed series
ampSum = E(fitInd,:) \ fidSum(fitInd);
for jj = find(ismember(names,pars.mainPoleOnly))'
    ip = find(poleSel(:,jj));
    if numel(ip)<2, continue, end
    [~,im] = max(abs(ampSum(ip)));
    ip(im) = [];
    poleSel(ip,jj) = false;
end
poleOther = ~any(poleSel,2);

poleSum = @(a,jj) poleArea(a,poleSel(:,jj),E,ppm);

sig   = nan(nFlip,nPeaks);
amps  = zeros(numel(lw),nFlip);
specFit = zeros(npts,nFlip);
for ii=1:nFlip
    amps(:,ii) = E(fitInd,:) \ fidAll(fitInd,ii);
    specFit(:,ii) = fftshift(fft(E*amps(:,ii),[],1),1);
    for jj=1:nPeaks
        sig(ii,jj) = poleSum(amps(:,ii),jj);
    end
end

% component decomposition of the summed series
specSum = sum(spec,2);
fitSum  = fftshift(fft(E*ampSum,[],1),1);
compSum = zeros(npts,nPeaks+1);
for jj=1:nPeaks
    if any(poleSel(:,jj))
        compSum(:,jj) = fftshift(fft(E(:,poleSel(:,jj))*ampSum(poleSel(:,jj)),[],1),1);
    end
end
if any(poleOther)
    compSum(:,nPeaks+1) = fftshift(fft(E(:,poleOther)*ampSum(poleOther),[],1),1);
end

fprintf('poles per peak: ');
for jj=1:nPeaks, fprintf('%s %d  ',names{jj},sum(poleSel(:,jj))); end
fprintf('unassigned %d\n',sum(poleOther));

%% fit each peak separately
% offset of each peak from the excitation centre; the readout centre (ppm
% axis) is independent of it
peakOffHz = nan(1,nPeaks);
if ~isempty(pars.excCtr)
    peakOffHz = (cell2mat(pars.peaks(:,2))' - pars.excCtr) * f0;
end
switch lower(pars.pulseModel)
    case 'sine'
        fitFunc = @(par,FAdeg,dfHz) par(1) * sind(par(2) * FAdeg);
    case 'rect'
        fitFunc = @(par,FAdeg,dfHz) par(1) * rectMxy(FAdeg,pars.pulseDur,dfHz,par(2));
    otherwise
        error('pars.pulseModel must be ''rect'' or ''sine'', not %s',pars.pulseModel)
end
opts    = optimset('Display','off','TolFun',1e-14,'TolX',1e-14,'MaxFunEvals',4000);

relB1Peak = nan(1,nPeaks); aPeak = nan(1,nPeaks); rmsPeak = nan(1,nPeaks);
for jj=1:nPeaks
    y = sig(:,jj);
    if ~any(isfinite(y)) || max(y)<=0, continue, end
    sc  = max(y);
    ff  = @(par,FAdeg) fitFunc(par,FAdeg,peakOffHz(jj));
    par = lsqcurvefit(ff,[1 0.9],fa,y/sc,[0 0.2],[10 2],opts);
    relB1Peak(jj)   = par(2);
    aPeak(jj)   = par(1)*sc;
    rmsPeak(jj) = sqrt(mean((y/sc - ff(par,fa)).^2));
end

vPeak  = TxVolt ./ relB1Peak;
true90 = 90 ./ relB1Peak;

switch lower(pars.combine)
    case 'median'
        relB1Use = median(relB1Peak,'omitnan');
    case 'mean'
        relB1Use = mean(relB1Peak,'omitnan');
    otherwise
        iu = find(strcmpi(names,pars.combine),1);
        assert(~isempty(iu),'combine peak %s is not in pars.peaks',pars.combine);
        relB1Use = relB1Peak(iu);
end
refVolt = TxVolt / relB1Use;

%% report
if strcmpi(pars.pulseModel,'rect')
    fprintf('\nexcitation: %g us rectangular, centre %g ppm\n',pars.pulseDur,pars.excCtr);
else
    fprintf('\nexcitation: pure sine model (profile ignored)\n');
end
if ~isempty(pars.isisBW) && ~isempty(pars.excCtr)
    edgeHz = max(abs(peakOffHz));
    fprintf('ISIS inversion bandwidth %g Hz (+-%g); peaks reach %.0f Hz off centre\n', ...
        pars.isisBW, pars.isisBW/2, edgeHz);
    if edgeHz > 0.5*pars.isisBW
        warning('%s is %.0f Hz off the excitation centre, outside the +-%g Hz ISIS inversion band', ...
            names{find(abs(peakOffHz)==edgeHz,1)}, edgeHz, pars.isisBW/2);
    end
end
fprintf('\n%-6s %10s %10s %10s %10s %10s %10s\n', ...
    'peak','ppm','off Hz','relB1','true 90','refVolt','fit RMS');
for jj=1:nPeaks
    fprintf('%-6s %10.2f %10.0f %10.3f %10.1f %10.1f %10.3f\n', ...
        names{jj}, pars.peaks{jj,2}, peakOffHz(jj), relB1Peak(jj), true90(jj), vPeak(jj), rmsPeak(jj));
end
fprintf('%-6s %10s %10s %10.3f %10.1f %10.1f\n', ...
    upper(pars.combine),'','',relB1Use,90/relB1Use,refVolt);
if nPeaks>1
    fprintf('\nper peak refVolt spread: %.1f to %.1f V (CV %.1f%%)\n', ...
        min(vPeak), max(vPeak), 100*std(vPeak,'omitnan')/mean(vPeak,'omitnan'));
end
fprintf('Prescribed ref voltage = %6.1f V\n', TxVolt);
fprintf('Calibrated ref voltage = %6.1f V\n', refVolt);

%% plots
if pars.plt
    figure('Name','ISIS flip angle series')
    plot(ppm,abs(spec)), set(gca,'xdir','reverse')
    xlim(sort(pars.hsvd.range))
    xlabel('ppm'), ylabel('|spectrum|')
    legend(arrayfun(@(x)sprintf('%g deg',x),fa,'UniformOutput',false))

    figure('Name','summed spectrum components')
    hold on
    co  = lines(nPeaks);
    h   = gobjects(nPeaks+3,1);
    lab = cell(nPeaks+3,1);
    h(1) = plot(ppm,real(specSum),'-k','LineWidth',1);
    lab{1} = 'data (summed series)';
    h(2) = plot(ppm,real(fitSum),'--','Color',[.4 .4 .4],'LineWidth',1);
    lab{2} = 'total fit';
    for jj=1:nPeaks
        h(2+jj) = plot(ppm,real(compSum(:,jj)),'-','Color',co(jj,:),'LineWidth',1.5);
        lab{2+jj} = sprintf('%s (%d poles)',names{jj},sum(poleSel(:,jj)));
    end
    h(nPeaks+3) = plot(ppm,real(compSum(:,nPeaks+1)),'-','Color',[.75 .75 .75],'LineWidth',1);
    lab{nPeaks+3} = sprintf('unassigned (%d poles)',sum(poleOther));
    plot(ppm,real(specSum)-real(fitSum),':g')
    set(gca,'xdir','reverse'), xlim(sort(pars.hsvd.range)), grid on
    xlabel('ppm'), ylabel('Real')
    legend(h,lab,'Location','best')
    title('HSVD component fit of the summed series')

    figure('Name','fitted spectra')
    nr = ceil(sqrt(nFlip)); nc = ceil(nFlip/nr);
    inRange = ppm>=min(pars.hsvd.range) & ppm<=max(pars.hsvd.range);
    ylims = [min([real(spec(inRange,:)); real(spec(inRange,:)-specFit(inRange,:))],[],'all') ...
             max(real(spec(inRange,:)),[],'all')];
    ylims = ylims + 0.05*diff(ylims)*[-1 1];
    for jj=1:nFlip
        subplot(nr,nc,jj), hold on
        plot(ppm,real(spec(:,jj)),'-k')
        plot(ppm,real(specFit(:,jj)),'-r')
        plot(ppm,real(spec(:,jj))-real(specFit(:,jj)),':g')
        set(gca,'xdir','reverse')
        xlim(sort(pars.hsvd.range)), ylim(ylims), grid on
        if jj>nFlip-nc, xlabel('ppm'), end
        title(sprintf('%g deg',fa(jj)))
        if jj==1, legend({'data','HSVD fit','residual'},'Location','best'), end
    end

    figure('Name','flip angle dependence by peak')
    nr = ceil(sqrt(nPeaks)); nc = ceil(nPeaks/nr);
    for jj=1:nPeaks
        subplot(nr,nc,jj), hold on
        FAinterp = 0:1:max([fa(:); true90(jj)],[],'omitnan');
        sc = max(sig(:,jj));
        plot(fa,sig(:,jj)/sc,'bo','MarkerFaceColor','b')
        plot(FAinterp,fitFunc([aPeak(jj)/sc relB1Peak(jj)],FAinterp,peakOffHz(jj)),'k--')
        yl = ylim;
        plot([true90(jj) true90(jj)],yl,'-r')
        axis tight
        ylabel('Signal'), grid on
        if jj>nPeaks-nc
            xlabel('Flip Angle')
        end
        title(sprintf('%s: relB1 = %.3f, refVolt = %.0f V',names{jj},relB1Peak(jj),vPeak(jj)))
    end

    % figure('Name','reference voltage')
    % plot(1:nP,vPeak,'bo','MarkerFaceColor','b'), hold on
    % plot([0 nP+1],[refVolt refVolt],'k--')
    % plot([0 nP+1],[TxVolt TxVolt],'r-')
    % set(gca,'xtick',1:nP,'xticklabel',names), xlim([0 nP+1])
    % ylabel('Reference Voltage (V)'), grid on
    % legend({'per peak',sprintf('%s = %.0f V',pars.combine,refVolt), ...
    %         sprintf('prescribed = %.0f V',TxVolt)},'Location','best')
    % title(sprintf('Calibrated Ref Voltage = %6.1f V',refVolt),'fontweight','bold','fontsize',14)
end

%% output
output.files    = {hdr.file};
output.nucleus  = pars.nucleus;
output.FA       = fa;
output.TxVolt   = TxVolt;
output.refVolt  = refVolt;
output.relB1    = relB1Use;
output.peakName = names;
output.peakPpm  = cell2mat(pars.peaks(:,2));
output.sig      = sig;
output.relB1Peak= relB1Peak;
output.peakOffHz= peakOffHz;
output.pulseModel= pars.pulseModel;
output.pulseDur = pars.pulseDur;
output.excCtr   = pars.excCtr;
output.isisBW   = pars.isisBW;
output.aPeak    = aPeak;
output.rmsPeak  = rmsPeak;
output.vPeak    = vPeak;
output.true90   = true90;
output.polePpm  = poleppm;
output.poleLW   = lw;
output.poleAmp  = amps;
output.poleSel  = poleSel;
output.spec     = spec;
output.specFit  = specFit;
output.specSum  = specSum;
output.fitSum   = fitSum;
output.compSum  = compSum;
output.ampSum   = ampSum;
output.hzshift  = hzshift;
output.ph0      = ph0;
output.ppm      = ppm;
output.time     = t;
output.pars     = pars;

fclose('all');

rmpath([thisPath filesep 'mapVBVD']);
rmpath([thisPath filesep 'utils']);

end

% =====================================================================
function v = poleArea(a,sel,E,ppm)
% Area under the magnitude spectrum of the coherently summed poles.
v = 0;
if ~any(sel), return, end
v = trapz( ppm, abs( fftshift(fft(E(:,sel)*a(sel),[],1),1) ) );
end

% =====================================================================
function m = rectMxy(thDeg,Tus,dfHz,rel)
% Transverse magnetization from M0 = z after a rectangular pulse of duration
% Tus at offset dfHz.  Rotation by w_eff*T about the tilted effective field.
T  = Tus*1e-6;
w1 = rel*deg2rad(thDeg(:))/T;
dw = 2*pi*dfHz;
we = hypot(w1,dw);
we(we==0) = eps;
al = we*T;
sb = w1./we;
cb = dw./we;
m  = sb .* sqrt( cb.^2 .* (1-cos(al)).^2 + sin(al).^2 );
m  = reshape(m,size(thDeg));
end

% =====================================================================
function c = readSpecTwix(fname,pars)
% One acquisition from a Siemens raw .dat file.
twix = mapVBVD(fname);
if iscell(twix), twix = twix{end}; end
twix.image.flagRemoveOS    = false;
twix.image.flagAverageSets = false;
h = twix.hdr;

if isfield(h.Config,'DwellTime')
    c.bw = 1e9/h.Config.DwellTime;
elseif isfield(h.Config,'DwellTimeSig')
    c.bw = 1e9/h.Config.DwellTimeSig;
else
    error('unknown dwell time in %s',fname)
end

c.source = 'twix';
c.npts = twix.image.NCol;
c.nCha = twix.image.NCha;
fid = reshape(twix.image{''},c.npts,c.nCha,[]);
c.nAve = size(fid,3);
c.fid = mean(fid,3);                        % all averages and reps -> one FID per coil

c.fa      = h.MeasYaps.adFlipAngleDegree{1};
c.TxVolt  = h.MeasYaps.sTXSPEC.asNucleusInfo{1}.flReferenceAmplitude;
c.nucleus = regexprep(h.MeasYaps.sTXSPEC.asNucleusInfo{1}.tNucleus,'[^0-9A-Za-z]','');
c.rfVolt  = excPulseVolt(h);
c.f0      = 1e-6*h.Config.Frequency;
if isfield(h.Phoenix.sSpecPara,'dDeltaFrequency')       % absent when 0
    c.adcshift = h.Phoenix.sSpecPara.dDeltaFrequency;
else
    c.adcshift = 0;
end
c.pulseDur = wipVal(h,'alFree',pars.wip.durIdx);
c.excCtr   = wipVal(h,'adFree',pars.wip.excCtrIdx);
c.isisBW   = wipVal(h,'adFree',pars.wip.isisBWIdx);
end

% =====================================================================
function c = readSpecDicom(fname,pars)
% One acquisition from a Siemens spectroscopy DICOM (.IMA/.dcm).  The FID is
% already coil combined and averaged.
info = dicominfo(fname);

dtag = '';
for f = {'Private_7fe1_1010','Private_7fe1_1110'}
    if isfield(info,f{1}), dtag = f{1}; break, end
end
assert(~isempty(dtag),'no spectroscopy data element (7FE1,1010) in %s',fname);
v = typecast(uint8(info.(dtag)(:))','single');
% the DICOM FID is the conjugate of the raw twix FID (verified against S077,
% correlation 1.0000 after conjugation); without this the spectrum is mirrored
% about the ADC centre
c.fid  = double(v(1:2:end)).' - 1i*double(v(2:2:end)).';
c.npts = numel(c.fid);
c.nCha = 1;

A = ascconvParse(info);
c.source   = 'dicom';
c.bw       = 1e9/aget(A,'sRXSPEC.alDwellTime[0]',fname);
c.fa       = aget(A,'adFlipAngleDegree[0]',fname);
c.TxVolt   = aget(A,'sTXSPEC.asNucleusInfo[0].flReferenceAmplitude',fname);
c.nucleus  = regexprep(aget(A,'sTXSPEC.asNucleusInfo[0].tNucleus',fname),'[^0-9A-Za-z]','');
c.f0       = 1e-6*aget(A,'sTXSPEC.asNucleusInfo[0].lFrequency',fname);
c.adcshift = aget(A,'sSpecPara.dDeltaFrequency',[]);         % absent when 0
if isempty(c.adcshift), c.adcshift = 0; end
c.pulseDur = aget(A,sprintf('sWipMemBlock.alFree[%d]',pars.wip.durIdx-1),[]);
c.excCtr   = aget(A,sprintf('sWipMemBlock.adFree[%d]',pars.wip.excCtrIdx-1),[]);
c.isisBW   = aget(A,sprintf('sWipMemBlock.adFree[%d]',pars.wip.isisBWIdx-1),[]);
if isempty(c.excCtr), c.excCtr = 0; end                      % empty WIP entries are 0
if isempty(c.isisBW), c.isisBW = 0; end

c.nAve = aget(A,'lAverages',[]);
if isempty(c.nAve), c.nAve = 1; end

% the excitation pulse amplitude, for the clipping check
c.rfVolt = NaN;
for kk = 0:63
    nm = aget(A,sprintf('sTXSPEC.aRFPULSE[%d].tName',kk),[]);
    if isempty(nm), continue, end
    if contains(nm,'exc','IgnoreCase',true)
        c.rfVolt = aget(A,sprintf('sTXSPEC.aRFPULSE[%d].flAmplitude',kk),[]);
        break
    end
end
end

% =====================================================================
function A = ascconvParse(info)
% key/value map of the ASCCONV protocol block in the Siemens private header.
pv = '';
for f = {'Private_0029_1120','Private_0029_1020','Private_0029_1220'}
    if isfield(info,f{1}), pv = char(info.(f{1}))'; break, end
end
assert(~isempty(pv),'no Siemens private protocol header in this DICOM');
i1 = strfind(pv,'### ASCCONV BEGIN');
i2 = strfind(pv,'### ASCCONV END');
assert(~isempty(i1) && ~isempty(i2),'no ASCCONV block in this DICOM');
ln = regexp(pv(i1(1):i2(1)),'\n','split');
A = struct('key',{{}},'val',{{}});
for k = 1:numel(ln)
    tk = regexp(ln{k},'^\s*([^\s=]+)\s*=\s*(.+?)\s*$','tokens','once');
    if isempty(tk) || contains(tk{1},'__attribute__'), continue, end
    A.key{end+1} = tk{1};
    if startsWith(tk{2},'"')
        A.val{end+1} = strrep(tk{2},'"','');
    else
        A.val{end+1} = str2double(tk{2});
    end
end
end

% =====================================================================
function v = aget(A,key,fname)
% One ASCCONV value.  fname empty means optional ([] if absent).
k = find(strcmp(A.key,key),1);
if isempty(k)
    if isempty(fname)
        v = [];
        return
    end
    error('%s not found in the DICOM protocol header of %s',key,fname)
end
v = A.val{k};
end

% =====================================================================
function v = wipVal(hdr,fld,idx)
% One WIP entry (sWipMemBlock.alFree / .adFree).  Empty entries are 0, as in
% dat2mat_svs31p; [] only if the array is absent or shorter than idx.
v = [];
if isempty(idx), return, end
try
    w = hdr.MeasYaps.sWipMemBlock.(fld);
catch
    return
end
w(cellfun('isempty',w)) = {0};
w = cell2mat(w);
if numel(w)>=idx
    v = w(idx);
end
end

% =====================================================================
function v = excPulseVolt(hdr)
% Amplitude (V) of the excitation pulse, for the clipping check.
v = NaN;
try
    rf = hdr.MeasYaps.sTXSPEC.aRFPULSE;
catch
    return
end
iexc = [];
for ii=1:numel(rf)
    if isempty(rf{ii}) || ~isfield(rf{ii},'tName'), continue, end
    if contains(rf{ii}.tName,'exc','IgnoreCase',true), iexc = ii; break, end
end
if isempty(iexc), iexc = 1; end
if numel(rf)>=iexc && isfield(rf{iexc},'flAmplitude')
    v = rf{iexc}.flAmplitude;
end
end
