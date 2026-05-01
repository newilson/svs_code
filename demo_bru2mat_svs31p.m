%% demo dat2mat_svs31p
clear,clc

%% Bruker data - 500 MHz Chemistry scanner 

scan = 0;
switch scan
    case 0
        bru_dir = 'C:\Users\CAMIPM-NW\Downloads\rbnmr\rbnmr\data\Human_Blood_1\7';
    case 1
        bru_dir = 'C:\Users\CAMIPM-NW\Downloads\rbnmr\rbnmr\data\Human_Blood_2\7';
end

pars.lb         = 10; % Hz
pars.peakAddLB = 10; % Hz
pars.initPhDeg = 180; % degrees
pars.initPPMShift = -3.0; % ppm
pars.PC         = 0;            % skip manual phase GUI

pars.hsvd.flag  = false;
pars.hsvd.sv    = 30;

% --- peak-based phase + frequency correction
%     With >=2 peaks, bru2mat_svs31p fits a linear phase ramp through them
%     and applies it (zero-order + first-order in one shot). bATP / gATP /
%     aATP are the most reliable, well-isolated 31P landmarks in blood.
%     Note that there is no appreciable PCr.
pars.peaks      = [-16.18, -2.53, -7.56];   % bATP, gATP, aATP (ppm)
pars.peakRanges = [0.5 0.5;            % bATP   ±0.5 ppm
                   0.4 0.4;            % gATP  ±0.4 ppm
                   0.4 0.4];           % aATP  ±0.3 ppm (avoid NADH at -8.11)
% Alternative: zero-order phase only at alpha-ATP (anchor near NAD/UDP fit
% window so residual phase across that window stays small; pivot is
% auto-set from pars.peaks(1) downstream).
% pars.peaks      = -7.6;
% pars.peakRanges = [0.2 0.2];

pars.dofit            = true;
pars.fit.mode         = 'basisvarpro';
pars.plt              = true;


% --- basisVarpro: full-spectrum fit (LCModel convention, ±19.5 ppm)
pars.fit.ppm_range    = [-19.5 19.5];
pars.fit.metabs       = { ...
    'PE','PC','DPG1','DPG2','Pi','GPE','GPC','BroadMP', ...
    'PCr','gATP','aATP','bATP', ...
    'NADH','NADplus', ...
    'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc', ...
    ...'Broad8_5',...
    'Broad10_5'};
% Alternative: narrow NAD/UDP region only (faster, fewer free parameters)
% pars.fit.ppm_range = [-10.5 -6.5];
% pars.fit.metabs    = {'aATP','NADH','NADplus', ...
%                       'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'};
pars.fit.phase1Enable = true;

% --- per-metabolite shift / lbL bounds (built from pars.fit.metabs)
%   Categories:
%     NAD/UDP : tight chemical-shift bounds (Ren ±0.05 ppm => ±6 Hz at 7T;
%               allow ~2x as a safety margin). UDP lbL capped tight so a
%               peg at the cap is a diagnostic that the basis is being
%               used as a baseline absorber rather than fitting a peak.
%     Broad   : wide shift freedom so each broad component can find its
%               own center; lbL upper bound stays at the default (basis
%               already carries lwHz=300 Hz baked in).
%     others  : default ±60 Hz shift, 0-100 Hz lbL.
nMet    = numel(pars.fit.metabs);
isUDP   = startsWith(pars.fit.metabs,'UDP');
isNAD   = ismember(pars.fit.metabs,{'NADH','NADplus'});
isBroad = startsWith(pars.fit.metabs,'Broad');

shB = 4*repmat([-60 60], nMet, 1);                          % default
shB(isUDP,:)   = 4*repmat([-6   6],   sum(isUDP),   1);     % ~±0.05 ppm (Ren 2021 prior)
shB(isNAD,:)   = 4*repmat([-12  12],  sum(isNAD),   1);
% shB(isBroad,:) = repmat([-300 300], sum(isBroad), 1);     % ~±2.5 ppm
% shB(strcmp(pars.fit.metabs,'Broad8_5'),:)  = [-60 60];   % ~±0.5 ppm
shB(strcmp(pars.fit.metabs,'Broad10_5'),:) = 4*[-300 60]; % keep wide

lbB = repmat([0 100], nMet, 1);                           % default
lbB(isUDP,:)   = repmat([15 50], sum(isUDP), 1);           % UDP narrow
lbB(isNAD,:)   = repmat([15 50], sum(isNAD), 1);           % NAD a bit wider
lbB(isBroad,:) = repmat([20 200], sum(isBroad), 1);        % allow more flexibility over broad peaks

% --- per-metabolite lbL initial values (must lie inside lbLBounds rows)
%     Default lbLInit=2 Hz violates the UDP/NAD floor (15 Hz) and Broad
%     floor (50 Hz). Seed each category mid-window so the optimizer starts
%     interior on every component.
lbInit = 50  * ones(nMet,1);    % default
lbInit(isUDP)   = 25;           % interior of [15 50]
lbInit(isNAD)   = 25;
lbInit(isBroad) = 100;          % interior of [50 200]

pars.fit.fitOpt.shiftBounds = shB;
pars.fit.fitOpt.lbLBounds   = lbB;
pars.fit.fitOpt.lbLInit     = lbInit;
pars.fit.fitOpt.phaseBounds = [-20 20];

% --- soft amplitude-ratio constraints (Ren NMR Biomed 2021 brain UDP fractions)
sc(1) = struct('idxA','UDPGal',    'idxB','UDPGlcNAc','ratio',0.53,'lambda',1e4);
sc(2) = struct('idxA','UDPGlc',    'idxB','UDPGlcNAc','ratio',1.00,'lambda',1e4);
sc(3) = struct('idxA','UDPGalNAc', 'idxB','UDPGlcNAc','ratio',0.80,'lambda',1e4);
pars.fit.softConstraints = sc;

% --- parameter ties (share nonlinear parameters across listed members)
%     The four UDP species are chemically near-identical and live in the same
%     shim, so their Lorentzian linewidth and frequency shift must be equal
%     to within ~Hz in vivo.  Tying lbL+shift across all four prevents UDP
%     bases from soaking up unrelated baseline/aATP/NAD residual via wildly
%     divergent linewidths.  Phase (phi0) is left independent -- it absorbs
%     small residual phase that may differ in the real data even if the
%     basis is shared.  ATPs are intentionally NOT tied.
tg(1) = struct('members', {{'UDPGlcNAc','UDPGalNAc','UDPGlc','UDPGal'}}, ...
               'params',  {{'lbL','shift'}});
% pars.fit.tieGroups = tg;
pars.fit.tieGroups = [];

%% process
[output,pars,fullHdr] = bru2mat_svs31p(pars,bru_dir);
