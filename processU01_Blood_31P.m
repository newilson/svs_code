function processU01_Blood_31P(SessionDir, FinalOutDir, exp31P, plt)
% PROCESSU01_BLOOD_31P  Absolute 31P quantification for one U01 blood
%                       session acquired on the Bruker high-field NMR.
%
%   processU01_Blood_31P(SessionDir, FinalOutDir, exp31P)
%   processU01_Blood_31P(SessionDir, FinalOutDir, exp31P, plt)
%
% Reads the 31P acquisition from <SessionDir>\<exp31P>\ (Bruker experiment
% subfolder containing fid/ser + acqus), runs the basis-Varpro fit from
% demo_bru2mat_svs31p.m via bru2mat_svs31p.m, then converts each fitted
% metabolite area to absolute [mM] using aATP as an internal reference with
% per-metabolite T1/T2 relaxation correction (Calf-style):
%
%   [M] = (S_M / S_aATP) * (R_aATP / R_M) * (nP_aATP / nP_M) * [aATP]
%
%   R_x = (1 - exp(-TR/T1_x)) * exp(-TE/T2_x)
%
% Blood has no PCr to use as a reference (unlike muscle); ATP is the most
% chemically defined 31P signal in whole blood, dominated by intra-RBC ATP
% (~1.8 mM in whole blood -- PLACEHOLDER, verify with site-specific value).
%
% Output: <FinalOutDir>\<scanID>_blood_31P_abs.mat
%   (output, parsOut, hdr, absQ, pars, metProps)

arguments
    SessionDir  (1,:) char
    FinalOutDir (1,:) char
    exp31P      (1,1) double {mustBeInteger,mustBePositive}
    plt         (1,1) logical = false
end

assert(isfolder(SessionDir), 'SessionDir not found: %s', SessionDir);
if ~isfolder(FinalOutDir), mkdir(FinalOutDir); end

bruDir = fullfile(SessionDir, num2str(exp31P));
assert(isfolder(bruDir), '31P exp folder not found: %s', bruDir);
fprintf('31P : %s\n', bruDir);

% =====================================================================
% Processing/Fitting parameters (from demo_bru2mat_svs31p.m -- blood preset)
% =====================================================================
pars = struct();
pars.lb         = 10;            % Hz
pars.peakAddLB  = 10;            % Hz
pars.initPhDeg  = 180;           % degrees
pars.initPPMShift = -3.0;        % ppm
pars.PC         = 0;             % skip manual phase GUI
pars.plt        = plt;

pars.hsvd.flag  = false;
pars.hsvd.sv    = 30;

% Peak-based phase + frequency correction.  bATP / gATP / aATP are the most
% reliable 31P landmarks in blood (no appreciable PCr).
pars.peaks      = [-16.18, -2.53, -7.56];   % bATP, gATP, aATP (ppm)
pars.peakRanges = [0.5 0.5;
                   0.4 0.4;
                   0.4 0.4];

pars.dofit            = true;
pars.fit.mode         = 'basisvarpro';

pars.fit.ppm_range    = [-19.5 19.5];
pars.fit.metabs       = { ...
    'PE','PC','DPG1','DPG2','Pi','GPE','GPC','BroadMP', ...
    'PCr','gATP','aATP','bATP', ...
    'NADH','NADplus', ...
    'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc', ...
    'Broad10_5'};
pars.fit.phase1Enable = true;

nMet    = numel(pars.fit.metabs);
isUDP   = startsWith(pars.fit.metabs,'UDP');
isNAD   = ismember(pars.fit.metabs,{'NADH','NADplus'});
isBroad = startsWith(pars.fit.metabs,'Broad');

shB = repmat([-30 30], nMet, 1);
shB(isUDP,:) = repmat([-6   6],  sum(isUDP), 1);
shB(isNAD,:) = repmat([-12  12], sum(isNAD), 1);
shB(strcmp(pars.fit.metabs,'Broad10_5'),:) = 4*[-300 60];

lbB = repmat([10 100], nMet, 1);
lbB(isUDP,:)   = repmat([15 50],  sum(isUDP),   1);
lbB(isNAD,:)   = repmat([15 50],  sum(isNAD),   1);
lbB(isBroad,:) = repmat([20 200], sum(isBroad), 1);

lbInit = 30 * ones(nMet,1);
lbInit(isUDP)   = 25;
lbInit(isNAD)   = 25;
lbInit(isBroad) = 100;

pars.fit.fitOpt.shiftBounds   = shB;
pars.fit.fitOpt.lbLBounds     = lbB;
pars.fit.fitOpt.lbLInit       = lbInit;
pars.fit.fitOpt.phaseBounds   = [-20 20];
pars.fit.baselineOpt.linearTol = 1e-8;

% Soft amplitude-ratio constraints (Ren NMR Biomed 2021 brain UDP fractions)
sc(1) = struct('idxA','UDPGal',    'idxB','UDPGlcNAc','ratio',0.53,'lambda',1e4);
sc(2) = struct('idxA','UDPGlc',    'idxB','UDPGlcNAc','ratio',1.00,'lambda',1e4);
sc(3) = struct('idxA','UDPGalNAc', 'idxB','UDPGlcNAc','ratio',0.80,'lambda',1e4);
pars.fit.softConstraints = sc;

pars.fit.tieGroups    = [];
pars.fit.lineShapeOpt = struct('enable',true,'nSide',2,'asymmetric',true,'maxSide',0.5);
pars.fit.verbose      = false;

% =====================================================================
% Per-metabolite T1 / T2 / number of equivalent 31P nuclei.
% PLACEHOLDER values at high field.  EDIT BEFORE PUBLICATION USE -- in
% particular, blood 31P relaxation is dominated by red-cell environment
% and differs substantially from brain/muscle literature values.
% =====================================================================
metProps = struct( ...
    'name',  {'PE','PC','DPG1','DPG2','Pi','GPE','GPC','BroadMP', ...
              'PCr','gATP','aATP','bATP', ...
              'NADH','NADplus', ...
              'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc', ...
              'Broad10_5'}, ...
    'T1_ms', {6000, 6000, 3000, 3000, 3500, 6000, 6000, 3000, ...
              3800, 1200, 1200, 1500, ...
               400,  400, ...
               500,  500,  500,  500, ...
              3000}, ...
    'T2_ms', {200, 200, 100, 100, 150, 200, 200, 150, ...
              180,  40,  40,  30, ...
               10,  10, ...
               10,  10,  10,  10, ...
              50}, ...
    'nPhos', {1,1,1,1,1,1,1,1,  1,1,1,1,  1,1,  1,1,1,1,  1});

refName    = 'aATP';
refConc_mM = 1.8;          % whole-blood [ATP] placeholder (RBC-dominated)
refIdx     = find(strcmp({metProps.name}, refName), 1);
assert(~isempty(refIdx), 'Reference %s missing from metProps', refName);

% =====================================================================
% Process / fit
% =====================================================================
[output, parsOut, hdr] = bru2mat_svs31p(pars, bruDir);

% bru2mat_svs31p rmpaths utils at exit; re-add for downstream helpers.
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

% sequence timing from Bruker acqus (D(2) is the relaxation delay in s;
% there is no TE field in 1D acqus, so TE is treated as 0 for absolute
% quant unless overridden later).
TR_ms = hdr.acqus.D(2) * 1e3;
TE_ms = 0;
fprintf('TR=%.0f ms  TE=%.2f ms\n', TR_ms, TE_ms);

% =====================================================================
% Absolute quant: aATP internal reference (single blood compartment)
% =====================================================================
% bru2mat_svs31p stores results per fit mode:
%   basisvarpro -> output.basisVarproResults.{names, amplScaled}
%                  (amplScaled = ampl ./ basis nVec, so it's already
%                   per-spin-normalized -> keep metProps.nPhos = 1)
%   amares      -> output.fit.{names, areas}  (stub; not yet implemented)
if isfield(output,'basisVarproResults')
    fitNames = output.basisVarproResults.names;
    fitAreas = output.basisVarproResults.amplScaled;
elseif isfield(output,'fit')
    fitNames = output.fit.names;
    fitAreas = output.fit.areas;
else
    error('processU01_Blood_31P: no fit results found in bru2mat_svs31p output');
end

[areaRef, T1ref, T2ref, nPref] = lookup_fit(refName, fitNames, fitAreas, metProps);
Rref = relax(TR_ms, TE_ms, T1ref, T2ref);

nFit = numel(fitNames);
absQ = struct();
absQ.reference = struct('name',refName,'conc_mM',refConc_mM, ...
                        'T1_ms',T1ref,'T2_ms',T2ref,'nPhos',nPref,'R',Rref);
absQ.TR_ms     = TR_ms;
absQ.TE_ms     = TE_ms;
absQ.metabolites = struct('name',{},'area',{},'conc_mM',{}, ...
                          'T1_ms',{},'T2_ms',{},'nPhos',{},'R',{},'ratio',{});

for k = 1:nFit
    nm_k = fitNames{k};
    pIdx = find(strcmp({metProps.name}, nm_k), 1);
    if isempty(pIdx)
        warning('No T1/T2/nPhos entry for %s -- skipping quant', nm_k);
        continue
    end
    T1_k = metProps(pIdx).T1_ms;
    T2_k = metProps(pIdx).T2_ms;
    nP_k = metProps(pIdx).nPhos;
    R_k  = relax(TR_ms, TE_ms, T1_k, T2_k);

    ratio   = fitAreas(k) / areaRef;
    conc_mM = ratio * (Rref / R_k) * (nPref / nP_k) * refConc_mM;

    absQ.metabolites(end+1) = struct( ...
        'name', nm_k, 'area', fitAreas(k), 'conc_mM', conc_mM, ...
        'T1_ms', T1_k, 'T2_ms', T2_k, 'nPhos', nP_k, ...
        'R', R_k, 'ratio', ratio); %#ok<AGROW>
end

for k = 1:numel(absQ.metabolites)
    fprintf('%-10s : %.4g mM\n', absQ.metabolites(k).name, absQ.metabolites(k).conc_mM);
end

% =====================================================================
% Save
% =====================================================================
[~, scanID] = fileparts(SessionDir);
outPath = fullfile(FinalOutDir, sprintf('%s_blood_31P_abs.mat', scanID));
save(outPath, 'output','parsOut','hdr','absQ','pars','metProps');
fprintf('Saved -> %s\n', outPath);

% remove path
rmpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

end

% =====================================================================
function R = relax(TR_ms, TE_ms, T1_ms, T2_ms)
R = (1 - exp(-TR_ms / T1_ms)) * exp(-TE_ms / T2_ms);
end

% =====================================================================
function [area, T1_ms, T2_ms, nPhos] = lookup_fit(name, fitNames, fitAreas, metProps)
iF = find(strcmp(fitNames, name), 1);
iP = find(strcmp({metProps.name}, name), 1);
assert(~isempty(iF), 'Metabolite %s not in fit names', name);
assert(~isempty(iP), 'Metabolite %s not in metProps',  name);
area  = fitAreas(iF);
T1_ms = metProps(iP).T1_ms;
T2_ms = metProps(iP).T2_ms;
nPhos = metProps(iP).nPhos;
end
