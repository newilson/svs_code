function [basisFIDs, basisInfo] = make31P_basisFIDs(t, hzpppm, varargin)
% MAKE31P_BASISFIDS  Simulate 31P brain metabolite basis FIDs for 7T fitting.
%
%   [basisFIDs, basisInfo] = make31P_basisFIDs(t, hzpppm)
%   [...] = make31P_basisFIDs(t, hzpppm, 'metabSubset', {...}, ...)
%
%   Returns a stack of complex FIDs (one column per metabolite) ready to
%   feed curvefitAuto_basisVarpro. Chemistry mirrors the AMARES
%   PK_7T_31P_Ren_NW prior knowledge (and make31P_ISIS_LCModel_basis.m):
%     PCr, ATP (alpha/beta/gamma), Pi, PE, PC, GPE, GPC, BroadMP, DPG (singlets/lines)
%     NADH (singlet)
%     NADplus (AB 2-spin system, J=20.14 Hz)
%     UDPGal, UDPGlc, UDPGalNAc, UDPGlcNAc (each as 2-spin alpha-P / beta-P)
%     Broad8_5, Broad10_5 (broad singlets at -8.5 and -10.5 ppm, modeling
%       broad envelope under aATP/NAD/UDP region)
%
%   Each FID is normalized so fids(t=0) = 1 (complex divide). With this
%   convention, the AB / multiplet structure is preserved, all metabolites
%   share the same t=0 magnitude scale, and per-component phi0 in the
%   downstream fit stays small.
%
%   INPUTS
%     t       time vector (s), length N. Expected t(1)=0.
%     hzpppm  scalar, Hz/ppm (e.g., 120.30 at 7T 31P).
%
%   NAME-VALUE OPTIONS
%     'metabSubset'  cell array of metabolite names to keep
%                    (default {} = all). Names are case-insensitive.
%     'adcshift'     ppm offset of the data ppm axis
%                    (e.g., fullHdr.Phoenix.sSpecPara.dDeltaFrequency).
%                    Each metabolite's ppm is shifted by -adcshift so the
%                    basis aligns with the data's ppm labeling. Default 0.
%     'lwHzDefault'  default Lorentzian FWHM applied in the basis (Hz).
%                    Default 0 (no broadening; lbL/lbG in the fit handle
%                    line broadening). BroadMP, Broad8_5 and Broad10_5 always start broad
%                    regardless of this value.
%
%   OUTPUTS
%     basisFIDs  [N x nMet] complex matrix
%     basisInfo  [nMet x 1] struct with fields:
%                  .name    metabolite label (char)
%                  .refPpm  nominal chemical shift (ppm, for reports)
%                  .group   integer group index
%                  .n       number of equivalent 31P nuclei contributing
%                           to this basis (e.g. 1 for ATP-alpha, 2 for
%                           NAD+, 2 for each UDP). Used downstream to
%                           convert fitted amplitudes into per-molecule
%                           concentrations: amplitude is proportional to
%                           C * n, so dividing by n gives a quantity
%                           proportional to molecular concentration.
%
%   See also: curvefitAuto_basisVarpro, make31P_ISIS_LCModel_basis

p = inputParser;
p.addParameter('metabSubset', {});
p.addParameter('adcshift', 0);
p.addParameter('lwHzDefault', 0);
p.parse(varargin{:});
metabSubset = p.Results.metabSubset;
adcshift    = p.Results.adcshift;
lwHzDefault = p.Results.lwHzDefault;

t = t(:);
N = length(t);

mets = getBasisParamStruct_31P7T_brain(lwHzDefault);

% filter
if ~isempty(metabSubset)
    keep = false(numel(mets),1);
    for k = 1:numel(mets)
        keep(k) = any(strcmpi(mets(k).name, metabSubset));
    end
    if ~any(keep)
        error('make31P_basisFIDs: no metabolites matched metabSubset.');
    end
    mets = mets(keep);
end

nMet = numel(mets);
basisFIDs = zeros(N, nMet);
basisInfo = repmat(struct('name','','refPpm',NaN,'group',NaN,'n',NaN,'lwHzBasis',0), nMet, 1);

for k = 1:nMet
    m = mets(k);

    switch lower(m.model)
        case 'singlet'
            fHz  = -(m.ppm - adcshift) * hzpppm;
            fids = m.amp * exp(-1i*2*pi*fHz*t) .* exp(-pi*m.lwHz*t);
            refPpm = m.ppm;

        case 'lines'
            [fHz_list, amp_list] = expandLineModel(m.lines, hzpppm, adcshift);
            fids = zeros(N,1);
            for j = 1:numel(fHz_list)
                fids = fids + amp_list(j) * exp(-1i*2*pi*fHz_list(j)*t);
            end
            fids = (m.amp * fids) .* exp(-pi*m.lwHz*t);
            refPpm = m.lines.ppm0;

        case 'pp_2spin'
            ppm_shifted = m.ppm_vec - adcshift;
            fids = simulate_isotropic_coupled_fid_2spin(ppm_shifted(:), m.J, hzpppm, t);
            fids = (m.amp * fids) .* exp(-pi*m.lwHz*t);
            refPpm = mean(m.ppm_vec);

        otherwise
            error('Unknown model "%s" for %s', m.model, m.name);
    end

    if abs(fids(1)) > 0
        fids = fids / fids(1);
    end

    basisFIDs(:,k)         = fids;
    basisInfo(k).name      = m.name;
    basisInfo(k).refPpm    = refPpm;
    basisInfo(k).group     = k;
    basisInfo(k).n         = m.n;
    basisInfo(k).lwHzBasis = m.lwHz;
end

end


% ========================================================================
% Metabolite prior knowledge (PCr = 0 ppm convention)
% ========================================================================
function mets = getBasisParamStruct_31P7T_brain(lwDefault)

mets = struct( ...
    'name',    {}, ...
    'model',   {}, ...
    'ppm',     {}, ...
    'ppm_vec', {}, ...
    'J',       {}, ...
    'lines',   {}, ...
    'amp',     {}, ...
    'lwHz',    {}, ...
    'n',       {} );

% ---- n = number of equivalent 31P nuclei contributing to the basis. The
%      simulated FID amplitude already accumulates the multiplet/AB lines
%      for each P, so for ratios within species sharing the same n, the
%      n factor cancels.  For cross-species comparisons (e.g., aATP vs
%      totalNAD), divide amplitude by n to get a quantity proportional to
%      molecular concentration.

% ---- Singlets (PCr-referenced)
mets(end+1) = newSinglet('PCr',  0.00, lwDefault, 1);
mets(end+1) = newSinglet('Pi',   4.84, lwDefault, 1);
mets(end+1) = newSinglet('DPG1', 5.71, lwDefault, 1);
mets(end+1) = newSinglet('DPG2', 5.23, lwDefault, 1);
mets(end+1) = newSinglet('BroadMP', 2.30, max(lwDefault,200), 1);   % broad
mets(end+1) = newSinglet('PC',   6.23, lwDefault, 1);
mets(end+1) = newSinglet('PE',   6.77, lwDefault, 1);
mets(end+1) = newSinglet('GPC',  2.94, lwDefault, 1);
mets(end+1) = newSinglet('GPE',  3.49, lwDefault, 1);

% ---- ATP (line models): each ATP basis is one P (alpha, gamma, or beta)
Jab = 16.3;
Jbg = 16.1;
mets(end+1) = newLines('aATP', lineModelDoublet(-7.56, Jab),       lwDefault, 1);
mets(end+1) = newLines('gATP', lineModelDoublet(-2.53, Jbg),       lwDefault, 1);
mets(end+1) = newLines('bATP', lineModelDD(-16.18, Jab, Jbg),      lwDefault, 1);

% ---- NAD(H) and UDP-sugars: pyrophosphate dinucleotides, 2 P each
mets(end+1) = newSinglet('NADH', -8.111, lwDefault, 2);
mets(end+1) = newPP2('NADplus',   [-8.131, -8.467], 20.14, lwDefault, 2);
mets(end+1) = newPP2('UDPGal',    [-8.022, -9.543], 20.3,  lwDefault, 2);
mets(end+1) = newPP2('UDPGlc',    [-8.162, -9.699], 20.3,  lwDefault, 2);
mets(end+1) = newPP2('UDPGalNAc', [-8.218, -9.782], 20.5,  lwDefault, 2);
mets(end+1) = newPP2('UDPGlcNAc', [-8.265, -9.963], 20.5,  lwDefault, 2);

% ---- Broad nuisance components (ad-hoc; no well-defined per-molecule n
%      -- treat as 1 for amplitude reporting). Two broad singlets that
%      can absorb the underlying broad envelope in the aATP/NAD/UDP region.
mets(end+1) = newSinglet('Broad8_5',  -8.5,  max(lwDefault,300), 1);
mets(end+1) = newSinglet('Broad10_5', -10.5, max(lwDefault,300), 1);

end


% ---- Struct constructors
function m = newSinglet(name, ppm, lwHz, n)
    m = struct('name',name,'model','singlet','ppm',ppm,'ppm_vec',[],'J',[], ...
               'lines',[],'amp',1,'lwHz',lwHz,'n',n);
end

function m = newLines(name, linesStruct, lwHz, n)
    m = struct('name',name,'model','lines','ppm',[],'ppm_vec',[],'J',[], ...
               'lines',linesStruct,'amp',1,'lwHz',lwHz,'n',n);
end

function m = newPP2(name, ppm_vec, J, lwHz, n)
    m = struct('name',name,'model','pp_2spin','ppm',[],'ppm_vec',ppm_vec,'J',J, ...
               'lines',[],'amp',1,'lwHz',lwHz,'n',n);
end


% ========================================================================
% Line-model helpers
% ========================================================================
function L = lineModelDoublet(ppm0, J1)
    L = struct('type','doublet','ppm0',ppm0,'J1',J1,'J2',[]);
end

function L = lineModelDD(ppm0, J1, J2)
    L = struct('type','dd','ppm0',ppm0,'J1',J1,'J2',J2);
end

function [fHz_list, amp_list] = expandLineModel(L, hzpppm, adcshift)
    if nargin < 3, adcshift = 0; end

    if isfield(L,'fHz') && ~isempty(L.fHz)
        fHz_list = L.fHz(:);
        amp_list = L.amp(:);
        return;
    end

    if ~isfield(L,'type') || ~isfield(L,'ppm0') || ~isfield(L,'J1')
        error('Line model needs either (fHz,amp) or (type,ppm0,J1[,J2]).');
    end

    f0 = -(L.ppm0 - adcshift) * hzpppm;

    switch lower(L.type)
        case 'singlet'
            fHz_list = f0;
            amp_list = 1;

        case 'doublet'
            J1 = L.J1;
            fHz_list = [f0 - J1/2; f0 + J1/2];
            amp_list = [0.5; 0.5];

        case {'triplet','t'}
            J1 = L.J1;
            fHz_list = [f0 - J1; f0; f0 + J1];
            amp_list = [0.25; 0.5; 0.25];

        case {'dd','doublet-of-doublets'}
            if ~isfield(L,'J2') || isempty(L.J2)
                error('dd line model requires J2.');
            end
            J1 = L.J1; J2 = L.J2;
            fHz_list = [f0 - J1/2 - J2/2;
                        f0 - J1/2 + J2/2;
                        f0 + J1/2 - J2/2;
                        f0 + J1/2 + J2/2];
            amp_list = 0.25 * ones(4,1);

        otherwise
            error('Unknown line model type "%s".', L.type);
    end
end


% ========================================================================
% 2-spin isotropic coupling simulation (exact AB handling)
% ========================================================================
function fids = simulate_isotropic_coupled_fid_2spin(ppm_vec, J, hzpppm, t)
% Hard 90x pulse on both spins, observe M = (Ix1+Ix2) + i(Iy1+Iy2).
% Frequency evolution uses FID-A sign convention: exp(-1i*2*pi*H*t).

    Ix = 0.5*[0 1; 1 0];
    Iy = 0.5*[0 -1i; 1i 0];
    Iz = 0.5*[1 0; 0 -1];

    I1x = kron(Ix, eye(2)); I2x = kron(eye(2), Ix);
    I1y = kron(Iy, eye(2)); I2y = kron(eye(2), Iy);
    I1z = kron(Iz, eye(2)); I2z = kron(eye(2), Iz);

    w1 = -ppm_vec(1) * hzpppm;
    w2 = -ppm_vec(2) * hzpppm;

    H = w1*I1z + w2*I2z + J*(I1x*I2x + I1y*I2y + I1z*I2z);

    rho0 = I1z + I2z;

    U = expm(-1i*(pi/2)*(I1x + I2x));
    rho = U * rho0 * U';

    M = (I1x + I2x) + 1i*(I1y + I2y);

    [V,D] = eig(H);
    evals = diag(D);

    rhoe = V' * rho * V;
    Me   = V' * M   * V;

    C = rhoe .* (Me.');
    F = evals - evals.';

    Cvec = C(:);
    Fvec = F(:);

    E = exp(-1i*2*pi*(t(:) * Fvec.'));
    fids = E * Cvec;

    % Sign-convention fix: align with the singlet/lines branch (and with
    % FID-A / Siemens FID convention used by dat2mat_svs31p.m). Without
    % this conjugation the AB peaks come out at +ppm instead of -ppm,
    % because raising-operator transitions in a positive-w Hamiltonian
    % give negative F = E_i-E_j. The original LCModel .RAW writer hid
    % this by writing RF(:,2) = -imag(fids) (= conj(fids)).
    fids = conj(fids);
end
