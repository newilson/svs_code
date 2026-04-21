function absQ = absoluteQuant_svs(data, seq, voxel)
% absQ = absoluteQuant_svs(data, seq, voxel)
%
% Absolute metabolite concentration from water-referenced 1H MRS, using the
% tissue-segmented formulation of Gasparovic et al., MRM 2006;55:1219-1226.
%
%   [M] = (S_M / S_W) * (N_W / N_M) * ...
%         ( sum_i f_i * rho_i * R_{W,i} ) / ( (1 - f_CSF) * R_M )
%
% with relaxation attenuation
%   R_{W,i} = (1 - exp(-TR_w/T1_{W,i})) * exp(-TE_w/T2_{W,i})
%   R_M     = (1 - exp(-TR/T1_M))      * exp(-TE/T2_M)
%
% The (1 - f_CSF) factor assumes the metabolite is confined to tissue
% parenchyma.  Setting CSF fraction to zero recovers the single-tissue /
% un-segmented case.
%
% -------------------------------------------------------------------------
% INPUT STRUCTURES
% -------------------------------------------------------------------------
% data (metabolite + water reference measurements)
%   .metabolites : struct array, one entry per metabolite peak
%       .name       : (optional) string label
%       .area       : peak area from fit (same arbitrary units as water)
%       .nProtons   : number of equivalent protons contributing to the peak
%       .T1_ms      : metabolite T1 in ms
%       .T2_ms      : metabolite T2 in ms
%   .water
%       .area       : water peak area (same units as metabolite areas)
%       .nAvg       : (optional) number of averages, if different from met.
%   .nAvg          : (optional) number of averages for the metabolite scan
%
% seq (sequence timing in ms)
%   .TR_ms          : TR for metabolite scan
%   .TE_ms          : TE for metabolite scan
%   .TR_ws_ms       : (optional) TR for water reference, defaults to TR_ms
%   .TE_ws_ms       : (optional) TE for water reference, defaults to TE_ms
%
% voxel (tissue composition)
%   .tissues        : struct array with one entry per tissue compartment
%       .name           : e.g. 'WM', 'GM', 'CSF'
%       .fraction       : volume fraction (0..1); fractions should sum ~ 1
%       .T1_ms          : water T1 in this tissue (ms)
%       .T2_ms          : water T2 in this tissue (ms)
%       .waterConc_M    : water concentration in this tissue (mol/L of tissue)
%                         Typical values: WM 36.08, GM 43.30, CSF 55.56
%       .flag_metSig    : (optional, logical) true for compartments that
%                         DO contain the metabolite (i.e. parenchyma).
%                         If absent, defaults to ~strcmpi(name,'CSF').
%
%   For a non-segmented / single-tissue voxel, pass a single tissue entry
%   with fraction = 1 and flag_metSig = true.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% absQ.metabolites(k)
%       .name        : metabolite label
%       .conc_mM     : absolute concentration in mmol per L of tissue
%                      (parenchyma, i.e. CSF excluded)
%       .R_M         : relaxation attenuation factor applied
%       .ratio       : raw S_M / S_W used
% absQ.water
%       .Rw_per_tissue : vector of R_{W,i} aligned with voxel.tissues
%       .weightedTerm  : sum_i f_i * rho_i * R_{W,i}   (mol/L)
%       .fCSF          : total CSF-like fraction
% absQ.inputs : echo of the effective inputs actually used
%
% -------------------------------------------------------------------------

    % ---- sequence timing ------------------------------------------------
    TR  = seq.TR_ms;
    TE  = seq.TE_ms;
    if isfield(seq,'TR_ws_ms') && ~isempty(seq.TR_ws_ms)
        TRw = seq.TR_ws_ms;
    else
        TRw = TR;
    end
    if isfield(seq,'TE_ws_ms') && ~isempty(seq.TE_ws_ms)
        TEw = seq.TE_ws_ms;
    else
        TEw = TE;
    end

    % ---- tissue table ---------------------------------------------------
    tissues = voxel.tissues;
    nT = numel(tissues);
    f       = zeros(nT,1);
    T1w     = zeros(nT,1);
    T2w     = zeros(nT,1);
    rho         = zeros(nT,1);
    flag_metSig = true(nT,1);
    for i = 1:nT
        f(i)   = tissues(i).fraction;
        T1w(i) = tissues(i).T1_ms;
        T2w(i) = tissues(i).T2_ms;
        rho(i) = tissues(i).waterConc_M;     % mol/L of that tissue
        if isfield(tissues(i),'flag_metSig') && ~isempty(tissues(i).flag_metSig)
            flag_metSig(i) = logical(tissues(i).flag_metSig);
        elseif isfield(tissues(i),'name') && ~isempty(tissues(i).name)
            flag_metSig(i) = ~strcmpi(tissues(i).name,'CSF');
        end
    end

    if abs(sum(f) - 1) > 1e-3
        warning('absoluteQuant_svs:fractions', ...
            'Tissue fractions sum to %.4f, not 1. Check segmentation inputs.', sum(f));
    end
    % fraction NOT containing metabolite = complement of flag_metSig tissues
    fCSF = sum(f(~flag_metSig));

    % ---- water-side relaxation attenuation per tissue -------------------
    Rw = (1 - exp(-TRw ./ T1w)) .* exp(-TEw ./ T2w);

    % sum_i f_i * rho_i * R_{W,i}   (mol/L)
    waterTerm = sum(f .* rho .* Rw);

    % ---- optional averages correction -----------------------------------
    % If water reference and metabolite scans used different averaging,
    % scale the water area so both areas represent "per-average" signal.
    Sw_raw = data.water.area;
    nAvg_m = getfield_or(data,         'nAvg', []);
    nAvg_w = getfield_or(data.water,   'nAvg', []);
    if ~isempty(nAvg_m) && ~isempty(nAvg_w) && nAvg_m ~= nAvg_w
        Sw = Sw_raw * (nAvg_m / nAvg_w);
    else
        Sw = Sw_raw;
    end

    % ---- per-metabolite concentration -----------------------------------
    N_W = 2;                                % two protons per water
    mets = data.metabolites;
    nM = numel(mets);
    outMet = repmat(struct('name','','conc_mM',NaN,'R_M',NaN,'ratio',NaN), nM, 1);

    denom_tissueOnly = max(1 - fCSF, eps);  % avoid /0 for pure-CSF voxel

    for k = 1:nM
        m = mets(k);
        Rm = (1 - exp(-TR/m.T1_ms)) * exp(-TE/m.T2_ms);
        ratio = m.area / Sw;
        conc_M = ratio * (N_W / m.nProtons) * waterTerm ...
                       / (denom_tissueOnly * Rm);

        if isfield(m,'name'), outMet(k).name = m.name; end
        outMet(k).conc_mM = conc_M * 1e3;   % mol/L -> mmol/L
        outMet(k).R_M     = Rm;
        outMet(k).ratio   = ratio;
    end

    % ---- package output -------------------------------------------------
    absQ.metabolites        = outMet;
    absQ.water.Rw_per_tissue = Rw(:)';
    absQ.water.weightedTerm  = waterTerm;
    absQ.water.fCSF          = fCSF;
    absQ.inputs = struct('TR_ms',TR,'TE_ms',TE, ...
                         'TR_ws_ms',TRw,'TE_ws_ms',TEw, ...
                         'fractions',f(:)','flag_metSig',flag_metSig(:)', ...
                         'Sw_used',Sw);
end

function v = getfield_or(s, fld, default)
    if isstruct(s) && isfield(s,fld) && ~isempty(s.(fld))
        v = s.(fld);
    else
        v = default;
    end
end
