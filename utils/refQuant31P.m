function out = refQuant31P(refName, refConc_mM, fitNames, fitAreas, metProps, TR_ms, TE_ms)
% REFQUANT31P  Absolute 31P quant against ONE internal reference metabolite.
%
%   out = refQuant31P(refName, refConc_mM, fitNames, fitAreas, metProps, ...
%                     TR_ms, TE_ms)
%
% Converts every fitted 31P amplitude to [mM] against the named internal
% reference, with per-metabolite T1/T2 relaxation correction:
%
%   [M] = (S_M / S_ref) * (R_ref / R_M) * (nP_ref / nP_M) * [ref]
%   R_x = (1 - exp(-TR/T1_x)) * exp(-TE/T2_x)
%
% Factored out of processU01_Brain_31P / processU01_Calf_31P so that a single
% scan can be reported against SEVERAL references (aATP and PCr) without
% duplicating the loop.  Switching reference rescales every concentration by
% one common factor -- it never changes a ratio between two metabolites in the
% same scan.  The point of reporting more than one is to expose how much that
% scale depends on which reference (and which assumed [ref] and T1/T2) you
% trust.
%
% No (1 - f_CSF) partial-volume term: with an internal reference that is
% itself CSF-absent, the CSF dilution cancels in the S_M/S_ref ratio.
%
% RefExt (the external phantom singlet) is skipped -- it is not a tissue
% metabolite and is quantified separately by extRefQuant31P.
%
% INPUTS
%   refName     e.g. 'aATP' or 'PCr'; must appear in BOTH fitNames and metProps
%   refConc_mM  assumed concentration of that reference in this tissue
%   fitNames    cellstr of fitted component names
%   fitAreas    per-spin-normalized amplitudes (amplScaled) matching fitNames
%   metProps    struct array with .name .T1_ms .T2_ms .nPhos
%   TR_ms,TE_ms sequence timing
%
% OUTPUT
%   out.reference    struct(name, conc_mM, T1_ms, T2_ms, nPhos, R)
%   out.ok           false when the reference was not fitted (all conc NaN)
%   out.metabolites  struct array(name, area, conc_mM, T1_ms, T2_ms, nPhos,
%                                 R, ratio)

out = struct();
out.ok = false;
out.metabolites = struct('name',{},'area',{},'conc_mM',{}, ...
                         'T1_ms',{},'T2_ms',{},'nPhos',{},'R',{},'ratio',{});

iF = find(strcmp(fitNames, refName), 1);
iP = find(strcmp({metProps.name}, refName), 1);
if isempty(iF) || isempty(iP)
    warning(['refQuant31P: reference %s not available (fit:%d metProps:%d); ' ...
             'this reference set is skipped.'], refName, ~isempty(iF), ~isempty(iP));
    out.reference = struct('name',refName,'conc_mM',refConc_mM, ...
                           'T1_ms',NaN,'T2_ms',NaN,'nPhos',NaN,'R',NaN);
    return
end

areaRef = fitAreas(iF);
T1ref   = metProps(iP).T1_ms;
T2ref   = metProps(iP).T2_ms;
nPref   = metProps(iP).nPhos;
Rref    = relax(TR_ms, TE_ms, T1ref, T2ref);

out.reference = struct('name',refName,'conc_mM',refConc_mM, ...
                       'T1_ms',T1ref,'T2_ms',T2ref,'nPhos',nPref,'R',Rref);

if ~(areaRef > 0)
    warning(['refQuant31P: reference %s amplitude <= 0 (%.3g); concentrations ' ...
             'against it are unreliable.'], refName, areaRef);
end

for k = 1:numel(fitNames)
    nm_k = fitNames{k};
    if strcmp(nm_k,'RefExt'), continue; end
    pIdx = find(strcmp({metProps.name}, nm_k), 1);
    if isempty(pIdx)
        warning('refQuant31P: no T1/T2/nPhos entry for %s -- skipping', nm_k);
        continue
    end
    T1_k = metProps(pIdx).T1_ms;
    T2_k = metProps(pIdx).T2_ms;
    nP_k = metProps(pIdx).nPhos;
    R_k  = relax(TR_ms, TE_ms, T1_k, T2_k);

    ratio   = fitAreas(k) / areaRef;
    conc_mM = ratio * (Rref / R_k) * (nPref / nP_k) * refConc_mM;

    out.metabolites(end+1) = struct( ...
        'name', nm_k, 'area', fitAreas(k), 'conc_mM', conc_mM, ...
        'T1_ms', T1_k, 'T2_ms', T2_k, 'nPhos', nP_k, ...
        'R', R_k, 'ratio', ratio);
end

out.ok = true;
end

% =====================================================================
function R = relax(TR_ms, TE_ms, T1_ms, T2_ms)
R = (1 - exp(-TR_ms / T1_ms)) * exp(-TE_ms / T2_ms);
end
