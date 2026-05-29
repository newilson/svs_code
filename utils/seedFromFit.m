function seeds = seedFromFit(prevFit, metabs, fitOpt)
% seedFromFit  Build per-component initial parameters from a prior basisVarpro fit.
%
% Usage:
%   seeds = seedFromFit(prevFit, metabs)
%   seeds = seedFromFit(prevFit, metabs, fitOpt)
%
% Inputs:
%   prevFit  - either the curvefitAuto_basisVarpro result struct (with
%              fields .pars.{shift,lbL,lbG,phase0,phase1}, .names,
%              .lineShapePars) OR a bru2mat_svs31p / dat2mat_svs31p
%              output struct that contains .basisVarproResults. The helper
%              auto-detects which level it was given.
%   metabs   - cell array of metab names for the *new* fit, in order
%              (typically pars.fit.metabs).
%   fitOpt   - (optional) struct with .shiftBounds, .lbLBounds, .lbGBounds,
%              .phaseBounds. If provided, seed values are clipped strictly
%              inside the bounds so the optimizer starts interior on every
%              component (avoids the warning + auto-clip in the wrapper).
%
% Output (struct):
%   .shiftInit, .lbLInit, .lbGInit, .phaseInit  - length(metabs) x 1 vectors.
%   .phase1Init                                 - scalar phi1 (deg/ppm).
%   .lineShapePars                              - prior lineshape side coefs
%                                                 (empty if not used).
%   .matched                                    - logical vector marking
%                                                 which components were
%                                                 found in prevFit.
%
% Components present in both prevFit and metabs get their converged values;
% components not in prevFit fall back to sensible interior defaults (0 for
% shift/phase0, 30 for lbL, 0 for lbG).
%
% Apply via:
%   pars.fit.fitOpt.shiftInit  = seeds.shiftInit;
%   pars.fit.fitOpt.lbLInit    = seeds.lbLInit;
%   pars.fit.fitOpt.lbGInit    = seeds.lbGInit;
%   pars.fit.fitOpt.phaseInit  = seeds.phaseInit;
%   pars.fit.fitOpt.phase1Init = seeds.phase1Init;   % only if wrapper supports it
%
% NW - 05/2026

if nargin < 3, fitOpt = struct(); end

% Accept either the bru2mat-level output (has .basisVarproResults) or the
% raw curvefitAuto_basisVarpro result.
if isfield(prevFit, 'basisVarproResults') && ~isempty(prevFit.basisVarproResults)
    src = prevFit.basisVarproResults;
else
    src = prevFit;
end

if ~isfield(src,'pars') || ~isfield(src,'names')
    error('seedFromFit: prevFit does not look like a basisVarpro result (missing .pars or .names).');
end

n = numel(metabs);
seeds.shiftInit = zeros(n,1);
seeds.lbLInit   = 30 * ones(n,1);
seeds.lbGInit   = zeros(n,1);
seeds.phaseInit = zeros(n,1);
seeds.matched   = false(n,1);

prevNames = src.names(:);
for k = 1:n
    idx = find(strcmp(prevNames, metabs{k}), 1);
    if isempty(idx), continue; end
    seeds.shiftInit(k) = src.pars.shift(idx);
    seeds.lbLInit(k)   = src.pars.lbL(idx);
    seeds.lbGInit(k)   = src.pars.lbG(idx);
    seeds.phaseInit(k) = src.pars.phase0(idx);
    seeds.matched(k)   = true;
end

% phi1 (scalar) and lineshape side coefs carry over straight.
if isfield(src.pars,'phase1') && ~isempty(src.pars.phase1)
    seeds.phase1Init = src.pars.phase1;
else
    seeds.phase1Init = 0;
end
if isfield(src,'lineShapePars')
    seeds.lineShapePars = src.lineShapePars;
else
    seeds.lineShapePars = [];
end

% Optional: clip to interior of provided bounds.
if ~isempty(fieldnames(fitOpt))
    seeds.shiftInit = clipInside(seeds.shiftInit, expandBounds(fitOpt,'shiftBounds',n,[-Inf Inf]));
    seeds.lbLInit   = clipInside(seeds.lbLInit,   expandBounds(fitOpt,'lbLBounds',  n,[0    Inf]));
    seeds.lbGInit   = clipInside(seeds.lbGInit,   expandBounds(fitOpt,'lbGBounds',  n,[0    Inf]));
    seeds.phaseInit = clipInside(seeds.phaseInit, expandBounds(fitOpt,'phaseBounds',n,[-Inf Inf]));
end

end

% -------------------------------------------------------------------------
function v = clipInside(v, B)
% Clip each row to its [lb ub] with a small interior margin so the
% optimizer starts strictly inside (lsqnonlin warns on equality with bound).
span   = B(:,2) - B(:,1);
margin = 1e-3 * max(span, 0);
margin(~isfinite(span)) = 0;
v = min(max(v, B(:,1) + margin), B(:,2) - margin);
end

% -------------------------------------------------------------------------
function B = expandBounds(fitOpt, name, n, fallback)
% Return n x 2 bounds. Accept scalar [lb ub], n x 2 matrix, or missing.
if isfield(fitOpt, name) && ~isempty(fitOpt.(name))
    B = fitOpt.(name);
    if size(B,1) == 1
        B = repmat(B, n, 1);
    end
else
    B = repmat(fallback, n, 1);
end
end
