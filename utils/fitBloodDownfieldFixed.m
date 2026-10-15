function [res, ref] = fitBloodDownfieldFixed(ppm, spec, f0, knot, ref, prange)
% FITBLOODDOWNFIELDFIXED  Fixed-lineshape downfield fit for blood relaxometry
% series: complex-Lorentzian {Broad9_4, NADH2, NADH6, Trp} over a narrow
% window with a spline baseline.  The reference pass fits pos+width+phase
% freely on the highest-SNR spectrum; every other spectrum in the series is
% then fit AMPLITUDE-ONLY with those pinned, so peaks cannot trade area with
% each other or the baseline across the series.
%
%   ppm, spec - processed spectrum (any ppm ordering; real part is fit)
%   f0        - 1H Larmor frequency (MHz)
%   knot      - spline knot spacing (ppm), default 0.5
%   ref       - [] -> reference pass (returns ref); struct (.pos/.width/.phase)
%               -> fixed-lineshape amplitude-only pass
%   prange    - fit window (ppm), default [8.9 10.3]
%
% res fields: ppm, spec, fit, baseline, comp, names, areas (=areasConv), pos,
%             width, f0
% ref fields: pos, width, phase, names

if nargin < 4 || isempty(knot),   knot = 0.5;          end
if nargin < 5,                    ref  = [];           end
if nargin < 6 || isempty(prange), prange = [8.9 10.3]; end

ppm  = ppm(:);
spec = spec(:);
inds = ppm > min(prange) & ppm < max(prange);
xw = ppm(inds);
yw = double(real(spec(inds)));

names = {'Broad9_4','NADH2','NADH6','Trp'};   % no NADH4 (out of narrow window)
n = numel(names);
minw = 10/f0;
ph_range = [-60 60];

bo = struct();
bo.enable      = true;
bo.style       = 'osprey';
bo.knotSpacing = knot;
bo.lambda      = 1;
bo.normalize   = true;
bo.TolFun      = 1e-7;

ls = struct('enable', false, 'areaPadFactor', 30);

if isempty(ref)
    % -------- reference pass: free lineshape (within physical windows) -----
    center0 = [9.60 9.25 9.06 10.02];
    width0  = [0.07 0.05 0.05 0.10];
    bo.centerBounds = [9.50 9.70;    % Broad9_4
                       9.17 9.33;    % NADH2
                       8.98 9.14;    % NADH6
                       9.88 10.16];  % Trp
    bo.widthBounds  = [0.04 0.15;    % Broad9_4
                       0.02 0.09;    % NADH2
                       0.02 0.09;    % NADH6
                       0.05 0.30];   % Trp (broad, wider cap)
    bo.amplUB = 1.5*max(abs(yw));
    out = curvefitAuto_varproBaseline(xw, yw, 3, ones(1,n), center0, width0, ph_range, minw, bo, ls);
    pos   = out.pars(  (1:n)+n   );
    width = out.pars(  (1:n)+2*n );
    phase = out.pars(  (1:n)+3*n );
    ref = struct('pos',pos,'width',width,'phase',phase,'names',{names});
else
    % -------- fixed-lineshape pass: amplitude only -----
    pos = ref.pos(:).'; width = ref.width(:).'; phase = ref.phase(:).';
    bo.centerBounds = [pos.'-0.004, pos.'+0.004];
    bo.widthBounds  = [max(width.'*0.95, minw), width.'*1.05];
    bo.phaseBounds  = [phase.'-3, phase.'+3];
    bo.phaseInit    = phase;
    bo.amplUB       = max(abs(yw)) * 2;
    out = curvefitAuto_varproBaseline(xw, yw, 3, ones(1,n), pos, width, ph_range, minw, bo, ls);
end

res = struct();
res.ppm=xw; res.spec=yw; res.fit=out.fit; res.baseline=out.baseline;
res.comp=out.comp; res.names=names; res.areas=out.areasConv(:).'; res.f0=f0;
res.pos = out.pars((1:n)+n); res.width = out.pars((1:n)+2*n);
end
