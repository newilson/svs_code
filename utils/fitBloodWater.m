function r = fitBloodWater(ppm, spec, watfit, f0)
% FITBLOODWATER  Fit one blood water spectrum with the shared blood water
% model (bloodWatfitPars).  Mirrors the watfit='fit' branch of
% bru2mat_svsNAD so a single-TE reference scan and a relaxometry series are
% fit identically.
%
%   ppm    - ppm axis
%   spec   - complex (or real) spectrum on that axis
%   watfit - bloodWatfitPars() struct
%   f0     - 1H Larmor frequency (MHz), for the minimum-linewidth bound
%
% r fields: ppm, spec, fit, baseline, comp, pars, parnames, areasConv,
%           area (= sum of areasConv), areaWin

ppm  = ppm(:);
spec = spec(:);
m = ppm > min(watfit.ppm_range) & ppm < max(watfit.ppm_range);
x = ppm(m);
y = double(real(spec(m)));

bo = watfit.baselineOpt;
bo.TolFun       = 1e-8;
bo.TolX         = 1e-8;
bo.weights      = [];
bo.centerBounds = watfit.centerBounds;
bo.knotSpacing  = max(bo.knotSpacing, diff(watfit.ppm_range));
bo.amplUB       = max(abs(y));

out = curvefitAuto_varproBaseline(x, y, watfit.mode, ones(1,watfit.nComp), ...
    watfit.centerInit, watfit.widthInit, watfit.ph_range, 10/f0, bo, watfit.lineShapeOpt);

r = struct();
r.ppm       = x;
r.spec      = y;
r.fit       = out.fit;
r.baseline  = out.baseline;
r.comp      = out.comp;
r.pars      = out.pars;
r.parnames  = {out.parnames};
r.areasConv = out.areasConv;
r.area      = sum(out.areasConv);
r.areaWin   = trapz(x, out.fit(:) - out.baseline(:));
end
