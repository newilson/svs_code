function r = fitBloodTSP(ppm, spec, opt, f0)
% FITBLOODTSP  Locate and fit the TSP singlet in a processed blood 0 ppm
% spectrum: one Voigt + linear baseline, on the raw (unreferenced) axis.
% Shared by processU01_Blood_1H, the TSP T2 series and the T1 series.
%
%   ppm  - ppm axis
%   spec - complex (or real) spectrum on that axis
%   opt  - pars.tspfit struct (search, window, mode, centerHalf, widthInit,
%          widthBounds, ph_range, baseline)
%   f0   - 1H Larmor frequency (MHz), for the minimum-linewidth bound
%
% r fields: ppm, spec, fit, baseline, pars, parnames, area, areaWin,
%           center, width, located, window

ppm  = ppm(:);
spec = spec(:);

s  = ppm > opt.search(1) & ppm < opt.search(2);
xs = ppm(s);
[~, i] = max(real(spec(s)));
c   = xs(i);
win = c + opt.window;
m   = ppm > win(1) & ppm < win(2);
x   = ppm(m);
y   = double(real(spec(m)));

bo = struct('enable', opt.baseline, 'knotSpacing', diff(win), 'lambda', 1, ...
            'TolFun', 1e-8, 'TolX', 1e-8, ...
            'centerBounds', c + opt.centerHalf*[-1 1], ...
            'widthBounds', opt.widthBounds, ...
            'amplUB', max(abs(y(abs(x - c) < 0.1))));
ls  = struct('enable', false, 'areaPadFactor', 100);
out = curvefitAuto_varproBaseline(x, y, opt.mode, 1, c, opt.widthInit, ...
                                  opt.ph_range, 5/f0, bo, ls);

r = struct();
r.ppm      = x;
r.spec     = y;
r.fit      = out.fit;
r.baseline = out.baseline;
r.pars     = out.pars;
r.parnames = {out.parnames};
r.area     = out.areasConv;
r.areaWin  = trapz(x, out.fit(:) - out.baseline(:));
r.center   = out.pars(strcmp(out.parnames, 'pos'));
r.width    = out.pars(strcmp(out.parnames, 'width'));
r.located  = c;
r.window   = win;
end
