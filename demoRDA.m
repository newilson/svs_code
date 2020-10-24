% fname = 'D:\CODE\MATLAB\neil\svs_code\20200914_DYER_MR_ES_CC017_19.rda';
fname = fullfile(pwd,'20200914_DYER_MR_ES_CC017_19.rda');

pars.fitmode = 3;

pars.lb = 3;

wsopts.type = 'hsvd';
wsopts.hsvd.bounds = [-0.025 0.025]; % normalized bounds for water
wsopts.hsvd.nsin = 25; % number of decaying sinusoids
wsopts.wavelet.zf = 0; % zero filling
wsopts.wavelet.scale = 7;
wsopts.wavelet.type = 'Daubechies';
wsopts.wavelet.par = 10;
wsopts.plt = true;
pars.wsopts = wsopts;

baseopts.method = 'pchip';
baseopts.stepsize = 0.5;
baseopts.windowsize = 0.5;
pars.baseopts = baseopts;

pars.den = 'wav';

pars.peaks = [];

output = postprocessSVS(fname,pars);