%% demo dat2mat_svsNAD - frequency domain fitting
clear,clc

mainDir = 'Z:\Data\Spectroscopy\4Nanga\857606_CEW\857606_S01_06022025\datfiles';

ws_dat = '857606_S01_06022025_meas_MID00070_FID09895_svssel_latest021224_NAD_512avg.dat';
nws_dat = '857606_S01_06022025_meas_MID00069_FID09894_svssel_latest021224_WatRef_16avg.dat';

use_nws = true;

pars.lb = 3;
pars.eccopt = 0; % -1: Brown, 0: Klose ECC
pars.PC = []; % []: manual phase correction (make sure to 'SAVE' and close figure when done), 1x1 value: constant phase correction (radians), 2x1: [ph0 pc1] use pars.PCpivot to control pivot point
pars.base = 'fit'; % nan: no baseline correction, 'full': semi-manual correction on full spectrum, 'fit': semi-manual correction on fitted range only (make sure to 'SAVE' and close figure when done)
pars.plt = true;

% USE ONLY MANUAL FITTING FOR NOW
pars.dofit = 'man'; % false: no fit, 'man': semi-manual fit, 'auto': auto fit, 'lcmodel': writes RDA files for LCmodel fitting
pars.fit.mode = 5; % 1: gaussian, 2: lorentzian, 3: complex Lorentzian, 4: complex Gaussian, 5: complex Voigt
pars.fit.peaks = {'NADH2','NADH6','NADH4'};
pars.fit.ph_range = [-35 35];
pars.fit.ppm_range = [8.7 9.5];

if use_nws
    [output,pars,hdr] = dat2mat_svsNAD(pars,[mainDir filesep ws_dat],[mainDir filesep nws_dat]);
else
    [output,pars,hdr] = dat2mat_svsNAD(pars,[mainDir filesep ws_dat],[]);
end
% figure, plot(twixOut.ppm,real(twixOut.met.spec),'-r'), set(gca,'xdir','reverse'), title('Twix')
% axis tight

% saving
outName = [ws_dat(1:end-4) '.mat'];
save(outName,'output','pars','hdr');

% csv file
csvName = [ws_dat(1:end-4) 'fit_results.csv'];
T = table([output.wat.fit.pars(2); output.wat.fit.areas],[output.met.fit.pars(4:6); output.met.fit.areas],'RowNames',{'ppm','area'},'VariableNames',{'Wat','NAD'});
T.Properties.DimensionNames = {'Var', 'Variables'};
writetable(T,csvName,'WriteRowNames',true,'WriteVariableNames',true);