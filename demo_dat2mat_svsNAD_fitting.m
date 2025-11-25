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

%% 
clear,clc

pars.lb = 3;
pars.eccopt = -1; % -1: Brown, 0: Klose ECC
pars.PC = 0; % []: manual phase correction (make sure to 'SAVE' and close figure when done), 1x1 value: constant phase correction (radians), 2x1: [ph0 pc1] use pars.PCpivot to control pivot point
pars.base = nan; % nan: no baseline correction, 'full': semi-manual correction on full spectrum, 'fit': semi-manual correction on fitted range only (make sure to 'SAVE' and close figure when done)
pars.plt = true;
pars.dofit = false;


scan = 20.2;
switch scan 
    case 10.1
        ws_dat = 'C:\Users\CAMIPM-NW\Downloads\datain\datain\854430_S010_02082024\datfiles\854430_S010_02082024_meas_MID01655_FID95761_svssel_NAD_TR1000ms_512avg.dat';
        nws_dat = 'C:\Users\CAMIPM-NW\Downloads\datain\datain\854430_S010_02082024\datfiles\854430_S010_02082024_meas_MID01654_FID95760_svssel_WATEREF.dat';

    case 10.2
        ws_dat = 'C:\Users\CAMIPM-NW\Downloads\datain\datain\854430_S010_02132024\datfiles\854430_S010_02132024_meas_MID02144_FID96293_svssel_NAD_TR1000ms_512avg.dat';
        nws_dat = 'C:\Users\CAMIPM-NW\Downloads\datain\datain\854430_S010_02132024\datfiles\854430_S010_02132024_meas_MID02143_FID96292_svssel_WATEREF.dat';

    case 20.2
        ws_dat = 'C:\Users\CAMIPM-NW\Downloads\datain\datain\854430_S020_07262024\datfiles\854430_S020_07262024_meas_MID02075_FID14271_svssel_NAD_TR1000ms_512avg.dat';
        nws_dat = 'C:\Users\CAMIPM-NW\Downloads\datain\datain\854430_S020_07262024\datfiles\854430_S020_07262024_meas_MID02074_FID14270_svssel_WATEREF.dat';

end


[output,pars,hdr] = dat2mat_svsNAD(pars,ws_dat,nws_dat);
figure, plot(output.ppm, real(output.met.spec)); set(gca,'xdir','reverse'), title(num2str(scan))
