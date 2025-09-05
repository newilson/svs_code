%% demo dat2mat_svsNAD - frequency domain fitting
clear,clc

mainDir = 'Y:\users\rnanga\857606_CEW\857606_S01_06022025\datfiles';

ws_dat = '857606_S01_06022025_meas_MID00070_FID09895_svssel_latest021224_NAD_512avg.dat';
nws_dat = '857606_S01_06022025_meas_MID00069_FID09894_svssel_latest021224_WatRef_16avg.dat';

pars.lb = 3;
pars.eccopt = 0;

[ws,nws,ppm] = dat2mat_svsNAD([mainDir filesep ws_dat],[mainDir filesep nws_dat],pars);