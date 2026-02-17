%% demo_postprocessSVS_dynamic.m

clear,clc

%% 31-P data - BH
mainDir = 'C:\Users\CAMIPM-NW\Desktop\S139_01122026\E200';

recov1 = 'FID_TEST1D_31P_SLICE-SEL_TR5S_BASE1R_10X_0004';
pathname1 = [mainDir filesep recov1]; 

% time domain fitting
pars.denoiseSeries = 3;
% pars.fitmode = 10.0; % HSVD 0 Hz only
pars.nComp = 30;
pars.lb = 15; % Hz

% frequency domain fitting
pars.denoiseSeries = 3;
pars.fitmode = 7;
pars.peaklistPPM = 0;

output = postprocessSVS_dynamic(pathname1,pars);
