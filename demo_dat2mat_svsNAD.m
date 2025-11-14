%% demo dat2mat_svsNAD - frequency domain fitting
clear,clc

mainDir = 'Z:\Data\Spectroscopy\4Nanga\857606_CEW\857606_S01_06022025\datfiles';

ws_dat = '857606_S01_06022025_meas_MID00070_FID09895_svssel_latest021224_NAD_512avg.dat';
nws_dat = '857606_S01_06022025_meas_MID00069_FID09894_svssel_latest021224_WatRef_16avg.dat';

use_nws = true;

pars.lb = 3;
pars.eccopt = 0; % -1: Brown, 0: Klose ECC
pars.PC = 0; % []: manual phase correction (make sure to 'SAVE' and close figure when done), 1x1 value: constant phase correction (radians), 2x1: [ph0 pc1] use pars.PCpivot to control pivot point
pars.base = true; % nan: no baseline correction, true: semi-manual correction (make sure to 'SAVE' and close figure when done)


if use_nws
    [twixOut,pars,hdr] = dat2mat_svsNAD(pars,[mainDir filesep ws_dat],[mainDir filesep nws_dat]);
else
    [twixOut,pars,hdr] = dat2mat_svsNAD(pars,[mainDir filesep ws_dat],[]);
end
% figure, plot(twixOut.ppm,real(twixOut.met.spec),'-r'), set(gca,'xdir','reverse'), title('Twix')
% axis tight

dcmPars.lb = 3;
dcmPars.base = nan;
dcmPars.PC = 0;
dcmPars.peaks = [];
dcmPars.fitmode = 0;
dcmPars.wsopts.type = 'none';
dcm = 'Z:\Data\Spectroscopy\4Nanga\857606_CEW\857606_S01_06022025\E100\SVSSEL_LATEST021224_NAD_512AVG_0014\857606_S01.MR.CAMRIS_BACKEDUP_USETHIS_WIERS.0014.0001.2025.06.02.11.21.26.820775.405632560.IMA';
dcmOut = postprocessSVS(dcm,dcmPars);
dcmOut.ppm = dcmOut.ppm - mean(dcmOut.ppm) + 4.72;
figure, plot(dcmOut.ppm,real(dcmOut.spec)/max(abs(real(dcmOut.spec))),'b-',twixOut.ppm,real(twixOut.met.spec)/max(abs(real(twixOut.met.spec))),'r-'), set(gca,'xdir','reverse'), legend({'dicom','raw'})
axis tight

% check phase
addpath('utils')
NWman_phase(twixOut.met.spec/max(abs(twixOut.met.spec)),twixOut.ppm,'test',dcmOut.spec/max(abs(dcmOut.spec)),dcmOut.ppm);

%% small voxel NAD
clear,clc

mainDir = 'Y:\users\nwilson\data\NAD_smallVox\MATHUR_25_09_12\RESEARCH_NWILSON_20250912';

% dicom file
dcm = 'SVSSEL_LATEST_SL_NAD_0015\MATHUR.MR.RESEARCH_NWILSON.0015.0001.2025.09.12.12.01.37.732089.553313816.IMA';

% no hsvd
dcmPars.fitmode = 0;
dcmPars.lb = 5;
dcmPars.base = nan;
dcmPars.PC = 0;
dcmPars.peaks = [];
dcmPars.wsopts.type = 'filt';
dcmPars.wsopts.hsvd.bounds = [-0.015 0.015]; % normalized bounds for water
dcmPars.wsopts.hsvd.nsin = 25; % number of decaying sinusoids
dcmPars.wsopts.filt.type = 'average';
dcmPars.wsopts.filt.N = 11;
dcmPars.wsopts.plt = true;
[dcmOut1, dcmPars] = postprocessSVS(fullfile(mainDir,dcm),dcmPars);

% with hsvd
dcmPars.den = 'hsvd';
dcmOut2 = postprocessSVS(fullfile(mainDir,dcm),dcmPars);

% raw data file
twixWS = 'meas_MID00224_FID23453_svssel_latest_sl_NAD.dat';
twixNWS = 'meas_MID00212_FID23441_svssel_latest_sl_WAT.dat';

pars.lb = 5;
pars.base = nan;
pars.pc = 0;
pars.wsopts.type = 'filt';
pars.wsopts.filt.type = 'average';
pars.wsopts.filt.N = 11;
pars.wsopts.plt = true;
pars.doWatSupPre = true;





