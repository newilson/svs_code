%% demo dat2mat_svsNAD - frequency domain fitting
clear,clc


%%  Case 1
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

%% Case 2
clear,clc

pars.lb = 0;
pars.eccopt = -1; % -1: Brown, 0: Klose ECC
pars.PC = 0; % []: manual phase correction (make sure to 'SAVE' and close figure when done), 1x1 value: constant phase correction (radians), 2x1: [ph0 pc1] use pars.PCpivot to control pivot point
pars.base = nan; % nan: no baseline correction, 'full': semi-manual correction on full spectrum, 'fit': semi-manual correction on fitted range only (make sure to 'SAVE' and close figure when done)
pars.plt = true;
pars.dofit = false;
pars.ccopt.minsig_frac = 0;
pars.WatSupPre.opts.type = 'none';
pars.WatSupPost.opts.type = 'none';
pars.block_size = []; % block averaging seems not to do anything - CHECK THIS LATER


scan = 10.1;
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
figure, 
plot(output.ppm, real(output.met.spec)); set(gca,'xdir','reverse'), title(num2str(scan))
% temp1 = real(output.met.spec);

%% Case 3 - varpro with spline baselines and hsvd initialization
clear,clc

pars.lb = 3;
pars.eccopt = -1; % -1: Brown, 0: Klose ECC
pars.PC = 0;
pars.base = nan;
pars.plt = true;
pars.block_size = [];
pars.ccopt.minsig_frac = 0.05;

pars.WatSupPre.opts.type = 'none';
pars.WatSupPost.opts.type = 'filt';
pars.WatSupPost.opts.filt.type = 'average';
pars.WatSupPost.opts.filt.N = 11;

pars.dofit = 'varpro';
pars.fit.mode = 5; % Voigt
% pars.fit.ppm_range = [8.8 10.5];
pars.fit.ppm_range = [8.9 9.7];
pars.fit.ph_range = [-35 35];

pars.peaks.name  = {'NADH2','NADH6','NADH4','Trp'};
pars.peaks.range = [9.25 9.45; 9.05 9.25; 8.8 9.05; 10.0 10.3];

pars.fit.baselineOpt.enable = true;
pars.fit.baselineOpt.knotSpacing = 2;
pars.fit.baselineOpt.lambda = 1;
pars.fit.baselineOpt.lambdaAmpl = 0.1;
pars.fit.baselineOpt.TolFun = 1e-12;
pars.fit.baselineOpt.TolX   = 1e-12;

pars.fit.lineShapeOpt.enable = true;
pars.fit.lineShapeOpt.nSide = 2;
pars.fit.lineShapeOpt.asymmetric = true;
pars.fit.lineShapeOpt.maxSide = 0.50;

% water area estimation: 'integrate', 'fid', or 'fit'
pars.watfit.method = 'integrate';
pars.watfit.ppm_range = [4 5.5];

scan = 10.1;
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

% display individual peak components
figure
plot(output.met.fit.ppm, real(output.met.fit.spec), 'k'), hold on
plot(output.met.fit.ppm, output.met.fit.spec_fit, 'r', 'LineWidth', 1.5)
plot(output.met.fit.ppm, output.met.fit.baseline, 'b--')
for kk = 1:size(output.met.fit.comp, 2)
    plot(output.met.fit.ppm, output.met.fit.comp(:,kk), '--')
end
plot(output.met.fit.ppm, real(output.met.fit.spec) - output.met.fit.spec_fit, 'g')
legend([{'data','fit','baseline'}, pars.peaks.name, {'residual'}])
set(gca,'xdir','reverse')
title(['VARPRO fit - scan ' num2str(scan)])

%% Case 4 - LCModel vs Varpro comparison across subjects
clear, clc

dataRoot = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_RR';

% --- Common parameters ---
pars.lb = 3;
pars.eccopt = -1;
pars.PC = 0; % no longer needed when doing 2-step fit
pars.base = nan;
pars.plt = true;
pars.pltSpec = false;
pars.block_size = [];
pars.ccopt.minsig_frac = 0.05;

pars.WatSupPre.opts.type = 'none';
pars.WatSupPost.opts.type = 'none';
pars.WatSupPost.opts.filt.type = 'average';
pars.WatSupPost.opts.filt.N = 11;
pars.WatSupPost.opts.hsvd.bounds = [-1 1] * 300/8000; % +/- 1 ppm away from center assuming 7T (300 Hz/ppm) and 4000 Hz spectral bandwidth
pars.WatSupPost.opts.hsvd.nsin = 25;

% fit parameters
pars.fit.mode = 5;
pars.fit.ppm_range = [8.9 9.7];
% pars.fit.ph_range = [-35 35];
pars.fit.ph_range = [-60 60];
pars.peaks.name  = {'NADH2','NADH6','NADH4'};
pars.peaks.range = [9.25 9.45; 9.05 9.25; 8.8 9.05];
pars.fit.baselineOpt.enable = true;
pars.fit.baselineOpt.knotSpacing = 2;
pars.fit.baselineOpt.lambda = 1;
pars.fit.baselineOpt.lambdaAmpl = 0.1;
pars.fit.baselineOpt.TolFun = 1e-6;
pars.fit.lineShapeOpt.enable = true;
pars.fit.lineShapeOpt.nSide = 2;
pars.fit.lineShapeOpt.asymmetric = true;
pars.fit.lineShapeOpt.maxSide = 0.50;

% water area estimation: 'integrate', 'fid', or 'fit'
pars.watfit.method = 'fid';
pars.watfit.ppm_range = [4 5.5];

% =====================================================================
% Define scan groups
% Each group: same subject, same body part, dates within a few weeks
% =====================================================================
groupsToProcess = [5, 17]; % <-- set which groups to run

nGrp = 0;

% --- BRAIN SCANS ---

% Group 1: S008, brain, 07/21/2025
if ismember(1, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S008 brain 07/21/2025';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S008_07212025_31P-1H-HeadCoil', 'datfiles', 'S08_07212025_meas_MID00285_FID16340_svssel_latest_M0_NAD_9_1ppm_256avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S008_07212025_31P-1H-HeadCoil', 'datfiles', 'S08_07212025_meas_MID00283_FID16338_svssel_latest_M0_Water_1avg_tra.dat');
groups(nGrp).ws  = {ws_dat};
groups(nGrp).nws = {nws_dat};
end

% Group 2: S076, brain, repeat scans (02/05 - 03/10/2026)
if ismember(2, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S076 brain 02-03/2026';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S076_02052026_50mm_31P-1H-HeadCoil', 'datfiles', 'S76_02052026_50mm_meas_MID00223_FID07639_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S076_02052026_50mm_31P-1H-HeadCoil', 'datfiles', 'S76_02052026_50mm_meas_MID00220_FID07636_svssel_latest_M0_Water_4avg_tra.dat');
ws_dat2  = fullfile(dataRoot, 'S076_02122026_31P-1H-HeadCoil', 'datfiles', 'S076_02122026_meas_MID00932_FID08354_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat2 = fullfile(dataRoot, 'S076_02122026_31P-1H-HeadCoil', 'datfiles', 'S076_02122026_meas_MID00930_FID08352_svssel_latest_M0_Water_4avg_tra.dat');
ws_dat3  = fullfile(dataRoot, 'S076_02192026_31P-1H-HeadCoil', 'datfiles', 'S076_02192026_meas_MID00615_FID09102_svssel_latest_M0_NAD_9_7ppm_512avg_tra_voxel_sinc.dat');
nws_dat3 = fullfile(dataRoot, 'S076_02192026_31P-1H-HeadCoil', 'datfiles', 'S076_02192026_meas_MID00613_FID09100_svssel_latest_M0_Water_4avg_tra.dat');
ws_dat4  = fullfile(dataRoot, 'S076_03102026_31P-1H-HeadCoil', 'datfiles', 'S076_03102026_meas_MID00459_FID10844_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat4 = fullfile(dataRoot, 'S076_03102026_31P-1H-HeadCoil', 'datfiles', 'S076_03102026_meas_MID00458_FID10843_svssel_latest_M0_Water_4avg_tra.dat');
groups(nGrp).ws  = {ws_dat, ws_dat2, ws_dat3, ws_dat4};
groups(nGrp).nws = {nws_dat, nws_dat2, nws_dat3, nws_dat4};
end

% Group 3: S092, brain, surface coil, 05/21/2025
if ismember(3, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S092 brain surface 05/21/2025';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S092_05212025_31P-1H-SurfaceCoil', 'datfiles', 'S092_05212025_meas_MID00345_FID08602_svssel_latest_M0_NAD_9_7ppm_256avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S092_05212025_31P-1H-SurfaceCoil', 'datfiles', 'S092_05212025_meas_MID00342_FID08599_svssel_latest_M0_Water_1avg_tra.dat');
groups(nGrp).ws  = {ws_dat};
groups(nGrp).nws = {nws_dat};
end

% Group 4: S092, brain, head coil, 07/17/2025
if ismember(4, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S092 brain head 07/17/2025';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S092_07172025_31P-1H-HeadCoil', 'datfiles', 'S092_07172025_meas_MID00474_FID16046_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S092_07172025_31P-1H-HeadCoil', 'datfiles', 'S092_07172025_meas_MID00471_FID16043_svssel_latest_M0_Water_1avg_tra.dat');
groups(nGrp).ws  = {ws_dat};
groups(nGrp).nws = {nws_dat};
end

% Group 5: S092, brain, head coil, repeat scans (03/03 - 03/05/2026)
if ismember(5, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S092 brain head 03/2026';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S092_03032026_31P-1H-HeadCoil', 'datfiles', 'S092_03032026_meas_MID01621_FID10109_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S092_03032026_31P-1H-HeadCoil', 'datfiles', 'S092_03032026_meas_MID01620_FID10108_svssel_latest_M0_Water_4avg_tra.dat');
ws_dat2  = fullfile(dataRoot, 'S092_03052026_31P-1H-HeadCoil', 'datfiles', 'S092_03052026_meas_MID01808_FID10296_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat2 = fullfile(dataRoot, 'S092_03052026_31P-1H-HeadCoil', 'datfiles', 'S092_03052026_meas_MID01806_FID10294_svssel_latest_M0_Water_4avg_tra.dat');
groups(nGrp).ws  = {ws_dat, ws_dat2};
groups(nGrp).nws = {nws_dat, nws_dat2};
end

% Group 6: S095, brain, surface coil, 05/22/2025
if ismember(6, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S095 brain surface 05/22/2025';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S095_05222025_31P-1H-SurfaceCoil', 'datfiles', 'S095_05222025_meas_MID00038_FID08644_svssel_latest_M0_NAD_9_7ppm_512avg_cor.dat');
nws_dat = fullfile(dataRoot, 'S095_05222025_31P-1H-SurfaceCoil', 'datfiles', 'S095_05222025_meas_MID00037_FID08643_svssel_latest_M0_Water_1avg_cor.dat');
groups(nGrp).ws  = {ws_dat};
groups(nGrp).nws = {nws_dat};
end

% Group 7: S095, brain, head coil, 06/17/2025
if ismember(7, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S095 brain head 06/17/2025';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S095_06172025_31P-1H-HeadCoil', 'datfiles', 'S095_06172025_meas_MID01937_FID12107_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S095_06172025_31P-1H-HeadCoil', 'datfiles', 'S095_06172025_meas_MID01935_FID12105_svssel_latest_M0_Water_1avg_tra.dat');
groups(nGrp).ws  = {ws_dat};
groups(nGrp).nws = {nws_dat};
end

% Group 8: S095, brain, head coil, 07/23/2025
if ismember(8, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S095 brain head 07/23/2025';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S095_07232025_31P-1H-HeadCoil', 'datfiles', 'S095_07232025_meas_MID00633_FID16690_svssel_latest_M0_NAD_9_7ppm_256avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S095_07232025_31P-1H-HeadCoil', 'datfiles', 'S095_07232025_meas_MID00631_FID16688_svssel_latest_M0_Water_1avg_tra.dat');
groups(nGrp).ws  = {ws_dat};
groups(nGrp).nws = {nws_dat};
end

% Group 9: S095, brain, head coil, 11/25/2025
if ismember(9, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S095 brain head 11/25/2025';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S95_11252025_31P-1H-HeadCoil', 'datfiles', 'S95_11252025_meas_MID01286_FID02750_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S95_11252025_31P-1H-HeadCoil', 'datfiles', 'S95_11252025_meas_MID01285_FID02749_svssel_latest_M0_Water_1avg_tra.dat');
groups(nGrp).ws  = {ws_dat};
groups(nGrp).nws = {nws_dat};
end

% Group 10: S138, brain, head coil, 11/20/2025
if ismember(10, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S138 brain head 11/20/2025';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S138_11202025_31P-1H-HeadCoil', 'datfiles', '11202025_S138_meas_MID00797_FID02261_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S138_11202025_31P-1H-HeadCoil', 'datfiles', '11202025_S138_meas_MID00796_FID02260_svssel_latest_M0_Water_1avg_tra.dat');
groups(nGrp).ws  = {ws_dat};
groups(nGrp).nws = {nws_dat};
end

% Group 11: S138, brain, head coil, repeat scans (02/17 - 03/19/2026)
if ismember(11, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S138 brain head 02-03/2026';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S138_02172026_31P-1H-HeadCoil', 'datfiles', 'S138_02172026_meas_MID00379_FID08871_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S138_02172026_31P-1H-HeadCoil', 'datfiles', 'S138_02172026_meas_MID00378_FID08870_svssel_latest_M0_Water_4avg_tra.dat');
ws_dat2  = fullfile(dataRoot, 'S138_03192026_31P-1H-HeadCoil', 'datfiles', 'S138_03192026_meas_MID00036_FID11602_svssel_latest_M0_NAD_9_7ppm_1024avg_tra.dat');
nws_dat2 = fullfile(dataRoot, 'S138_03192026_31P-1H-HeadCoil', 'datfiles', 'S138_03192026_meas_MID00035_FID11601_svssel_latest_M0_Water_8avg_tra.dat');
groups(nGrp).ws  = {ws_dat, ws_dat2};
groups(nGrp).nws = {nws_dat, nws_dat2};
end

% Group 12: S139, brain, head coil, repeat scans (01/29 - 02/03/2026)
if ismember(12, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S139 brain head 01-02/2026';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S139_01292026_31P-1H-HeadCoil', 'datfiles', 'S139_01292026_meas_MID01306_FID07054_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S139_01292026_31P-1H-HeadCoil', 'datfiles', 'S139_01292026_meas_MID01305_FID07053_svssel_latest_M0_Water_4avg_tra.dat');
ws_dat2  = fullfile(dataRoot, 'S139_02032026_31P-1H-HeadCoil', 'datfiles', 'S139_02032026_meas_MID00158_FID07348_svssel_latest_M0_NAD_9_1ppm_512avg_tra.dat');
nws_dat2 = fullfile(dataRoot, 'S139_02032026_31P-1H-HeadCoil', 'datfiles', 'S139_02032026_meas_MID00157_FID07347_svssel_latest_M0_Water_4avg_tra.dat');
groups(nGrp).ws  = {ws_dat, ws_dat2};
groups(nGrp).nws = {nws_dat, nws_dat2};
end

% Group 13: S167, brain, head coil, 02/02/2026
if ismember(13, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S167 brain head 02/02/2026';
groups(nGrp).bodypart = 'brain';
ws_dat  = fullfile(dataRoot, 'S167_02022026_31P-1H-HeadCoil', 'datfiles', 'S167_02022026_meas_MID00137_FID07200_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S167_02022026_31P-1H-HeadCoil', 'datfiles', 'S167_02022026_meas_MID00135_FID07198_svssel_latest_M0_Water_4avg_tra.dat');
groups(nGrp).ws  = {ws_dat};
groups(nGrp).nws = {nws_dat};
end

% --- CALF SCANS ---

% Group 14: S077, calf, 07/15/2025
if ismember(14, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S077 calf 07/15/2025';
groups(nGrp).bodypart = 'calf';
ws_dat  = fullfile(dataRoot, 'S077_07152025_31P-1H-HeadCoil_calf', 'datfiles', 'S077_07152025_meas_MID00124_FID15699_svssel_latest_M0_NAD_9_7ppm_512avg_tra_calf.dat');
nws_dat = fullfile(dataRoot, 'S077_07152025_31P-1H-HeadCoil_calf', 'datfiles', 'S077_07152025_meas_MID00123_FID15698_svssel_latest_M0_Water_1avg_tra_calf.dat');
groups(nGrp).ws  = {ws_dat};
groups(nGrp).nws = {nws_dat};
end

% Group 15: S095, calf, knee coil, 09/16/2025
if ismember(15, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S095 calf knee 09/16/2025';
groups(nGrp).bodypart = 'calf';
ws_dat  = fullfile(dataRoot, 'S095_09162025_28Ch-KneeCoil', 'datfiles', 'S095_09162025_meas_MID00442_FID23671_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S095_09162025_28Ch-KneeCoil', 'datfiles', 'S095_09162025_meas_MID00440_FID23669_svssel_latest_M0_Water_1avg_tra.dat');
groups(nGrp).ws  = {ws_dat};
groups(nGrp).nws = {nws_dat};
end

% Group 16: S138, calf, 07/24/2025
if ismember(16, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S138 calf 07/24/2025';
groups(nGrp).bodypart = 'calf';
ws_dat  = fullfile(dataRoot, 'S138_07242025_31P-1H-HeadCoil_calf', 'datfiles', 'S138_07242025_meas_MID00829_FID16886_svssel_latest_M0_NAD_9_7ppm_512avg_tra_calf.dat');
nws_dat = fullfile(dataRoot, 'S138_07242025_31P-1H-HeadCoil_calf', 'datfiles', 'S138_07242025_meas_MID00828_FID16885_svssel_latest_M0_Water_8avg_tra_calf.dat');
groups(nGrp).ws  = {ws_dat};
groups(nGrp).nws = {nws_dat};
end

% Group 17: S139, calf, repeat scans (03/12 - 03/18/2026)
if ismember(17, groupsToProcess)
nGrp = nGrp + 1;
groups(nGrp).label = 'S139 calf 03/2026';
groups(nGrp).bodypart = 'calf';
ws_dat  = fullfile(dataRoot, 'S139_03122026_31P-1H-HeadCoil_calf', 'datfiles', 'S139_03122026_calf_meas_MID00686_FID11071_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat = fullfile(dataRoot, 'S139_03122026_31P-1H-HeadCoil_calf', 'datfiles', 'S139_03122026_calf_meas_MID00685_FID11070_svssel_latest_M0_Water_4avg_tra.dat');
ws_dat2  = fullfile(dataRoot, 'S139_03182026_31P-1H-HeadCoil_calf', 'datfiles', 'S139_03182026_calf_meas_MID00239_FID11512_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat');
nws_dat2 = fullfile(dataRoot, 'S139_03182026_31P-1H-HeadCoil_calf', 'datfiles', 'S139_03182026_calf_meas_MID00238_FID11511_svssel_latest_M0_Water_8avg_tra.dat');
groups(nGrp).ws  = {ws_dat, ws_dat2};
groups(nGrp).nws = {nws_dat, nws_dat2};
end

% =====================================================================
% Process all groups: outer loop over groups, inner loop over repeats
% =====================================================================
for iGrp = 1:length(groups)
    nScans = length(groups(iGrp).ws);
    for iScan = 1:nScans
        fprintf('=== Group %d/%d (%s), scan %d/%d ===\n', ...
            iGrp, length(groups), groups(iGrp).label, iScan, nScans);

        % --- Varpro fitting ---
        pars_vp = pars;
        pars_vp.dofit = 'varpro';
        [output_vp, ~, ~] = dat2mat_svsNAD(pars_vp, groups(iGrp).ws{iScan}, groups(iGrp).nws{iScan});

        results(iGrp).label    = groups(iGrp).label;
        results(iGrp).bodypart = groups(iGrp).bodypart;
        results(iGrp).varpro.NADH2_amp(iScan) = output_vp.met.fit.areas(1);
        results(iGrp).varpro.NADH6_amp(iScan) = output_vp.met.fit.areas(2);
        results(iGrp).varpro.NADH4_amp(iScan) = output_vp.met.fit.areas(3);
        results(iGrp).varpro.wat_amp(iScan)   = output_vp.wat.fit.areaTotal;

    end

    % --- Compute NAD/water ratios and repeatability stats ---
    if nScans > 1
        for method = {'varpro'}
            m = method{1};
            results(iGrp).(m).NADH2_ratio = results(iGrp).(m).NADH2_amp ./ results(iGrp).(m).wat_amp;
            results(iGrp).(m).NADH6_ratio = results(iGrp).(m).NADH6_amp ./ results(iGrp).(m).wat_amp;
            results(iGrp).(m).NADH4_ratio = results(iGrp).(m).NADH4_amp ./ results(iGrp).(m).wat_amp;

            for pk = {'NADH2','NADH6','NADH4'}
                r = results(iGrp).(m).([pk{1} '_ratio']);
                results(iGrp).(m).([pk{1} '_mean']) = mean(r);
                results(iGrp).(m).([pk{1} '_std'])  = std(r);
                results(iGrp).(m).([pk{1} '_cv'])   = std(r) / mean(r) * 100;
            end
        end
    end
end

% =====================================================================
% Write repeatability results to TSV (groups with >1 scan only)
% =====================================================================
repGrps = find(arrayfun(@(g) length(g.ws), groups) > 1);
if ~isempty(repGrps)
    maxScans = max(arrayfun(@(g) length(g.ws), groups(repGrps)));
    peakNames = {'NADH2','NADH6','NADH4'};

    fid = fopen('varpro_repeatability.tsv', 'w');
    sep = '\t';

    % --- Build header ---
    hdr = 'Group\tBodyPart\tMethod';
    for s = 1:maxScans
        hdr = [hdr sprintf('\tWat_amp_s%d', s)];
    end
    for iPk = 1:length(peakNames)
        for s = 1:maxScans
            hdr = [hdr sprintf('\t%s_amp_s%d', peakNames{iPk}, s)];
        end
    end
    for iPk = 1:length(peakNames)
        for s = 1:maxScans
            hdr = [hdr sprintf('\t%s_ratio_s%d', peakNames{iPk}, s)];
        end
    end
    for iPk = 1:length(peakNames)
        hdr = [hdr sprintf('\t%s_mean\t%s_std\t%s_cv', peakNames{iPk}, peakNames{iPk}, peakNames{iPk})];
    end
    fprintf(fid, '%s\n', hdr);

    % --- Write rows ---
    for ii = 1:length(repGrps)
        iGrp = repGrps(ii);
        nScans = length(groups(iGrp).ws);

        for method = {'varpro'}
            m = method{1};
            row = sprintf('%s\t%s\t%s', results(iGrp).label, results(iGrp).bodypart, m);

            % water amplitudes
            for s = 1:maxScans
                if s <= nScans
                    row = [row sprintf('\t%.6g', results(iGrp).(m).wat_amp(s))];
                else
                    row = [row sprintf('\t')];
                end
            end

            % NAD+ peak amplitudes
            for iPk = 1:length(peakNames)
                ampField = [peakNames{iPk} '_amp'];
                for s = 1:maxScans
                    if s <= nScans
                        row = [row sprintf('\t%.6g', results(iGrp).(m).(ampField)(s))];
                    else
                        row = [row sprintf('\t')];
                    end
                end
            end

            % NAD+ / water ratios
            for iPk = 1:length(peakNames)
                ratioField = [peakNames{iPk} '_ratio'];
                for s = 1:maxScans
                    if s <= nScans
                        row = [row sprintf('\t%.6g', results(iGrp).(m).(ratioField)(s))];
                    else
                        row = [row sprintf('\t')];
                    end
                end
            end

            % mean, std, cv
            for iPk = 1:length(peakNames)
                row = [row sprintf('\t%.6g\t%.6g\t%.2f', ...
                    results(iGrp).(m).([peakNames{iPk} '_mean']), ...
                    results(iGrp).(m).([peakNames{iPk} '_std']), ...
                    results(iGrp).(m).([peakNames{iPk} '_cv']))];
                disp([peakNames{iPk} ' CV : ' num2str(results(iGrp).(m).([peakNames{iPk} '_cv']))]);
            end

            fprintf(fid, '%s\n', row);
        end
    end

    fclose(fid);
    fprintf('Repeatability results written to varpro_repeatability.tsv\n');
end