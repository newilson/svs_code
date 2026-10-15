%% characterize_data_differences.m
% Data-only characterization of sinc vs U01 MRS scans.
% NO fit-based metrics: all numbers come from raw spectra/FIDs.
%
% Per scan:
%   ---- Water reference ----
%   wat_A_int      : trapz [4 5.5] of real(spec) minus linear BL
%   wat_FWHM_Hz    : magnitude half-max
%   wat_SNR        : peak mag / std(real(spec), |ppm|∉[4,5.5])
%   wat_T2star_ms  : fit to log|FID| in first 80 ms
%   wat_noiseStd   : std(real(spec) outside water region)
%   wat_fid1       : |fid(1)| (total water magnetization)
%   wat_refVolt    : flTransRefAmpl
%   wat_vecSize    : lVectorSize
%
%   ---- Metabolite scan ----
%   met_fid1       : |fid(1)|
%   met_A_H2       : trapz [9.25 9.45] real(spec) - linBL
%   met_A_H6       : trapz [9.05 9.25] real(spec) - linBL
%   met_A_H4       : trapz [8.80 9.05] real(spec) - linBL
%   met_A_nuisance : trapz [7.8 8.5] real(spec) - linBL
%   met_A_resWater : max real(spec) in [4 5.5] (residual water peak height)
%   met_FWHM_H2_Hz : FWHM of H2 region (mag half-max within [9.25 9.45])
%   met_SNR_H2     : H2 peak / std(real) in [1 3] ppm
%   met_noiseStd   : std(real(spec)) in [1 3] ppm (met-side noise)
%   met_T2star_ms  : fit to log|FID| first 80 ms (envelope)
%   met_refVolt, met_vecSize, met_nAvg, met_TR_ms, met_TE_ms
%
%   ---- Core quantification ratios (data, not model) ----
%   ratio_H2_wat   = met_A_H2 / wat_A_int
%   ratio_H6_wat   = met_A_H6 / wat_A_int
%   ratio_H4_wat   = met_A_H4 / wat_A_int
%   ratio_nuis_wat = met_A_nuisance / wat_A_int

clear; clc
addpath('Z:\code\svs_code');
addpath('Z:\code\svs_code\utils');

outDir  = 'Z:\Data\Spectroscopy\Pilot-Sinc-DFMRS\PROCESSED';
cacheMat = fullfile(outDir,'datadiff_cache.mat');

% use the recommended processing
pars.lb_met = 3;
pars.lb_wat = 3;
pars.eccopt = 0;
pars.PC = 0;
pars.base = nan;
pars.plt = false; pars.pltSpec = false;
pars.block_align.size = 64; pars.block_align.method = 'fd';
pars.block_align.Nfit = 256; pars.block_align.ppmRange = [8.5 10];
pars.ccopt.minsig_frac = 0.05;
pars.WatSupPre.opts.type = 'none';
pars.WatSupPost.opts.type = 'none';
pars.WatSupPost.opts.filt.type = 'average'; pars.WatSupPost.opts.filt.N = 11;
pars.WatSupPost.opts.hsvd.bounds = [-1 1]*300/8000; pars.WatSupPost.opts.hsvd.nsin = 25;
pars.dofit = 'varpro';
pars.fit.mode = 5;
pars.fit.ppm_range = [8.7 10.0];
pars.fit.ph_range = [-60 60];
pars.peaks.name = {'NADH2','NADH6','NADH4'};
pars.peaks.range = [9.25 9.45; 9.05 9.25; 8.8 9.05];
pars.fit.baselineOpt.enable = true;
pars.fit.baselineOpt.knotSpacing = 2;
pars.fit.baselineOpt.lambda = 1;
pars.fit.baselineOpt.lambdaAmpl = 0;
pars.fit.baselineOpt.TolFun = 1e-6;
pars.fit.lineShapeOpt.enable = true;
pars.fit.lineShapeOpt.nSide = 4;
pars.fit.lineShapeOpt.asymmetric = true;
pars.fit.lineShapeOpt.maxSide = 0.50;
pars.fit.hsvdClean.enable = false;
pars.watfit.method = 'integrate';
pars.watfit.ppm_range = [4 5.5];

f0_MHz = 298.06;

% groups (same as sweep)
sincRoot = 'Z:\Data\Spectroscopy\Pilot-Sinc-DFMRS';
u01Root  = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_RR';
nGrp = 0;
nGrp=nGrp+1; groups(nGrp).label='S004 sinc 9.7ppm';
groups(nGrp).ws={fullfile(sincRoot,'S004_01092025','datfiles','S004_01092025_meas_MID01677_FID03067_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
                 fullfile(sincRoot,'S004_01142025','datfiles','S004_01142025_meas_MID02134_FID03526_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')};
groups(nGrp).nws={fullfile(sincRoot,'S004_01092025','datfiles','S004_01092025_meas_MID01675_FID03065_svssel_WATERREF_2x4_7ppm.dat')
                  fullfile(sincRoot,'S004_01142025','datfiles','S004_01142025_meas_MID02132_FID03524_svssel_WATERREF_2x4_7ppm.dat')};
nGrp=nGrp+1; groups(nGrp).label='S092 sinc 9.7ppm';
groups(nGrp).ws={fullfile(sincRoot,'S092_01072025','datfiles','S092_01072025_meas_MID01374_FID02760_svssel_sinc2x2_2X9_7ppm_TR1000_256avg_TE13ms.dat')
                 fullfile(sincRoot,'S092_01162025','datfiles','S092_01162025_meas_MID02350_FID03744_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')};
groups(nGrp).nws={fullfile(sincRoot,'S092_01072025','datfiles','S092_01072025_meas_MID01373_FID02759_svssel_WATERREF_2x4_7ppm_TE13ms.dat')
                  fullfile(sincRoot,'S092_01162025','datfiles','S092_01162025_meas_MID02348_FID03742_svssel_WATERREF_2x4_7ppm.dat')};
nGrp=nGrp+1; groups(nGrp).label='S092 U01 9.7ppm';
groups(nGrp).ws={fullfile(u01Root,'S092_03032026_31P-1H-HeadCoil','datfiles','S092_03032026_meas_MID01621_FID10109_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
                 fullfile(u01Root,'S092_03052026_31P-1H-HeadCoil','datfiles','S092_03052026_meas_MID01808_FID10296_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')};
groups(nGrp).nws={fullfile(u01Root,'S092_03032026_31P-1H-HeadCoil','datfiles','S092_03032026_meas_MID01620_FID10108_svssel_latest_M0_Water_4avg_tra.dat')
                  fullfile(u01Root,'S092_03052026_31P-1H-HeadCoil','datfiles','S092_03052026_meas_MID01806_FID10294_svssel_latest_M0_Water_4avg_tra.dat')};
nGrp=nGrp+1; groups(nGrp).label='S139 sinc 9.7ppm';
groups(nGrp).ws={fullfile(sincRoot,'S139_01092025','datfiles','S139_01092025_meas_MID01648_FID03038_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
                 fullfile(sincRoot,'S139_01162025','datfiles','S139_01162025_meas_MID02381_FID03775_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')};
groups(nGrp).nws={fullfile(sincRoot,'S139_01092025','datfiles','S139_01092025_meas_MID01646_FID03036_svssel_WATERREF_2x4_7ppm.dat')
                  fullfile(sincRoot,'S139_01162025','datfiles','S139_01162025_meas_MID02379_FID03773_svssel_WATERREF_2x4_7ppm.dat')};
nGrp=nGrp+1; groups(nGrp).label='S145 sinc 9.7ppm';
groups(nGrp).ws={fullfile(sincRoot,'S145_01082025','datfiles','S145_01082025_meas_MID01618_FID03008_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
                 fullfile(sincRoot,'S145_02042025','datfiles','S145_02042025_meas_MID01442_FID05524_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')};
groups(nGrp).nws={fullfile(sincRoot,'S145_01082025','datfiles','S145_01082025_meas_MID01616_FID03006_svssel_WATERREF_2x4_7ppm.dat')
                  fullfile(sincRoot,'S145_02042025','datfiles','S145_02042025_meas_MID01440_FID05522_svssel_WATERREF_2x4_7ppm.dat')};
nGrp=nGrp+1; groups(nGrp).label='S076 U01 9.7ppm';
groups(nGrp).ws={fullfile(u01Root,'S076_02052026_50mm_31P-1H-HeadCoil','datfiles','S76_02052026_50mm_meas_MID00223_FID07639_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
                 fullfile(u01Root,'S076_02122026_31P-1H-HeadCoil','datfiles','S076_02122026_meas_MID00932_FID08354_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
                 fullfile(u01Root,'S076_03102026_31P-1H-HeadCoil','datfiles','S076_03102026_meas_MID00459_FID10844_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')};
groups(nGrp).nws={fullfile(u01Root,'S076_02052026_50mm_31P-1H-HeadCoil','datfiles','S76_02052026_50mm_meas_MID00220_FID07636_svssel_latest_M0_Water_4avg_tra.dat')
                  fullfile(u01Root,'S076_02122026_31P-1H-HeadCoil','datfiles','S076_02122026_meas_MID00930_FID08352_svssel_latest_M0_Water_4avg_tra.dat')
                  fullfile(u01Root,'S076_03102026_31P-1H-HeadCoil','datfiles','S076_03102026_meas_MID00458_FID10843_svssel_latest_M0_Water_4avg_tra.dat')};
nGrp=nGrp+1; groups(nGrp).label='S138 U01 9.7ppm';
groups(nGrp).ws={fullfile(u01Root,'S138_02172026_31P-1H-HeadCoil','datfiles','S138_02172026_meas_MID00379_FID08871_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
                 fullfile(u01Root,'S138_03192026_31P-1H-HeadCoil','datfiles','S138_03192026_meas_MID00036_FID11602_svssel_latest_M0_NAD_9_7ppm_1024avg_tra.dat')
                 fullfile(u01Root,'S138_04092026_31P-1H-HeadCoil','datfiles','S138_04092026_meas_MID02333_FID13984_svssel_latest_M0_NAD_9_7ppm_1024avg_tra.dat')
                 fullfile(u01Root,'S138_04102026_31P-1H-HeadCoil','datfiles','S138_04102026_meas_MID02435_FID14086_svssel_latest_M0_NAD_9_7ppm_1024avg_tra.dat')};
groups(nGrp).nws={fullfile(u01Root,'S138_02172026_31P-1H-HeadCoil','datfiles','S138_02172026_meas_MID00378_FID08870_svssel_latest_M0_Water_4avg_tra.dat')
                  fullfile(u01Root,'S138_03192026_31P-1H-HeadCoil','datfiles','S138_03192026_meas_MID00035_FID11601_svssel_latest_M0_Water_8avg_tra.dat')
                  fullfile(u01Root,'S138_04092026_31P-1H-HeadCoil','datfiles','S138_04092026_meas_MID02332_FID13983_svssel_latest_M0_Water_8avg_tra.dat')
                  fullfile(u01Root,'S138_04102026_31P-1H-HeadCoil','datfiles','S138_04102026_meas_MID02434_FID14085_svssel_latest_M0_Water_8avg_tra.dat')};

% ---------- build/load cache ----------
if exist(cacheMat,'file')
    fprintf('Loading cache %s\n', cacheMat);
    S = load(cacheMat); cache = S.cache;
else
    cache = struct([]); idx = 0;
    for iGrp = 1:length(groups)
        nS = length(groups(iGrp).ws);
        for iScan = 1:nS
            idx = idx + 1;
            fprintf('[%d] %s scan %d/%d\n', idx, groups(iGrp).label, iScan, nS);
            try
                [output,~,~] = dat2mat_svsNAD(pars, groups(iGrp).ws{iScan}, groups(iGrp).nws{iScan});
                cache(idx).iGrp = iGrp; cache(idx).iScan = iScan;
                cache(idx).label = groups(iGrp).label;
                cache(idx).wat_ppm = output.wat.ppm;
                cache(idx).wat_fid = output.wat.fid;
                cache(idx).wat_spec = output.wat.spec;
                cache(idx).wat_refvolt = output.wat.refvolt;
                cache(idx).wat_TR_ms = output.wat.TR_ms;
                cache(idx).wat_TE_ms = output.wat.TE_ms;
                cache(idx).met_ppm = output.met.ppm;
                cache(idx).met_fid = output.met.fid;
                cache(idx).met_spec = output.met.spec;
                cache(idx).met_refvolt = output.met.refvolt;
                cache(idx).met_TR_ms = output.met.TR_ms;
                cache(idx).met_TE_ms = output.met.TE_ms;
                cache(idx).metAreas = output.met.fit.areas(1:3)';
            catch ME
                fprintf('  FAILED: %s\n', ME.message);
            end
        end
    end
    save(cacheMat,'cache','groups','-v7.3');
end

nC = numel(cache);
proto = strings(nC,1);
for i=1:nC
    if contains(cache(i).label,'sinc'), proto(i)="sinc";
    elseif contains(cache(i).label,'U01'), proto(i)="U01"; end
end

% ---------- metrics ----------
M = struct();
for i = 1:nC
    c = cache(i);
    % WATER
    ppm_w = c.wat_ppm(:); sp_w = c.wat_spec(:); fid_w = c.wat_fid(:);
    Nw = length(fid_w);
    span_Hz = abs(ppm_w(end)-ppm_w(1))*f0_MHz;
    dt_w = 1/span_Hz;

    sel = ppm_w>4 & ppm_w<5.5;
    xp = ppm_w(sel); yp = double(real(sp_w(sel)));
    bl = linspace(yp(1),yp(end),length(yp))';
    M(i).wat_A_int = trapz(xp, yp - bl);

    mag = abs(sp_w); [pk, imax] = max(mag); half = pk/2;
    Lidx = find(mag(1:imax)<=half,1,'last');
    Ridx = find(mag(imax:end)<=half,1,'first');
    if ~isempty(Lidx)&&~isempty(Ridx)
        M(i).wat_FWHM_Hz = abs(ppm_w(imax+Ridx-1)-ppm_w(Lidx))*f0_MHz;
    else, M(i).wat_FWHM_Hz = NaN;
    end
    noiseSel = (ppm_w<3)|(ppm_w>7);
    M(i).wat_noiseStd = std(double(real(sp_w(noiseSel))));
    M(i).wat_SNR = pk / M(i).wat_noiseStd;

    t = (0:Nw-1)'*dt_w;
    nFit = min(Nw, round(0.08/dt_w));
    y = abs(fid_w(1:nFit))/abs(fid_w(1));
    p = polyfit(t(1:nFit), log(max(y,1e-6)), 1);
    M(i).wat_T2star_ms = -1000/p(1);
    M(i).wat_fid1 = abs(fid_w(1));
    M(i).wat_refVolt = c.wat_refvolt;
    M(i).wat_N = Nw;
    M(i).wat_TR_ms = c.wat_TR_ms;
    M(i).wat_TE_ms = c.wat_TE_ms;

    % METABOLITE
    ppm_m = c.met_ppm(:); sp_m = c.met_spec(:); fid_m = c.met_fid(:);
    Nm = length(fid_m);
    span_m_Hz = abs(ppm_m(end)-ppm_m(1))*f0_MHz;
    dt_m = 1/span_m_Hz;

    M(i).met_A_H2 = integ_bl(ppm_m, sp_m, 9.25, 9.45);
    M(i).met_A_H6 = integ_bl(ppm_m, sp_m, 9.05, 9.25);
    M(i).met_A_H4 = integ_bl(ppm_m, sp_m, 8.80, 9.05);
    M(i).met_A_nuisance = integ_bl(ppm_m, sp_m, 7.8, 8.5);
    M(i).met_A_0to2      = integ_bl(ppm_m, sp_m, 0.5, 2.0);   % aliphatic region
    M(i).met_A_resWater  = max(double(real(sp_m(ppm_m>4 & ppm_m<5.5))));

    % H2 FWHM (mag, within [9.2 9.5])
    selH2 = ppm_m>9.20 & ppm_m<9.50;
    mag_m = abs(sp_m); xH2 = ppm_m(selH2); yH2 = mag_m(selH2);
    [pkH2, iH2] = max(yH2); half = pkH2/2;
    Lidx = find(yH2(1:iH2)<=half,1,'last'); Ridx = find(yH2(iH2:end)<=half,1,'first');
    if ~isempty(Lidx)&&~isempty(Ridx)
        M(i).met_FWHM_H2_Hz = abs(xH2(iH2+Ridx-1) - xH2(Lidx))*f0_MHz;
    else, M(i).met_FWHM_H2_Hz = NaN;
    end

    noiseSel_m = ppm_m>1 & ppm_m<3;
    M(i).met_noiseStd = std(double(real(sp_m(noiseSel_m))));
    M(i).met_SNR_H2 = pkH2 / M(i).met_noiseStd;

    t_m = (0:Nm-1)'*dt_m;
    nFit_m = min(Nm, round(0.08/dt_m));
    ym = abs(fid_m(1:nFit_m))/abs(fid_m(1));
    pm = polyfit(t_m(1:nFit_m), log(max(ym,1e-6)), 1);
    M(i).met_T2star_ms = -1000/pm(1);
    M(i).met_fid1 = abs(fid_m(1));
    M(i).met_refVolt = c.met_refvolt;
    M(i).met_N = Nm;
    M(i).met_TR_ms = c.met_TR_ms;
    M(i).met_TE_ms = c.met_TE_ms;

    % Key ratios
    M(i).ratio_H2_wat   = M(i).met_A_H2 / M(i).wat_A_int;
    M(i).ratio_H6_wat   = M(i).met_A_H6 / M(i).wat_A_int;
    M(i).ratio_H4_wat   = M(i).met_A_H4 / M(i).wat_A_int;
    M(i).ratio_nuis_wat = M(i).met_A_nuisance / M(i).wat_A_int;
    M(i).ratio_fid_wat  = M(i).met_fid1 / M(i).wat_fid1;
    M(i).ratio_H2_fid   = M(i).met_A_H2 / M(i).met_fid1;

    M(i).proto = proto(i);
    M(i).label = c.label;
    M(i).iGrp  = c.iGrp;
end

% ---------- print table ----------
fprintf('\n#   label                     proto  vecM  vecW   TR_m   TR_w  avg_m  avg_w   refV_m refV_w\n');
for i=1:nC
    fprintf('%2d  %-22s  %-5s  %4d  %4d  %5.0f  %5.0f   %5d  %5d   %5.0f %5.0f\n', ...
        i, M(i).label, M(i).proto, M(i).met_N, M(i).wat_N, ...
        M(i).met_TR_ms, M(i).wat_TR_ms, NaN, NaN, ...
        M(i).met_refVolt, M(i).wat_refVolt);
end

fprintf('\n#   label                     WAT:    A_int        fid1    FWHM_Hz    SNR   noise   T2*ms\n');
for i=1:nC
    fprintf('%2d  %-22s   %10.2f  %10.2f  %7.1f  %6.0f  %7.3f  %5.1f\n', ...
        i, M(i).label, M(i).wat_A_int, M(i).wat_fid1, M(i).wat_FWHM_Hz, M(i).wat_SNR, M(i).wat_noiseStd, M(i).wat_T2star_ms);
end

fprintf('\n#   label                     MET:    fid1       H2       H6       H4     nuisc   resW    FWHM_H2  SNR_H2  noise   T2*ms\n');
for i=1:nC
    fprintf('%2d  %-22s   %10.3f  %7.3f  %7.3f  %7.3f  %7.3f  %5.2f  %6.1f  %5.0f  %6.3f  %5.1f\n', ...
        i, M(i).label, M(i).met_fid1, M(i).met_A_H2, M(i).met_A_H6, M(i).met_A_H4, ...
        M(i).met_A_nuisance, M(i).met_A_resWater, M(i).met_FWHM_H2_Hz, M(i).met_SNR_H2, M(i).met_noiseStd, M(i).met_T2star_ms);
end

fprintf('\n#   label                     RATIOS:  H2/wat     H6/wat    H4/wat    nuis/wat   fid_m/fid_w    H2/fid_m\n');
for i=1:nC
    fprintf('%2d  %-22s  %10.3e %10.3e %10.3e %10.3e  %10.3e  %10.3e\n', ...
        i, M(i).label, M(i).ratio_H2_wat, M(i).ratio_H6_wat, M(i).ratio_H4_wat, ...
        M(i).ratio_nuis_wat, M(i).ratio_fid_wat, M(i).ratio_H2_fid);
end

% ---------- medians by protocol ----------
fldList = {'wat_A_int','wat_fid1','wat_FWHM_Hz','wat_SNR','wat_noiseStd','wat_T2star_ms',...
           'met_fid1','met_A_H2','met_A_H6','met_A_H4','met_A_nuisance','met_A_resWater',...
           'met_FWHM_H2_Hz','met_SNR_H2','met_noiseStd','met_T2star_ms',...
           'ratio_H2_wat','ratio_H6_wat','ratio_H4_wat','ratio_nuis_wat','ratio_fid_wat','ratio_H2_fid',...
           'met_refVolt','wat_refVolt','met_N','wat_N','met_TR_ms'};

fprintf('\n==== MEDIAN by protocol  (U01/sinc ratio in parens) ====\n');
fprintf('%-20s  %12s  %12s  %10s\n', 'metric','sinc','U01','U01/sinc');
for f = fldList
    vs = []; vu = [];
    for i=1:nC
        if proto(i)=="sinc", vs(end+1) = M(i).(f{1}); end
        if proto(i)=="U01",  vu(end+1) = M(i).(f{1}); end
    end
    ms = median(vs,'omitnan'); mu = median(vu,'omitnan');
    fprintf('%-20s  %12.4g  %12.4g  %10.3f\n', f{1}, ms, mu, mu/ms);
end

% ---------- S092 direct sinc vs U01 (same subject) ----------
fprintf('\n==== S092 same-subject: sinc vs U01 ====\n');
iSs = find(strcmp({M.label},'S092 sinc 9.7ppm'));
iSu = find(strcmp({M.label},'S092 U01 9.7ppm'));
fprintf('%-20s  %12s  %12s  %10s\n','metric','sinc','U01','U01/sinc');
for f = fldList
    vs = median([M(iSs).(f{1})],'omitnan');
    vu = median([M(iSu).(f{1})],'omitnan');
    fprintf('%-20s  %12.4g  %12.4g  %10.3f\n', f{1}, vs, vu, vu/vs);
end

save(fullfile(outDir,'datadiff_metrics.mat'),'M','proto','-v7.3');

% ---------- figures ----------
% (A) Overlaid full metabolite spectra, normalized per-protocol
fig = figure('Position',[50 50 1400 800],'Color','w');

subplot(2,2,1); hold on; grid on
colS = [0.1 0.35 0.8]; colU = [0.85 0.25 0.1];
for i=1:nC
    sp = real(cache(i).met_spec); ppm = cache(i).met_ppm(:);
    % normalize by MEDIAN peak height in met region for visual compare
    mxH2 = max(sp(ppm>9.20 & ppm<9.50));
    col = colS; if proto(i)=="U01", col=colU; end
    plot(ppm, sp/mxH2, 'Color',[col 0.5], 'LineWidth',1.0);
end
set(gca,'xdir','reverse'); xlim([0.5 10]); ylim([-0.5 2])
xlabel('ppm'); ylabel('real(spec) / max(H2 region)')
title('All metabolite spectra, normalized to H2 peak')

subplot(2,2,2); hold on; grid on
for i=1:nC
    sp = real(cache(i).met_spec); ppm = cache(i).met_ppm(:);
    mxH2 = max(sp(ppm>9.20 & ppm<9.50));
    col = colS; if proto(i)=="U01", col=colU; end
    plot(ppm, sp/mxH2, 'Color',[col 0.5], 'LineWidth',1.0);
end
set(gca,'xdir','reverse'); xlim([7.5 10]); ylim([-0.5 1.5])
xlabel('ppm'); ylabel('normalized'); title('Zoom: NAD + nuisance [7.5 10]')

subplot(2,2,3); hold on; grid on
for i=1:nC
    sp = real(cache(i).met_spec); ppm = cache(i).met_ppm(:);
    col = colS; if proto(i)=="U01", col=colU; end
    plot(ppm, sp, 'Color',[col 0.5], 'LineWidth',1.0);
end
set(gca,'xdir','reverse'); xlim([3 6]);
xlabel('ppm'); ylabel('real(spec) (absolute)')
title('Residual water region (raw scale)')

subplot(2,2,4); hold on; grid on
ms = [M.met_A_H2]; mw = [M.wat_A_int]; mr = ms./mw;
ix_s = proto=="sinc"; ix_u = proto=="U01";
scatter(1:sum(ix_s), mr(ix_s),60,colS,'filled','DisplayName','sinc');
scatter(sum(ix_s)+(1:sum(ix_u)), mr(ix_u),60,colU,'filled','DisplayName','U01');
ylabel('H2 / water_{int}'); title('Core ratio: H2-area / water-area')
xlabel('scan #'); legend('Location','best'); grid on
set(gca,'YScale','log')

sgtitle('Data-only characterization: sinc vs U01','FontWeight','bold')
saveas(fig, fullfile(outDir,'datadiff_overlays.png'));
fprintf('\nFigure -> %s\n', fullfile(outDir,'datadiff_overlays.png'));

function A = integ_bl(x,y,lo,hi)
    sel_ = x>lo & x<hi;
    if nnz(sel_) < 2, A = NaN; return; end
    xs = x(sel_); ys = double(real(y(sel_)));
    bl_ = linspace(ys(1), ys(end), length(ys))';
    A = trapz(xs, ys-bl_);
end
