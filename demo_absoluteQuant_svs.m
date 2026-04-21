%% demo_absoluteQuant_svs.m
% Demo: absolute NAD+ concentration from all repeatability scans.
%
% Runs dat2mat_svsNAD.m with the optimal varpro/Voigt parameters used in
% processRepeatability.m, then calls absoluteQuant_svs.m to convert each
% fit into absolute mmol/L using a 3-compartment (WM/GM/CSF) tissue model.
%
% Fixed tissue composition (as requested): 50% WM, 40% GM, 10% CSF.
% Water T1/T2/concentration and metabolite T1/T2 use literature values
% (7T human brain); NAD+ relaxation times from
%   https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/nbm.5324
%
% Output mimics processRepeatability.m but reports absolute [M] (mM) and
% repeatability CVs on the concentration itself rather than area ratios.

clear; clc

addpath('Z:\code\svs_code');
addpath('Z:\code\svs_code\utils');

outDir = 'Z:\Data\Spectroscopy\Pilot-Sinc-DFMRS\PROCESSED';
outName = 'AbsoluteQuant_9_7.tsv';
if ~exist(outDir,'dir'), mkdir(outDir); end

% =====================================================================
% Fit parameters (optimal settings from processRepeatability.m)
% =====================================================================
pars.lb_met = 3;
pars.lb_wat = 3;
pars.eccopt = 0;
pars.PC = 0;
pars.base = nan;
pars.plt = false;
pars.pltSpec = false;

pars.block_align.size = 64;
pars.block_align.method = 'fd';
pars.block_align.Nfit = 256;
pars.block_align.ppmRange = [8.5 10];
pars.ccopt.minsig_frac = 0.05;

pars.WatSupPre.opts.type = 'none';
pars.WatSupPost.opts.type = 'none';
pars.WatSupPost.opts.filt.type = 'average';
pars.WatSupPost.opts.filt.N = 11;
pars.WatSupPost.opts.hsvd.bounds = [-1 1] * 300/8000;
pars.WatSupPost.opts.hsvd.nsin = 25;

pars.dofit = 'varpro';
pars.fit.mode = 5; % Voigt
pars.fit.ppm_range = [8.7 10.0];
pars.fit.ph_range = [-60 60];
pars.peaks.name  = {'NADH2','NADH6','NADH4'};
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

pars.watfit.method = 'fit';
pars.watfit.nComp = 1;
pars.watfit.ppm_range = [4 5.5];

% =====================================================================
% Tissue composition and literature relaxation/concentration values
% =====================================================================
% Water T1/T2 and water content [W] per tissue taken from Table 1 of
% Swago et al., NMR Biomed 2025 (doi:10.1002/nbm.5324), the same paper
% from which the NAD+ T1/T2 values below are drawn. [W] is in mol/L.
% Voxel fractions fixed here at 50% WM / 40% GM / 10% CSF (user spec).

voxel.tissues(1).name        = 'WM';
voxel.tissues(1).fraction    = 0.50;
voxel.tissues(1).T1_ms       = 1130;
voxel.tissues(1).T2_ms       = 40;
voxel.tissues(1).waterConc_M = 38.5;
voxel.tissues(1).flag_metSig = true;

voxel.tissues(2).name        = 'GM';
voxel.tissues(2).fraction    = 0.40;
voxel.tissues(2).T1_ms       = 1939;
voxel.tissues(2).T2_ms       = 40;
voxel.tissues(2).waterConc_M = 44;
voxel.tissues(2).flag_metSig = true;

voxel.tissues(3).name        = 'CSF';
voxel.tissues(3).fraction    = 0.10;
voxel.tissues(3).T1_ms       = 4400;
voxel.tissues(3).T2_ms       = 300;
voxel.tissues(3).waterConc_M = 55;
voxel.tissues(3).flag_metSig = false;   % CSF excluded from metabolite volume

% --- NAD+ T1/T2 per peak from Swago et al., NBM 2025 (Table 2 means) ---
% Single aromatic proton per peak (H2, H6, H4 of nicotinamide).
met_T1_ms   = [205.6, 211.6, 237.3];   % NADH2, NADH6, NADH4
met_T2_ms   = [ 33.6,  29.1,  42.3];
met_nProt   = [    1,     1,     1];
peakNames   = {'NADH2','NADH6','NADH4'};
nPeaks = numel(peakNames);

% =====================================================================
% Groups (copied from processRepeatability.m)
% =====================================================================
sincRoot = 'Z:\Data\Spectroscopy\Pilot-Sinc-DFMRS';
u01Root  = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_RR';
nGrp = 0;

nGrp=nGrp+1;
groups(nGrp).label='S004 sinc 9.7ppm'; groups(nGrp).subject='S004';
groups(nGrp).centerFreq='9.7 ppm'; groups(nGrp).dataset='Pilot-Sinc';
groups(nGrp).ws = {
    fullfile(sincRoot,'S004_01092025','datfiles','S004_01092025_meas_MID01677_FID03067_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
    fullfile(sincRoot,'S004_01142025','datfiles','S004_01142025_meas_MID02134_FID03526_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')};
groups(nGrp).nws = {
    fullfile(sincRoot,'S004_01092025','datfiles','S004_01092025_meas_MID01675_FID03065_svssel_WATERREF_2x4_7ppm.dat')
    fullfile(sincRoot,'S004_01142025','datfiles','S004_01142025_meas_MID02132_FID03524_svssel_WATERREF_2x4_7ppm.dat')};

nGrp=nGrp+1;
groups(nGrp).label='S092 sinc 9.7ppm'; groups(nGrp).subject='S092';
groups(nGrp).centerFreq='9.7 ppm'; groups(nGrp).dataset='Pilot-Sinc';
groups(nGrp).ws = {
    fullfile(sincRoot,'S092_01072025','datfiles','S092_01072025_meas_MID01374_FID02760_svssel_sinc2x2_2X9_7ppm_TR1000_256avg_TE13ms.dat')
    fullfile(sincRoot,'S092_01162025','datfiles','S092_01162025_meas_MID02350_FID03744_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')};
groups(nGrp).nws = {
    fullfile(sincRoot,'S092_01072025','datfiles','S092_01072025_meas_MID01373_FID02759_svssel_WATERREF_2x4_7ppm_TE13ms.dat')
    fullfile(sincRoot,'S092_01162025','datfiles','S092_01162025_meas_MID02348_FID03742_svssel_WATERREF_2x4_7ppm.dat')};

nGrp=nGrp+1;
groups(nGrp).label='S092 U01 9.7ppm'; groups(nGrp).subject='S092';
groups(nGrp).centerFreq='9.7 ppm'; groups(nGrp).dataset='U01_RR';
groups(nGrp).ws = {
    fullfile(u01Root,'S092_03032026_31P-1H-HeadCoil','datfiles','S092_03032026_meas_MID01621_FID10109_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
    fullfile(u01Root,'S092_03052026_31P-1H-HeadCoil','datfiles','S092_03052026_meas_MID01808_FID10296_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')};
groups(nGrp).nws = {
    fullfile(u01Root,'S092_03032026_31P-1H-HeadCoil','datfiles','S092_03032026_meas_MID01620_FID10108_svssel_latest_M0_Water_4avg_tra.dat')
    fullfile(u01Root,'S092_03052026_31P-1H-HeadCoil','datfiles','S092_03052026_meas_MID01806_FID10294_svssel_latest_M0_Water_4avg_tra.dat')};

nGrp=nGrp+1;
groups(nGrp).label='S139 sinc 9.7ppm'; groups(nGrp).subject='S139';
groups(nGrp).centerFreq='9.7 ppm'; groups(nGrp).dataset='Pilot-Sinc';
groups(nGrp).ws = {
    fullfile(sincRoot,'S139_01092025','datfiles','S139_01092025_meas_MID01648_FID03038_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
    fullfile(sincRoot,'S139_01162025','datfiles','S139_01162025_meas_MID02381_FID03775_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')};
groups(nGrp).nws = {
    fullfile(sincRoot,'S139_01092025','datfiles','S139_01092025_meas_MID01646_FID03036_svssel_WATERREF_2x4_7ppm.dat')
    fullfile(sincRoot,'S139_01162025','datfiles','S139_01162025_meas_MID02379_FID03773_svssel_WATERREF_2x4_7ppm.dat')};

nGrp=nGrp+1;
groups(nGrp).label='S145 sinc 9.7ppm'; groups(nGrp).subject='S145';
groups(nGrp).centerFreq='9.7 ppm'; groups(nGrp).dataset='Pilot-Sinc';
groups(nGrp).ws = {
    fullfile(sincRoot,'S145_01082025','datfiles','S145_01082025_meas_MID01618_FID03008_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
    fullfile(sincRoot,'S145_02042025','datfiles','S145_02042025_meas_MID01442_FID05524_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')};
groups(nGrp).nws = {
    fullfile(sincRoot,'S145_01082025','datfiles','S145_01082025_meas_MID01616_FID03006_svssel_WATERREF_2x4_7ppm.dat')
    fullfile(sincRoot,'S145_02042025','datfiles','S145_02042025_meas_MID01440_FID05522_svssel_WATERREF_2x4_7ppm.dat')};

nGrp=nGrp+1;
groups(nGrp).label='S076 U01 9.7ppm'; groups(nGrp).subject='S076';
groups(nGrp).centerFreq='9.7 ppm'; groups(nGrp).dataset='U01_RR';
groups(nGrp).ws = {
    fullfile(u01Root,'S076_02052026_50mm_31P-1H-HeadCoil','datfiles','S76_02052026_50mm_meas_MID00223_FID07639_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
    fullfile(u01Root,'S076_02122026_31P-1H-HeadCoil','datfiles','S076_02122026_meas_MID00932_FID08354_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
    fullfile(u01Root,'S076_03102026_31P-1H-HeadCoil','datfiles','S076_03102026_meas_MID00459_FID10844_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')};
groups(nGrp).nws = {
    fullfile(u01Root,'S076_02052026_50mm_31P-1H-HeadCoil','datfiles','S76_02052026_50mm_meas_MID00220_FID07636_svssel_latest_M0_Water_4avg_tra.dat')
    fullfile(u01Root,'S076_02122026_31P-1H-HeadCoil','datfiles','S076_02122026_meas_MID00930_FID08352_svssel_latest_M0_Water_4avg_tra.dat')
    fullfile(u01Root,'S076_03102026_31P-1H-HeadCoil','datfiles','S076_03102026_meas_MID00458_FID10843_svssel_latest_M0_Water_4avg_tra.dat')};

nGrp=nGrp+1;
groups(nGrp).label='S138 U01 9.7ppm'; groups(nGrp).subject='S138';
groups(nGrp).centerFreq='9.7 ppm'; groups(nGrp).dataset='U01_RR';
groups(nGrp).ws = {
    fullfile(u01Root,'S138_02172026_31P-1H-HeadCoil','datfiles','S138_02172026_meas_MID00379_FID08871_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
    fullfile(u01Root,'S138_03192026_31P-1H-HeadCoil','datfiles','S138_03192026_meas_MID00036_FID11602_svssel_latest_M0_NAD_9_7ppm_1024avg_tra.dat')
    fullfile(u01Root,'S138_04092026_31P-1H-HeadCoil','datfiles','S138_04092026_meas_MID02333_FID13984_svssel_latest_M0_NAD_9_7ppm_1024avg_tra.dat')
    fullfile(u01Root,'S138_04102026_31P-1H-HeadCoil','datfiles','S138_04102026_meas_MID02435_FID14086_svssel_latest_M0_NAD_9_7ppm_1024avg_tra.dat')};
groups(nGrp).nws = {
    fullfile(u01Root,'S138_02172026_31P-1H-HeadCoil','datfiles','S138_02172026_meas_MID00378_FID08870_svssel_latest_M0_Water_4avg_tra.dat')
    fullfile(u01Root,'S138_03192026_31P-1H-HeadCoil','datfiles','S138_03192026_meas_MID00035_FID11601_svssel_latest_M0_Water_8avg_tra.dat')
    fullfile(u01Root,'S138_04092026_31P-1H-HeadCoil','datfiles','S138_04092026_meas_MID02332_FID13983_svssel_latest_M0_Water_8avg_tra.dat')
    fullfile(u01Root,'S138_04102026_31P-1H-HeadCoil','datfiles','S138_04102026_meas_MID02434_FID14085_svssel_latest_M0_Water_8avg_tra.dat')};

% =====================================================================
% Process all groups -> absolute concentration per scan
% =====================================================================
for iGrp = 1:length(groups)
    nScans = length(groups(iGrp).ws);
    conc_uM = zeros(nScans, nPeaks);   % [M] per peak per scan (uM)
    watAreas = zeros(nScans,1);
    metAreas = zeros(nScans, nPeaks);

    fprintf('\n=== Group %d/%d: %s (%d scans) ===\n', ...
        iGrp, length(groups), groups(iGrp).label, nScans);

    for iScan = 1:nScans
        fprintf('  Scan %d/%d\n', iScan, nScans);
        [output,~,~] = dat2mat_svsNAD(pars, ...
                        groups(iGrp).ws{iScan}, groups(iGrp).nws{iScan});

        % --- Assemble inputs to absoluteQuant_svs ---
        data = struct();
        data.water.area = output.wat.fit.areaTotal;
        for k = 1:nPeaks
            data.metabolites(k).name     = peakNames{k};
            data.metabolites(k).area     = output.met.fit.areas(k);
            data.metabolites(k).nProtons = met_nProt(k);
            data.metabolites(k).T1_ms    = met_T1_ms(k);
            data.metabolites(k).T2_ms    = met_T2_ms(k);
        end

        seq = struct();
        seq.TR_ms     = output.met.TR_ms;
        seq.TE_ms     = output.met.TE_ms;
        seq.TR_ws_ms  = output.wat.TR_ms;
        seq.TE_ws_ms  = output.wat.TE_ms;

        absQ = absoluteQuant_svs(data, seq, voxel);

        conc_uM(iScan,:)  = [absQ.metabolites.conc_mM] * 1e3;  % mM -> uM
        watAreas(iScan)   = output.wat.fit.areaTotal;
        metAreas(iScan,:) = output.met.fit.areas(1:nPeaks)';
    end

    % --- Repeatability stats on concentration ---
    if nScans > 1
        cvs = std(conc_uM,0,1) ./ mean(conc_uM,1) * 100;
    else
        cvs = nan(1,nPeaks);
    end

    groups(iGrp).conc_uM   = conc_uM;
    groups(iGrp).watAreas  = watAreas;
    groups(iGrp).metAreas  = metAreas;
    groups(iGrp).cvs       = cvs;
    groups(iGrp).nScans    = nScans;
    groups(iGrp).age       = getAgeForGroup(groups(iGrp).ws);
end

% =====================================================================
% Per-peak concentration table (all groups)
% =====================================================================
fprintf('\n\n');
maxScans = max([groups.nScans]);
for iPk = 1:nPeaks
    fprintf('=== %s absolute concentration (uM) ===\n', peakNames{iPk});
    fprintf('%-25s  %6s', 'Group', 'nScans');
    for iS = 1:maxScans, fprintf('  %10s', sprintf('Scan%d',iS)); end
    fprintf('  %8s  %8s  %8s\n','Mean','Std','CV(%)');
    fprintf('%s\n', repmat('-',1, 35 + maxScans*12 + 30));
    for iGrp = 1:length(groups)
        g = groups(iGrp);
        fprintf('%-25s  %6d', g.label, g.nScans);
        for iS = 1:maxScans
            if iS <= g.nScans
                fprintf('  %10.1f', g.conc_uM(iS, iPk));
            else
                fprintf('  %10s','');
            end
        end
        if g.nScans > 1
            fprintf('  %8.1f  %8.1f  %7.1f%%\n', ...
                mean(g.conc_uM(:,iPk)), std(g.conc_uM(:,iPk)), g.cvs(iPk));
        else
            fprintf('  %8.1f  %8s  %8s\n', g.conc_uM(1,iPk),'','');
        end
    end
    fprintf('\n');
end

% =====================================================================
% Mean(H2+H6) and mean(H2+H6+H4) summary
% =====================================================================
fprintf('=== Mean concentrations across peaks (uM) ===\n');
fprintf('%-25s  %6s  %10s  %10s  %10s  %10s\n', ...
    'Group','nScans','H2+H6_mean','H2+H6_CV%','H2+H6+H4_m','H2+H6+H4_CV%');
fprintf('%s\n', repmat('-',1,85));
for iGrp = 1:length(groups)
    g = groups(iGrp);
    m26  = mean(g.conc_uM(:,1:2), 2);
    m246 = mean(g.conc_uM(:,1:3), 2);
    if g.nScans > 1
        cv26  = std(m26) /mean(m26)  * 100;
        cv246 = std(m246)/mean(m246) * 100;
    else
        cv26 = nan; cv246 = nan;
    end
    fprintf('%-25s  %6d  %10.1f  %9.1f%%  %10.1f  %9.1f%%\n', ...
        g.label, g.nScans, mean(m26), cv26, mean(m246), cv246);
end

% =====================================================================
% Write TSV
% =====================================================================
fid = fopen(fullfile(outDir,outName),'w');
hdr = 'Group\tSubject\tDataset\tCenterFreq\tPeak\tnScans';
for iS = 1:maxScans
    hdr = [hdr sprintf('\tconc_uM_s%d',iS)];
end
hdr = [hdr '\tMean_uM\tStd_uM\tCV_%%'];
fprintf(fid,'%s\n',hdr);
for iGrp = 1:length(groups)
    g = groups(iGrp);
    for iPk = 1:nPeaks
        row = sprintf('%s\t%s\t%s\t%s\t%s\t%d', ...
            g.label,g.subject,g.dataset,g.centerFreq,peakNames{iPk},g.nScans);
        for iS = 1:maxScans
            if iS <= g.nScans
                row = [row sprintf('\t%.6g', g.conc_uM(iS,iPk))];
            else
                row = [row '\t'];
            end
        end
        if g.nScans > 1
            row = [row sprintf('\t%.6g\t%.6g\t%.2f', ...
                mean(g.conc_uM(:,iPk)), std(g.conc_uM(:,iPk)), g.cvs(iPk))];
        else
            row = [row sprintf('\t%.6g\t\t', g.conc_uM(1,iPk))];
        end
        fprintf(fid,'%s\n',row);
    end
end
fclose(fid);
fprintf('\nAbsolute concentrations written to %s\n', fullfile(outDir,outName));

save(fullfile(outDir,'groups_absQuant.mat'),'groups','voxel','-v7.3');
fprintf('Groups struct saved to %s\n', fullfile(outDir,'groups_absQuant.mat'));

% =====================================================================
% NAD H2+H6 concentration vs age (sinc vs U01 separate)
% =====================================================================
isSinc = strcmp({groups.dataset},'Pilot-Sinc');
isU01  = strcmp({groups.dataset},'U01_RR');

fig = figure('Position',[100 100 900 500],'Color','w'); hold on; grid on
colS = [0.1 0.35 0.8]; colU = [0.85 0.25 0.1];

for iGrp = 1:length(groups)
    g = groups(iGrp);
    if isnan(g.age), continue; end
    m26 = mean(g.conc_uM(:,1:2),2);       % per-scan H2+H6 mean (uM)
    mu  = mean(m26);
    sd  = std(m26);
    if strcmp(g.dataset,'Pilot-Sinc'), col = colS; else, col = colU; end
    if g.nScans > 1
        errorbar(g.age, mu, sd, 'o', 'Color',col, 'MarkerFaceColor',col, ...
            'MarkerSize',8, 'LineWidth',1.2, 'CapSize',6);
    else
        plot(g.age, mu, 'o', 'Color',col, 'MarkerFaceColor',col, 'MarkerSize',8);
    end
    text(g.age+0.3, mu, g.subject, 'Color',col, 'FontSize',8);
end

% dummy handles for legend
hS = plot(nan,nan,'o','Color',colS,'MarkerFaceColor',colS,'MarkerSize',8);
hU = plot(nan,nan,'o','Color',colU,'MarkerFaceColor',colU,'MarkerSize',8);

% linear fits per group (if enough points)
ageS = arrayfun(@(g)g.age, groups(isSinc));
cS   = arrayfun(@(g)mean(mean(g.conc_uM(:,1:2),2)), groups(isSinc));
ageU = arrayfun(@(g)g.age, groups(isU01));
cU   = arrayfun(@(g)mean(mean(g.conc_uM(:,1:2),2)), groups(isU01));
okS = ~isnan(ageS) & ~isnan(cS);
okU = ~isnan(ageU) & ~isnan(cU);
if nnz(okS) >= 2
    p = polyfit(ageS(okS), cS(okS), 1);
    xr = linspace(min(ageS(okS)), max(ageS(okS)), 50);
    plot(xr, polyval(p,xr), '--', 'Color',colS, 'LineWidth',1.2);
end
if nnz(okU) >= 2
    p = polyfit(ageU(okU), cU(okU), 1);
    xr = linspace(min(ageU(okU)), max(ageU(okU)), 50);
    plot(xr, polyval(p,xr), '--', 'Color',colU, 'LineWidth',1.2);
end

xlabel('Age (years)'); ylabel('NAD+ H2+H6 concentration (\muM)')
title('NAD+ (H2+H6 mean) vs age')
legend([hS hU], {'Pilot-Sinc','U01\_RR'}, 'Location','best')
saveas(fig, fullfile(outDir,'NAD_H2H6_vs_age.png'));
fprintf('Figure -> %s\n', fullfile(outDir,'NAD_H2H6_vs_age.png'));

% =====================================================================
% Local: age at first scan; reads birth year from any session's rda
% =====================================================================
function age = getAgeForGroup(datPaths)
    age = NaN;
    birthYear = NaN;
    for iP = 1:numel(datPaths)
        [datDir,~,~] = fileparts(datPaths{iP});
        rdaDir = fullfile(fileparts(datDir),'rdafiles');
        if ~isfolder(rdaDir), continue; end
        rdaList = dir(fullfile(rdaDir,'*.rda'));
        for k = 1:numel(rdaList)
            try
                txt = fileread(fullfile(rdaDir, rdaList(k).name));
                tok = regexp(txt,'PatientBirthDate:\s*(\d{4})','tokens','once');
                if ~isempty(tok)
                    birthYear = str2double(tok{1});
                    break
                end
            catch, continue
            end
        end
        if ~isnan(birthYear), break; end
    end
    if isnan(birthYear), return; end

    % scan date = first 8-digit MMDDYYYY token in first scan's filename
    [~,fn,~] = fileparts(datPaths{1});
    tok = regexp(fn,'_(\d{8})(?=_)','tokens','once');
    if isempty(tok), return; end
    d = tok{1};
    scanMonth = str2double(d(1:2));
    scanYear  = str2double(d(5:8));
    age = scanYear - birthYear;
    if scanMonth < 7, age = age - 0.5; end
end
