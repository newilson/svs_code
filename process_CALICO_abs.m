%% process_CALICO_abs - Absolute NAD+ quantification for CALICO brain data
%
% Reruns the CALICO brain fits with the optimized fitting parameters from
% demo_absoluteQuant_svs.m and converts each fit to absolute [mM] using the
% 3-compartment (WM/GM/CSF) water-reference model (Gasparovic 2006) via
% absoluteQuant_svs.m.
%
% Per-scan tissue composition and subject age come from
%   Z:\Data\Spectroscopy\CALICO\datain\PROCESSED\CALICO_brain_20May2024.xlsx
% (columns DataName / PatAge / CSFfraction / GMfraction / WMfraction and the
% reference concentration [NAD] (mM)).  DataName has form
% 'SC7T_YYYYMMDD_854430_SXXX' whereas folder dates use MMDDYYYY.

clear; clc

addpath('Z:\code\svs_code');
addpath('Z:\code\svs_code\utils');

dataRoot = 'Z:\Data\Spectroscopy\CALICO\datain';
procDir  = fullfile(dataRoot, 'PROCESSED');
xlsxFile = fullfile(procDir, 'CALICO_brain_20May2024.xlsx');

% =====================================================================
% Fit parameters (from demo_absoluteQuant_svs.m)
% =====================================================================
pars.lb_met = 3;
pars.lb_wat = 3;
pars.eccopt = 0;
pars.PC     = 0;
pars.base   = nan;
pars.plt    = false;
pars.pltSpec = false;

pars.block_align.size     = 64;
pars.block_align.method   = 'fd';
pars.block_align.Nfit     = 256;
pars.block_align.ppmRange = [8.5 10];
pars.ccopt.minsig_frac    = 0.05;

pars.WatSupPre.opts.type  = 'none';
pars.WatSupPost.opts.type = 'none';
pars.WatSupPost.opts.filt.type = 'average';
pars.WatSupPost.opts.filt.N    = 11;
pars.WatSupPost.opts.hsvd.bounds = [-1 1] * 300/8000;
pars.WatSupPost.opts.hsvd.nsin   = 25;

pars.dofit = 'varpro';
pars.fit.mode       = 5;                    % Voigt
pars.fit.ppm_range  = [8.7 10.0];
pars.fit.ph_range   = [-60 60];
pars.peaks.name  = {'NADH2','NADH6','NADH4'};
pars.peaks.range = [9.25 9.45; 9.05 9.25; 8.8 9.05];

pars.fit.baselineOpt.enable      = true;
pars.fit.baselineOpt.knotSpacing = 2;
pars.fit.baselineOpt.lambda      = 1;
pars.fit.baselineOpt.lambdaAmpl  = 0;
pars.fit.baselineOpt.TolFun      = 1e-6;

pars.fit.lineShapeOpt.enable     = true;
pars.fit.lineShapeOpt.nSide      = 4;
pars.fit.lineShapeOpt.asymmetric = true;
pars.fit.lineShapeOpt.maxSide    = 0.50;

pars.fit.hsvdClean.enable = false;

pars.watfit.method    = 'fit';
pars.watfit.nComp     = 1;
pars.watfit.ppm_range = [4 5.5];

% =====================================================================
% Literature relaxation / water content (from demo_absoluteQuant_svs.m)
% Swago et al., NMR Biomed 2025 (doi:10.1002/nbm.5324)
% =====================================================================
WM  = struct('name','WM', 'T1_ms',1130, 'T2_ms',40, 'waterConc_M',38.5, 'flag_metSig',true);
GM  = struct('name','GM', 'T1_ms',1939, 'T2_ms',40, 'waterConc_M',44,   'flag_metSig',true);
CSF = struct('name','CSF','T1_ms',4400, 'T2_ms',300,'waterConc_M',55,   'flag_metSig',false);

% From Bagga paper
% WM  = struct('name','WM', 'T1_ms',1500, 'T2_ms',60, 'waterConc_M',38.5, 'flag_metSig',true);
% GM  = struct('name','GM', 'T1_ms',1500, 'T2_ms',60, 'waterConc_M',44,   'flag_metSig',true);
% CSF = struct('name','CSF','T1_ms',1500, 'T2_ms',60,'waterConc_M',55,   'flag_metSig',true);

met_T1_ms = [205.6, 211.6, 237.3];   % NADH2, NADH6, NADH4
met_T2_ms = [ 33.6,  29.1,  42.3];
met_nProt = [    1,     1,     1];
peakNames = {'NADH2','NADH6','NADH4'};
nPeaks    = numel(peakNames);

% =====================================================================
% Read xlsx
% =====================================================================
T = readtable(xlsxFile, 'VariableNamingRule','preserve');
dataNames = string(T.DataName);

% =====================================================================
% Subjects and session dates (brain only; dates in MMDDYYYY)
% =====================================================================
subjects = {'S001','S002','S003','S004','S005','S006','S007','S008','S009','S010', ...
            'S011','S012','S013','S014','S015','S016','S017','S018','S019','S020'};

dates = { ...
    {'12082023','12192023'}, ... % S001
    {'12112023','01022024'}, ... % S002
    {'12142023','01022024'}, ... % S003
    {'12142023','01022024'}, ... % S004
    {'01042024','01112024'}, ... % S005
    {'01082024','01292024'}, ... % S006
    {'01092024','01312024'}, ... % S007
    {'01162024','02052024'}, ... % S008
    {'01232024','02062024'}, ... % S009
    {'02082024','02132024'}, ... % S010
    {'02132024','04022024'}, ... % S011
    {'02152024','04012024'}, ... % S012
    {'04012024','04152024'}, ... % S013
    {'04182024','04222024'}, ... % S014
    {'04232024','05102024'}, ... % S015
    {'05202024','05242024'}, ... % S016
    {'05292024','06052024'}, ... % S017
    {'07032024','07192024'}, ... % S018
    {'07112024','07242024'}, ... % S019
    {'07232024','07262024'}, ... % S020
};

nSubj = length(subjects);
nSess = 2;

% =====================================================================
% Build brain file paths
% =====================================================================
scans = struct('ws',{cell(1,nSess)},'nws',{cell(1,nSess)});
scans = repmat(scans, 1, nSubj);
for iSubj = 1:nSubj
    subj = subjects{iSubj};
    for iSess = 1:nSess
        d = dates{iSubj}{iSess};
        datDir = fullfile(dataRoot, sprintf('854430_%s_%s', subj, d), 'datfiles');
        allFiles = dir(fullfile(datDir, '*.dat'));
        for iF = 1:length(allFiles)
            fn = allFiles(iF).name;
            if contains(fn, 'eja_svs_slaser') || contains(fn, '256points') || contains(fn, '_Rep')
                continue
            end
            if contains(fn, 'calf'), continue; end
            if contains(fn, 'WATEREF')
                scans(iSubj).nws{iSess} = fullfile(datDir, fn);
            elseif contains(fn, 'NAD')
                scans(iSubj).ws{iSess} = fullfile(datDir, fn);
            end
        end
    end
end

% =====================================================================
% Match xlsx rows, pull per-scan tissue fractions, age, reference [NAD]
% =====================================================================
ages       = nan(nSubj,1);
sessFrac   = cell(nSubj, nSess);          % [CSF GM WM] per scan
xlsxNAD_mM = nan(nSubj, nSess);           % reference [NAD] in mM

for iSubj = 1:nSubj
    subj = subjects{iSubj};
    patAges = [];
    for iSess = 1:nSess
        d = dates{iSubj}{iSess};                        % MMDDYYYY
        yyyymmdd = [d(5:8) d(1:2) d(3:4)];
        key = sprintf('SC7T_%s_854430_%s', yyyymmdd, subj);
        idx = find(strcmp(dataNames, key), 1);
        if isempty(idx)
            warning('No xlsx row for %s / %s (%s)', subj, d, key);
            continue
        end
        sessFrac{iSubj,iSess}  = [T.CSFfraction(idx), T.GMfraction(idx), T.WMfraction(idx)];
        xlsxNAD_mM(iSubj,iSess) = T.('[NAD] (mM)')(idx);
        if T.PatAge(idx) > 0
            patAges(end+1) = T.PatAge(idx); %#ok<AGROW>
        end
    end

    if ~isempty(patAges)
        ages(iSubj) = max(patAges);     % anonymized rows have PatAge=0
    else
        % fall back to RDA birth year if both sessions were anonymized
        d1 = dates{iSubj}{1};
        rdaFile = fullfile(dataRoot, sprintf('854430_%s_%s', subj, d1), ...
                            'rdafiles', sprintf('854430_%s_%s_NAD.rda', subj, d1));
        if iSubj == 19 && strcmp(d1, '07112024')
            rdaFile = fullfile(dataRoot, sprintf('854430_%s_%s', subj, d1), ...
                               'rdafiles', sprintf('8554430_%s_%s_NAD.rda', subj, d1));
        end
        if isfile(rdaFile)
            rdaTxt = fileread(rdaFile);
            tok = regexp(rdaTxt, 'PatientBirthDate:\s*(\d{4})', 'tokens', 'once');
            if ~isempty(tok)
                birthYear = str2double(tok{1});
                scanYear  = str2double(d1(5:8));
                scanMonth = str2double(d1(1:2));
                ages(iSubj) = scanYear - birthYear - (scanMonth < 7) * 0.5;
            end
        end
    end
end

% =====================================================================
% Process all brain scans -> absolute concentration
% =====================================================================
conc_uM  = nan(nSubj, nSess, nPeaks);
metAreas = nan(nSubj, nSess, nPeaks);
watAreas = nan(nSubj, nSess);

for iSubj = 1:nSubj
    subj = subjects{iSubj};
    fprintf('\n=== %s (age %.1f) ===\n', subj, ages(iSubj));
    for iSess = 1:nSess
        if isempty(sessFrac{iSubj,iSess})
            fprintf('  Session %d: skipped (no xlsx match)\n', iSess);
            continue
        end
        if isempty(scans(iSubj).ws{iSess}) || isempty(scans(iSubj).nws{iSess})
            fprintf('  Session %d: skipped (missing dat files)\n', iSess);
            continue
        end
        fprintf('  Session %d/%d\n', iSess, nSess);

        try
            [output, parsOut, hdr] = dat2mat_svsNAD(pars, ...
                scans(iSubj).ws{iSess}, scans(iSubj).nws{iSess});
        catch ME
            warning('dat2mat_svsNAD failed for %s s%d: %s', subj, iSess, ME.message);
            continue
        end

        % --- Voxel with per-scan tissue fractions ---
        frac = sessFrac{iSubj,iSess};
        voxel = struct();
        voxel.tissues = [CSF GM WM];
        voxel.tissues(1).fraction = frac(1);
        voxel.tissues(2).fraction = frac(2);
        voxel.tissues(3).fraction = frac(3);

        % --- absoluteQuant_svs inputs ---
        data = struct();
        data.water.area = output.wat.fit.areaTotal;
        for k = 1:nPeaks
            data.metabolites(k).name     = peakNames{k};
            data.metabolites(k).area     = output.met.fit.areas(k);
            data.metabolites(k).nProtons = met_nProt(k);
            data.metabolites(k).T1_ms    = met_T1_ms(k);
            data.metabolites(k).T2_ms    = met_T2_ms(k);
        end
        seq = struct('TR_ms',   output.met.TR_ms, 'TE_ms',   output.met.TE_ms, ...
                     'TR_ws_ms',output.wat.TR_ms, 'TE_ws_ms',output.wat.TE_ms);

        absQ = absoluteQuant_svs(data, seq, voxel);

        conc_uM(iSubj,iSess,:)  = [absQ.metabolites.conc_mM] * 1e3;
        watAreas(iSubj,iSess)   = output.wat.fit.areaTotal;
        metAreas(iSubj,iSess,:) = output.met.fit.areas(1:nPeaks);

        matName = sprintf('%s_brain_abs_s%d', subj, iSess);
        save(fullfile(procDir, [matName '.mat']), 'output','parsOut','hdr','absQ','voxel');
    end
end

save(fullfile(procDir, 'CALICO_abs_brain.mat'), ...
     'subjects','dates','ages','sessFrac','xlsxNAD_mM', ...
     'conc_uM','metAreas','watAreas','peakNames','pars','-v7.3');

% =====================================================================
% Write results to TSV
% =====================================================================
tsvPath = fullfile(procDir, 'CALICO_abs_brain.tsv');
fid = fopen(tsvPath, 'w');
hdr = 'Subject\tAge\tCSF_s1\tGM_s1\tWM_s1\tCSF_s2\tGM_s2\tWM_s2';
for iPk = 1:nPeaks
    hdr = [hdr sprintf('\t%s_uM_s1\t%s_uM_s2', peakNames{iPk}, peakNames{iPk})];
end
hdr = [hdr '\tH2H6_uM_s1\tH2H6_uM_s2\txlsxNAD_uM_s1\txlsxNAD_uM_s2'];
fprintf(fid, '%s\n', hdr);

for iSubj = 1:nSubj
    row = sprintf('%s\t%.1f', subjects{iSubj}, ages(iSubj));
    for iSess = 1:nSess
        f = sessFrac{iSubj,iSess};
        if isempty(f), f = nan(1,3); end
        row = [row sprintf('\t%.3f\t%.3f\t%.3f', f(1), f(2), f(3))];
    end
    for iPk = 1:nPeaks
        row = [row sprintf('\t%.4g\t%.4g', conc_uM(iSubj,1,iPk), conc_uM(iSubj,2,iPk))];
    end
    h26_s1 = mean(squeeze(conc_uM(iSubj,1,1:2)));
    h26_s2 = mean(squeeze(conc_uM(iSubj,2,1:2)));
    row = [row sprintf('\t%.4g\t%.4g\t%.4g\t%.4g', h26_s1, h26_s2, ...
            xlsxNAD_mM(iSubj,1)*1e3, xlsxNAD_mM(iSubj,2)*1e3)];
    fprintf(fid, '%s\n', row);
end
fclose(fid);
fprintf('\nResults -> %s\n', tsvPath);

% =====================================================================
% Plot: current-fit NADH2 [NAD+] vs age compared with xlsx [NAD](mM)*1000
% The xlsx reference is H2-only (columns: amp(NAD+ H2), lw(NAD+ H2),
% [NAD] (mM)), so only the NADH2 peak from the current fit is a valid
% comparison.
% =====================================================================
allAges = []; curr = []; refX = []; lbl = {};
for iSubj = 1:nSubj
    for iSess = 1:nSess
        c_curr = conc_uM(iSubj,iSess,1);                    % NADH2 only (uM)
        c_ref  = xlsxNAD_mM(iSubj,iSess) * 1e3;             % mM -> uM
        allAges(end+1) = ages(iSubj); %#ok<AGROW>
        curr(end+1)    = c_curr;      %#ok<AGROW>
        refX(end+1)    = c_ref;       %#ok<AGROW>
        lbl{end+1}     = sprintf('%s s%d', subjects{iSubj}, iSess); %#ok<AGROW>
    end
end

colC = [0.10 0.35 0.80];
colR = [0.85 0.25 0.10];

fig = figure('Position',[100 100 1400 500],'Color','w');

% --- overlay: current-fit NADH2 vs xlsx H2 ---
subplot(1,2,1); hold on; grid on
h1 = plot(allAges, curr, 'o','MarkerFaceColor',colC,'MarkerEdgeColor',colC,'MarkerSize',7);
% h2 = plot(allAges, refX, 's','MarkerFaceColor',colR,'MarkerEdgeColor',colR,'MarkerSize',7);
okC = ~isnan(curr) & ~isnan(allAges);
okR = ~isnan(refX) & ~isnan(allAges);
if nnz(okC)>=2
    p = polyfit(allAges(okC), curr(okC), 1);
    xr = linspace(min(allAges(okC)), max(allAges(okC)), 50);
    plot(xr, polyval(p,xr), '-', 'Color',colC, 'LineWidth',1.3);
end
% if nnz(okR)>=2
%     p = polyfit(allAges(okR), refX(okR), 1);
%     xr = linspace(min(allAges(okR)), max(allAges(okR)), 50);
%     plot(xr, polyval(p,xr), '-', 'Color',colR, 'LineWidth',1.3);
% end
xlabel('Age (years)'); ylabel('NAD+ H2 concentration (\muM)')
title('NADH2 vs age: current fit vs xlsx reference (H2-only)')
% legend([h1 h2], {'Current fit NADH2','xlsx [NAD] \times1000'}, 'Location','best')

% --- parity: current-fit NADH2 vs xlsx H2 ---
subplot(1,2,2); hold on; grid on
scatter(refX, curr, 40, 'k','filled')
text(refX, curr, lbl, 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',7)
lims = [min([refX curr],[],'omitnan') max([refX curr],[],'omitnan')];
if all(isfinite(lims))
    pad = 0.05 * diff(lims);
    lims = [lims(1)-pad lims(2)+pad];
    plot(lims, lims, 'r--','LineWidth',1);
    xlim(lims); ylim(lims); axis square
end
xlabel('xlsx [NAD] \times1000 (\muM)')
ylabel('Current fit NADH2 (\muM)')
title('Parity: current NADH2 vs xlsx reference')

figPath = fullfile(procDir, 'CALICO_abs_NADvsAge.png');
print(fig, figPath, '-dpng','-r300');
close(fig)
fprintf('Figure -> %s\n', figPath);

% --- per-peak NAD+ vs age (current fit only) ---
fig = figure('Position',[100 100 1500 420],'Color','w');
sgtitle('CALICO brain: NAD+ concentration vs age (per peak)')
ageMat = repmat(ages, 1, nSess);
for iPk = 1:nPeaks
    subplot(1, nPeaks, iPk); hold on; grid on
    cp = squeeze(conc_uM(:,:,iPk));
    scatter(ageMat(:), cp(:), 40, colC,'filled');
    okp = ~isnan(cp(:)) & ~isnan(ageMat(:));
    if nnz(okp)>=2
        p = polyfit(ageMat(okp), cp(okp), 1);
        xr = linspace(min(ageMat(okp)), max(ageMat(okp)), 50);
        plot(xr, polyval(p,xr), '-', 'Color',colC, 'LineWidth',1.3);
    end
    xlabel('Age (years)'); ylabel(sprintf('%s (\\muM)', peakNames{iPk}))
    title(peakNames{iPk})
end
figPath = fullfile(procDir, 'CALICO_abs_NADvsAge_perPeak.png');
print(fig, figPath, '-dpng','-r300');
close(fig)
fprintf('Figure -> %s\n', figPath);
