%% processRepeatability.m
% NAD+ repeatability analysis: Pilot-Sinc-DFMRS (9.1 & 9.7 ppm) + U01_RR brain (9.7 ppm)
% Runner-override hook: if a caller has set OVR_* variables in the workspace,
% use them and skip the default 'clear'. Otherwise behave exactly as before.
if exist('OVR_outName','var') || exist('OVR_outDir','var') || exist('OVR_block_align','var') || exist('OVR_lineshape','var') || exist('OVR_hsvd','var') || exist('OVR_baseline','var') || exist('OVR_ppm_range','var') || exist('OVR_peaks','var')
    % override mode: do NOT clear — the caller manages workspace state.
else
    clear
end
clc

addpath('Z:\code\svs_code');
addpath('Z:\code\svs_code\utils');

if exist('OVR_outDir','var') && ~isempty(OVR_outDir)
    outDir = OVR_outDir;
else
    outDir = 'Z:\Data\Spectroscopy\Pilot-Sinc-DFMRS\PROCESSED';
end
if exist('OVR_outName','var') && ~isempty(OVR_outName)
    outName = OVR_outName;
else
    outName = 'RepeatabilityResults_9_7.tsv';
end
if ~exist(outDir, 'dir'), mkdir(outDir); end

% --- Common parameters (from Case 4 of demo_dat2mat_svsNAD_fitting.m) ---
pars.lb = 3;
pars.eccopt = 0;
pars.PC = 0;
pars.base = nan;
pars.plt = false;
pars.pltSpec = false;
if exist('OVR_block_align','var') && ~isempty(OVR_block_align)
    pars.block_align = OVR_block_align;
else
    pars.block_align.size = 64;
    pars.block_align.method = 'fd';
    pars.block_align.Nfit = 256; % for td only
    pars.block_align.ppmRange = [8.5 10];
end
pars.ccopt.minsig_frac = 0.05;

pars.WatSupPre.opts.type = 'none';
pars.WatSupPost.opts.type = 'none';
pars.WatSupPost.opts.filt.type = 'average';
pars.WatSupPost.opts.filt.N = 11;
pars.WatSupPost.opts.hsvd.bounds = [-1 1] * 300/8000;
pars.WatSupPost.opts.hsvd.nsin = 25;

pars.dofit = 'varpro';
pars.fit.mode = 5; % Voigt
if exist('OVR_ppm_range','var') && ~isempty(OVR_ppm_range)
    pars.fit.ppm_range = OVR_ppm_range;
else
    pars.fit.ppm_range = [8.7 10.0];
end
pars.fit.ph_range = [-60 60];
if exist('OVR_peaks','var') && ~isempty(OVR_peaks)
    pars.peaks.name  = OVR_peaks.name;
    pars.peaks.range = OVR_peaks.range;
else
    pars.peaks.name  = {'NADH2','NADH6','NADH4'};
    pars.peaks.range = [9.25 9.45; 9.05 9.25; 8.8 9.05];
end
if exist('OVR_baseline','var') && ~isempty(OVR_baseline)
    pars.fit.baselineOpt = OVR_baseline;
else
    pars.fit.baselineOpt.enable = true;
    pars.fit.baselineOpt.knotSpacing = 2;
    pars.fit.baselineOpt.lambda = 1;
    pars.fit.baselineOpt.lambdaAmpl = 0;
    pars.fit.baselineOpt.TolFun = 1e-6;
end
if exist('OVR_lineshape','var') && ~isempty(OVR_lineshape)
    pars.fit.lineShapeOpt = OVR_lineshape;
else
    pars.fit.lineShapeOpt.enable = true;
    pars.fit.lineShapeOpt.nSide = 4;
    pars.fit.lineShapeOpt.asymmetric = true;
    pars.fit.lineShapeOpt.maxSide = 0.50;
end
if exist('OVR_hsvd','var') && ~isempty(OVR_hsvd)
    pars.fit.hsvdClean = OVR_hsvd;
else
    pars.fit.hsvdClean.enable = false;
    pars.fit.hsvdClean.K = 60;
    pars.fit.hsvdClean.range = [8.75 12];
end

pars.watfit.method = 'fid';
% pars.watfit.nComp = 7;
% pars.watfit.ppm_range = [4 5.5];

% =====================================================================
% Define all groups
% =====================================================================
sincRoot = 'Z:\Data\Spectroscopy\Pilot-Sinc-DFMRS';
u01Root  = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_RR';

nGrp = 0;

% --- Pilot-Sinc-DFMRS: S004, 9.7 ppm ---
nGrp = nGrp + 1;
groups(nGrp).label = 'S004 sinc 9.7ppm';
groups(nGrp).subject = 'S004';
groups(nGrp).centerFreq = '9.7 ppm';
groups(nGrp).dataset = 'Pilot-Sinc';
groups(nGrp).ws = {
    fullfile(sincRoot, 'S004_01092025', 'datfiles', 'S004_01092025_meas_MID01677_FID03067_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
    fullfile(sincRoot, 'S004_01142025', 'datfiles', 'S004_01142025_meas_MID02134_FID03526_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
};
groups(nGrp).nws = {
    fullfile(sincRoot, 'S004_01092025', 'datfiles', 'S004_01092025_meas_MID01675_FID03065_svssel_WATERREF_2x4_7ppm.dat')
    fullfile(sincRoot, 'S004_01142025', 'datfiles', 'S004_01142025_meas_MID02132_FID03524_svssel_WATERREF_2x4_7ppm.dat')
};

% --- S092 Pilot-Sinc, 9.7 ppm (2 scans) ---
nGrp = nGrp + 1;
groups(nGrp).label = 'S092 sinc 9.7ppm';
groups(nGrp).subject = 'S092';
groups(nGrp).centerFreq = '9.7 ppm';
groups(nGrp).dataset = 'Pilot-Sinc';
groups(nGrp).ws = {
    fullfile(sincRoot, 'S092_01072025', 'datfiles', 'S092_01072025_meas_MID01374_FID02760_svssel_sinc2x2_2X9_7ppm_TR1000_256avg_TE13ms.dat')
    fullfile(sincRoot, 'S092_01162025', 'datfiles', 'S092_01162025_meas_MID02350_FID03744_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
};
groups(nGrp).nws = {
    fullfile(sincRoot, 'S092_01072025', 'datfiles', 'S092_01072025_meas_MID01373_FID02759_svssel_WATERREF_2x4_7ppm_TE13ms.dat')
    fullfile(sincRoot, 'S092_01162025', 'datfiles', 'S092_01162025_meas_MID02348_FID03742_svssel_WATERREF_2x4_7ppm.dat')
};

% --- S092 U01_RR, 9.7 ppm (2 scans) ---
nGrp = nGrp + 1;
groups(nGrp).label = 'S092 U01 9.7ppm';
groups(nGrp).subject = 'S092';
groups(nGrp).centerFreq = '9.7 ppm';
groups(nGrp).dataset = 'U01_RR';
groups(nGrp).ws = {
    fullfile(u01Root, 'S092_03032026_31P-1H-HeadCoil', 'datfiles', 'S092_03032026_meas_MID01621_FID10109_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
    fullfile(u01Root, 'S092_03052026_31P-1H-HeadCoil', 'datfiles', 'S092_03052026_meas_MID01808_FID10296_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
};
groups(nGrp).nws = {
    fullfile(u01Root, 'S092_03032026_31P-1H-HeadCoil', 'datfiles', 'S092_03032026_meas_MID01620_FID10108_svssel_latest_M0_Water_4avg_tra.dat')
    fullfile(u01Root, 'S092_03052026_31P-1H-HeadCoil', 'datfiles', 'S092_03052026_meas_MID01806_FID10294_svssel_latest_M0_Water_4avg_tra.dat')
};

% --- Pilot-Sinc-DFMRS: S092, 9.1 ppm (session 2 only has 9.1 ppm; session 1 does not) ---
% Skipped: S092 session 1 (01072025) does not have a 2X9_1ppm scan

% --- Pilot-Sinc-DFMRS: S139, 9.7 ppm ---
nGrp = nGrp + 1;
groups(nGrp).label = 'S139 sinc 9.7ppm';
groups(nGrp).subject = 'S139';
groups(nGrp).centerFreq = '9.7 ppm';
groups(nGrp).dataset = 'Pilot-Sinc';
groups(nGrp).ws = {
    fullfile(sincRoot, 'S139_01092025', 'datfiles', 'S139_01092025_meas_MID01648_FID03038_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
    fullfile(sincRoot, 'S139_01162025', 'datfiles', 'S139_01162025_meas_MID02381_FID03775_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
};
groups(nGrp).nws = {
    fullfile(sincRoot, 'S139_01092025', 'datfiles', 'S139_01092025_meas_MID01646_FID03036_svssel_WATERREF_2x4_7ppm.dat')
    fullfile(sincRoot, 'S139_01162025', 'datfiles', 'S139_01162025_meas_MID02379_FID03773_svssel_WATERREF_2x4_7ppm.dat')
};

% --- Pilot-Sinc-DFMRS: S145, 9.7 ppm ---
nGrp = nGrp + 1;
groups(nGrp).label = 'S145 sinc 9.7ppm';
groups(nGrp).subject = 'S145';
groups(nGrp).centerFreq = '9.7 ppm';
groups(nGrp).dataset = 'Pilot-Sinc';
groups(nGrp).ws = {
    fullfile(sincRoot, 'S145_01082025', 'datfiles', 'S145_01082025_meas_MID01618_FID03008_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
    fullfile(sincRoot, 'S145_02042025', 'datfiles', 'S145_02042025_meas_MID01442_FID05524_svssel_sinc2x2_2X9_7ppm_TR1000_256avg.dat')
};
groups(nGrp).nws = {
    fullfile(sincRoot, 'S145_01082025', 'datfiles', 'S145_01082025_meas_MID01616_FID03006_svssel_WATERREF_2x4_7ppm.dat')
    fullfile(sincRoot, 'S145_02042025', 'datfiles', 'S145_02042025_meas_MID01440_FID05522_svssel_WATERREF_2x4_7ppm.dat')
};

% --- U01_RR Group 2: S076, brain, 9.7 ppm, 4 repeat scans ---
nGrp = nGrp + 1;
groups(nGrp).label = 'S076 U01 9.7ppm';
groups(nGrp).subject = 'S076';
groups(nGrp).centerFreq = '9.7 ppm';
groups(nGrp).dataset = 'U01_RR';
groups(nGrp).ws = {
    fullfile(u01Root, 'S076_02052026_50mm_31P-1H-HeadCoil', 'datfiles', 'S76_02052026_50mm_meas_MID00223_FID07639_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
    fullfile(u01Root, 'S076_02122026_31P-1H-HeadCoil', 'datfiles', 'S076_02122026_meas_MID00932_FID08354_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
    fullfile(u01Root, 'S076_03102026_31P-1H-HeadCoil', 'datfiles', 'S076_03102026_meas_MID00459_FID10844_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
};
groups(nGrp).nws = {
    fullfile(u01Root, 'S076_02052026_50mm_31P-1H-HeadCoil', 'datfiles', 'S76_02052026_50mm_meas_MID00220_FID07636_svssel_latest_M0_Water_4avg_tra.dat')
    fullfile(u01Root, 'S076_02122026_31P-1H-HeadCoil', 'datfiles', 'S076_02122026_meas_MID00930_FID08352_svssel_latest_M0_Water_4avg_tra.dat')
    fullfile(u01Root, 'S076_03102026_31P-1H-HeadCoil', 'datfiles', 'S076_03102026_meas_MID00458_FID10843_svssel_latest_M0_Water_4avg_tra.dat')
};

% --- U01_RR Group 11: S138, brain, 9.7 ppm, 2 repeat scans ---
nGrp = nGrp + 1;
groups(nGrp).label = 'S138 U01 9.7ppm';
groups(nGrp).subject = 'S138';
groups(nGrp).centerFreq = '9.7 ppm';
groups(nGrp).dataset = 'U01_RR';
groups(nGrp).ws = {
    fullfile(u01Root, 'S138_02172026_31P-1H-HeadCoil', 'datfiles', 'S138_02172026_meas_MID00379_FID08871_svssel_latest_M0_NAD_9_7ppm_512avg_tra.dat')
    fullfile(u01Root, 'S138_03192026_31P-1H-HeadCoil', 'datfiles', 'S138_03192026_meas_MID00036_FID11602_svssel_latest_M0_NAD_9_7ppm_1024avg_tra.dat')
    fullfile(u01Root, 'S138_04092026_31P-1H-HeadCoil', 'datfiles', 'S138_04092026_meas_MID02333_FID13984_svssel_latest_M0_NAD_9_7ppm_1024avg_tra.dat')
    fullfile(u01Root, 'S138_04102026_31P-1H-HeadCoil', 'datfiles', 'S138_04102026_meas_MID02435_FID14086_svssel_latest_M0_NAD_9_7ppm_1024avg_tra.dat')
};
groups(nGrp).nws = {
    fullfile(u01Root, 'S138_02172026_31P-1H-HeadCoil', 'datfiles', 'S138_02172026_meas_MID00378_FID08870_svssel_latest_M0_Water_4avg_tra.dat')
    fullfile(u01Root, 'S138_03192026_31P-1H-HeadCoil', 'datfiles', 'S138_03192026_meas_MID00035_FID11601_svssel_latest_M0_Water_8avg_tra.dat')
    fullfile(u01Root, 'S138_04092026_31P-1H-HeadCoil', 'datfiles', 'S138_04092026_meas_MID02332_FID13983_svssel_latest_M0_Water_8avg_tra.dat')
    fullfile(u01Root, 'S138_04102026_31P-1H-HeadCoil', 'datfiles', 'S138_04102026_meas_MID02434_FID14085_svssel_latest_M0_Water_8avg_tra.dat')
};

% =====================================================================
% Process all groups
% =====================================================================
peakNames = {'NADH2', 'NADH6', 'NADH4'};
nPeaks = length(peakNames);

for iGrp = 1:length(groups)
    nScans = length(groups(iGrp).ws);
    areas    = zeros(nScans, nPeaks);
    watAreas = zeros(nScans, 1);
    outputs  = cell(nScans, 1);

    fprintf('\n\n=== Group %d/%d: %s (%d scans) ===\n', ...
        iGrp, length(groups), groups(iGrp).label, nScans);

    scanDiag = struct([]);
    for iScan = 1:nScans
        fprintf('  Scan %d/%d\n', iScan, nScans);
        [output, ~, ~] = dat2mat_svsNAD(pars, groups(iGrp).ws{iScan}, groups(iGrp).nws{iScan});
        outputs{iScan} = output;
        tmpA = output.met.fit.areas(:)';
        areas(iScan, :) = tmpA(1:nPeaks);  % only track NADH2,H6,H4 (extra peaks if any are ignored here)
        watAreas(iScan) = output.wat.fit.areaTotal;

        % ----- Diagnostic metrics -----
        d = struct();
        d.scanDate = output.met.scandate;
        d.TR_ms    = output.met.TR_ms;
        d.TE_ms    = output.met.TE_ms;
        d.nAvg     = output.met.nAvg;
        d.refvolt  = output.met.refvolt;
        d.seqname  = output.met.seqname;

        % Voigt fit parameters (mode=5): amp(1:n), pos(n+1:2n), widthL(2n+1:3n),
        % ph(3n+1:4n), widthG(4n+1:5n), width(5n+1:6n)
        n = nPeaks;
        p = output.met.fit.pars;
        d.amp     = p(1:n);
        d.pos_ppm = p(n+1 : 2*n);
        d.widthL_Hz = p(2*n+1 : 3*n);
        d.widthG_Hz = p(4*n+1 : 5*n);
        d.width_Hz  = p(5*n+1 : 6*n); % total Voigt FWHM per peak
        d.areas     = output.met.fit.areas(:)';

        % Noise std from a quiet downfield region (far from metabolites)
        noiseMask = output.met.ppm >= 11.5 & output.met.ppm <= 13;
        if nnz(noiseMask) < 30
            noiseMask = output.met.ppm >= 12 & output.met.ppm <= 14;
        end
        if nnz(noiseMask) < 30
            % fallback: first 10% of time-domain FID
            noiseMask = false(size(output.met.ppm));
            noiseMask(1:round(numel(noiseMask)*0.05)) = true;
        end
        d.noiseStd = std(real(output.met.spec(noiseMask)));
        d.snr = d.amp ./ d.noiseStd;

        % Fit residual RMS (in fit window)
        d.residRMS = sqrt(mean((real(output.met.fit.spec(:)) - output.met.fit.spec_fit(:)).^2));

        % Water FWHM from spectrum (magnitude-based, robust)
        watMag = abs(output.wat.spec);
        [wMax, iMax] = max(watMag);
        halfMax = wMax / 2;
        above = watMag >= halfMax;
        % find contiguous run containing iMax
        iL = iMax; while iL > 1 && above(iL-1), iL = iL - 1; end
        iR = iMax; while iR < numel(above) && above(iR+1), iR = iR + 1; end
        d.watFWHM_Hz  = abs(output.wat.hz(iR)  - output.wat.hz(iL));
        d.watFWHM_ppm = abs(output.wat.ppm(iR) - output.wat.ppm(iL));
        d.watPeakAmp  = wMax;
        d.watArea     = output.wat.fit.areaTotal;

        % Keep the spectra needed for post-hoc overlay comparisons
        d.fit_ppm      = output.met.fit.ppm(:);
        d.fit_spec     = real(output.met.fit.spec(:));
        d.fit_spec_fit = output.met.fit.spec_fit(:);
        d.fit_baseline = output.met.fit.baseline(:);
        d.fit_comp     = output.met.fit.comp;
        d.wat_ppm      = output.wat.ppm(:);
        d.wat_spec     = real(output.wat.spec(:));

        if iScan == 1
            scanDiag = d;
        else
            scanDiag(iScan,1) = d; %#ok<AGROW>
        end
    end
    groups(iGrp).diag = scanDiag;

    ratios = areas ./ watAreas;
    cvs = std(ratios, 0, 1) ./ mean(ratios, 1) * 100;

    % Store results
    groups(iGrp).areas    = areas;
    groups(iGrp).watAreas = watAreas;
    groups(iGrp).ratios   = ratios;
    groups(iGrp).cvs      = cvs;
    groups(iGrp).nScans   = nScans;

    % --- Create Nx3 figure ---
    fig = figure('Position', [100 100 1600 250*nScans], 'Visible', 'off');
    for iScan = 1:nScans
        out = outputs{iScan};

        % Column 1: full metabolite spectrum
        subplot(nScans, 3, (iScan-1)*3 + 1)
        plot(out.met.ppm, real(out.met.spec), 'k')
        set(gca, 'xdir', 'reverse')
        xlabel('ppm'), ylabel('Signal')
        title(sprintf('Scan %d - Full Spectrum', iScan))
        axis tight

        % Column 2: zoomed NAD+ region with fit
        subplot(nScans, 3, (iScan-1)*3 + 2)
        plot(out.met.fit.ppm, real(out.met.fit.spec), 'k'), hold on
        plot(out.met.fit.ppm, out.met.fit.spec_fit, 'r', 'LineWidth', 1.5)
        plot(out.met.fit.ppm, out.met.fit.baseline, 'b--')
        for kk = 1:size(out.met.fit.comp, 2)
            plot(out.met.fit.ppm, out.met.fit.comp(:,kk), '--')
        end
        plot(out.met.fit.ppm, real(out.met.fit.spec) - out.met.fit.spec_fit, 'g')
        set(gca, 'xdir', 'reverse')
        xlabel('ppm'), ylabel('Signal')
        title(sprintf('Scan %d - NAD+ Fit', iScan))
        if iScan == 1
            legend([{'data','fit','baseline'}, peakNames, {'residual'}], ...
                'Location', 'best', 'FontSize', 7)
        end
        axis tight

        % Column 3: full water spectrum
        subplot(nScans, 3, (iScan-1)*3 + 3)
        plot(out.wat.ppm, real(out.wat.spec), 'k')
        set(gca, 'xdir', 'reverse')
        xlabel('ppm'), ylabel('Signal')
        title(sprintf('Scan %d - Water', iScan))
        axis tight
    end

    cvStr = '';
    for iPk = 1:nPeaks
        if iPk > 1, cvStr = [cvStr ', ']; end
        cvStr = [cvStr sprintf('%s CV=%.1f%%', peakNames{iPk}, cvs(iPk))];
    end
    sgtitle(sprintf('%s\n%s', groups(iGrp).label, cvStr))

    figName = strrep(groups(iGrp).label, ' ', '_');
    saveas(fig, fullfile(outDir, [figName '.png']));
    saveas(fig, fullfile(outDir, [figName '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', figName);
end

% =====================================================================
% Summary tables
% =====================================================================

%% --- Per-peak CV table for all 9.7 ppm groups ---
fprintf('\n\n');
idx97 = find(strcmp({groups.centerFreq}, '9.7 ppm'));

for iPk = 1:nPeaks
    fprintf('=== %s Varpro CV (9.7 ppm, all groups) ===\n', peakNames{iPk});
    fprintf('%-25s  %6s', 'Group', 'nScans');
    % Dynamic scan columns
    maxScans = max([groups(idx97).nScans]);
    for iS = 1:maxScans
        fprintf('  %12s', sprintf('Scan%d', iS));
    end
    fprintf('  %8s\n', 'CV');
    fprintf('%s\n', repmat('-', 1, 35 + maxScans*14 + 10));
    for jj = 1:length(idx97)
        g = groups(idx97(jj));
        fprintf('%-25s  %6d', g.label, g.nScans);
        for iS = 1:maxScans
            if iS <= g.nScans
                fprintf('  %12.4e', g.ratios(iS, iPk));
            else
                fprintf('  %12s', '');
            end
        end
        fprintf('  %7.1f%%\n', g.cvs(iPk));
    end
    fprintf('\n');
end

%% --- CV of mean(H2+H6+H4) ratio ---
% fprintf('\n=== CV of mean(H2+H6+H4) ratio (9.7 ppm, all groups) ===\n');
% fprintf('%-25s  %6s', 'Group', 'nScans');
% for iS = 1:maxScans
%     fprintf('  %12s', sprintf('Scan%d', iS));
% end
% fprintf('  %8s\n', 'CV');
% fprintf('%s\n', repmat('-', 1, 35 + maxScans*14 + 10));
% for jj = 1:length(idx97)
%     g = groups(idx97(jj));
%     meanR = mean(g.ratios(:, 1:3), 2); % mean of 3 peaks per scan
%     cv = std(meanR) / mean(meanR) * 100;
%     fprintf('%-25s  %6d', g.label, g.nScans);
%     for iS = 1:maxScans
%         if iS <= g.nScans
%             fprintf('  %12.4e', meanR(iS));
%         else
%             fprintf('  %12s', '');
%         end
%     end
%     fprintf('  %7.1f%%\n', cv);
% end

%% --- CV of mean(H2+H6) ratio ---
fprintf('\n=== CV of mean(H2+H6) ratio (9.7 ppm, all groups) ===\n');
fprintf('%-25s  %6s', 'Group', 'nScans');
for iS = 1:maxScans
    fprintf('  %12s', sprintf('Scan%d', iS));
end
fprintf('  %8s\n', 'CV');
fprintf('%s\n', repmat('-', 1, 35 + maxScans*14 + 10));
for jj = 1:length(idx97)
    g = groups(idx97(jj));
    meanR = mean(g.ratios(:, 1:2), 2); % mean of H2+H6 per scan
    cv = std(meanR) / mean(meanR) * 100;
    fprintf('%-25s  %6d', g.label, g.nScans);
    for iS = 1:maxScans
        if iS <= g.nScans
            fprintf('  %12.4e', meanR(iS));
        else
            fprintf('  %12s', '');
        end
    end
    fprintf('  %7.1f%%\n', cv);
end

%% --- HSVD comparison (Pilot-Sinc 9.7 ppm only) ---
% hsvd.subjects = {'S092', 'S139', 'S004', 'S145'};
% hsvd.scan1    = [0.390, 0.447, 0.215, 0.340]; % [NAD-H2] (mM), 2x9.7 ppm
% hsvd.scan2    = [0.450, 0.418, 0.227, 0.341];
% for ii = 1:length(hsvd.subjects)
%     m = mean([hsvd.scan1(ii), hsvd.scan2(ii)]);
%     s = std([hsvd.scan1(ii), hsvd.scan2(ii)]);
%     hsvd.cv(ii) = s / m * 100;
% end
% 
% fprintf('\n\n=== NAD+ H2 CV Comparison: HSVD vs Varpro (Pilot-Sinc 2x9.7 ppm) ===\n');
% fprintf('%-8s  %12s %12s  %12s %12s  %8s %8s\n', ...
%     'Subject', 'HSVD Scan1', 'HSVD Scan2', 'VP Scan1', 'VP Scan2', 'HSVD CV', 'VP CV');
% fprintf('%s\n', repmat('-', 1, 80));
% for ii = 1:length(hsvd.subjects)
%     % Find matching Pilot-Sinc 9.7 ppm group
%     idx = find(strcmp({groups.subject}, hsvd.subjects{ii}) & ...
%                strcmp({groups.centerFreq}, '9.7 ppm') & ...
%                strcmp({groups.dataset}, 'Pilot-Sinc'));
%     if ~isempty(idx)
%         g = groups(idx(1));
%         vp_s1 = g.ratios(1, 1);
%         vp_s2 = g.ratios(2, 1);
%     else
%         vp_s1 = NaN; vp_s2 = NaN;
%     end
%     fprintf('%-8s  %12.4f %12.4f  %12.4e %12.4e  %7.1f%% %7.1f%%\n', ...
%         hsvd.subjects{ii}, ...
%         hsvd.scan1(ii), hsvd.scan2(ii), ...
%         vp_s1, vp_s2, ...
%         hsvd.cv(ii), g.cvs(1));
% end

%% --- Write summary TSV ---
fid = fopen(fullfile(outDir, outName), 'w');
% Header
hdr = 'Group\tSubject\tDataset\tCenterFreq\tPeak\tnScans';
maxScans = max([groups.nScans]);
for iS = 1:maxScans
    hdr = [hdr sprintf('\tarea_s%d\twatArea_s%d\tratio_s%d', iS, iS, iS)];
end
hdr = [hdr '\tMean_ratio\tStd_ratio\tCV'];
fprintf(fid, '%s\n', hdr);
for iGrp = 1:length(groups)
    g = groups(iGrp);
    for iPk = 1:nPeaks
        row = sprintf('%s\t%s\t%s\t%s\t%s\t%d', ...
            g.label, g.subject, g.dataset, g.centerFreq, peakNames{iPk}, g.nScans);
        for iS = 1:maxScans
            if iS <= g.nScans
                row = [row sprintf('\t%.6g\t%.6g\t%.6g', ...
                    g.areas(iS, iPk), g.watAreas(iS), g.ratios(iS, iPk))];
            else
                row = [row sprintf('\t\t\t')];
            end
        end
        row = [row sprintf('\t%.6g\t%.6g\t%.2f', ...
            mean(g.ratios(:, iPk)), std(g.ratios(:, iPk)), g.cvs(iPk))];
        fprintf(fid, '%s\n', row);
    end
end
fclose(fid);
fprintf('\nResults written to %s\n', fullfile(outDir, outName));

%% --- Save full groups struct for later re-analysis ---
save(fullfile(outDir, 'groups.mat'), 'groups', '-v7.3');
fprintf('Groups struct saved to %s\n', fullfile(outDir, 'groups.mat'));

%% --- Summary printout ---
fprintf('\n=== All Groups Summary ===\n');
for iGrp = 1:length(groups)
    g = groups(iGrp);
    fprintf('%-25s (%d scans): ', g.label, g.nScans);
    for iPk = 1:nPeaks
        fprintf('%s CV=%.1f%%  ', peakNames{iPk}, g.cvs(iPk));
    end
    meanR26  = mean(g.ratios(:, 1:2), 2);
    meanR246 = mean(g.ratios(:, 1:3), 2);
    cv26  = std(meanR26)  / mean(meanR26)  * 100;
    cv246 = std(meanR246) / mean(meanR246) * 100;
    fprintf('H2+H6 CV=%.1f%%  H2+H6+H4 CV=%.1f%%', cv26, cv246);
    fprintf('\n');
end

% =====================================================================
% Per-scan diagnostics (focus on what drives CV)
% =====================================================================

%% --- CV of met area vs CV of water-normalized ratio (to isolate water's role) ---
fprintf('\n\n=== CV Decomposition: raw met area vs water area vs H2+H6 ratio (9.7 ppm) ===\n');
fprintf('%-22s  %6s  %10s  %10s  %10s\n', 'Group', 'nScans', 'CV_met(H2+H6)', 'CV_water', 'CV_ratio(H2+H6)');
fprintf('%s\n', repmat('-', 1, 70));
for jj = 1:length(idx97)
    g = groups(idx97(jj));
    meanMet26 = mean(g.areas(:, 1:2), 2);
    cvMet26   = std(meanMet26) / mean(meanMet26) * 100;
    cvWat     = std(g.watAreas) / mean(g.watAreas) * 100;
    meanR26   = mean(g.ratios(:, 1:2), 2);
    cvR26     = std(meanR26) / mean(meanR26) * 100;
    fprintf('%-22s  %6d  %12.1f%%  %8.1f%%  %12.1f%%\n', ...
        g.label, g.nScans, cvMet26, cvWat, cvR26);
end
fprintf('Interpretation: if CV_ratio >> CV_met, the water reference is adding variability.\n');
fprintf('                if CV_ratio ~< CV_met, water normalization is stabilizing the metric.\n');

%% --- Per-scan diagnostics table for each group ---
fprintf('\n\n=== Per-scan diagnostics (9.7 ppm groups) ===\n');
diagFile = fullfile(outDir, 'repeatability_diagnostics.tsv');
fd = fopen(diagFile, 'w');
fprintf(fd, ['Group\tSubject\tScan\tScanDate\tTR_ms\tTE_ms\tnAvg\tRefVolt\t' ...
             'watArea\twatFWHM_Hz\tnoiseStd\tresidRMS\t' ...
             'H2_amp\tH2_SNR\tH2_widthHz\tH2_area\tH2_ratio\t' ...
             'H6_amp\tH6_SNR\tH6_widthHz\tH6_area\tH6_ratio\t' ...
             'H4_amp\tH4_SNR\tH4_widthHz\tH4_area\tH4_ratio\t' ...
             'mean_H26_ratio\tdev_from_grpMean_%%\n']);

for jj = 1:length(idx97)
    g = groups(idx97(jj));
    meanR26 = mean(g.ratios(:, 1:2), 2);
    grpMean = mean(meanR26);

    fprintf('\n--- %s (%d scans) ---\n', g.label, g.nScans);
    fprintf('%-12s  %-10s  %5s %5s %5s  %10s %10s  %8s  %8s  %8s  %8s  %10s  %7s\n', ...
        'Scan', 'Date', 'TR', 'TE', 'nAvg', 'WatArea', 'WatFWHM_Hz', ...
        'H2_SNR', 'H6_SNR', 'H2_LW', 'H6_LW', 'H26_ratio', 'Dev%');
    fprintf('%s\n', repmat('-', 1, 130));
    for iS = 1:g.nScans
        d = g.diag(iS);
        dateStr = '';
        if ~isnat(d.scanDate), dateStr = datestr(d.scanDate, 'mm-dd-yyyy'); end
        devPct = (meanR26(iS) - grpMean) / grpMean * 100;
        fprintf('%-12s  %-10s  %5.0f %5.0f %5d  %10.3f %10.2f  %8.1f  %8.1f  %8.2f  %8.2f  %10.3e  %+6.1f%%\n', ...
            sprintf('Scan%d', iS), dateStr, d.TR_ms, d.TE_ms, d.nAvg, ...
            d.watArea, d.watFWHM_Hz, ...
            d.snr(1), d.snr(2), d.width_Hz(1), d.width_Hz(2), ...
            meanR26(iS), devPct);

        % Write to TSV
        fprintf(fd, '%s\t%s\t%d\t%s\t%.1f\t%.1f\t%d\t%.1f\t%.4g\t%.3f\t%.4g\t%.4g', ...
            g.label, g.subject, iS, dateStr, d.TR_ms, d.TE_ms, d.nAvg, d.refvolt, ...
            d.watArea, d.watFWHM_Hz, d.noiseStd, d.residRMS);
        for k = 1:nPeaks
            fprintf(fd, '\t%.4g\t%.2f\t%.3f\t%.4g\t%.4g', ...
                d.amp(k), d.snr(k), d.width_Hz(k), d.areas(k), g.ratios(iS, k));
        end
        fprintf(fd, '\t%.4g\t%+.2f\n', meanR26(iS), devPct);
    end
end
fclose(fd);
fprintf('\nDiagnostics TSV written to %s\n', diagFile);

%% --- Focus: 4-scan groups, check for temporal trends ---
idx4 = idx97(arrayfun(@(g) g.nScans >= 4, groups(idx97)));
if ~isempty(idx4)
    fprintf('\n\n=== 4-scan groups: temporal/date analysis ===\n');
    for jj = 1:length(idx4)
        g = groups(idx4(jj));
        dates = [g.diag.scanDate];
        meanR26 = mean(g.ratios(:, 1:2), 2);

        % Sort by date
        [sortedDates, iSrt] = sort(dates);
        sortedRatios = meanR26(iSrt);

        fprintf('\n%s:\n', g.label);
        fprintf('  %-12s  %-12s  %10s  %7s\n', 'Scan#', 'Date', 'H26_ratio', 'dDays');
        for kk = 1:length(sortedDates)
            if kk == 1
                dd = 0;
            else
                dd = days(sortedDates(kk) - sortedDates(1));
            end
            fprintf('  %-12d  %-12s  %10.3e  %7d\n', iSrt(kk), ...
                datestr(sortedDates(kk),'mm-dd-yyyy'), sortedRatios(kk), dd);
        end

        % Correlation of ratio vs date (days)
        dD = days(sortedDates - sortedDates(1));
        if std(sortedRatios) > 0 && std(dD) > 0
            x = dD(:) - mean(dD); y = sortedRatios(:) - mean(sortedRatios);
            r = sum(x.*y) / sqrt(sum(x.^2)*sum(y.^2));
            fprintf('  Pearson r (ratio vs days) = %+.3f\n', r);
        end

        % Which scan drives the CV? (remove each scan, recompute CV)
        fprintf('  Leave-one-out CV(H2+H6 ratio):\n');
        cvFull = std(meanR26) / mean(meanR26) * 100;
        fprintf('    Full (n=%d):       CV = %5.1f%%\n', g.nScans, cvFull);
        for kk = 1:g.nScans
            keep = true(g.nScans,1); keep(kk) = false;
            sub = meanR26(keep);
            cvSub = std(sub) / mean(sub) * 100;
            fprintf('    Remove scan %d:    CV = %5.1f%%  (delta = %+5.1f pp)\n', ...
                kk, cvSub, cvSub - cvFull);
        end

        % --- Diagnostic figure: per-scan metrics ---
        fig = figure('Position',[100 100 1400 700],'Visible','off');
        nS = g.nScans;
        xLabels = arrayfun(@(iS) sprintf('S%d\n%s', iS, ...
            datestr(g.diag(iS).scanDate,'mm/dd')), 1:nS, 'UniformOutput', false);
        xIdx = 1:nS;

        subplot(2,3,1)
        bar(xIdx, meanR26, 'FaceColor',[0.3 0.5 0.8])
        yline(mean(meanR26),'--k','mean')
        set(gca,'XTick',xIdx,'XTickLabel',xLabels)
        ylabel('H2+H6 ratio'), title(sprintf('Ratio per scan (CV=%.1f%%)',cvFull))
        grid on

        subplot(2,3,2)
        bar(xIdx, g.watAreas, 'FaceColor',[0.2 0.7 0.5])
        cvW = std(g.watAreas)/mean(g.watAreas)*100;
        yline(mean(g.watAreas),'--k','mean')
        set(gca,'XTick',xIdx,'XTickLabel',xLabels)
        ylabel('Water area'), title(sprintf('Water area (CV=%.1f%%)',cvW))
        grid on

        subplot(2,3,3)
        watFW = [g.diag.watFWHM_Hz];
        bar(xIdx, watFW, 'FaceColor',[0.8 0.5 0.2])
        set(gca,'XTick',xIdx,'XTickLabel',xLabels)
        ylabel('Water FWHM (Hz)'), title('Water linewidth')
        grid on

        subplot(2,3,4)
        snrMat = vertcat(g.diag.snr);
        bar(xIdx, snrMat(:,1:2))
        legend('H2','H6','Location','best')
        set(gca,'XTick',xIdx,'XTickLabel',xLabels)
        ylabel('SNR'), title('NAD+ SNR')
        grid on

        subplot(2,3,5)
        lwMat = vertcat(g.diag.width_Hz);
        bar(xIdx, lwMat(:,1:2))
        legend('H2','H6','Location','best')
        set(gca,'XTick',xIdx,'XTickLabel',xLabels)
        ylabel('FWHM (Hz)'), title('NAD+ linewidth')
        grid on

        subplot(2,3,6)
        residAll = [g.diag.residRMS];
        bar(xIdx, residAll, 'FaceColor',[0.6 0.3 0.6])
        set(gca,'XTick',xIdx,'XTickLabel',xLabels)
        ylabel('Residual RMS'), title('Fit residual RMS')
        grid on

        sgtitle(sprintf('%s — per-scan diagnostics', g.label))

        dName = [strrep(g.label,' ','_') '_diag'];
        saveas(fig, fullfile(outDir, [dName '.png']));
        saveas(fig, fullfile(outDir, [dName '.fig']));
        close(fig);

        % --- Side-by-side spectrum comparison (overlay + small multiples) ---
        fig2 = figure('Position',[100 100 1600 900],'Visible','off');
        nS = g.nScans;
        cmap = lines(nS);
        legEntries = arrayfun(@(iS) sprintf('Scan %d (%s)  r=%.3e', iS, ...
            datestr(g.diag(iS).scanDate,'mm/dd/yy'), meanR26(iS)), ...
            1:nS, 'UniformOutput', false);

        % Panel 1 (top-left): NAD+ fit region overlay — data
        subplot(3,3,[1 2])
        hold on
        for iS = 1:nS
            plot(g.diag(iS).fit_ppm, g.diag(iS).fit_spec, 'Color', cmap(iS,:), 'LineWidth', 1)
        end
        set(gca,'xdir','reverse'), grid on, axis tight
        xlabel('ppm'), ylabel('Signal')
        title('NAD+ region data overlay'), legend(legEntries,'Location','bestoutside','FontSize',7)

        % Panel 2 (top-right): fit overlay
        subplot(3,3,3)
        hold on
        for iS = 1:nS
            plot(g.diag(iS).fit_ppm, g.diag(iS).fit_spec_fit, 'Color', cmap(iS,:), 'LineWidth', 1.2)
        end
        set(gca,'xdir','reverse'), grid on, axis tight
        xlabel('ppm'), title('Total fit overlay')

        % Panel 3 (middle-left): residuals overlay
        subplot(3,3,[4 5])
        hold on
        for iS = 1:nS
            resid = g.diag(iS).fit_spec - g.diag(iS).fit_spec_fit;
            plot(g.diag(iS).fit_ppm, resid, 'Color', cmap(iS,:), 'LineWidth', 1)
        end
        yline(0,'k:')
        set(gca,'xdir','reverse'), grid on, axis tight
        xlabel('ppm'), ylabel('Residual')
        title('Fit residuals overlay')

        % Panel 4 (middle-right): water spectrum overlay (zoomed around peak)
        subplot(3,3,6)
        hold on
        watXLim = [4.2 5.2];
        for iS = 1:nS
            plot(g.diag(iS).wat_ppm, g.diag(iS).wat_spec, 'Color', cmap(iS,:), 'LineWidth', 1)
        end
        set(gca,'xdir','reverse'), xlim(watXLim), grid on
        xlabel('ppm'), title('Water spectrum overlay')

        % Panel 5-8 (bottom row): per-scan NAD+ fit panels with components
        for iS = 1:nS
            subplot(3, nS, 2*nS + iS)
            d = g.diag(iS);
            plot(d.fit_ppm, d.fit_spec, 'k'), hold on
            plot(d.fit_ppm, d.fit_spec_fit, 'r', 'LineWidth', 1.5)
            plot(d.fit_ppm, d.fit_baseline, 'b--')
            for kk = 1:size(d.fit_comp,2)
                plot(d.fit_ppm, d.fit_comp(:,kk), '--')
            end
            plot(d.fit_ppm, d.fit_spec - d.fit_spec_fit, 'g')
            set(gca,'xdir','reverse'), grid on, axis tight
            xlabel('ppm')
            title(sprintf('Scan %d: r=%.2e  dev=%+.1f%%', iS, meanR26(iS), ...
                (meanR26(iS)-mean(meanR26))/mean(meanR26)*100))
        end

        sgtitle(sprintf('%s — side-by-side scan comparison', g.label))

        cName = [strrep(g.label,' ','_') '_compare'];
        saveas(fig2, fullfile(outDir, [cName '.png']));
        saveas(fig2, fullfile(outDir, [cName '.fig']));
        close(fig2);
    end
end
