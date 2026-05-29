%% process_CALICO - Varpro fitting for CALICO study (brain + calf)
clear, clc

%%
addpath('Z:\code\svs_code'); % for dat2mat_svsNAD.m

dataRoot = 'Z:\Data\Spectroscopy\CALICO\datain';

% --- Common parameters (matching demo Case 4) ---
parsIn.lb = 3;
parsIn.eccopt = -1;
parsIn.PC = 0;
parsIn.base = nan;
parsIn.plt = false;
parsIn.pltSpec = false;
parsIn.block_size = [];
parsIn.ccopt.minsig_frac = 0.05;

parsIn.WatSupPre.opts.type = 'none';
parsIn.WatSupPost.opts.type = 'none';
parsIn.WatSupPost.opts.filt.type = 'average';
parsIn.WatSupPost.opts.filt.N = 11;
parsIn.WatSupPost.opts.hsvd.bounds = [-1 1] * 300/8000;
parsIn.WatSupPost.opts.hsvd.nsin = 25;

% fit parameters
parsIn.dofit = 'varpro';
parsIn.fit.mode = 5;
parsIn.fit.ppm_range_coarse = [7.0 9.7];
parsIn.fit.ppm_range = [8.9 9.7];          % focus range for Step 2
parsIn.fit.focusWeight = 5;
parsIn.fit.ph_range = [-60 60];
parsIn.peaks.name  = {'NADH2','NADH6','NADH4','NAA_NH', ...
                      'S1','S2','S3','S4','S5','S6','S7','S8'};
parsIn.peaks.range = [9.28 9.38; 9.09 9.19; 8.78 8.88; 7.82 7.88; ...
                      7.24 7.30; 7.48 7.54; 7.90 7.96; 8.01 8.07; ...
                      8.13 8.19; 8.24 8.30; 8.36 8.42; 8.47 8.53];
parsIn.fit.baselineOpt.enable = true;
parsIn.fit.baselineOpt.knotSpacing = 2;
parsIn.fit.baselineOpt.lambda = 1;
parsIn.fit.baselineOpt.lambdaAmpl = 0.1;
parsIn.fit.baselineOpt.TolFun = 1e-6;
parsIn.fit.lineShapeOpt.enable = true;
parsIn.fit.lineShapeOpt.nSide = 2;
parsIn.fit.lineShapeOpt.asymmetric = true;
parsIn.fit.lineShapeOpt.maxSide = 0.50;

% water area estimation: 'integrate', 'fid', or 'fit'
parsIn.watfit.method = 'fid';
parsIn.watfit.ppm_range = [4 5.5];

% =====================================================================
% Define subjects and sessions
% =====================================================================
subjects = {'S001','S002','S003','S004','S005','S006','S007','S008','S009','S010', ...
            'S011','S012','S013','S014','S015','S016','S017','S018','S019','S020'};

% Session dates for each subject (chronological order): {date1, date2}
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

% =====================================================================
% Build file paths and read birth years
% =====================================================================
nSubj = length(subjects);
ages = zeros(nSubj, 1);

for iSubj = 1:nSubj
    subj = subjects{iSubj};
    d1 = dates{iSubj}{1};
    folder1 = fullfile(dataRoot, sprintf('854430_%s_%s', subj, d1));

    % Read birth year from first session's RDA file
    rdaFile = fullfile(folder1, 'rdafiles', sprintf('854430_%s_%s_NAD.rda', subj, d1));
    if iSubj == 19 && strcmp(d1, '07112024')
        % S019 session 1 has typo in rda filename
        rdaFile = fullfile(folder1, 'rdafiles', sprintf('8554430_%s_%s_NAD.rda', subj, d1));
    end
    rdaTxt = fileread(rdaFile);
    tok = regexp(rdaTxt, 'PatientBirthDate:\s*(\d{4})', 'tokens', 'once');
    birthYear = str2double(tok{1});

    % Compute age at first scan
    scanMonth = str2double(d1(1:2));
    scanYear  = str2double(d1(5:8));
    ages(iSubj) = scanYear - birthYear;
    if scanMonth < 7
        ages(iSubj) = ages(iSubj) - 0.5; % approximate
    end

    for iSess = 1:2
        d = dates{iSubj}{iSess};
        folder = fullfile(dataRoot, sprintf('854430_%s_%s', subj, d));
        datDir = fullfile(folder, 'datfiles');

        % Find brain files (no 'calf' in name)
        allFiles = dir(fullfile(datDir, '*.dat'));
        for iF = 1:length(allFiles)
            fn = allFiles(iF).name;
            if contains(fn, 'eja_svs_slaser') || contains(fn, '256points') || contains(fn, '_Rep')
                continue
            end
            if ~contains(fn, 'calf')
                if contains(fn, 'WATEREF')
                    scans(iSubj).brain.nws{iSess} = fullfile(datDir, fn);
                elseif contains(fn, 'NAD')
                    scans(iSubj).brain.ws{iSess} = fullfile(datDir, fn);
                end
            else
                if contains(fn, 'WATEREF')
                    scans(iSubj).calf.nws{iSess} = fullfile(datDir, fn);
                elseif contains(fn, 'NAD')
                    scans(iSubj).calf.ws{iSess} = fullfile(datDir, fn);
                end
            end
        end
    end
end

% =====================================================================
% Process all subjects: brain and calf, 2 sessions each
% =====================================================================
peakNames = {'NADH2','NADH6','NADH4'};

for iSubj = 1:nSubj
    for bp = {'brain','calf'}
        bpStr = bp{1};
        fprintf('=== %s (age ~%d) %s ===\n', subjects{iSubj}, round(ages(iSubj)), bpStr);

        for iSess = 1:2
            fprintf('  Session %d/2\n', iSess);
            [output, pars, hdr] = dat2mat_svsNAD(parsIn, scans(iSubj).(bpStr).ws{iSess}, scans(iSubj).(bpStr).nws{iSess});

            results(iSubj).(bpStr).NADH2_amp(iSess) = output.met.fit.areas(1);
            results(iSubj).(bpStr).NADH6_amp(iSess) = output.met.fit.areas(2);
            results(iSubj).(bpStr).NADH4_amp(iSess) = output.met.fit.areas(3);
            results(iSubj).(bpStr).wat_amp(iSess)   = output.wat.fit.areaTotal;

            % Save output to .mat file
            matName = sprintf('%s_%s_s%d', subjects{iSubj}, bpStr, iSess);
            save([dataRoot filesep 'PROCESSED' filesep matName '.mat'], 'output', 'pars', 'hdr');
        end

        % NAD/water ratios
        for iPk = 1:length(peakNames)
            pk = peakNames{iPk};
            r = results(iSubj).(bpStr).([pk '_amp']) ./ results(iSubj).(bpStr).wat_amp;
            results(iSubj).(bpStr).([pk '_ratio']) = r;
            results(iSubj).(bpStr).([pk '_mean'])  = mean(r);
            results(iSubj).(bpStr).([pk '_std'])   = std(r);
            results(iSubj).(bpStr).([pk '_cv'])    = std(r) / mean(r) * 100;
        end
    end
end

%% =====================================================================
% Write results to TSV
% =====================================================================
nSess = 2;
fid = fopen([dataRoot filesep 'PROCESSED' filesep 'CALICO_varpro_repeatability.tsv'], 'w');

% Header
hdr = sprintf('Subject\tAge\tBodyPart');
for s = 1:nSess
    hdr = [hdr sprintf('\tWat_amp_s%d', s)];
end
for iPk = 1:length(peakNames)
    for s = 1:nSess
        hdr = [hdr sprintf('\t%s_amp_s%d', peakNames{iPk}, s)];
    end
end
for iPk = 1:length(peakNames)
    for s = 1:nSess
        hdr = [hdr sprintf('\t%s_ratio_s%d', peakNames{iPk}, s)];
    end
end
for iPk = 1:length(peakNames)
    hdr = [hdr sprintf('\t%s_mean\t%s_std\t%s_cv', peakNames{iPk}, peakNames{iPk}, peakNames{iPk})];
end
fprintf(fid, '%s\n', hdr);

% Rows
for iSubj = 1:nSubj
    for bp = {'brain','calf'}
        bpStr = bp{1};
        row = sprintf('%s\t%d\t%s', subjects{iSubj}, round(ages(iSubj)), bpStr);

        % water amplitudes
        for s = 1:nSess
            row = [row sprintf('\t%.6g', results(iSubj).(bpStr).wat_amp(s))];
        end

        % NAD+ peak amplitudes
        for iPk = 1:length(peakNames)
            for s = 1:nSess
                row = [row sprintf('\t%.6g', results(iSubj).(bpStr).([peakNames{iPk} '_amp'])(s))];
            end
        end

        % NAD+ / water ratios
        for iPk = 1:length(peakNames)
            for s = 1:nSess
                row = [row sprintf('\t%.6g', results(iSubj).(bpStr).([peakNames{iPk} '_ratio'])(s))];
            end
        end

        % mean, std, cv
        for iPk = 1:length(peakNames)
            row = [row sprintf('\t%.6g\t%.6g\t%.2f', ...
                results(iSubj).(bpStr).([peakNames{iPk} '_mean']), ...
                results(iSubj).(bpStr).([peakNames{iPk} '_std']), ...
                results(iSubj).(bpStr).([peakNames{iPk} '_cv']))];
        end

        fprintf(fid, '%s\n', row);
    end
end

fclose(fid);
fprintf('Results written to CALICO_varpro_repeatability.tsv\n');

%% NAD+ ratio vs age plots
for bp = {'brain','calf'}
    bpStr = bp{1};
    fig = figure('Position', [100 100 1400 400], 'Visible', 'off');
    sgtitle(sprintf('NAD+ / Water vs Age (%s)', bpStr))

    for iPk = 1:length(peakNames)
        subplot(1, 3, iPk)
        ratios = zeros(nSubj, 1);
        for iSubj = 1:nSubj
            ratios(iSubj) = results(iSubj).(bpStr).([peakNames{iPk} '_mean']);
        end
        scatter(ages, ratios, 50, 'k', 'filled'), hold on
        p = polyfit(ages, ratios, 1);
        xline = linspace(min(ages), max(ages), 100);
        plot(xline, polyval(p, xline), 'r-', 'LineWidth', 1.5)
        xlabel('Age (years)'), ylabel('NAD+ / Water')
        title(peakNames{iPk})
    end

    print(fig, [dataRoot filesep 'PROCESSED' filesep 'NAD_vs_age_' bpStr '.png'], '-dpng', '-r300')
    close(fig)
    fprintf('Saved NAD_vs_age_%s.png\n', bpStr)
end

%% Session 1 vs Session 2 repeatability scatter plots
for bp = {'brain','calf'}
    bpStr = bp{1};
    fig = figure('Position', [100 100 1400 400], 'Visible', 'off');
    sgtitle(sprintf('Session 1 vs Session 2 NAD+/Water (%s)', bpStr))

    for iPk = 1:length(peakNames)
        subplot(1, 3, iPk)
        s1 = zeros(nSubj, 1);
        s2 = zeros(nSubj, 1);
        for iSubj = 1:nSubj
            s1(iSubj) = results(iSubj).(bpStr).([peakNames{iPk} '_ratio'])(1);
            s2(iSubj) = results(iSubj).(bpStr).([peakNames{iPk} '_ratio'])(2);
        end
        scatter(s1, s2, 50, 'k', 'filled'), hold on
        text(s1, s2, subjects, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 7)
        % identity line
        ax = [min([s1; s2]) max([s1; s2])];
        margin = 0.05 * diff(ax);
        ax = [ax(1)-margin ax(2)+margin];
        plot(ax, ax, 'r--', 'LineWidth', 1)
        xlim(ax), ylim(ax), axis square
        xlabel('Session 1 NAD+/Water'), ylabel('Session 2 NAD+/Water')
        title(peakNames{iPk})
    end

    print(fig, [dataRoot filesep 'PROCESSED' filesep 'NAD_s1_vs_s2_' bpStr '.png'], '-dpng', '-r300')
    close(fig)
    fprintf('Saved NAD_s1_vs_s2_%s.png\n', bpStr)
end

%% Peak-vs-peak scatter plots (H2 vs H6, H4 vs H6)
peakPairs = {{'NADH2','NADH6'}, {'NADH4','NADH6'}};
for bp = {'brain','calf'}
    bpStr = bp{1};
    fig = figure('Position', [100 100 900 400], 'Visible', 'off');
    sgtitle(sprintf('Peak-vs-Peak NAD+/Water (%s)', bpStr))

    for iPair = 1:length(peakPairs)
        pkX = peakPairs{iPair}{1};
        pkY = peakPairs{iPair}{2};
        subplot(1, 2, iPair)
        xAll = zeros(nSubj * 2, 1);
        yAll = zeros(nSubj * 2, 1);
        labels = cell(nSubj * 2, 1);
        for iSubj = 1:nSubj
            for iSess = 1:2
                idx = (iSubj - 1) * 2 + iSess;
                xAll(idx) = results(iSubj).(bpStr).([pkX '_ratio'])(iSess);
                yAll(idx) = results(iSubj).(bpStr).([pkY '_ratio'])(iSess);
                labels{idx} = sprintf('%s s%d', subjects{iSubj}, iSess);
            end
        end
        scatter(xAll, yAll, 50, 'k', 'filled'), hold on
        text(xAll, yAll, labels, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 6)
        ax = [min([xAll; yAll]) max([xAll; yAll])];
        margin = 0.05 * diff(ax);
        ax = [ax(1)-margin ax(2)+margin];
        plot(ax, ax, 'r--', 'LineWidth', 1)
        xlabel([pkX ' / Water']), ylabel([pkY ' / Water'])
        title(sprintf('%s vs %s', pkX, pkY))
    end

    print(fig, [dataRoot filesep 'PROCESSED' filesep 'NAD_peak_vs_peak_' bpStr '.png'], '-dpng', '-r300')
    close(fig)
    fprintf('Saved NAD_peak_vs_peak_%s.png\n', bpStr)
end

%% Plot saved results
% Load both sessions, create a 2x4 figure (left=s1, right=s2), save as PNG

for iSubj = 1:length(subjects)
    for bp = {'brain','calf'}
        bpStr = bp{1};
        figName = sprintf('%s_%s', subjects{iSubj}, bpStr);

        fig = figure('Position', [50 100 2000 900], 'Visible', 'off');
        sgtitle(strrep(figName, '_', ' '))

        for iSess = 1:2
            matName = sprintf('%s_%s_s%d', subjects{iSubj}, bpStr, iSess);
            matPath = [dataRoot filesep 'PROCESSED' filesep matName '.mat'];
            if ~isfile(matPath), continue; end
            load(matPath, 'output', 'pars');

            col = (iSess - 1) * 2;  % column offset: 0 for s1, 2 for s2

            % Row 1, Col 1/3: full metabolite spectrum
            subplot(2, 4, col + 1)
            plot(output.met.ppm, real(output.met.spec), 'k')
            set(gca, 'XDir', 'reverse')
            xlabel('ppm'), title(sprintf('Full spectrum (s%d)', iSess))

            % Row 1, Col 2/4: phased full spectrum
            subplot(2, 4, col + 2)
            plot(output.met.ppm, real(output.met.fit.spec_full), 'k')
            set(gca, 'XDir', 'reverse')
            xlabel('ppm'), title(sprintf('Phased spectrum (s%d)', iSess))

            % Row 2, Col 1/3: fit over full fitting range
            subplot(2, 4, col + 5)
            x = output.met.fit.ppm;
            plot(x, real(output.met.fit.spec), 'k'), hold on
            plot(x, output.met.fit.spec_fit, 'r', 'LineWidth', 1.5)
            plot(x, output.met.fit.baseline, 'b--')
            plot(x, real(output.met.fit.spec) - output.met.fit.spec_fit, 'g')
            legend({'data','fit','baseline','residual'}, 'Location', 'best', 'FontSize', 7)
            set(gca, 'XDir', 'reverse')
            xlabel('ppm'), title(sprintf('VARPRO fit full (s%d)', iSess))

            % Row 2, Col 2/4: NAD+ region zoom
            subplot(2, 4, col + 6)
            x = output.met.fit.ppm;
            zoomInds = x > 8.8 & x < 9.7;
            xz = x(zoomInds);
            plot(xz, real(output.met.fit.spec(zoomInds)), 'k'), hold on
            plot(xz, output.met.fit.spec_fit(zoomInds), 'r', 'LineWidth', 1.5)
            plot(xz, output.met.fit.baseline(zoomInds), 'b--')
            plot(xz, real(output.met.fit.spec(zoomInds)) - output.met.fit.spec_fit(zoomInds), 'g')
            legend({'data','fit','baseline','residual'}, 'Location', 'best', 'FontSize', 7)
            set(gca, 'XDir', 'reverse')
            xlabel('ppm'), title(sprintf('NAD+ region (s%d)', iSess))
        end

        print(fig, [dataRoot filesep 'PROCESSED' filesep figName '.png'], '-dpng', '-r300')
        close(fig)
        fprintf('Saved %s.png\n', figName)
    end
end
