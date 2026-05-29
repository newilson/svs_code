%% Plot AMARES fit results in the NAD/UDP region
% Run after dat2mat_svs31p has produced 'Results' (the Stage 2 output)
% and 'amaresStruct' in the workspace.
%
% Usage:
%   plot_amares_nadudp          % uses defaults
%   plot_amares_nadudp(Results) % pass Results explicitly

function plot_amares_nadudp(Results)

%% Extract fit internals
fs = Results.fitStatus;

% Reconstruct per-line FIDs (Npts x Nlines)
[modelFid, ~, modelFids] = AMARES.makeModelFidAndJacobianReIm( ...
    fs.xFit, fs.constraintsCellArray, ...
    fs.exptParams.beginTime, fs.exptParams.dwellTime, ...
    fs.exptParams.imagingFrequency, fs.exptParams.samples, ...
    'complexOutput', true);

%% Build line-to-compound mapping from prior knowledge
pk = fs.pkWithLinLsq;
nCompounds = numel(pk.initialValues);
compoundNames = cell(nCompounds, 1);
lineToCompound = [];  % which compound each spectral line belongs to

for i = 1:nCompounds
    pn = pk.initialValues(i).peakName;
    if iscell(pn)
        nLines = numel(pn);
        compoundNames{i} = pn{1}(1:end-1);  % strip trailing digit
    else
        nLines = 1;
        compoundNames{i} = pn;
    end
    lineToCompound = [lineToCompound; repmat(i, nLines, 1)]; %#ok<AGROW>
end

%% Group compounds into molecules for plotting
% Auto-detect whether we have combined UDPS or four separate UDP species
molecules = struct('name', {}, 'compoundIdx', {}, 'color', {});

molecules(end+1).name = 'NADH';
molecules(end).compoundIdx = find(strcmp(compoundNames, 'NADH'));
molecules(end).color = [0.8 0.2 0.2];  % red

molecules(end+1).name = 'NAD+';
molecules(end).compoundIdx = [find(strcmp(compoundNames, 'NADp_A')), ...
                               find(strcmp(compoundNames, 'NADp_B'))];
molecules(end).color = [0.0 0.5 0.8];  % blue

% Check for individual UDP species vs combined UDPS
if ~isempty(find(strcmp(compoundNames, 'UDPGal_a'), 1))
    % Four separate UDP species
    molecules(end+1).name = 'UDP-Gal';
    molecules(end).compoundIdx = [find(strcmp(compoundNames, 'UDPGal_a')), ...
                                   find(strcmp(compoundNames, 'UDPGal_b'))];
    molecules(end).color = [0.2 0.7 0.3];  % green

    molecules(end+1).name = 'UDP-Glc';
    molecules(end).compoundIdx = [find(strcmp(compoundNames, 'UDPGlc_a')), ...
                                   find(strcmp(compoundNames, 'UDPGlc_b'))];
    molecules(end).color = [0.9 0.6 0.0];  % orange

    molecules(end+1).name = 'UDP-GalNAc';
    molecules(end).compoundIdx = [find(strcmp(compoundNames, 'UDPGaNAc_a')), ...
                                   find(strcmp(compoundNames, 'UDPGaNAc_b'))];
    molecules(end).color = [0.6 0.2 0.8];  % purple

    molecules(end+1).name = 'UDP-GlcNAc';
    molecules(end).compoundIdx = [find(strcmp(compoundNames, 'UDPGcNAc_a')), ...
                                   find(strcmp(compoundNames, 'UDPGcNAc_b'))];
    molecules(end).color = [0.8 0.4 0.6];  % pink

elseif ~isempty(find(strcmp(compoundNames, 'UDPS_a'), 1))
    % Combined UDP species
    molecules(end+1).name = 'UDP-sugars';
    molecules(end).compoundIdx = [find(strcmp(compoundNames, 'UDPS_a')), ...
                                   find(strcmp(compoundNames, 'UDPS_b'))];
    molecules(end).color = [0.2 0.7 0.3];  % green
end

% Include broad components if present
broadIdx10 = find(strcmp(compoundNames, 'Broad10_5'));
if ~isempty(broadIdx10)
    molecules(end+1).name = 'Broad10_5';
    molecules(end).compoundIdx = broadIdx10;
    molecules(end).color = [0.7 0.7 0.7];  % light grey
end
broadIdx8 = find(strcmp(compoundNames, 'Broad8_5'));
if ~isempty(broadIdx8)
    molecules(end+1).name = 'Broad8_5';
    molecules(end).compoundIdx = broadIdx8;
    molecules(end).color = [0.55 0.55 0.55];  % grey
end

% Also include ATP-alpha since it overlaps this region
molecules(end+1).name = 'ATP-\alpha';
molecules(end).compoundIdx = find(strcmp(compoundNames, 'ATPa'));
molecules(end).color = [0.5 0.5 0.5];  % dark grey

%% Compute spectra
ppmAxis = fs.exptParams.ppmAxis;

% No phase correction or additional apodization — plot on original data
dataSpec  = specFft(fs.inputFid, 1);
totalFit  = specFft(modelFid);
residual  = real(dataSpec) - real(totalFit);

%% Sum per-line spectra into per-molecule spectra
molSpec = zeros(fs.exptParams.samples, numel(molecules));
for m = 1:numel(molecules)
    for c = molecules(m).compoundIdx(:)'
        lines = find(lineToCompound == c);
        for k = lines(:)'
            molSpec(:, m) = molSpec(:, m) + ...
                real(specFft(modelFids(:,k)));
        end
    end
end

%% Plot — NAD/UDP region
ppmRange = [-11.5 -6.5];
inRange = ppmAxis >= min(ppmRange) & ppmAxis <= max(ppmRange);

figure('Name', 'AMARES Fit: NAD/UDP Region', 'Position', [100 100 900 700])

% --- Subplot 1: Data + total fit + individual molecules ---
ax1 = subplot(3,1,1);
plot(ppmAxis(inRange), real(dataSpec(inRange)), 'k', 'LineWidth', 1)
hold on
plot(ppmAxis(inRange), real(totalFit(inRange)), 'r', 'LineWidth', 1.5)
set(gca, 'XDir', 'reverse')
legend('Data', 'Total fit', 'Location', 'best')
ylabel('Signal (a.u.)')
title('NAD / UDP Region')
box off

% --- Subplot 2: Individual molecule contributions ---
ax2 = subplot(3,1,2);
plot(ppmAxis(inRange), real(dataSpec(inRange)), 'k', 'LineWidth', 0.5)
hold on
legendEntries = {'Data'};
for m = 1:numel(molecules)
    plot(ppmAxis(inRange), molSpec(inRange, m), ...
        'Color', molecules(m).color, 'LineWidth', 1.5)
    legendEntries{end+1} = molecules(m).name; %#ok<AGROW>
end
set(gca, 'XDir', 'reverse')
legend(legendEntries, 'Location', 'bestoutside', 'FontSize', 8)
ylabel('Signal (a.u.)')
box off

% --- Subplot 3: Residual ---
ax3 = subplot(3,1,3);
plot(ppmAxis(inRange), residual(inRange), 'k', 'LineWidth', 0.5)
hold on
plot(ppmAxis(inRange), zeros(sum(inRange),1), 'r--')
noiseStd = std(residual(inRange));
plot(ppmAxis(inRange), +noiseStd * ones(sum(inRange),1), 'r:')
plot(ppmAxis(inRange), -noiseStd * ones(sum(inRange),1), 'r:')
set(gca, 'XDir', 'reverse')
ylabel('Residual')
xlabel('\delta / ppm')
box off

linkaxes([ax1 ax2 ax3], 'x')
xlim(ax1, sort(ppmRange, 'ascend'))

%% Print fitted parameters — ALL compounds
fprintf('\n%-20s %12s %12s %10s %10s   %s\n', ...
    'Compound', 'TotalAmpl', 'CRB(Ampl)', 'LW (Hz)', 'Phase', 'Line positions (ppm)')
fprintf('%s\n', repmat('-', 1, 95))

for i = 1:nCompounds
    lines = find(lineToCompound == i);
    totalAmp = sum(Results.Amplitudes(lines));
    totalCRB = sqrt(sum(Results.Standard_deviation_of_Amplitudes(lines).^2));
    lw = Results.Linewidths(lines(1));
    ph = Results.Phases(lines(1));

    fprintf('%-20s %12.6f %12.6f %10.1f %10.1f   ', ...
        compoundNames{i}, totalAmp, totalCRB, lw, ph)
    fprintf('%.3f ', Results.ChemicalShifts(lines))
    fprintf('\n')
end

%% Print molecule-level summary (grouped)
fprintf('\n%-20s %12s %12s %10s   %s\n', ...
    'Molecule', 'TotalAmpl', 'CRB(Ampl)', 'LW (Hz)', 'Line positions (ppm)')
fprintf('%s\n', repmat('-', 1, 85))

for m = 1:numel(molecules)
    molAmp = 0;
    molCRBsq = 0;
    molLW = [];
    molCS = [];

    for c = molecules(m).compoundIdx(:)'
        lines = find(lineToCompound == c);
        molAmp = molAmp + sum(Results.Amplitudes(lines));
        molCRBsq = molCRBsq + sum(Results.Standard_deviation_of_Amplitudes(lines).^2);
        molLW(end+1) = Results.Linewidths(lines(1)); %#ok<AGROW>
        molCS = [molCS Results.ChemicalShifts(lines)']; %#ok<AGROW>
    end

    fprintf('%-20s %12.6f %12.6f %10.1f   ', ...
        molecules(m).name, molAmp, sqrt(molCRBsq), molLW(1))
    fprintf('%.3f ', molCS)
    fprintf('\n')
end
fprintf('\n')
