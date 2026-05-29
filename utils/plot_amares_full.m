%% Plot AMARES fit results — full spectrum
% Run after any AMARES stage to visualise all fitted peaks.
%
% Usage:
%   plot_amares_full(Results)

function plot_amares_full(Results)

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

%% Auto-group compounds into molecules
% Known multi-compound molecules; everything else is its own molecule
groups = {
    'NAD+',        {'NADp_A', 'NADp_B'},       [0.0 0.5 0.8]    % blue
    'UDP-sugars',  {'UDPS_a', 'UDPS_b'},        [0.2 0.7 0.3]    % green
    'UDP-Gal',     {'UDPGal_a', 'UDPGal_b'},    [0.2 0.7 0.3]    % green
    'UDP-Glc',     {'UDPGlc_a', 'UDPGlc_b'},    [0.9 0.6 0.0]    % orange
    'UDP-GalNAc',  {'UDPGaNAc_a', 'UDPGaNAc_b'},[0.6 0.2 0.8]    % purple
    'UDP-GlcNAc',  {'UDPGcNAc_a', 'UDPGcNAc_b'},[0.8 0.4 0.6]    % pink
};

molecules = struct('name', {}, 'compoundIdx', {}, 'color', {});
assigned = false(nCompounds, 1);

% Try each known group
for g = 1:size(groups, 1)
    idx = [];
    for k = 1:numel(groups{g, 2})
        found = find(strcmp(compoundNames, groups{g, 2}{k}));
        idx = [idx, found(:)']; %#ok<AGROW>
    end
    if ~isempty(idx)
        molecules(end+1).name = groups{g, 1}; %#ok<AGROW>
        molecules(end).compoundIdx = idx;
        molecules(end).color = groups{g, 3};
        assigned(idx) = true;
    end
end

% Display names for known singlet compounds
displayNames = containers.Map(...
    {'PCr','ATPb','ATPa','ATPg','Pi','PE','PC','GPE','GPC','BroadMP','DPG','NADH','Broad10_5','Broad8_5'}, ...
    {'PCr','ATP-\beta','ATP-\alpha','ATP-\gamma','Pi','PE','PC','GPE','GPC','BroadMP','DPG','NADH','Broad10_5','Broad8_5'});

% Assign colours by spectral region for unassigned compounds
defaultColors = lines(nCompounds);
for i = 1:nCompounds
    if ~assigned(i)
        molecules(end+1).name = compoundNames{i}; %#ok<AGROW>
        if displayNames.isKey(compoundNames{i})
            molecules(end).name = displayNames(compoundNames{i});
        end
        molecules(end).compoundIdx = i;
        molecules(end).color = defaultColors(i, :);
    end
end

%% Compute spectra — no phase correction, no additional apodization
ppmAxis = fs.exptParams.ppmAxis;

dataSpec  = specFft(fs.inputFid, 1);
totalFit  = specFft(modelFid);
residual  = real(dataSpec) - real(totalFit);

%% Sum per-line spectra into per-molecule spectra
nMol = numel(molecules);
molSpec = zeros(fs.exptParams.samples, nMol);
for m = 1:nMol
    for c = molecules(m).compoundIdx(:)'
        lns = find(lineToCompound == c);
        for k = lns(:)'
            molSpec(:, m) = molSpec(:, m) + real(specFft(modelFids(:, k)));
        end
    end
end

%% Plot — full spectrum
ppmRange = [floor(min(ppmAxis)) ceil(max(ppmAxis))];
if max(ppmAxis) > 15 && min(ppmAxis) < -20
    ppmRange = [-25 15];  % standard 31P range
end
inRange = ppmAxis >= min(ppmRange) & ppmAxis <= max(ppmRange);

figure('Name', 'AMARES Fit: Full Spectrum', 'Position', [100 100 1100 750])

% --- Subplot 1: Data + total fit ---
ax1 = subplot(3,1,1);
plot(ppmAxis(inRange), real(dataSpec(inRange)), 'k', 'LineWidth', 1)
hold on
plot(ppmAxis(inRange), real(totalFit(inRange)), 'r', 'LineWidth', 1.5)
set(gca, 'XDir', 'reverse')
legend('Data', 'Total fit', 'Location', 'best')
ylabel('Signal (a.u.)')
title('AMARES Fit — Full Spectrum')
box off

% --- Subplot 2: Individual molecule contributions ---
ax2 = subplot(3,1,2);
plot(ppmAxis(inRange), real(dataSpec(inRange)), 'k', 'LineWidth', 0.5)
hold on
legendEntries = {'Data'};
for m = 1:nMol
    plot(ppmAxis(inRange), molSpec(inRange, m), ...
        'Color', molecules(m).color, 'LineWidth', 1.5)
    legendEntries{end+1} = molecules(m).name; %#ok<AGROW>
end
set(gca, 'XDir', 'reverse')
legend(legendEntries, 'Location', 'bestoutside', 'FontSize', 7)
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

%% Print fitted parameters — all compounds
fprintf('\n%-20s %12s %12s %10s %10s   %s\n', ...
    'Compound', 'TotalAmpl', 'CRB(Ampl)', 'LW (Hz)', 'Phase', 'Line positions (ppm)')
fprintf('%s\n', repmat('-', 1, 95))

for i = 1:nCompounds
    lns = find(lineToCompound == i);
    totalAmp = sum(Results.Amplitudes(lns));
    totalCRB = sqrt(sum(Results.Standard_deviation_of_Amplitudes(lns).^2));
    lw = Results.Linewidths(lns(1));
    ph = Results.Phases(lns(1));

    fprintf('%-20s %12.6f %12.6f %10.1f %10.1f   ', ...
        compoundNames{i}, totalAmp, totalCRB, lw, ph)
    fprintf('%.3f ', Results.ChemicalShifts(lns))
    fprintf('\n')
end

%% Print molecule-level summary (grouped)
fprintf('\n%-20s %12s %12s %10s   %s\n', ...
    'Molecule', 'TotalAmpl', 'CRB(Ampl)', 'LW (Hz)', 'Line positions (ppm)')
fprintf('%s\n', repmat('-', 1, 85))

for m = 1:nMol
    molAmp = 0;
    molCRBsq = 0;
    molLW = [];
    molCS = [];

    for c = molecules(m).compoundIdx(:)'
        lns = find(lineToCompound == c);
        molAmp = molAmp + sum(Results.Amplitudes(lns));
        molCRBsq = molCRBsq + sum(Results.Standard_deviation_of_Amplitudes(lns).^2);
        molLW(end+1) = Results.Linewidths(lns(1)); %#ok<AGROW>
        molCS = [molCS Results.ChemicalShifts(lns)']; %#ok<AGROW>
    end

    fprintf('%-20s %12.6f %12.6f %10.1f   ', ...
        molecules(m).name, molAmp, sqrt(molCRBsq), molLW(1))
    fprintf('%.3f ', molCS)
    fprintf('\n')
end
fprintf('\n')
