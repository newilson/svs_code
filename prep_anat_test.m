function prep_anat_test(sessionDir, outDir, force)
% PREP_ANAT_TEST  Stage 1 of testBrainSeg pipeline.
%   For each GRE3D* and MPRAGE* DICOM series under sessionDir
%   (excluding MPRAGE_P0*, which has no matching svsmask):
%     1) DICOM -> NIfTI (dicom2nifti) into outDir
%     2) SPM bias correction (me_spm_biascor) -> <stem>_SPM8bias.nii
%     3) Match the bias-corrected scan's grid (dims + sform) against the
%        svsmask*.nii files in outDir; pick the one whose grid matches
%        exactly. Skip the scan if no mask matches.
%   Writes outDir/manifest.csv with columns scanName, biasNii, svsMask.

if nargin < 3 || isempty(force), force = false; end
if ~exist(outDir,'dir'), mkdir(outDir); end

% --- inventory candidate masks and cache their grids ---
maskFiles = dir(fullfile(outDir,'svsmask*.nii'));
maskFiles = arrayfun(@(f) fullfile(f.folder,f.name), maskFiles, 'UniformOutput', false);
assert(~isempty(maskFiles), 'No svsmask*.nii in %s', outDir);
maskGrids = cellfun(@gridSig, maskFiles, 'UniformOutput', false);

dirs = [dir(fullfile(sessionDir,'GRE3D*')); dir(fullfile(sessionDir,'MPRAGE*'))];
dirs = dirs([dirs.isdir]);
% --- skip MPRAGE_P0* (no matching mask) ---
keep = ~startsWith({dirs.name}, 'MPRAGE_P0');
if any(~keep)
    fprintf('Skipping (no matching mask): %s\n', strjoin({dirs(~keep).name}, ', '));
end
dirs = dirs(keep);
assert(~isempty(dirs), 'No GRE3D*/MPRAGE* series under %s after filtering', sessionDir);

rows = cell(0,3);
for k = 1:numel(dirs)
    d   = dirs(k);
    raw = fullfile(d.folder, d.name);

    % Strip Siemens trailing series number (e.g. MPRAGE_P2_1MM_0002 -> MPRAGE_P2_1MM)
    stem = regexprep(d.name, '_\d{4}$', '');

    niiFile  = fullfile(outDir, [stem '.nii']);
    biasFile = fullfile(outDir, [stem '_SPM8bias.nii']);

    % --- DICOM -> NIfTI ---
    if force || ~isfile(niiFile)
        files = dir(fullfile(raw,'*'));
        files = files(~[files.isdir]);
        files = arrayfun(@(f) fullfile(f.folder,f.name), files, 'UniformOutput', false);
        assert(~isempty(files), 'No DICOM files in %s', raw);
        fprintf('[%d/%d] dicom2nifti: %s\n', k, numel(dirs), stem);
        dicom2nifti(files, niiFile);
        assert(isfile(niiFile), 'dicom2nifti did not produce %s', niiFile);
    else
        fprintf('[%d/%d] skip dicom2nifti (exists): %s\n', k, numel(dirs), stem);
    end

    % --- SPM bias correction ---
    %   me_spm_biascor writes <inpath>/<outroot>_SPM8bias.nii. Pass bare
    %   stem so output lands next to niiFile (which is in outDir).
    if force || ~isfile(biasFile)
        fprintf('[%d/%d] me_spm_biascor: %s\n', k, numel(dirs), stem);
        me_spm_biascor(niiFile, stem);
        assert(isfile(biasFile), 'me_spm_biascor did not produce %s', biasFile);
    else
        fprintf('[%d/%d] skip bias (exists): %s\n', k, numel(dirs), stem);
    end

    % --- Match svsmask by exact grid (dims + sform) ---
    g = gridSig(biasFile);
    matchIdx = find(cellfun(@(m) gridsEqual(g,m), maskGrids), 1);
    if isempty(matchIdx)
        fprintf('         no svsmask matches grid (dim=[%d %d %d] vox=[%.3f %.3f %.3f]); skipping %s\n', ...
            g.dim(1), g.dim(2), g.dim(3), g.vox(1), g.vox(2), g.vox(3), stem);
        continue
    end
    sel = maskFiles{matchIdx};
    [~,maskName] = fileparts(sel);
    fprintf('         grid dim=[%d %d %d] vox=[%.3f %.3f %.3f] -> %s\n', ...
        g.dim(1), g.dim(2), g.dim(3), g.vox(1), g.vox(2), g.vox(3), maskName);

    rows(end+1,:) = {stem, biasFile, sel}; %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', {'scanName','biasNii','svsMask'});
manifestPath = fullfile(outDir,'manifest.csv');
writetable(T, manifestPath);
fprintf('Wrote %s (%d scans)\n', manifestPath, height(T));

end

% ----- helpers -----
function g = gridSig(niiFile)
% Return a struct holding dims + sform for grid comparison.
info  = niftiinfo(niiFile);
g.dim = info.ImageSize(1:3);
g.vox = info.PixelDimensions(1:3);
g.aff = info.Transform.T;   % 4x4, MATLAB convention (row-major)
end

function tf = gridsEqual(a, b)
% Same dims and same affine within float tolerance.
tf = isequal(a.dim, b.dim) && norm(a.aff(:) - b.aff(:), Inf) < 1e-3;
end
