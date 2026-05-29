function prep_anatomic_u01(sessionDir, outDir)
% U01_PREP_ANAT  Stage 1 of U01 SVS pipeline.
%   Locates T1w DICOMs and a 1H SVS DICOM under sessionDir on Z:\,
%   converts T1w to NIfTI, SPM-bias-corrects it, and runs imscribe to
%   produce svsmask.nii. All outputs land in outDir (the staging dir
%   inside the VM-shared folder).
%
%   Each stage is skipped if its output already exists in outDir, so
%   reruns only recompute what's missing.

if ~exist(outDir,'dir'), mkdir(outDir); end

% ---- Expected output paths -----------------------------------------
t1Nii      = fullfile(outDir, 't1w.nii');
[~, t1Stem] = fileparts(t1Nii);
biasNii8   = fullfile(outDir, [t1Stem '_SPM8bias.nii']);
biasNii12  = fullfile(outDir, [t1Stem '_SPM12bias.nii']);
svsmaskNii = fullfile(outDir, 'svsmask.nii');

needT1   = ~isfile(t1Nii);
needBias = ~(isfile(biasNii8) || isfile(biasNii12));
needMask = ~isfile(svsmaskNii);

if ~(needT1 || needBias || needMask)
    fprintf('prep_anatomic: all outputs already present in %s -- skipping.\n', outDir);
    return
end

% ---- Locate exam dir (needed for any DICOM step) -------------------
if needT1 || needMask
    examDirs = dir(fullfile(sessionDir,'E*'));
    examDirs = examDirs([examDirs.isdir]);
    assert(numel(examDirs) >= 1, 'No E* exam folder under %s', sessionDir);
    examDir  = fullfile(examDirs(1).folder, examDirs(1).name);
end

% ---- 1) DICOM -> NIfTI ---------------------------------------------
if needT1
    % T1w (MPRAGE_*) -- assume one such series per exam
    t1Dirs = dir(fullfile(examDir,'MPRAGE*'));
    t1Dirs = t1Dirs([t1Dirs.isdir]);
    assert(~isempty(t1Dirs), 'No MPRAGE* folder in %s', examDir);
    t1Dir  = fullfile(t1Dirs(1).folder, t1Dirs(1).name);

    t1Files = dir(fullfile(t1Dir,'*'));
    t1Files = t1Files(~[t1Files.isdir]);
    t1Files = arrayfun(@(d) fullfile(d.folder,d.name), t1Files, 'UniformOutput', false);
    assert(~isempty(t1Files), 'No DICOM files in %s', t1Dir);

    [~,~,t1Nii] = dicom2nifti(t1Files, t1Nii);
    assert(isfile(t1Nii), 'dicom2nifti did not produce %s', t1Nii);
else
    fprintf('prep_anatomic: %s exists -- skipping DICOM->NIfTI\n', t1Nii);
end

% ---- 2) SPM bias correction ----------------------------------------
% me_spm_biascor writes <inpath>/<stem>_<SPMver>bias.nii (inpath is
% taken from the input file's folder; second arg is just the stem).
if needBias
    [~, stem] = fileparts(t1Nii);
    me_spm_biascor(t1Nii, stem);
else
    existing = biasNii8;
    if ~isfile(existing), existing = biasNii12; end
    fprintf('prep_anatomic: %s exists -- skipping bias correction\n', existing);
end

% ---- 3) imscribe SVS mask -----------------------------------------
% Mirrors imscribe_svsmask_batchdemo.m: NIfTI structural + DICOM SVS,
% scratchpath = outDir (where 'svsmask.nii' will be written).
if needMask
    % 1H SVS NAD series -- used to define the voxel mask
    svsDirs = dir(fullfile(examDir,'SVSSEL*NAD*'));
    svsDirs = svsDirs([svsDirs.isdir]);
    assert(~isempty(svsDirs), 'No SVSSEL*NAD* folder in %s', examDir);
    svsDir  = fullfile(svsDirs(1).folder, svsDirs(1).name);

    svsFiles = dir(fullfile(svsDir,'*'));
    svsFiles = svsFiles(~[svsFiles.isdir]);
    assert(~isempty(svsFiles), 'No DICOM files in %s', svsDir);
    svsDcm   = fullfile(svsFiles(1).folder, svsFiles(1).name);

    opt = struct();
    opt.datatypes   = {'NIFTI','DICOM',''};
    opt.scratchpath = [outDir filesep];
    opt.svsmask     = 1;
    imscribe_engine({t1Nii}, {svsDcm}, {[]}, opt);

    assert(isfile(svsmaskNii), ...
        'imscribe_engine did not produce svsmask.nii in %s', outDir);
else
    fprintf('prep_anatomic: %s exists -- skipping imscribe\n', svsmaskNii);
end

end
