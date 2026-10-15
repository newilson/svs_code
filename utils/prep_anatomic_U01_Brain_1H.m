function prep_anatomic_u01(sessionDir, outDir, opts)
% U01_PREP_ANAT  Stage 1 of U01 SVS pipeline.
%   Locates T1w DICOMs and a 1H SVS DICOM under sessionDir on Z:\,
%   converts T1w to NIfTI, SPM-bias-corrects it, and runs imscribe to
%   produce svsmask.nii. All outputs land in outDir (the staging dir
%   inside the VM-shared folder).
%
%   prep_anatomic_U01_Brain_1H(sessionDir, outDir, 'ExamPattern', p, ...
%       'SvsPattern', q, 'MaskName', m)
%     ExamPattern  exam folder under sessionDir     (default 'E*')
%     SvsPattern   SVS DICOM series in the exam     (default 'SVSSEL*NAD*')
%     MaskName     mask file written to outDir      (default 'svsmask.nii')
%
%   Each stage is skipped if its output already exists in outDir, so
%   reruns only recompute what's missing.

arguments
    sessionDir       (1,:) char
    outDir           (1,:) char
    opts.ExamPattern (1,:) char = 'E*'
    opts.SvsPattern  (1,:) char = 'SVSSEL*NAD*'
    opts.MaskName    (1,:) char = 'svsmask.nii'
end

if ~exist(outDir,'dir'), mkdir(outDir); end

% ---- Expected output paths -----------------------------------------
t1Nii      = fullfile(outDir, 't1w.nii');
[~, t1Stem] = fileparts(t1Nii);
biasNii8   = fullfile(outDir, [t1Stem '_SPM8bias.nii']);
biasNii12  = fullfile(outDir, [t1Stem '_SPM12bias.nii']);
svsmaskNii = fullfile(outDir, opts.MaskName);

needT1   = ~isfile(t1Nii);
needBias = ~(isfile(biasNii8) || isfile(biasNii12));
needMask = ~isfile(svsmaskNii);

if ~(needT1 || needBias || needMask)
    fprintf('prep_anatomic: all outputs already present in %s -- skipping.\n', outDir);
    return
end

% ---- Locate exam dir (needed for any DICOM step) -------------------
if needT1 || needMask
    examDirs = dir(fullfile(sessionDir,opts.ExamPattern));
    examDirs = examDirs([examDirs.isdir]);
    assert(numel(examDirs) >= 1, 'No %s exam folder under %s', opts.ExamPattern, sessionDir);
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
    svsDirs = dir(fullfile(examDir,opts.SvsPattern));
    svsDirs = svsDirs([svsDirs.isdir]);
    assert(~isempty(svsDirs), 'No %s folder in %s', opts.SvsPattern, examDir);
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

    engineMask = fullfile(outDir, 'svsmask.nii');
    assert(isfile(engineMask), ...
        'imscribe_engine did not produce svsmask.nii in %s', outDir);
    if ~strcmp(engineMask, svsmaskNii)
        movefile(engineMask, svsmaskNii);
    end
else
    fprintf('prep_anatomic: %s exists -- skipping imscribe\n', svsmaskNii);
end

end
