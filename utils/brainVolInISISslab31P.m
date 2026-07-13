function [V_M_mL, info] = brainVolInISISslab31P(brainMaskNii, t1Nii, isisDcm, scratchDir)
% BRAINVOLINISISSLAB31P  Brain volume within the ISIS spectroscopy slab.
%
%   [V_M_mL, info] = brainVolInISISslab31P(brainMaskNii, t1Nii, isisDcm, scratchDir)
%
% Computes the brain-tissue volume inside the excited 1D-ISIS slab, entirely
% in the T1 anatomical grid (no cross-frame resampling):
%   - imscribe_engine rasterizes the ISIS SVS DICOM VOI into the t1w grid,
%   - the HD-BET brain mask (skull-stripped T1, brain>0) is intersected with it.
%
% This was cross-validated against a direct B1-slab resample: the two agreed
% to <0.1% on the S095 test dataset.  Requires the CAMIPM toolbox (imscribe).
%
% INPUTS
%   brainMaskNii  HD-BET output (e.g. *_hdbet.nii.gz or *_bet.nii.gz); brain>0
%   t1Nii         the T1 NIfTI the mask was derived from (same grid)
%   isisDcm       one ISIS_SVS DICOM file (defines the VOI geometry)
%   scratchDir    (optional) writable dir for the imscribe svsmask; default tempdir
%
% OUTPUTS
%   V_M_mL  brain-in-slab volume (mL)
%   info    struct: V_slab_mL, V_brain_total_mL, voxVol_mm3, svsmaskPath

if nargin < 4 || isempty(scratchDir), scratchDir = tempname; end
if ~isfolder(scratchDir), mkdir(scratchDir); end

assert(isfile(brainMaskNii), 'brainVolInISISslab31P: brain mask not found: %s', brainMaskNii);
assert(isfile(t1Nii),        'brainVolInISISslab31P: t1 NIfTI not found: %s', t1Nii);
assert(isfile(isisDcm),      'brainVolInISISslab31P: ISIS DICOM not found: %s', isisDcm);
assert(~isempty(which('imscribe_engine')), ...
    'brainVolInISISslab31P: imscribe_engine not on path (add CAMIPM toolbox).');

t1info = niftiinfo(t1Nii);
voxVol = prod(t1info.PixelDimensions);          % mm^3
brain  = niftiread(brainMaskNii) > 0;

opt = struct('datatypes', {{'NIFTI','DICOM',''}}, ...
             'scratchpath', [scratchDir filesep], 'svsmask', 1);
imscribe_engine({t1Nii}, {isisDcm}, {[]}, opt);

svsmaskPath = fullfile(scratchDir,'svsmask.nii');
assert(isfile(svsmaskPath), 'brainVolInISISslab31P: imscribe did not write svsmask.nii');
slab = niftiread(svsmaskPath) > 0;
assert(isequal(size(slab), size(brain)), 'slab/brain grid mismatch');

V_M_mL = nnz(brain & slab) * voxVol / 1000;
info = struct('V_slab_mL', nnz(slab)*voxVol/1000, ...
              'V_brain_total_mL', nnz(brain)*voxVol/1000, ...
              'voxVol_mm3', voxVol, 'svsmaskPath', svsmaskPath);
end
