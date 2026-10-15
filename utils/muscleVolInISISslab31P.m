function [V_M_mL, info] = muscleVolInISISslab31P(dixonMagDir, isisDcm, opt)
% MUSCLEVOLINISISSLAB31P  Muscle volume within the ISIS spectroscopy slab.
%
%   [V_M_mL, info] = muscleVolInISISslab31P(dixonMagDir, isisDcm, opt)
%
% The calf analog of brainVolInISISslab31P, and it ends by calling that
% function -- the slab intersection is identical once a mask and a grid
% reference exist.  What is calf-specific is producing them:
%
%   1) 3-point Dixon water/fat separation (readB0FWgre, Z:\code\dixonProcessing)
%   2) muscle = signal-bearing voxels with fat fraction below ffThresh, so
%      subcutaneous fat and marrow -- which carry no PCr or ATP -- are excluded
%   3) both written as NIfTI on the Dixon grid
%
% This replaces an intensity threshold on the B1 magnitude image, which had
% three defects: it included fat and marrow, it ran at 20 mm slice thickness,
% and it covered only 20-40 mm of a 60 mm ISIS slab so the volume had to be
% extrapolated.  The Dixon is 1.25 x 1.25 x 3 mm over ~165 mm and covers the
% slab outright.
%
% INPUTS
%   dixonMagDir  GRE_*DIXON* MAGNITUDE series dir.  The phase series is found
%                by readB0FWgre incrementing the trailing _NNNN.
%   isisDcm      one ISIS_SVS DICOM (defines the VOI geometry)
%   opt          optional struct:
%                  .ffThresh   fat-fraction cutoff for muscle   (default 0.50)
%                  .sigFrac    body-mask threshold, frac of max (default 0.10)
%                  .dixonOpt   readB0FWgre algorithm            (default 3,
%                              Berglund 2017 native MATLAB graph-cut)
%                  .scratchDir writable dir                     (default tempname)
%
% OUTPUTS
%   V_M_mL  muscle-in-slab volume (mL)
%   info    struct with the intermediate volumes and the QC fields below
%
% ORIENTATION.  readB0FWgre returns arrays in dicomread layout; dicom2nifti
% writes the transpose-and-flip of that.  Verified empirically to machine
% precision on S139_06052026:
%       dicomread_stack == flip(permute(niftiread, [2 1 3]), 1)
% so the inverse, applied here to get Dixon output onto the NIfTI grid, is
%       nii = permute(flip(A, 1), [2 1 3])
% Do not "simplify" this without re-checking against the DICOMs.

if nargin < 3 || isempty(opt), opt = struct(); end
if ~isfield(opt,'ffThresh')   || isempty(opt.ffThresh),   opt.ffThresh   = 0.50;   end
if ~isfield(opt,'sigFrac')    || isempty(opt.sigFrac),    opt.sigFrac    = 0.10;   end
if ~isfield(opt,'dixonOpt')   || isempty(opt.dixonOpt),   opt.dixonOpt   = 3;      end
if ~isfield(opt,'scratchDir') || isempty(opt.scratchDir), opt.scratchDir = tempname; end

assert(isfolder(dixonMagDir), 'muscleVolInISISslab31P: Dixon dir not found: %s', dixonMagDir);
assert(isfile(isisDcm),       'muscleVolInISISslab31P: ISIS DICOM not found: %s', isisDcm);
assert(~isempty(which('readB0FWgre')), ...
    'muscleVolInISISslab31P: readB0FWgre not on path (add Z:\\code\\dixonProcessing).');
assert(~isempty(which('dicom2nifti')), ...
    'muscleVolInISISslab31P: dicom2nifti not on path (add the CAMIPM toolbox).');

sc = opt.scratchDir;
if ~isfolder(sc), mkdir(sc); end

% ---- 1) water / fat separation ----------------------------------------
[~, wat, fat] = readB0FWgre(dixonMagDir, [], opt.dixonOpt, false);
wat = abs(squeeze(wat(:,:,:,1)));
fat = abs(squeeze(fat(:,:,:,1)));

% ---- 2) muscle mask ----------------------------------------------------
S    = wat + fat;
body = S > opt.sigFrac * max(S(:));
body = imfill(keepLargestCC3(body), 'holes');       % the leg, fat included
FF   = fat ./ max(S, eps);
muscle = body & (FF < opt.ffThresh);

% ---- 3) NIfTI on the Dixon grid ---------------------------------------
% dicom2nifti on the first-echo magnitude defines the grid and the affine that
% imscribe expects; the Dixon output is mapped onto it by the inverse of the
% verified layout relation above.
e1  = firstEchoFiles(dixonMagDir, size(wat,3));
ref = fullfile(sc, 'dixon_mag.nii');
[~,~,ref] = dicom2nifti(e1, ref);
assert(isfile(ref), 'muscleVolInISISslab31P: dicom2nifti produced no %s', ref);

ni = niftiinfo(ref);
toNii = @(A) permute(flip(A, 1), [2 1 3]);
assert(isequal(size(toNii(muscle)), ni.ImageSize), ...
    'muscleVolInISISslab31P: Dixon %s does not match NIfTI grid %s', ...
    mat2str(size(toNii(muscle))), mat2str(ni.ImageSize));

maskNii = fullfile(sc, 'dixon_muscle');            % niftiwrite appends .nii
nm = ni; nm.Datatype = 'uint8'; nm.BitsPerPixel = 8;
niftiwrite(uint8(toNii(muscle)), maskNii, nm);
maskNii = [maskNii '.nii'];

% ---- 4) slab intersection (shared with the brain path) ----------------
[V_M_mL, sinfo] = brainVolInISISslab31P(maskNii, ref, isisDcm, sc);

% ---- QC ----------------------------------------------------------------
voxVol = prod(ni.PixelDimensions);
info = sinfo;
info.V_muscle_total_mL = nnz(muscle) * voxVol / 1000;
info.V_leg_total_mL    = nnz(body)   * voxVol / 1000;
info.muscleFrac        = nnz(muscle) / max(nnz(body),1);
info.ffThresh          = opt.ffThresh;
info.dixonOpt          = opt.dixonOpt;
info.dixonMagDir       = dixonMagDir;
info.maskNii           = maskNii;
info.refNii            = ref;

% The rasterized slab should be as thick as the prescription.  If the affine
% were wrong this is what would catch it, so check rather than assume.
slab = niftiread(sinfo.svsmaskPath) > 0;
zs   = find(squeeze(any(any(slab,1),2)));
info.slabExtent_mm = numel(zs) * ni.PixelDimensions(3);
if isempty(zs)
    warning('muscleVolInISISslab31P: ISIS slab does not intersect the Dixon volume at all.');
elseif info.muscleFrac < 0.2 || info.muscleFrac > 0.95
    warning(['muscleVolInISISslab31P: muscle is %.0f%% of the leg cross-section ' ...
             '(expected roughly 50-80%%); check ffThresh=%.2f.'], ...
             100*info.muscleFrac, opt.ffThresh);
end
end

% =====================================================================
function f = firstEchoFiles(d, nz)
% The first nz DICOMs in dir() order, which is the first echo: readdicomfiles2d
% consumes files in this order and readB0FWgre reshapes as [nx ny nz nte], so
% echo blocks are contiguous and slices ascend within each.  Verified on the
% U01_2026 calf data (3 echoes x 56 slices, InstanceNumber monotonic).
q = dir(d); f = {};
for k = 1:numel(q)
    a = q(k).name;
    if a(1) ~= '.' && a(1) ~= '$' && ~q(k).isdir && isdicom(fullfile(d,a))
        f{end+1} = fullfile(d, a); %#ok<AGROW>
    end
end
assert(numel(f) >= nz, 'firstEchoFiles: %d DICOMs found, need at least %d', numel(f), nz);
f = f(1:nz);
end

% =====================================================================
function out = keepLargestCC3(bw)
cc = bwconncomp(bw, 6);
if cc.NumObjects <= 1, out = bw; return; end
n = cellfun(@numel, cc.PixelIdxList);
[~, i] = max(n);
out = false(size(bw));
out(cc.PixelIdxList{i}) = true;
end
