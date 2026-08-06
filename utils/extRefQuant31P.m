function absExt = extRefQuant31P(fitNames, fitAreas, metProps, TR_ms, TE_ms, extRef, internalAbsQ, outDir, scanID)
% EXTREFQUANT31P  Absolute 31P quant from an EXTERNAL reference phantom.
%
%   absExt = extRefQuant31P(fitNames, fitAreas, metProps, TR_ms, TE_ms, ...
%                           extRef, internalAbsQ, outDir, scanID)
%
% Converts each fitted metabolite amplitude to absolute [mM] against an
% external phantom of known concentration, correcting for:
%   - reference vs metabolite signal ratio            (S_M / S_ref)
%   - phantom-vs-tissue volume within the excited slab (V_ref / V_M)
%   - T1/T2 relaxation weighting                        (R_ref / R_M)
%   - equivalent-31P count                              (n_ref / n_M)
%
%   [M] = [ref] * (S_M/S_ref) * (V_ref/V_M) * (R_ref/R_M) * (n_ref/n_M)
%
% *** B1 CORRECTION IS CURRENTLY DISABLED (extRef.applyB1 = false). ***
% The B1+ excitation-efficiency factor (eta_ref/eta_M) and the reciprocity
% receive factor (relB1_ref/relB1_M) are both computed and reported for QC
% but are NOT applied to the concentrations.  Set extRef.applyB1 = true to
% restore them; see the applyB1 block below.  Rationale: the external
% phantom is faint (D2O / PBS) so its relB1 is either masked out by
% calc_B1map_new or must be extrapolated, and neither is trustworthy yet.
%
% The isis_svs sequence is 1D-ISIS (slab select via adiabatic pulses + hard
% pulse excitation), so a relative B1+ of relB1 would give an actual flip of
% relB1*flip and excitation efficiency sin(relB1*flip).
%
% ANATOMY: extRef.anatomy selects the tissue segmentation.
%   'brain' - HD-BET brain mask INTERSECT ISIS slab when available, else an
%             intensity-threshold head mask (includes skull/scalp).
%   'calf'  - intensity-threshold leg mask (largest bright component).  NOTE
%             this includes subcutaneous fat and marrow, which carry no PCr /
%             ATP, so V_M is an OVERestimate and the external-reference calf
%             concentrations are correspondingly biased LOW.  There is no
%             HD-BET analogue for the leg.
%
% extRef fields (see extRefDefaults in processU01_Brain_31P / _Calf_31P):
%   b1Dir, refPpm, refConc_mM, nP, flipDeg, T1_ms, T2_ms, anatomy, applyB1,
%   dropEndSlices, brainThr, phantThr, phantMinVox, phantomMask (optional
%   override), reciprocityReceive, plt.
%
% fitAreas are the per-spin-normalized amplitudes (amplScaled), so n cancels
% when nP=1 for all species; the n_ref/n_M factor is kept explicit for clarity.

absExt = struct();
absExt.ok = false;

% ---------------------------------------------------------------
% Reference peak amplitude
% ---------------------------------------------------------------
iRef = find(strcmp(fitNames, 'RefExt'), 1);
if isempty(iRef)
    warning('extRefQuant31P: RefExt not in fit results; skipping external quant.');
    absExt.reason = 'RefExt not in fit results';
    return
end
areaRef = fitAreas(iRef);
absExt.areaRef = areaRef;   % recorded even if the quant below cannot run
if ~(areaRef > 0)
    warning('extRefQuant31P: RefExt amplitude <= 0 (%.3g); external quant unreliable.', areaRef);
end

% ---------------------------------------------------------------
% B1 map: relative-B1 volume + magnitude images for segmentation
% ---------------------------------------------------------------

if isempty(extRef.b1Dir) || ~isfolder(extRef.b1Dir)
    if isempty(extRef.b1Dir)
        msg = 'no B1 map series found in session';
    else
        msg = sprintf('B1 map dir not found: %s', extRef.b1Dir);
    end
    warning(['extRefQuant31P: %s. RefExt was fitted (amplitude %.4g) but cannot ' ...
             'be converted to concentrations without the phantom/tissue volumes. ' ...
             'Continuing with internal reference only.'], msg, areaRef);
    absExt.reason = msg;
    return
end

try
    [B1map, ~, b1alpha, images] = readB1(extRef.b1Dir);
catch ME
    warning(['extRefQuant31P: could not read B1 map (%s: %s). RefExt was fitted ' ...
             '(amplitude %.4g) but external quant is skipped. Continuing.'], ...
            extRef.b1Dir, ME.message, areaRef);
    absExt.reason = sprintf('readB1 failed: %s', ME.message);
    return
end

% images is [nx ny nz nFlip]; use the mean over flip angles for a stable
% anatomical magnitude for segmentation.
if ndims(images) == 4
    magVol = mean(images, 4);
else
    magVol = images;
end
nz = size(B1map, 3);

% voxel dimensions from DICOM
[px, slthk] = b1VoxelDims(extRef.b1Dir);
voxVol_mm3  = px(1) * px(2) * slthk;

% ---------------------------------------------------------------
% Slices spanning the excited slab
% ---------------------------------------------------------------
% The B1 map shares the ISIS slab orientation, so the slab is the whole B1
% stack minus the first and last slice (slice-profile edge effects).  A 2D
% (single-slice) B1 map has nothing to drop -- keep it as-is.
drop = extRef.dropEndSlices;
if nz > 2*drop
    sel = (1 + drop) : (nz - drop);
else
    sel = 1 : nz;
    if nz > 1
        warning(['extRefQuant31P: only %d B1 slice(s); cannot drop %d from each ' ...
                 'end -- using all slices.'], nz, drop);
    end
end
mag  = magVol(:,:,sel);
B1c  = B1map(:,:,sel);
gmax = max(mag(:));

% ---------------------------------------------------------------
% Segment tissue (brain or calf) and phantom
% ---------------------------------------------------------------
% Tissue mask.  brain: caller-supplied B1-space mask -> HD-BET NIfTI
% (resampled, when wired) -> intensity-threshold head fallback.
% calf: caller-supplied mask -> intensity-threshold leg (largest bright CC).
[tissueMask, brainSource] = computeTissueMask(extRef, mag, gmax, sel);

% Phantom: bright, compact object OUTSIDE the subject. Either supplied by the
% caller or found automatically as the largest CC outside the dilated subject.
% Use an intensity-threshold WHOLE-SUBJECT mask for the exclusion (independent
% of the tissueMask, so for brain the skull/scalp rim is not mistaken for
% phantom, and for calf the fat rim is not either).
headForExcl = mag > extRef.brainThr * gmax;
headForExcl = imfill(keepLargestCC(headForExcl), 'holes');
if isfield(extRef,'phantomMask') && ~isempty(extRef.phantomMask)
    pm = logical(extRef.phantomMask);
    if size(pm,3) == nz, pm = pm(:,:,sel); end
    phantomMask = pm;
else
    outside = (mag > extRef.phantThr * gmax) & ~imdilate3(headForExcl, 2);
    outside = bwareaopen(outside, extRef.phantMinVox);
    phantomMask = keepLargestCC(outside);
end

nVoxRef = nnz(phantomMask);
nVoxM   = nnz(tissueMask);
% An empty phantom or tissue mask would silently drive every external
% concentration to 0 or Inf, which is worse than reporting nothing.  Bail out
% the same way a missing B1 scan does.
if nVoxRef == 0 || nVoxM == 0
    if nVoxRef == 0
        msg = 'no phantom voxels segmented (check b1Dir / phantThr / phantMinVox)';
    else
        msg = 'no tissue voxels segmented (check brainThr)';
    end
    warning(['extRefQuant31P: %s. RefExt amplitude %.4g kept, but external ' ...
             'concentrations are not computed. Continuing.'], msg, areaRef);
    absExt.reason  = msg;
    absExt.nVoxRef = nVoxRef;
    absExt.nVoxTissue = nVoxM;
    return
end

V_ref_mm3 = nVoxRef * voxVol_mm3;

% Tissue volume V_M.  Priority:
%   1) explicit extRef.brainVol_mL
%   2) BRAIN ONLY: HD-BET brain mask INTERSECT ISIS slab via imscribe
%      (validated Way 1).  There is no calf analogue -- no skull-stripping
%      equivalent exists for the leg -- so calf always uses (3).
%   3) segmented tissueMask voxel count over the slab slices.  Because the B1
%      map shares the ISIS slab orientation and `sel` already drops the end
%      slices, this count IS the in-slab tissue volume.
V_M_mm3    = nVoxM * voxVol_mm3;
V_M_source = sprintf('%s (%d B1 voxels in slab)', brainSource, nVoxM);
if isfield(extRef,'brainVol_mL') && ~isempty(extRef.brainVol_mL)
    V_M_mm3    = extRef.brainVol_mL * 1000;
    V_M_source = sprintf('caller brainVol_mL = %.1f mL', extRef.brainVol_mL);
elseif strcmpi(extRef.anatomy,'brain') ...
        && all(isfield(extRef,{'brainMaskNii','t1Nii','isisDcm'})) ...
        && ~isempty(extRef.brainMaskNii) && ~isempty(extRef.t1Nii) && ~isempty(extRef.isisDcm) ...
        && isfile(extRef.brainMaskNii) && isfile(extRef.t1Nii) && isfile(extRef.isisDcm)
    try
        vml        = brainVolInISISslab31P(extRef.brainMaskNii, extRef.t1Nii, extRef.isisDcm);
        V_M_mm3    = vml * 1000;
        V_M_source = sprintf('HD-BET brain INT ISIS slab (imscribe) = %.1f mL', vml);
    catch ME
        warning('extRefQuant31P: brainVolInISISslab31P failed (%s); using threshold V_M.', ME.message);
    end
end

% ---------------------------------------------------------------
% B1+ excitation efficiency, eta = sin(relB1 * flip)   [appropriate model for adiabatic ISIS]
% ---------------------------------------------------------------
% Everything in this block is diagnostic while extRef.applyB1 is false: the
% values are stored and printed, but the concentration loop below multiplies
% by 1 instead. 
relB1_M = NaN; relB1_ref = NaN; b1method = 'not computed';
try
    % Tissue is bright -> the processed (masked, smoothed) B1 map is valid there.
    relB1_M = safeMean(B1c(tissueMask & B1c > 0));

    % Phantom B1+.  The phantom is faint (D2O/PBS), so calc_B1map_new masks it
    % to 0 AND a direct double-angle measurement there is unreliable/too high.
    % Default is to EXTRAPOLATE the smooth, valid tissue B1+ field to the
    % phantom location by a low-order polynomial fit (B1+ varies slowly).
    switch lower(extRef.phantomB1Method)
        case 'value'
            relB1_ref  = extRef.phantomB1Value;
            b1method   = 'fixed value';
        case 'measure'
            if ndims(images) == 4 && size(images,4) >= 3
                relVol   = calcRelB1_unmasked(images(:,:,sel,1), images(:,:,sel,2), images(:,:,sel,3), b1alpha);
                relB1_ref = safeMean(relVol(phantomMask & isfinite(relVol) & relVol > 0 & relVol < 1.8));
            else
                relB1_ref = safeMean(B1c(phantomMask & B1c > 0));
            end
            b1method = 'direct double-angle (raw images)';
        otherwise  % 'extrap'
            relB1_ref = extrapB1ToMask(B1c, tissueMask, phantomMask, extRef.b1ExtrapOrder);
            b1method  = sprintf('order-%d polynomial extrapolation of tissue B1+', extRef.b1ExtrapOrder);
    end
catch ME
    warning('extRefQuant31P: B1 estimation failed (%s); B1 diagnostics left NaN.', ME.message);
    b1method = sprintf('FAILED: %s', ME.message);
end
% actual flip is relB1 * (prescribed flip), so eta = sin(relB1 * flip).
flipRad = extRef.flipDeg * pi/180;
eta_ref = sin(relB1_ref * flipRad);
eta_M   = sin(relB1_M   * flipRad);

% ---------------------------------------------------------------
% Reference relaxation weighting + properties
% ---------------------------------------------------------------
R = @(T1,T2) (1 - exp(-TR_ms/T1)) * exp(-TE_ms/T2);
Rref = R(extRef.T1_ms, extRef.T2_ms);
nPref = extRef.nP;

% ---------------------------------------------------------------
% Per-metabolite external-reference concentration
% ---------------------------------------------------------------
absExt.ok         = true;
absExt.refPpm     = extRef.refPpm;
absExt.refConc_mM = extRef.refConc_mM;
absExt.areaRef    = areaRef;
absExt.T1ref_ms   = extRef.T1_ms;
absExt.T2ref_ms   = extRef.T2_ms;
absExt.nPref      = nPref;
absExt.Rref       = Rref;
absExt.V_ref_mL   = V_ref_mm3 / 1000;
absExt.V_M_mL     = V_M_mm3   / 1000;
absExt.nVoxRef    = nVoxRef;
absExt.nVoxTissue = nVoxM;
absExt.voxVol_mm3 = voxVol_mm3;
absExt.relB1_ref  = relB1_ref;
absExt.relB1_M    = relB1_M;
absExt.eta_ref    = eta_ref;
absExt.eta_M      = eta_M;
absExt.slicesUsed = sel;
absExt.TR_ms      = TR_ms;
absExt.TE_ms      = TE_ms;
absExt.metabolites = struct('name',{},'area',{},'conc_mM',{}, ...
                            'T1_ms',{},'T2_ms',{},'nPhos',{},'R',{}, ...
                            'ratio',{},'conc_internal_mM',{},'ext_over_int',{});

volFac = (V_ref_mm3 / V_M_mm3);

% ---------------------------------------------------------------
% B1 CORRECTION FACTORS -- DISABLED by default (extRef.applyB1 = false)
% ---------------------------------------------------------------
% etaFac  : B1+ transmit efficiency ratio, sin(relB1_ref*flip)/sin(relB1_M*flip)
% recvFac : B1- receive ratio under reciprocity (B1- ~ B1+), relB1_ref/relB1_M
% Both are reported below as etaFac_qc / recvFac_qc regardless.
etaFac_qc = (eta_ref / eta_M);
if isfield(extRef,'reciprocityReceive') && extRef.reciprocityReceive
    recvFac_qc = relB1_ref / relB1_M;
else
    recvFac_qc = 1;
end

if extRef.applyB1
    etaFac  = etaFac_qc;
    recvFac = recvFac_qc;
else
    etaFac  = 1;
    recvFac = 1;
end

for k = 1:numel(fitNames)
    nm = fitNames{k};
    if strcmp(nm,'RefExt'), continue; end
    pIdx = find(strcmp({metProps.name}, nm), 1);
    if isempty(pIdx), continue; end
    T1k = metProps(pIdx).T1_ms;
    T2k = metProps(pIdx).T2_ms;
    nPk = metProps(pIdx).nPhos;
    Rk  = R(T1k, T2k);

    ratio   = fitAreas(k) / areaRef;
    conc_mM = ratio * volFac * etaFac * recvFac * (Rref / Rk) * (nPref / nPk) * extRef.refConc_mM;

    % internal-reference value for the same metabolite (for comparison)
    ci = NaN;
    if ~isempty(internalAbsQ) && isfield(internalAbsQ,'metabolites')
        j = find(strcmp({internalAbsQ.metabolites.name}, nm), 1);
        if ~isempty(j), ci = internalAbsQ.metabolites(j).conc_mM; end
    end

    absExt.metabolites(end+1) = struct( ...
        'name', nm, 'area', fitAreas(k), 'conc_mM', conc_mM, ...
        'T1_ms', T1k, 'T2_ms', T2k, 'nPhos', nPk, 'R', Rk, ...
        'ratio', ratio, 'conc_internal_mM', ci, ...
        'ext_over_int', conc_mM / ci);
end

% ---------------------------------------------------------------
% Report
% ---------------------------------------------------------------
absExt.anatomy     = extRef.anatomy;
absExt.applyB1     = extRef.applyB1;
absExt.recvFac     = recvFac;       % as APPLIED (1 when applyB1 is false)
absExt.etaFac      = etaFac;        % as APPLIED (1 when applyB1 is false)
absExt.recvFac_qc  = recvFac_qc;    % as COMPUTED, for QC
absExt.etaFac_qc   = etaFac_qc;     % as COMPUTED, for QC
absExt.reciprocity = isfield(extRef,'reciprocityReceive') && extRef.reciprocityReceive;
absExt.b1method    = b1method;
absExt.volFac      = volFac;
absExt.brainSource = brainSource;
absExt.V_M_source  = V_M_source;

fprintf('\n=== External-reference 31P quant (%s) ===\n', extRef.anatomy);
fprintf('  ref %.2f ppm, [ref]=%.0f mM, nP=%g, T1/T2=%.0f/%.0f ms (placeholder)\n', ...
    extRef.refPpm, extRef.refConc_mM, nPref, extRef.T1_ms, extRef.T2_ms);
fprintf('  slab slices %s of %d  (voxel %.3g mm^3)\n', mat2str(sel), nz, voxVol_mm3);
fprintf('  tissue seg: %s\n', brainSource);
fprintf('  V_M source: %s\n', V_M_source);
fprintf('  V_phantom=%.2f mL (%d vox)   V_tissue=%.2f mL   V_ref/V_M=%.4g\n', ...
    absExt.V_ref_mL, nVoxRef, absExt.V_M_mL, volFac);
if extRef.applyB1
    fprintf('  phantom B1+ via %s\n', b1method);
    fprintf('  prescribed flip=%.1f deg -> eta = sin(relB1*flip)\n', extRef.flipDeg);
    fprintf('  relB1 phantom=%.3f (eta=%.3f)   tissue=%.3f (eta=%.3f)   eta_ref/eta_M=%.4g [APPLIED]\n', ...
        relB1_ref, eta_ref, relB1_M, eta_M, etaFac);
    if absExt.reciprocity
        fprintf('  reciprocity receive factor relB1_ref/relB1_M = %.4g [APPLIED]\n', recvFac);
    end
else
    fprintf('  *** B1 CORRECTION DISABLED (applyB1=false): etaFac and recvFac forced to 1 ***\n');
    fprintf('  QC only -- phantom B1+ via %s\n', b1method);
    fprintf('  QC only -- relB1 phantom=%.3f  tissue=%.3f  eta_ref/eta_M=%.4g  recv=%.4g\n', ...
        relB1_ref, relB1_M, etaFac_qc, recvFac_qc);
end
fprintf('  %-10s %12s %12s %10s\n','metabolite','ext[mM]','int[mM]','ext/int');
for k = 1:numel(absExt.metabolites)
    m = absExt.metabolites(k);
    fprintf('  %-10s %12.4g %12.4g %10.3f\n', m.name, m.conc_mM, m.conc_internal_mM, m.ext_over_int);
end

% ---------------------------------------------------------------
% Optional QC figure: segmentation overlay
% ---------------------------------------------------------------
if isfield(extRef,'plt') && extRef.plt
    try
        f = figure('Name','extRef segmentation','Position',[50 50 1300 700]);
        tiledlayout(f,'flow','Padding','compact','TileSpacing','compact');
        for kk = 1:numel(sel)
            nexttile;
            rgb = repmat(mat2gray(mag(:,:,kk)), [1 1 3]);
            tk = tissueMask(:,:,kk); pk = phantomMask(:,:,kk);
            r = rgb(:,:,1); g = rgb(:,:,2); b = rgb(:,:,3);
            g(tk) = min(1, g(tk) + 0.35);            % tissue -> green tint
            r(pk) = 1; g(pk) = 0.2; b(pk) = 0.2;      % phantom -> red
            rgb = cat(3,r,g,b);
            image(rgb); axis image off; title(sprintf('sl %d',sel(kk)),'FontSize',8);
        end
        if nargin >= 8 && ~isempty(outDir)
            if nargin < 9 || isempty(scanID), scanID = 'extRef'; end
            print(f, fullfile(outDir, sprintf('%s_extRef_seg.png', scanID)), '-dpng','-r100');
        end
    catch ME
        warning('extRefQuant31P: QC figure failed: %s', ME.message);
    end
end

end

% ====================================================================
function [px, slthk] = b1VoxelDims(b1Dir)
d = dir(fullfile(b1Dir,'*.IMA'));
if isempty(d), d = dir(fullfile(b1Dir,'*.dcm')); end
assert(~isempty(d), 'b1VoxelDims: no DICOM in %s', b1Dir);
di = dicominfo(fullfile(b1Dir, d(1).name));
px = double(di.PixelSpacing(:))';
if isfield(di,'SliceThickness') && ~isempty(di.SliceThickness)
    slthk = double(di.SliceThickness);
elseif isfield(di,'SpacingBetweenSlices')
    slthk = double(di.SpacingBetweenSlices);
else
    slthk = px(1);
end
end

% ====================================================================
function BW = keepLargestCC(BW)
CC = bwconncomp(BW);
if CC.NumObjects <= 1, return; end
n = cellfun(@numel, CC.PixelIdxList);
[~,imax] = max(n);
out = false(size(BW));
out(CC.PixelIdxList{imax}) = true;
BW = out;
end

% ====================================================================
function BW = imdilate3(BW, r)
% simple isotropic 3D dilation by r voxels (cube structuring element)
se = true(2*r+1, 2*r+1, 2*r+1);
BW = imdilate(BW, se);
end

% ====================================================================
function [mask, src] = computeTissueMask(extRef, mag, gmax, sel)
% Tissue mask in B1-map (slab-slice) space.  Brain and calf share the
% intensity-threshold path; only brain has the HD-BET branch.
if isfield(extRef,'brainMask') && ~isempty(extRef.brainMask)
    bm = logical(extRef.brainMask);
    if size(bm,3) ~= size(mag,3) && size(bm,3) >= max(sel)
        bm = bm(:,:,sel);
    end
    assert(isequal(size(bm), size(mag)), 'brainMask size %s ~= B1 slab-slice size %s', ...
        mat2str(size(bm)), mat2str(size(mag)));
    mask = bm; src = 'caller-supplied tissue mask (B1 space)';
    return
end
if strcmpi(extRef.anatomy,'brain') && isfield(extRef,'brainMaskNii') && ~isempty(extRef.brainMaskNii)
    warning(['extRefQuant31P: brainMaskNii set but NIfTI->B1 resampling is not ' ...
             'yet validated; using intensity-threshold head for the B1 average.']);
end
% Largest bright connected component.
%   brain -> head (includes skull/scalp)
%   calf  -> leg  (includes subcutaneous fat and marrow, which carry no
%                  PCr/ATP, so V_M is an overestimate -- see header)
body = mag > extRef.brainThr * gmax;
body = keepLargestCC(body);
body = imfill(body, 'holes');
mask = body;
if strcmpi(extRef.anatomy,'calf')
    src = 'intensity-threshold leg, incl. fat/marrow [calf]';
else
    src = 'intensity-threshold head, incl. skull/scalp [FALLBACK]';
end
end

% ====================================================================
function val = extrapB1ToMask(B1c, fitMask, tgtMask, order)
% Fit a low-order 3D polynomial to the valid tissue B1+ samples and evaluate
% it (mean) over the target (phantom) mask.  Coords are centered/scaled by
% the tissue-sample statistics for conditioning and controlled extrapolation.
sz = size(B1c, 1:3);
[I,J,K] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
valid = fitMask & B1c > 0 & B1c < 1.8;
xf = I(valid); yf = J(valid); zf = K(valid); bf = B1c(valid);
assert(numel(bf) > 20, 'extrapB1ToMask: too few valid brain B1 samples (%d)', numel(bf));
mu = [mean(xf) mean(yf) mean(zf)];
sc = [std(xf)  std(yf)  std(zf)];  sc(sc == 0) = 1;
A  = polyDesign((xf-mu(1))/sc(1), (yf-mu(2))/sc(2), (zf-mu(3))/sc(3), order);
coef = A \ bf;
xt = I(tgtMask); yt = J(tgtMask); zt = K(tgtMask);
At = polyDesign((xt-mu(1))/sc(1), (yt-mu(2))/sc(2), (zt-mu(3))/sc(3), order);
val = mean(At * coef);
end

% ====================================================================
function A = polyDesign(x, y, z, order)
x = x(:); y = y(:); z = z(:);
A = [];
for dx = 0:order
    for dy = 0:order-dx
        for dz = 0:order-dx-dy
            A = [A, (x.^dx).*(y.^dy).*(z.^dz)]; %#ok<AGROW>
        end
    end
end
end

% ====================================================================
function r1 = calcRelB1_unmasked(i1, i2, i3, alpha)
i1 = double(i1); i2 = double(i2); i3 = double(i3);
r  = i2 ./ (1.0 + i1);
r(i2 > i1) = -r(i2 > i1);
x  = 0.25 * (r + sqrt(r.*r + 8.0));
r1 = acos(min(max(x,-1),1)) / alpha;
r  = i3 ./ (1.0 + i2);
r(i3 > i2) = -r(i3 > i2);
x  = 0.25 * (r + sqrt(r.*r + 8.0));
r2 = acos(min(max(x,-1),1)) / (2*alpha);
r1(r1 < 1.0) = r2(r1 < 1.0);
end

% ====================================================================
function m = safeMean(v)
v = v(:); v = v(isfinite(v));
if isempty(v), m = NaN; else, m = mean(v); end
end
