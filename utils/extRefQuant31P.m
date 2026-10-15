function absExt = extRefQuant31P(fitNames, fitAreas, metProps, TR_ms, TE_ms, extRef, internalAbsQ)
% EXTREFQUANT31P  Absolute 31P quant from an EXTERNAL reference phantom.
%
%   absExt = extRefQuant31P(fitNames, fitAreas, metProps, TR_ms, TE_ms, ...
%                           extRef, internalAbsQ)
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
% Every input has exactly one source.  There are no fallbacks: a missing input
% stops the quant with a message naming it, rather than substituting something
% that yields a plausible wrong number.
%
%   relB1  the 31P ISIS B1CALIB flip series, fit per peak by relB1FromCalib31P
%          -- tissue from PCr, phantom from RefExt.  Located automatically
%          under extRef.sessionDir\E*.
%   V_ref  refArea_mm2 * slabThickness_mm, geometric.
%   V_M    explicit brainVol_mL, or HD-BET (brain) / Dixon (calf) mask
%          INTERSECT the ISIS slab.
%
% The 1H TurboFLASH B1 map is GONE, along with its threshold segmentation.  It
% was the wrong nucleus for relB1, its masks under-read the faint D2O/PBS
% phantom by ~3x for V_ref, and its voxel-count V_M silently returned a
% skull-and-scalp head mask (brain) or under-read the excited leg by 1.5-3x
% (calf).  Do not reintroduce it.
%
% B1 correction is applied only when extRef.applyB1 is true; it is still off by
% default, but the numbers behind it are now measured rather than extrapolated.
%
% The isis_svs sequence is 1D-ISIS (slab select via adiabatic pulses + hard
% pulse excitation).  The hard pulse is not spatially selective, so there is
% no chemical-shift displacement from it, but it does have a SPECTRAL profile:
% peaks at different shifts are excited differently.  With the acquisition's
% pulse duration and center available (extRef.acqPulseDur / .acqExcCtr /
% .hzPerPpm) the factor is rectMxy(flip, T, (ppm-ctr)*hz, relB1), evaluated
% per metabolite; without them it degrades to the on-resonance sine,
% sin(relB1*flip), which is the previous behavior.  Note the calibration's
% pulse and center are deliberately NOT the acquisition's -- the calibration
% uses a long pulse to reach flips past 90 deg, which narrows its profile --
% so the factor is recomputed here rather than carried over from the fit.
%
% ANATOMY: extRef.anatomy selects which V_M route applies.
%   'brain' - HD-BET brain mask INTERSECT the ISIS slab (imscribe).
%   'calf'  - Dixon water/fat -> muscle mask INTERSECT the ISIS slab, which
%             excludes subcutaneous fat and marrow (they carry no PCr / ATP,
%             so counting them would bias the calf concentrations LOW).
%
% extRef fields (see extRefDefaults in processU01_Brain_31P / _Calf_31P):
%   sessionDir, refPpm, refConc_mM, refArea_mm2, slabThickness_mm, nP, flipDeg,
%   T1_ms, T2_ms, anatomy, applyB1, reciprocityReceive, b1CalibFiles,
%   acqPulseDur / acqExcCtr / hzPerPpm (acquisition excitation geometry),
%   basisInfo, line (locateRefExtLine31P result; not found stops the quant),
%   brainVol_mL | brainMaskNii+t1Nii+isisDcm | dixonMagDir+isisDcm.
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
if isfield(extRef,'line') && ~isempty(extRef.line)
    absExt.line = extRef.line;
    if ~extRef.line.found
        msg = sprintf('no phantom line in the spectrum (best %.2f ppm, SNR %.1f)', ...
                      extRef.line.ppmFound, extRef.line.snr);
        warning('extRefQuant31P: %s. RefExt amplitude %.4g kept; no concentrations.', msg, areaRef);
        absExt.reason = msg;
        return
    end
end

% Locate and fit the calibration here rather than in the drivers: it is only
% ever used for absolute quantitation, so it belongs with the code that
% consumes it.  A caller may still pre-supply extRef.b1calib (already fitted)
% or extRef.b1CalibFiles (explicit DICOMs), both of which win over the search.
if strcmpi(getfielddef(extRef,'b1Source','b1calib'),'b1calib') && ...
        ~isfield(extRef,'b1calib')
    sd = getfielddef(extRef,'sessionDir','');
    if isempty(sd)
        extRef.b1calib = struct('ok',false, ...
            'reason','extRef.sessionDir not set, cannot locate the B1CALIB series');
    else
        extRef.b1calib = relB1FromCalib31P(sd, extRef.refPpm, ...
                            struct('files',{getfielddef(extRef,'b1CalibFiles',{})}));
    end
    if extRef.b1calib.ok
        fprintf(['extRef: 31P B1CALIB relB1  tissue(PCr) %.4f  phantom(RefExt) %.4f' ...
                 '  [refVolt %.1f V]\n'], extRef.b1calib.relB1_M, ...
                extRef.b1calib.relB1_ref, extRef.b1calib.refVolt);
    else
        fprintf('extRef: 31P B1CALIB unavailable (%s)\n', extRef.b1calib.reason);
    end
end
% Recorded now so the relB1 numbers survive the early returns below when
% the volume path cannot run.
absExt.b1calib = getfielddef(extRef,'b1calib',struct('ok',false));


% ---------------------------------------------------------------
% Phantom volume V_ref
% ---------------------------------------------------------------
% The phantom is a 15 mL tube standing across the slab, so the EXCITED
% reference volume is its cross-section times the slab thickness -- geometric,
% nothing to segment.
%
% refArea_mm2 is measured off the images (200 mm2, observed 180-220, so ~+-10%
% straight onto every external concentration).  slabThickness_mm is read from
% the header per scan and is NOT constant (40 mm brain, 60 mm calf), so it must
% never be hard-coded.  Both are required: there is no fallback.
haveArea = isfield(extRef,'refArea_mm2')      && ~isempty(extRef.refArea_mm2)      && extRef.refArea_mm2 > 0;
haveSlab = isfield(extRef,'slabThickness_mm') && ~isempty(extRef.slabThickness_mm) && extRef.slabThickness_mm > 0;
if ~haveArea || ~haveSlab
    if ~haveArea
        msg = 'extRef.refArea_mm2 not set';
    else
        msg = 'ISIS slab thickness not readable from the header';
    end
    warning(['extRefQuant31P: %s, so V_ref cannot be computed. RefExt amplitude ' ...
             '%.4g kept; no concentrations.'], msg, areaRef);
    absExt.reason = msg;
    return
end
V_ref_mm3    = extRef.refArea_mm2 * extRef.slabThickness_mm;
V_ref_source = sprintf('cross-section %.0f mm2 * %.0f mm slab (hdr)', ...
                       extRef.refArea_mm2, extRef.slabThickness_mm);

% ---------------------------------------------------------------
% Tissue volume V_M
% ---------------------------------------------------------------
% Priority:
%   1) explicit extRef.brainVol_mL
%   2) BRAIN: HD-BET brain mask INTERSECT ISIS slab via imscribe (validated
%      Way 1)
%   3) CALF: Dixon water/fat -> muscle mask INTERSECT ISIS slab
% There is no fallback.  The old route -- counting bright voxels in the 1H B1
% map and scaling that count from the B1 slab to the ISIS slab -- is gone with
% the rest of the B1 mapping: for brain it silently returned a head mask
% including skull and scalp, and for calf it under-read the excited leg by
% 1.5-3x depending on B1 coverage.  A missing anatomic volume now stops the
% quant instead of producing a plausible wrong number.
V_M_mm3     = NaN;
V_M_source  = 'not determined';
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
        warning('extRefQuant31P: brainVolInISISslab31P failed (%s); V_M unavailable.', ME.message);
    end
elseif strcmpi(extRef.anatomy,'calf') ...
        && isfield(extRef,'dixonMagDir') && ~isempty(extRef.dixonMagDir) ...
        && isfolder(extRef.dixonMagDir) ...
        && isfield(extRef,'isisDcm') && ~isempty(extRef.isisDcm) && isfile(extRef.isisDcm)
    % Calf analog: Dixon water/fat -> muscle mask INTERSECT ISIS slab.
    % Excludes subcutaneous fat and marrow, and needs no slab extrapolation --
    % the Dixon spans ~165 mm against a 60 mm slab.
    try
        [vml, dxi]  = muscleVolInISISslab31P(extRef.dixonMagDir, extRef.isisDcm, ...
                          struct('ffThresh', extRef.ffThresh));
        V_M_mm3     = vml * 1000;
        V_M_source  = sprintf(['Dixon muscle INT ISIS slab (imscribe) = %.1f mL ' ...
                               '(FF<%.2f, muscle %.0f%% of leg)'], ...
                               vml, dxi.ffThresh, 100*dxi.muscleFrac);
        absExt.dixon = dxi;
    catch ME
        warning('extRefQuant31P: muscleVolInISISslab31P failed (%s); V_M unavailable.', ME.message);
    end
end

% ---------------------------------------------------------------
% Both volumes settled -- this is where a missing input becomes fatal
% ---------------------------------------------------------------
% Checked on the VOLUMES rather than on the mask voxel counts, because the
% masks are only one of several ways to get here and are no longer the
% expected one.  A zero or NaN volume would drive every concentration to 0,
% Inf or NaN, which is worse than reporting nothing.
if ~(isfinite(V_ref_mm3) && V_ref_mm3 > 0)
    msg = ['V_ref could not be determined: set extRef.refArea_mm2 and make the ' ...
           'ISIS slab thickness readable from the header'];
    warning('extRefQuant31P: %s. RefExt amplitude %.4g kept; no concentrations.', msg, areaRef);
    absExt.reason = msg;
    return
end
if ~(isfinite(V_M_mm3) && V_M_mm3 > 0)
    msg = ['V_M could not be determined: supply extRef.brainVol_mL, or the ' ...
           'HD-BET mask + T1 + ISIS DICOM (brain), or the Dixon series + ISIS ' ...
           'DICOM (calf)'];
    warning('extRefQuant31P: %s. RefExt amplitude %.4g kept; no concentrations.', msg, areaRef);
    absExt.reason = msg;
    return
end

% ---------------------------------------------------------------
% relB1, from the 31P B1CALIB flip series -- the only source
% ---------------------------------------------------------------
relB1_M = NaN; relB1_ref = NaN; b1method = 'unavailable'; voltScale = NaN;
bc = getfielddef(extRef,'b1calib',struct('ok',false));
acqV = getfielddef(extRef,'acqRefVolt',NaN);
calV = getfielddef(bc,'TxVolt',NaN);
if isfield(bc,'ok') && bc.ok && ~(isfinite(acqV) && acqV > 0 && isfinite(calV) && calV > 0)
    bc.ok = false;
    bc.reason = sprintf(['reference voltage missing (acquisition %g V, calibration %g V), ' ...
                         'so relB1 cannot be carried to the acquisition'], acqV, calV);
end
if isfield(bc,'ok') && bc.ok
    voltScale = acqV / calV;
    relB1_M   = bc.relB1_M   * voltScale;
    relB1_ref = bc.relB1_ref * voltScale;
    % The calibration's pulse and center are recorded because they are NOT the
    % acquisition's and must never be reused as if they were (see
    % relB1FromCalib31P).  Only relB1 crosses over, rescaled to the
    % acquisition's reference voltage.
    b1method  = sprintf(['31P ISIS B1CALIB flip series (PCr / RefExt), ' ...
                         'calib %g us pulse at %g ppm, relB1 x %g/%g V'], ...
                         bc.pulseDur, bc.excCtr, acqV, calV);
else
    reason = 'not available';
    if isfield(bc,'reason') && ~isempty(bc.reason), reason = bc.reason; end
    warning('extRefQuant31P: 31P B1CALIB unusable (%s) -- there is no relB1.', reason);
end

% With applyB1 on, a missing relB1 would multiply every concentration by NaN.
% Stop instead, naming what to fix.
if extRef.applyB1 && ~(isfinite(relB1_M) && relB1_M > 0 && isfinite(relB1_ref) && relB1_ref > 0)
    msg = sprintf(['applyB1 is on but relB1 is unavailable (tissue %.4g, phantom %.4g; ' ...
                   'source: %s)'], relB1_M, relB1_ref, b1method);
    warning('extRefQuant31P: %s. RefExt amplitude %.4g kept; no concentrations.', msg, areaRef);
    absExt.reason = msg;
    return
end

% ---------------------------------------------------------------
% Excitation factor at the ACQUISITION's own geometry
% ---------------------------------------------------------------
% relB1 is the SPATIAL transmit scaling; the excitation pulse also has a
% SPECTRAL profile, so two peaks at different shifts are excited differently
% even at the same relB1.  With a hard pulse both are in rectMxy.  The
% calibration's pulse and center are NOT the acquisition's (the calibration
% deliberately uses a long pulse to reach flips past 90 deg, which narrows its
% profile), so this is recomputed here rather than reused from the fit.
%
% Without the acquisition geometry this degrades to the on-resonance sine,
% eta = sin(relB1 * flip), i.e. no profile term -- the previous behavior.
acqTus = getfielddef(extRef,'acqPulseDur',[]);
acqCtr = getfielddef(extRef,'acqExcCtr',[]);
hzPpm  = getfielddef(extRef,'hzPerPpm',[]);
useProfile = ~isempty(acqTus) && ~isempty(acqCtr) && ~isempty(hzPpm) && ...
             isfinite(acqTus) && isfinite(acqCtr) && isfinite(hzPpm);

if useProfile
    etaOf   = @(ppm,rel) rectMxy(extRef.flipDeg, acqTus, (ppm-acqCtr)*hzPpm, rel);
    excModel = sprintf('rect %g us at %g ppm, flip %g deg', acqTus, acqCtr, extRef.flipDeg);
else
    flipRad = extRef.flipDeg * pi/180;
    etaOf   = @(~,rel) sin(rel * flipRad);       % ppm ignored: no profile term
    excModel = sprintf('on-resonance sine, flip %g deg (no excitation profile)', extRef.flipDeg);
end
eta_ref = etaOf(extRef.refPpm, relB1_ref);
eta_M   = etaOf(0,             relB1_M);   % tissue at PCr; per-metabolite below

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
absExt.V_ref_mL     = V_ref_mm3 / 1000;
absExt.V_ref_source = V_ref_source;
absExt.V_M_mL     = V_M_mm3   / 1000;
absExt.relB1_ref  = relB1_ref;
absExt.relB1_M    = relB1_M;
absExt.relB1_ref_calib = getfielddef(bc,'relB1_ref',NaN);
absExt.relB1_M_calib   = getfielddef(bc,'relB1_M',NaN);
absExt.relB1VoltScale  = voltScale;
absExt.excModel   = excModel;
absExt.eta_ref    = eta_ref;
absExt.eta_M      = eta_M;
absExt.TR_ms      = TR_ms;
absExt.TE_ms      = TE_ms;
absExt.metabolites = struct('name',{},'area',{},'conc_mM',{}, ...
                            'T1_ms',{},'T2_ms',{},'nPhos',{},'R',{}, ...
                            'ratio',{},'conc_internal_mM',{},'ext_over_int',{}, ...
                            'ppm',{},'eta',{},'etaFac',{});

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

    % Excitation factor at THIS metabolite's shift.  All tissue metabolites
    % share the slab and so share relB1_M, but they sit at different offsets
    % from the excitation center and are excited differently.  ppmOfMet
    % returns [] when the shift is unknown, which falls back to the PCr value.
    ppmK   = ppmOfMet(extRef, nm);
    if useProfile && ~isempty(ppmK)
        etaK = etaOf(ppmK, relB1_M);
    else
        etaK = eta_M;
    end
    etaFac_k = etaFac;
    if extRef.applyB1 && etaK > 0
        etaFac_k = eta_ref / etaK;
    end

    ratio   = fitAreas(k) / areaRef;
    conc_mM = ratio * volFac * etaFac_k * recvFac * (Rref / Rk) * (nPref / nPk) * extRef.refConc_mM;

    % internal-reference value for the same metabolite (for comparison)
    ci = NaN;
    if ~isempty(internalAbsQ) && isfield(internalAbsQ,'metabolites')
        j = find(strcmp({internalAbsQ.metabolites.name}, nm), 1);
        if ~isempty(j), ci = internalAbsQ.metabolites(j).conc_mM; end
    end
    if isempty(ppmK), ppmK = NaN; end

    absExt.metabolites(end+1) = struct( ...
        'name', nm, 'area', fitAreas(k), 'conc_mM', conc_mM, ...
        'T1_ms', T1k, 'T2_ms', T2k, 'nPhos', nPk, 'R', Rk, ...
        'ratio', ratio, 'conc_internal_mM', ci, ...
        'ext_over_int', conc_mM / ci, ...
        'ppm', ppmK, 'eta', etaK, 'etaFac', etaFac_k);
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
absExt.V_M_source  = V_M_source;

fprintf('\n=== External-reference 31P quant (%s) ===\n', extRef.anatomy);
fprintf('  ref %.2f ppm, [ref]=%.0f mM, nP=%g, T1/T2=%.0f/%.0f ms (placeholder)\n', ...
    extRef.refPpm, extRef.refConc_mM, nPref, extRef.T1_ms, extRef.T2_ms);
fprintf('  V_M source: %s\n', V_M_source);
fprintf('  V_phantom=%.3f mL   V_tissue=%.1f mL   V_ref/V_M=%.4g\n', ...
    absExt.V_ref_mL, absExt.V_M_mL, volFac);
if extRef.applyB1
    fprintf('  phantom B1+ via %s\n', b1method);
    fprintf('  excitation model: %s\n', excModel);
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
fprintf('  %-10s %8s %8s %8s %12s %12s %10s\n','metabolite','ppm','eta','etaFac','ext[mM]','int[mM]','ext/int');
fprintf('  %-10s %8.2f %8.4f %8s\n', 'RefExt', extRef.refPpm, eta_ref, '-');
for k = 1:numel(absExt.metabolites)
    m = absExt.metabolites(k);
    fprintf('  %-10s %8.2f %8.4f %8.4f %12.4g %12.4g %10.3f\n', m.name, m.ppm, m.eta, m.etaFac, ...
        m.conc_mM, m.conc_internal_mM, m.ext_over_int);
end


end

% =====================================================================
function v = getfielddef(s, f, dflt)
% Struct field with a default, so older extRef structs still work.
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

% =====================================================================
function ppm = ppmOfMet(extRef, name)
% Chemical shift of one fitted component, from the basis table the driver
% passes through as extRef.basisInfo (basisInfo(k).name / .refPpm from
% make31P_basisFIDs).  Empty when unavailable -- callers fall back to the
% tissue reference shift rather than guessing.
ppm = [];
if ~isfield(extRef,'basisInfo') || isempty(extRef.basisInfo), return, end
bi = extRef.basisInfo;
j  = find(strcmp({bi.name}, name), 1);
if isempty(j) || ~isfield(bi,'refPpm') || ~isfinite(bi(j).refPpm), return, end
ppm = bi(j).refPpm;
end
