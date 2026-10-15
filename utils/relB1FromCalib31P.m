function b1 = relB1FromCalib31P(SessionDir, refPpm, opts)
% RELB1FROMCALIB31P  Per-location relB1 from the on-scanner 31P B1 calibration.
%
%   b1 = relB1FromCalib31P(SessionDir, refPpm)
%   b1 = relB1FromCalib31P(SessionDir, refPpm, opts)
%
% Locates the ISIS_SVS *B1CALIB* flip-angle series under <SessionDir>\E* and
% fits relB1 = achieved/prescribed separately for the tissue peak (PCr) and
% the external phantom singlet (RefExt), by handing both to isis_svs_refVolt
% as a two-peak list.  PCr and every other tissue metabolite share the ISIS
% slab, so one tissue relB1 covers all of them; only the phantom, at a
% different place in the coil, needs its own.
%
% This replaces the 1H TurboFLASH B1 map as the source of the transmit
% correction for the external-reference quant.  The map was the right SHAPE
% but not the right MAGNITUDE (wrong nucleus), which is why applyB1 defaulted
% to false.  These numbers are 31P, on the same coil, in the same session, at
% the actual peak positions.
%
% *** THE CALIBRATION'S SEQUENCE PARAMETERS ARE NOT THE ACQUISITION'S. ***
% They are not expected to match and there is nothing wrong when they do not
% -- the two scans are tuned for different jobs:
%
%   calibration   long excitation pulse, so prescribed flips can go well past
%                 90 deg and the sine curve is sampled through its maximum.
%                 A long hard pulse has a NARROW spectral profile (first null
%                 at 1/T), which is the accepted cost.
%   acquisition   short pulse, so the whole 31P spectrum is excited uniformly.
%                 Flip range does not matter here; profile flatness does.
%
% TR, flip list, excitation center and reference amplitude routinely differ
% too.  So: take from this function only what is a PROPERTY OF THE COIL AND
% THE LOCATION -- relB1, which is relative to .TxVolt and must be scaled by
% (acquisition / calibration reference voltage) -- and recompute every
% sequence-dependent term
% (excitation factor, relaxation weighting) with the acquisition's own header
% values.  b1.pulseDur and b1.excCtr are returned for QC and provenance, NOT
% for reuse in the quant.  Never assume the two match; read both headers.
%
% INPUTS
%   SessionDir  session folder holding E100/E200/...
%   refPpm      external-reference chemical shift (PCr = 0), e.g. 8.0
%   opts        optional struct:
%     .files       cellstr of DICOMs, overriding auto-location
%     .tissuePpm   tissue peak shift  (default 0, PCr)
%     .tissueWin   assignment half-width ppm (default 0.6)
%     .refWin      external-ref half-width ppm (default 0.5)
%     .sigMode     'hsvd' (default) or 'integrate', see isis_svs_refVolt
%     .plt         false (default)
%
% OUTPUT b1
%   .ok          false when the calibration is missing or unusable
%   .reason      why, when ok is false
%   .relB1_M     tissue (PCr) relB1
%   .relB1_ref   phantom (RefExt) relB1
%   .rms_M/.rms_ref   per-peak sine-fit RMS; a ref RMS far above the tissue
%                     one means the phantom sat too near an excitation-profile
%                     null for its relB1 to be trusted
%   .refVolt     true 90 deg voltage from the tissue peak
%   .TxVolt      reference amplitude the calibration ran at
%   .pulseDur, .excCtr   CALIBRATION geometry (us, ppm) -- NOT the
%                acquisition's; the excitation factor for the quant must be
%                recomputed with the acquisition's own values
%   .files, .FA, .peakOffHz

if nargin < 3 || isempty(opts), opts = struct(); end
d = struct('files',{{}}, 'tissuePpm',0, 'tissueWin',0.6, 'refWin',0.5, ...
           'sigMode','hsvd', 'plt',false);
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(opts,fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = d.(fn{i}); end
end

b1 = struct('ok',false,'reason','','relB1_M',NaN,'relB1_ref',NaN, ...
            'rms_M',NaN,'rms_ref',NaN,'refVolt',NaN,'TxVolt',NaN, ...
            'pulseDur',[],'excCtr',[],'files',{{}},'FA',[],'peakOffHz',[]);

files = opts.files;
if isempty(files)
    files = findB1CalibFiles(SessionDir);
end
if isempty(files)
    b1.reason = sprintf('no *B1CALIB* series with DICOMs under %s\\E*', SessionDir);
    return
end
if numel(files) < 3
    b1.reason = sprintf('only %d B1CALIB acquisition(s); need at least 3 flip angles', numel(files));
    return
end
b1.files = files;

pars = struct();
pars.peaks   = {'PCr', opts.tissuePpm, opts.tissueWin; 'RefExt', refPpm, opts.refWin};
pars.combine = 'PCr';
pars.sigMode = opts.sigMode;
pars.plt     = opts.plt;

% isis_svs_refVolt addpaths utils + mapVBVD on entry and rmpaths BOTH on exit,
% so calling it mid-pipeline would strip utils from under the caller and the
% next utils function would not resolve.  Restore whatever the caller had.
pathSnapshot = path;
cleanupPath  = onCleanup(@() path(pathSnapshot));

try
    [refVolt, out] = isis_svs_refVolt(files, pars);
catch ME
    % A protocol change partway through the series (excitation center or pulse
    % moved) trips isis_svs_refVolt's consistency check.  Report it rather
    % than silently fitting whichever subset happens to agree -- the halves
    % measure different things and must not be pooled.
    b1.reason = sprintf('isis_svs_refVolt failed: %s', ME.message);
    return
end

iM = find(strcmp(out.peakName,'PCr'),   1);
iR = find(strcmp(out.peakName,'RefExt'),1);
b1.relB1_M   = out.relB1Peak(iM);
b1.relB1_ref = out.relB1Peak(iR);
b1.rms_M     = out.rmsPeak(iM);
b1.rms_ref   = out.rmsPeak(iR);
b1.refVolt   = refVolt;
b1.TxVolt    = out.TxVolt;
b1.pulseDur  = out.pulseDur;
b1.excCtr    = out.excCtr;
b1.FA        = out.FA;
b1.peakOffHz = out.peakOffHz;

if ~(isfinite(b1.relB1_M) && b1.relB1_M > 0)
    b1.reason = 'tissue (PCr) relB1 did not fit';
    return
end
if ~(isfinite(b1.relB1_ref) && b1.relB1_ref > 0)
    b1.reason = 'phantom (RefExt) relB1 did not fit -- peak too weak, or too near an excitation-profile null';
    return
end

% The phantom sees less B1 than the head, so its achieved flip is lower and
% the series may never reach the sine maximum, leaving relB1 and amplitude
% partly degenerate.  Warn when the highest achieved flip stays well short of
% 90 deg, which is where the estimate stops being well determined.
achMax = max(out.FA) * b1.relB1_ref;
if achMax < 85
    warning(['relB1FromCalib31P: highest ACHIEVED flip at the phantom is %.0f deg ' ...
             '(prescribed %.0f x relB1 %.3f).  The flip series does not reach the ' ...
             'sine maximum there, so relB1_ref is weakly determined.'], ...
            achMax, max(out.FA), b1.relB1_ref);
end
if isfinite(b1.rms_M) && isfinite(b1.rms_ref) && b1.rms_ref > 5*b1.rms_M
    warning(['relB1FromCalib31P: RefExt sine-fit RMS (%.4g) is %.0fx the PCr RMS ' ...
             '(%.4g).  The phantom is likely near an excitation-profile null at ' ...
             'this excitation center (%g ppm, %g us pulse).'], ...
            b1.rms_ref, b1.rms_ref/b1.rms_M, b1.rms_M, b1.excCtr, b1.pulseDur);
end

b1.ok = true;
end

% =====================================================================
function files = findB1CalibFiles(SessionDir)
% One DICOM per *B1CALIB* series under SessionDir\E*, series order preserved.
% Empty when the session has no calibration -- not an error here; the caller
% decides whether to fall back or give up.
files = {};
exams = dir(fullfile(SessionDir,'E*'));
exams = exams([exams.isdir]);
cand  = [];
for k = 1:numel(exams)
    c = dir(fullfile(exams(k).folder, exams(k).name, '*B1CALIB*'));
    c = c([c.isdir]);
    if ~isempty(c), cand = [cand; c(:)]; end %#ok<AGROW>
end
if isempty(cand), return, end

[~, ord] = sort({cand.name});
cand = cand(ord);
for k = 1:numel(cand)
    p = fullfile(cand(k).folder, cand(k).name);
    f = dir(fullfile(p,'*.IMA'));
    if isempty(f), f = dir(fullfile(p,'*.dcm')); end
    f = f(~[f.isdir]);
    if isempty(f), continue, end          % series exported with no files
    files{end+1} = fullfile(f(1).folder, f(1).name); %#ok<AGROW>
end
end
