function collectU01_2026_TSV(dataRoot, outDir)
% COLLECTU01_2026_TSV  Gather the 2026 U01 reprocessing results into TSVs.
%
%   collectU01_2026_TSV()
%   collectU01_2026_TSV(dataRoot, outDir)
%
% Walks every session under dataRoot, reads whatever
% <session>\processed\*_abs.mat files exist, and writes TWO tab-separated
% tables:
%
% A session folder is one holding datfiles\, processed\ or an E<n> exam
% folder.  Any other folder is searched up to 3 levels down, so dataRoot can
% be a flat folder of sessions (U01_2026\2026, with Repeatability\ nested), or
% one level up from per-anatomy folders (859875_U01\Brain, ...\Calf).
%
%   <outDir>\U01_2026_brain.tsv
%   <outDir>\U01_2026_calf.tsv
%
% Sessions with no processed output are omitted.  Any individual quantity that
% is missing comes out as NaN rather than dropping the row, so a partially
% processed session is still visible.
%
% COLUMNS
%   session, subject, scandate, anatomy
%   1H  [uM] : NADH2, NADH6, NADH4, NAD_1H_mean, Trp
%              NAD_1H_mean is the mean of the H2/H4/H6 resonances, which are
%              three views of the SAME NAD+ pool -- they are reported
%              individually as well because they do not agree (H4 in
%              particular sits near the fit-window edge).  Pick deliberately;
%              do not assume the mean is the right summary.
%   31P [mM] : NADplus, NADH, UDP (= UDPGal+UDPGlc+UDPGalNAc+UDPGlcNAc),
%              aATP, PCr -- reported THREE times, once per reference:
%                aATPref_* : internal aATP  (brain 2.7 mM, calf 8.2 mM)
%                PCrref_*  : internal PCr   (brain 3.5 mM, calf 33 mM)
%                extref_*  : external phantom (500 mM), only for sessions run
%                            with -ExtRef and a usable B1 map
%              The three differ by ONE common scale factor per scan -- no
%              ratio between metabolites changes.  Comparing them shows how
%              much the absolute scale depends on the reference choice.
%   RefExt_area, extRef_ok, extRef_reason, extRef_refConc_mM,
%   extRef_V_ref_mL, extRef_V_M_mL, extRef_applyB1
%   ref_aATP_conc_mM, ref_PCr_conc_mM : the assumed [ref] behind each set
%   brain only: f_GM, f_WM, f_CSF  (voxel tissue fractions from brainseg.sh)
%   TR_ms, TE_ms, flip_deg for both nuclei, for traceability
%   prot_1H, seq_1H, prot_31P, seq_31P, noe_31P : acquisition protocol, read
%              back out of the Siemens header that each _abs.mat stores.  The
%              protocol was revised mid-study, so these are what identifies
%              which family a session belongs to -- folder names do NOT (e.g.
%              S095_05212026 still carries a "...hann..." protocol name but
%              was acquired with the newer settings).  The 31P revision moved
%              TE 0.3 -> 0.1 ms, excitation flip 90 -> 70 deg and the NOE
%              pulse from hann to rect (noe_31P, sWipMemBlock.alFree(5:6),
%              goes 9/500 -> 3/200); TE_31P_ms alone separates the two.
%
% NOTE: the internal-reference concentrations depend on the PLACEHOLDER T1/T2
% table in processU01_*_31P.m, and the external ones scale linearly with the
% assumed phantom concentration.  Neither is publication-ready as-is.

if nargin < 1 || isempty(dataRoot)
    dataRoot = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026\2026';
end
if nargin < 2 || isempty(outDir)
    outDir = fileparts(dataRoot);
end
assert(isfolder(dataRoot), 'dataRoot not found: %s', dataRoot);
if ~isfolder(outDir), mkdir(outDir); end

[sessions, cohorts] = findSessions(dataRoot);
fprintf('Found %d session folder(s) under %s\n', numel(sessions), dataRoot);

rowsBrain = {};
rowsCalf  = {};

for k = 1:numel(sessions)
    s      = sessions{k};
    [~, scanID] = fileparts(s);
    procDir = fullfile(s, 'processed');
    if ~isfolder(procDir), continue; end

    for anat = ["brain" "calf"]
        f1H  = fullfile(procDir, sprintf('%s_%s_1H_abs.mat',  scanID, anat));
        f31P = fullfile(procDir, sprintf('%s_%s_31P_abs.mat', scanID, anat));
        if ~isfile(f1H) && ~isfile(f31P), continue; end

        r = buildRow(scanID, anat, f1H, f31P);
        r.cohort = cohorts{k};
        if anat == "brain", rowsBrain{end+1} = r; else, rowsCalf{end+1} = r; end %#ok<AGROW>
        fprintf('  %-50s %-5s  1H:%d 31P:%d\n', scanID, anat, isfile(f1H), isfile(f31P));
    end
end

writeTsv(rowsBrain, fullfile(outDir,'U01_2026_brain.tsv'), true);
writeTsv(rowsCalf,  fullfile(outDir,'U01_2026_calf.tsv'),  false);

end

% =====================================================================
function [sessions, cohorts] = findSessions(root, cohort, depthLeft)
% Session folders under root (see isSessionDir); other folders are searched
% down to depthLeft levels.  cohorts{k} is 'repeat' for sessions living under
% a Repeatability\ folder (the test-retest visits) and 'main' otherwise, so
% the repeatability analysis can select them without re-deriving membership
% from folder names.
if nargin < 2, cohort = 'main'; end
if nargin < 3, depthLeft = 3; end
sessions = {}; cohorts = {};
d = subdirs(root);
for k = 1:numel(d)
    p = fullfile(d(k).folder, d(k).name);
    c = cohort;
    if strcmpi(d(k).name, 'Repeatability'), c = 'repeat'; end
    if isSessionDir(p)
        sessions{end+1} = p; %#ok<AGROW>
        cohorts{end+1}  = c; %#ok<AGROW>
    elseif depthLeft > 1
        [s, cc] = findSessions(p, c, depthLeft - 1);
        sessions = [sessions s];  %#ok<AGROW>
        cohorts  = [cohorts cc];  %#ok<AGROW>
    end
end
end

function tf = isSessionDir(p)
d  = subdirs(p);
n  = {d.name};
tf = any(strcmpi(n, 'datfiles')) || any(strcmpi(n, 'processed')) || ...
     any(~cellfun(@isempty, regexp(n, '^E\d+$', 'once')));
end

function d = subdirs(p)
d = dir(p);
d = d([d.isdir] & ~startsWith({d.name}, '.'));
end

% =====================================================================
function r = buildRow(scanID, anat, f1H, f31P)
r = struct();
r.session  = scanID;
[r.subject, r.scandate] = parseScanID(scanID);
r.anatomy  = char(anat);

% ---- defaults (everything NaN until proven otherwise) ----
% NR is the 9.7 ppm NR/NMN resonance -- fitted in the calf basis only, so it
% comes out NaN for brain.  Peak10_3 is deliberately not carried here.
% NAD_1H_mean averages H2/H4/H6; NAD_H2H6_mean averages only the two
% resonances that behave (H4 sits at the fit-window edge -- see the note above).
oneH  = {'NADH2','NADH6','NADH4','NAD_1H_mean','NAD_H2H6_mean','Trp','NR'};
p31   = {'NADplus','NADH','UDP','aATP','PCr'};
sets  = {'aATPref','PCrref','extref'};
for i = 1:numel(oneH), r.(['uM_' oneH{i}]) = NaN; end
for s = 1:numel(sets)
    for i = 1:numel(p31), r.([sets{s} '_mM_' p31{i}]) = NaN; end
end
r.f_GM = NaN; r.f_WM = NaN; r.f_CSF = NaN;
r.TR_1H_ms = NaN; r.TE_1H_ms = NaN; r.TR_31P_ms = NaN; r.TE_31P_ms = NaN;
r.flip_1H_deg = NaN; r.flip_31P_deg = NaN;
r.prot_1H = ''; r.seq_1H = ''; r.prot_31P = ''; r.seq_31P = ''; r.noe_31P = '';
r.primary_ref_name = ''; r.ref_aATP_conc_mM = NaN; r.ref_PCr_conc_mM = NaN;
r.RefExt_area = NaN; r.extRef_ok = false; r.extRef_reason = '';
r.extRef_refConc_mM = NaN; r.extRef_V_ref_mL = NaN; r.extRef_V_M_mL = NaN;
r.extRef_applyB1 = false;

% ---- 1H ----
if isfile(f1H)
    try
        m = loadPresent(f1H, {'absQ','segFracs'});
        c = concMap(m.absQ);                      % name -> conc_mM
        for i = 1:3
            if isKey(c, oneH{i}), r.(['uM_' oneH{i}]) = c(oneH{i}) * 1e3; end
        end
        if isKey(c,'Trp'), r.uM_Trp = c('Trp') * 1e3; end
        if isKey(c,'NR'),  r.uM_NR  = c('NR')  * 1e3; end
        trio = [r.uM_NADH2 r.uM_NADH6 r.uM_NADH4];
        if any(~isnan(trio)), r.uM_NAD_1H_mean = mean(trio, 'omitnan'); end
        duo = [r.uM_NADH2 r.uM_NADH6];
        if any(~isnan(duo)), r.uM_NAD_H2H6_mean = mean(duo, 'omitnan'); end

        % tissue fractions live in the 1H mat as segFracs = [CSF GM WM]
        if isfield(m,'segFracs') && numel(m.segFracs) == 3
            r.f_CSF = m.segFracs(1); r.f_GM = m.segFracs(2); r.f_WM = m.segFracs(3);
        end
        [r.TR_1H_ms, r.TE_1H_ms] = seqTiming(m.absQ);
        % absQ carries no timing for 1H, but the raw header rides along in the
        % .mat -- take protocol identity and timing from there.
        if isfield(m,'hdr')
            p = protoInfo(m.hdr);
            r.prot_1H = p.protocol;  r.seq_1H = p.seqfile;
            r.flip_1H_deg = p.flipDeg;
            if isnan(r.TR_1H_ms), r.TR_1H_ms = p.TR_ms; end
            if isnan(r.TE_1H_ms), r.TE_1H_ms = p.TE_ms; end
        end
    catch ME
        warning('1H read failed for %s: %s', scanID, ME.message);
    end
end

% ---- 31P ----
if isfile(f31P)
    try
        m = loadPresent(f31P, {'absQ','segFracs'});

        % ---- the two internal-reference sets ----
        % absQ.byRef.<refName> is written by the current pipeline.  Older
        % .mat files only have absQ.metabolites (the primary reference), so
        % fall back to that and file it under whichever reference it used.
        if isfield(m.absQ,'byRef')
            if isfield(m.absQ.byRef,'aATP')
                r = fillSet(r, 'aATPref', m.absQ.byRef.aATP);
                r.ref_aATP_conc_mM = m.absQ.byRef.aATP.reference.conc_mM;
            end
            if isfield(m.absQ.byRef,'PCr')
                r = fillSet(r, 'PCrref', m.absQ.byRef.PCr);
                r.ref_PCr_conc_mM = m.absQ.byRef.PCr.reference.conc_mM;
            end
        elseif isfield(m.absQ,'reference')
            switch m.absQ.reference.name
                case 'aATP'
                    r = fillSet(r, 'aATPref', m.absQ);
                    r.ref_aATP_conc_mM = m.absQ.reference.conc_mM;
                case 'PCr'
                    r = fillSet(r, 'PCrref', m.absQ);
                    r.ref_PCr_conc_mM = m.absQ.reference.conc_mM;
            end
        end
        if isfield(m.absQ,'reference')
            r.primary_ref_name = m.absQ.reference.name;
        end
        if isfield(m.absQ,'TR_ms'), r.TR_31P_ms = m.absQ.TR_ms; end
        if isfield(m.absQ,'TE_ms'), r.TE_31P_ms = m.absQ.TE_ms; end
        if isfield(m,'hdr')
            p = protoInfo(m.hdr);
            r.prot_31P = p.protocol;  r.seq_31P = p.seqfile;
            r.flip_31P_deg = p.flipDeg;  r.noe_31P = p.noe;
            if isnan(r.TR_31P_ms), r.TR_31P_ms = p.TR_ms; end
            if isnan(r.TE_31P_ms), r.TE_31P_ms = p.TE_ms; end
        end

        % brain 31P also stores segFracs; only use it if 1H did not supply them
        if isnan(r.f_GM) && isfield(m,'segFracs') && numel(m.segFracs) == 3 ...
                && ~all(isnan(m.segFracs))
            r.f_CSF = m.segFracs(1); r.f_GM = m.segFracs(2); r.f_WM = m.segFracs(3);
        end

        % ---- external reference ----
        if isfield(m.absQ,'external')
            e = m.absQ.external;
            if isfield(e,'areaRef'),    r.RefExt_area       = e.areaRef;    end
            if isfield(e,'refConc_mM'), r.extRef_refConc_mM = e.refConc_mM; end
            if isfield(e,'reason'),     r.extRef_reason     = e.reason;     end
            if isfield(e,'error'),      r.extRef_reason     = e.error;      end
            if isfield(e,'applyB1'),    r.extRef_applyB1    = logical(e.applyB1); end
            r.extRef_ok = isfield(e,'ok') && e.ok;
            if r.extRef_ok
                r = fillSet(r, 'extref', e);
                if isfield(e,'V_ref_mL'), r.extRef_V_ref_mL = e.V_ref_mL; end
                if isfield(e,'V_M_mL'),   r.extRef_V_M_mL   = e.V_M_mL;   end
            end
        end
    catch ME
        warning('31P read failed for %s: %s', scanID, ME.message);
    end
end
end

% =====================================================================
function m = loadPresent(matFile, ~)
% Load the .mat in ONE open and let the caller isfield() its way around.
%
% Do NOT reintroduce whos('-file',...) + a selective load() here: that is two
% opens per file, and on the \\flip\nwilson SMB share the handles are not
% released promptly, so a full 21-session pass exhausts them and MATLAB dies
% with 'Too many open files'.  A selective load() alone is no good either --
% naming a variable the file lacks (segFracs, which the calf pipeline never
% writes) warns once per file and drowns the real output.  Loading everything
% costs ~0.8 MB per 31P file and sidesteps both problems.
try
    m = load(matFile);
catch ME
    warning('loadPresent: %s (%s)', matFile, ME.message);
    m = struct();
end
% Downstream code reads m.absQ unconditionally; an empty struct makes every
% isfield/concMap check below fall through to NaN instead of erroring.
if ~isfield(m,'absQ'), m.absQ = struct(); end
end

% =====================================================================
function r = fillSet(r, prefix, q)
% Copy the five reported 31P species out of one quant struct (anything with a
% .metabolites array) into the r.<prefix>_mM_* columns.
c = concMap(q);
r.([prefix '_mM_NADplus']) = getOr(c,'NADplus');
r.([prefix '_mM_NADH'])    = getOr(c,'NADH');
r.([prefix '_mM_aATP'])    = getOr(c,'aATP');
r.([prefix '_mM_PCr'])     = getOr(c,'PCr');
r.([prefix '_mM_UDP'])     = sumOr(c, {'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'});
end

% =====================================================================
function c = concMap(absQ)
% name -> conc_mM for whichever metabolite list this struct carries.
c = containers.Map('KeyType','char','ValueType','double');
if ~isstruct(absQ) || ~isfield(absQ,'metabolites'), return; end
M = absQ.metabolites;
for k = 1:numel(M)
    if isfield(M(k),'conc_mM') && isscalar(M(k).conc_mM)
        c(M(k).name) = double(M(k).conc_mM);
    end
end
end

function v = getOr(c, name)
if isKey(c, name), v = c(name); else, v = NaN; end
end

function v = sumOr(c, names)
% Sum of the UDP species; NaN only if none of them were fitted.
vals = [];
for k = 1:numel(names)
    if isKey(c, names{k}), vals(end+1) = c(names{k}); end %#ok<AGROW>
end
if isempty(vals), v = NaN; else, v = sum(vals); end
end

function p = protoInfo(h)
% Protocol identity + timing straight out of the stored Siemens header.
% Everything is best-effort: a missing field just leaves its default.
p = struct('protocol','', 'seqfile','', 'TR_ms',NaN, 'TE_ms',NaN, ...
           'flipDeg',NaN, 'noe','');
if ~isstruct(h), return; end
p.protocol = hstr(h, 'Config', 'ProtocolName');
p.seqfile  = hstr(h, 'Config', 'SequenceFileName');
try, p.TR_ms   = double(h.MeasYaps.alTR{1})/1000; catch, end
try, p.TE_ms   = double(h.MeasYaps.alTE{1})/1000; catch, end
try, p.flipDeg = double(h.MeasYaps.adFlipAngleDegree{1}); catch, end
% alFree(5:6) are the NOE pulse settings; they track the hann->rect change.
try
    a = h.MeasYaps.sWipMemBlock.alFree;
    p.noe = sprintf('%g/%g', double(a{5}), double(a{6}));
catch
end
end

function s = hstr(h, a, b)
s = '';
try
    v = h.(a).(b);
    if ischar(v), s = strtrim(v); end
catch
end
end

% =====================================================================
function [tr, te] = seqTiming(absQ)
tr = NaN; te = NaN;
if isfield(absQ,'sequence')
    if isfield(absQ.sequence,'TR_ms'), tr = absQ.sequence.TR_ms; end
    if isfield(absQ.sequence,'TE_ms'), te = absQ.sequence.TE_ms; end
end
if isfield(absQ,'TR_ms'), tr = absQ.TR_ms; end
if isfield(absQ,'TE_ms'), te = absQ.TE_ms; end
end

% =====================================================================
function [subj, dstr] = parseScanID(scanID)
% "S095_07242026_31P-1H-HeadCoil_external_ref" -> "S095",       "2026-07-24"
% "859875_S01_09012026_brain"                  -> "859875_S01", "2026-09-01"
subj = ''; dstr = '';
tok = regexp(scanID, '^((?:\d+_)?S\d+)_(\d{2})(\d{2})(\d{4})', 'tokens', 'once');
if ~isempty(tok)
    subj = tok{1};
    dstr = sprintf('%s-%s-%s', tok{4}, tok{2}, tok{3});   % yyyy-mm-dd
else
    p = strsplit(scanID, '_');
    subj = p{1};
end
end

% =====================================================================
function writeTsv(rows, outPath, isBrain)
if isempty(rows)
    fprintf('No rows for %s -- not written.\n', outPath);
    return
end
T = struct2table([rows{:}], 'AsArray', true);

% column order
lead = {'session','subject','scandate','anatomy','cohort'};
oneH = {'uM_NADH2','uM_NADH6','uM_NADH4','uM_NAD_1H_mean','uM_NAD_H2H6_mean', ...
        'uM_Trp','uM_NR'};
% three reference sets, kept as three adjacent blocks so they read together
spec = {'NADplus','NADH','UDP','aATP','PCr'};
qAA  = strcat('aATPref_mM_', spec);
qPCr = strcat('PCrref_mM_',  spec);
qExt = strcat('extref_mM_',  spec);
refs = {'primary_ref_name','ref_aATP_conc_mM','ref_PCr_conc_mM'};
ext  = {'RefExt_area','extRef_ok','extRef_refConc_mM','extRef_V_ref_mL', ...
        'extRef_V_M_mL','extRef_applyB1','extRef_reason'};
seg  = {'f_GM','f_WM','f_CSF'};
tail = {'TR_1H_ms','TE_1H_ms','flip_1H_deg','TR_31P_ms','TE_31P_ms', ...
        'flip_31P_deg','noe_31P','prot_1H','seq_1H','prot_31P','seq_31P'};

if isBrain
    order = [lead oneH qAA qPCr qExt refs ext seg tail];
else
    order = [lead oneH qAA qPCr qExt refs ext tail];
end
order = order(ismember(order, T.Properties.VariableNames));
T = T(:, order);
T = sortrows(T, {'subject','scandate'});

writetable(T, outPath, 'FileType','text', 'Delimiter','\t');
fprintf('Wrote %s  (%d rows, %d cols)\n', outPath, height(T), width(T));
end
