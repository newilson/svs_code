function R = compare_hsvd_varpro(opt)
%COMPARE_HSVD_VARPRO  HSVD vs VARPRO fit comparison over every 1H data set.
%
%   compare_hsvd_varpro
%   compare_hsvd_varpro('cohorts',{'u01','calf'})
%   compare_hsvd_varpro('run',false)             % report from existing outputs
%   R = compare_hsvd_varpro(...)                 % per-cohort matrices
%
% Runs each session through its production driver (processU01_Brain_1H /
% processU01_Calf_1H), then compares, per metabolite peak and for the water
% reference, the coarse HSVD fit saved in output.met.fit.hsvd /
% output.wat.fit.hsvd against the VARPRO fit in the same output.
%
% Both sides of the concentration are HSVD: the metabolite area is the HSVD
% peak area and the water area is the HSVD water area.  dat2mat_svsNAD does not
% compute HSVD concentrations and is not asked to.  absoluteQuant_svs is linear
% in (metabolite area / water area) and the segmentation, T1, T2 and sequence
% timing are identical for the two fits, so
%
%   conc_HSVD = conc_VARPRO * (areaHSVD_met / areaVARPRO_met)
%                           / (areaHSVD_wat / areaVARPRO_wat)
%
% Options
%   outDir   results folder (default <repo>\_hsvd_varpro)
%   cohorts  any of 'u01', 'pilotsinc', 'dfs', 'calf'
%   run      process sessions whose output is missing or predates the HSVD
%            fields (default true)
%   force    reprocess every session (default false)
%   plot     write one PNG per cohort (default true)
%
% Writes <outDir>\hsvd_vs_varpro_peaks.tsv, _water.tsv, _report.txt

arguments
    opt.outDir  (1,:) char    = fullfile(fileparts(mfilename('fullpath')), '_hsvd_varpro')
    opt.cohorts (1,:) cell    = {'u01','pilotsinc','dfs','calf'}
    opt.run     (1,1) logical = true
    opt.force   (1,1) logical = false
    opt.plot    (1,1) logical = true
end

REPO = fileparts(mfilename('fullpath'));
addpath(REPO);
if ~isfolder(opt.outDir), mkdir(opt.outDir); end

repFile = fullfile(opt.outDir,'hsvd_vs_varpro_report.txt');
say();                                   % reset the report buffer

say('HSVD vs VARPRO comparison -- %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
say('repo    : %s\n', REPO);
say('results : %s\n', opt.outDir);

R = struct([]);
for c = 1:numel(opt.cohorts)
    spec = cohortSpec(opt.cohorts{c}, opt.outDir, REPO);
    if isempty(spec)
        say('\nSKIP %s: data root not found\n', opt.cohorts{c});
        continue
    end
    say('\n%s: %d sessions\n', spec.label, numel(spec.dirs));

    cohDir = fullfile(opt.outDir, lower(opt.cohorts{c}));
    if ~isfolder(cohDir), mkdir(cohDir); end

    S = {};
    for k = 1:numel(spec.dirs)
        tag = spec.tags{k};
        mat = fullfile(cohDir, sprintf('%s_%s_1H_abs.mat', tag, spec.kind));
        if opt.run && (opt.force || ~hasHsvd(mat))
            fprintf('\n[%s %d/%d] %s\n', spec.label, k, numel(spec.dirs), tag);
            try
                if strcmp(spec.kind,'calf')
                    spec.driver(spec.dirs{k}, cohDir, false);
                else
                    spec.driver(spec.dirs{k}, cohDir, spec.seg(spec.dirs{k}), false, tag);
                end
            catch ME
                say('  FAILED %s: %s\n', tag, ME.message);
                continue
            end
        end
        if ~isfile(mat)
            say('  no output for %s\n', tag);
            continue
        end
        if ~hasHsvd(mat)
            say('  %s has no HSVD fields (reprocess with run=true)\n', tag);
            continue
        end
        S{end+1} = readSession(mat, spec.label, tag); %#ok<AGROW>
    end
    if isempty(S), continue; end

    D = stack(S);
    D.label = spec.label;
    D.pairs = spec.pairs;
    reportCohort(D);
    if opt.plot, plotCohort(D, fullfile(opt.outDir, sprintf('hsvd_vs_varpro_%s.png', lower(opt.cohorts{c})))); end
    if isempty(R), R = D; else, R(end+1) = D; end %#ok<AGROW>
end

if ~isempty(R)
    writePeakTsv(fullfile(opt.outDir,'hsvd_vs_varpro_peaks.tsv'), R);
    writeWaterTsv(fullfile(opt.outDir,'hsvd_vs_varpro_water.tsv'), R);
    say('\nTSV  -> %s\n', fullfile(opt.outDir,'hsvd_vs_varpro_peaks.tsv'));
    say('TSV  -> %s\n', fullfile(opt.outDir,'hsvd_vs_varpro_water.tsv'));
end
say('report -> %s\n', repFile);
f = fopen(repFile,'w');
assert(f > 0, 'cannot write %s', repFile);
fwrite(f, say());
fclose(f);
end

% =====================================================================
% cohorts
% =====================================================================
function spec = cohortSpec(name, outDir, repo)
spec = [];
switch lower(name)
    case 'u01'
        % S007_07222026 excluded (bad scan)
        root = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026\2026';
        rep  = fullfile(root,'Repeatability');
        if ~isfolder(root), return; end
        spec.label  = 'U01 brain';
        spec.kind   = 'brain';
        spec.dirs   = { ...
            fullfile(root,'S007_05182026_31P-1H-HeadCoil'), ...
            fullfile(root,'S092_07082026_31P-1H-HeadCoil_external_d2o'), ...
            fullfile(root,'S092_07242026_31P-1H-HeadCoil_external_ref'), ...
            fullfile(root,'S095_06092026_31P-1H-HeadCoil_external_d2o'), ...
            fullfile(root,'S095_07012026_31P-1H-HeadCoil_external_d2o'), ...
            fullfile(root,'S095_07242026_31P-1H-HeadCoil_external_ref'), ...
            fullfile(rep, 'S092_03032026_31P-1H-HeadCoil'), ...
            fullfile(rep, 'S092_03052026_31P-1H-HeadCoil'), ...
            fullfile(rep, 'S138_04092026_31P-1H-HeadCoil'), ...
            fullfile(rep, 'S138_04102026_31P-1H-HeadCoil')};
        spec.pairs  = {'S092_03032026_31P-1H-HeadCoil','S092_03052026_31P-1H-HeadCoil'
                       'S138_04092026_31P-1H-HeadCoil','S138_04102026_31P-1H-HeadCoil'};
        spec.driver = @processU01_Brain_1H;
        spec.seg    = @u01Seg;

    case 'pilotsinc'
        root = 'Z:\Data\Spectroscopy\Pilot-Sinc-DFMRS';
        if ~isfolder(root), return; end
        tags = {'S004_01092025','S004_01142025','S092_01072025','S092_01162025', ...
                'S139_01092025','S139_01162025','S145_01082025','S145_02042025'};
        spec.label  = 'Pilot-Sinc';
        spec.kind   = 'brain';
        spec.dirs   = fullfile(root, tags);
        spec.pairs  = {'S004_01092025','S004_01142025'
                       'S092_01072025','S092_01162025'
                       'S139_01092025','S139_01162025'
                       'S145_01082025','S145_02042025'};
        spec.driver = str2func(pilotSincDriver(outDir, repo));
        segTsv      = fixedSeg(outDir);
        spec.seg    = @(~) segTsv;

    case 'dfs'
        root = 'C:\Users\CAMIPM-NW\Downloads\DFS_CONTROL_STUDY (1)\DFS_CONTROL_STUDY';
        if ~isfolder(root), return; end
        tags = {'S002','S008','S009','S075','S092','S095','S138','S145', ...
                'S167','S169','S170','S179'};
        spec.label  = 'DFS control';
        spec.kind   = 'brain';
        spec.dirs   = fullfile(root, tags);
        spec.pairs  = {};
        spec.driver = @processU01_Brain_1H;
        spec.seg    = @(d) dfsSeg(d, outDir);

    case 'calf'
        root = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026\2026';
        if ~isfolder(root), return; end
        f = dir(fullfile(root,'**','processed','*_calf_1H_abs.mat'));
        if isempty(f), return; end
        spec.label  = 'U01 calf';
        spec.kind   = 'calf';
        spec.dirs   = unique(cellfun(@fileparts, {f.folder}, 'UniformOutput', false));
        spec.pairs  = {};
        spec.driver = @processU01_Calf_1H;
        spec.seg    = [];

    otherwise
        error('unknown cohort: %s', name);
end
spec.tags = cell(1, numel(spec.dirs));
for k = 1:numel(spec.dirs)
    [~, spec.tags{k}] = fileparts(spec.dirs{k});
end
end

function st = u01Seg(sdir)
d = dir(fullfile(sdir,'processed','*_seg_svsmask.tsv'));
assert(isscalar(d), 'expected one seg TSV in %s', sdir);
st = fullfile(d(1).folder, d(1).name);
end

function st = fixedSeg(outDir)
% Pilot-Sinc has no segmentation; WM 50 / GM 40 / CSF 10 as in the ridge and
% lambdaWidth runs.
sd = fullfile(outDir,'seg');
if ~isfolder(sd), mkdir(sd); end
st = fullfile(sd,'fixed_WM50_GM40_CSF10_seg.tsv');
f = fopen(st,'w');
fprintf(f,'labelval\tvoxels\tvolume(mm^3)\tlabelname\n');
fprintf(f,'1\t100000\t400000.00\tCSF\n2\t400000\t1600000.00\tGM\n3\t500000\t2000000.00\tWM\n');
fclose(f);
end

function st = dfsSeg(sdir, outDir)
% Tissue fractions come from the session's existing baseline output.
[~, tag] = fileparts(sdir);
fb = dir(fullfile(sdir,'*_brain_1H_abs.mat'));
fb = fb(~ismember(lower({fb.name}), ...
     lower({'second_brain_1H_abs.mat','second5_brain_1H_abs.mat','second6_brain_1H_abs.mat'})));
assert(isscalar(fb), '%s: expected one baseline mat, found %d', tag, numel(fb));
L = load(fullfile(sdir, fb(1).name), 'segFracs');
vx = round(double(L.segFracs(:)).' * 1e9);
sd = fullfile(outDir,'seg');
if ~isfolder(sd), mkdir(sd); end
st = fullfile(sd, [tag '_dfs_seg.tsv']);
f = fopen(st,'w');
fprintf(f,'labelval\tvoxels\tvolume(mm^3)\tlabelname\n');
nm = {'CSF','GM','WM'};
for q = 1:3, fprintf(f,'%d\t%d\t%.2f\t%s\n', q, vx(q), vx(q)*4, nm{q}); end
fclose(f);
end

function fn = pilotSincDriver(outDir, repo)
% processU01_Brain_1H with the Pilot-Sinc .dat pair pinned: its picker reads
% '*WATERREF_2x4_7ppm*' as water-suppressed and takes the highest MID, which is
% the 3X10 voxel.  Fit settings are the driver's own.
fn  = 'hv_drv_pilotsinc';
shd = fullfile(outDir,'shadow');
if ~isfolder(shd), mkdir(shd); end
line = '[wsDat, watDat] = pick_latest_pair(datDir);';
pick = [ ...
    'dws = dir(fullfile(datDir, ''*2X9_7ppm*.dat''));'                      newline ...
    'dwr = dir(fullfile(datDir, ''*WATERREF_2x4_7ppm*.dat''));'             newline ...
    'assert(~isempty(dws), ''no *2X9_7ppm* .dat in %s'', datDir);'          newline ...
    'assert(~isempty(dwr), ''no *WATERREF_2x4_7ppm* .dat in %s'', datDir);' newline ...
    'wsDat  = fullfile(datDir, dws(1).name);'                               newline ...
    'watDat = fullfile(datDir, dwr(1).name);'];
src = fileread(fullfile(repo,'processU01_Brain_1H.m'));
t = regexprep(src, 'function\s+processU01_Brain_1H\(', ['function ' fn '('], 'once');
t = strrep(t, "fileparts(mfilename('fullpath'))", ['''' repo '''']);
assert(~contains(t, "mfilename('fullpath')"), 'path substitution incomplete');
assert(contains(t, line), 'picker line not found in processU01_Brain_1H');
t = strrep(t, line, pick);
f = fopen(fullfile(shd,[fn '.m']),'w');  fwrite(f, t);  fclose(f);
addpath(shd);
rehash path;
end

% =====================================================================
% read one processed session
% =====================================================================
function tf = hasHsvd(mat)
tf = false;
if ~isfile(mat), return; end
try
    L = load(mat,'output');
catch
    return
end
tf = isfield(L.output.met.fit,'hsvd') && isfield(L.output.met.fit.hsvd,'area') && ...
     isfield(L.output.wat.fit,'hsvd');
end

function S = readSession(mat, cohort, tag)
L = load(mat,'output','absQ');
o = L.output;  F = o.met.fit;  W = o.wat;
S.cohort = cohort;  S.tag = tag;  S.names = F.hsvd.names(:)';
np = numel(S.names);
f0 = abs((o.met.hz(2) - o.met.hz(1)) / (o.met.ppm(2) - o.met.ppm(1)));

% ---- VARPRO (mode 5 layout: amp pos widthL ph widthG width) ----------
S.areaV = F.areasConv(:)';
[S.ctrV, S.wV, S.etaV, S.phV] = deal(nan(1,np));
if isfield(F,'region')
    for r = 1:numel(F.region)
        Rg  = F.region(r);  sel = Rg.peakSelect;  npr = numel(sel);
        p   = Rg.fit.pars;
        if numel(p) < 6*npr, continue; end
        wt  = p(5*npr + (1:npr)) * f0;
        wg  = abs(p(4*npr + (1:npr))) * f0;
        S.ctrV(sel) = p(npr + (1:npr));
        S.phV(sel)  = p(3*npr + (1:npr));
        S.wV(sel)   = wt;
        S.etaV(sel) = wg ./ max(wt, eps);
    end
end

% ---- coarse-phase QC (absent in outputs written before it existed) ----
S.qcEnabled = false;  S.qcFailed = false;  S.qcNegFrac = NaN;
if isfield(F,'coarseQC')
    S.qcEnabled = logical(F.coarseQC.enabled);
    S.qcFailed  = logical(F.coarseQC.failed);
    S.qcNegFrac = max(F.coarseQC.negFrac);
end

% ---- HSVD ------------------------------------------------------------
H = F.hsvd;
S.foundH = logical(H.found(:)');
S.nCompH = H.nComp(:)';
S.areaH  = H.area(:)';
S.ctrH   = H.ppm(:)';
S.wH     = H.widthHz(:)';
S.phH    = H.phaseDegFit(:)';
S.lamH   = H.lambda(:)';

% ---- water -----------------------------------------------------------
hw   = W.fit.hsvd;
f0w  = abs((W.hz(2) - W.hz(1)) / (W.ppm(2) - W.ppm(1)));
S.watAreaV = W.fit.areaTotal;
S.watAreaH = hw.area;
S.watRatio = S.watAreaH / S.watAreaV;
S.watNComp = hw.nComp;
S.watPpmH  = hw.ppm;
S.watWH    = hw.widthHz;
S.watLam   = hw.lambda;
S.watK     = hw.K;
a = hw.compAmp(logical(hw.assigned));
S.watCoh = abs(sum(a)) / max(sum(abs(a)), eps);   % 1 = components add, 0 = cancel
[S.watPpmV, S.watWV] = deal(NaN);
if isfield(W.fit,'pars') && isfield(W.fit,'names') && isscalar(W.fit.names) && numel(W.fit.pars) >= 6
    S.watPpmV = W.fit.pars(2);
    S.watWV   = W.fit.pars(6) * f0w;
end

% ---- concentrations --------------------------------------------------
S.concV = nan(1,np);
for q = 1:numel(L.absQ.metabolites)
    j = find(strcmp(S.names, L.absQ.metabolites(q).name), 1);
    if ~isempty(j), S.concV(j) = L.absQ.metabolites(q).conc_mM * 1e3; end
end
S.ratio = S.areaH ./ S.areaV;
S.concH = S.concV .* S.ratio / S.watRatio;
end

function D = stack(S)
D.names = S{1}.names;
D.tags  = cellfun(@(s) s.tag, S, 'UniformOutput', false);
n  = numel(S);
np = numel(D.names);
per = {'areaV','areaH','ratio','ctrV','ctrH','wV','wH','etaV','phV','phH', ...
       'concV','concH','nCompH','foundH','lamH'};
for f = per
    M = nan(n, np);
    for i = 1:n
        v = S{i}.(f{1});
        if numel(v) == np, M(i,:) = double(v); end
    end
    D.(f{1}) = M;
end
sc = {'watAreaV','watAreaH','watRatio','watNComp','watCoh','watPpmV','watPpmH', ...
      'watWV','watWH','watLam','watK','qcEnabled','qcFailed','qcNegFrac'};
for f = sc
    v = nan(n,1);
    for i = 1:n, v(i) = double(S{i}.(f{1})); end
    D.(f{1}) = v;
end
end

% =====================================================================
% report
% =====================================================================
function reportCohort(D)
nm = D.names;  n = numel(D.tags);  np = numel(nm);
say('\n\n#################### %s  (n = %d) ####################\n', D.label, n);

say('\n--- A. HSVD component assignment ---\n');
head(nm, 'peak');
row('peaks found (%)',             100*mean(D.foundH == 1, 1), '%9.0f');
row('components per found peak',   sum(D.nCompH,1) ./ max(sum(D.foundH == 1,1),1), '%9.2f');
say('  ridge lambda (metabolites): median %.5f    water: K %d, lambda median %.5f\n', ...
    median(D.lamH(:),'omitnan'), round(median(D.watK)), median(D.watLam,'omitnan'));

say('\n--- B. peak area, HSVD / VARPRO (areasConv) ---\n');
head(nm, 'peak');
row('median', median(D.ratio,1,'omitnan'), '%9.3f');
row('  IQR low',  prctile(D.ratio,25,1), '%9.3f');
row('  IQR high', prctile(D.ratio,75,1), '%9.3f');
row('group CV% of the ratio', cv_(D.ratio), '%9.1f');

say('\n--- C. line position, width, phase ---\n');
head(nm, 'peak');
row('median HSVD - VARPRO center (ppm)', median(D.ctrH - D.ctrV,1,'omitnan'), '%9.4f');
row('median HSVD width (Hz)',            median(D.wH,1,'omitnan'), '%9.1f');
row('median VARPRO width (Hz)',          median(D.wV,1,'omitnan'), '%9.1f');
row('median HSVD / VARPRO width',        median(D.wH./D.wV,1,'omitnan'), '%9.2f');
row('median VARPRO eta (wG/wT)',         median(D.etaV,1,'omitnan'), '%9.3f');
row('median |HSVD - VARPRO phase| (deg)',median(abs(angdiff_(D.phH, D.phV)),1,'omitnan'), '%9.1f');

say('\n--- D. water reference ---\n');
say('  HSVD / VARPRO area : median %.3f  [%.3f %.3f]\n', median(D.watRatio,'omitnan'), min(D.watRatio), max(D.watRatio));
say('  components         : median %d  [%d %d]\n', round(median(D.watNComp)), min(D.watNComp), max(D.watNComp));
say('  component coherence: median %.2f  [%.2f %.2f]\n', median(D.watCoh,'omitnan'), min(D.watCoh), max(D.watCoh));
say('  HSVD - VARPRO ppm  : median %+.4f\n', median(D.watPpmH - D.watPpmV,'omitnan'));
say('  width HSVD / VARPRO: %.1f / %.1f Hz\n', median(D.watWH,'omitnan'), median(D.watWV,'omitnan'));

say('\n--- E. concentration (uM) ---\n');
head(nm, 'peak');
row('VARPRO median', median(D.concV,1,'omitnan'), '%9.1f');
row('HSVD median',   median(D.concH,1,'omitnan'), '%9.1f');
row('HSVD vs VARPRO, median % diff', median(100*(D.concH./D.concV - 1),1,'omitnan'), '%9.1f');
row('HSVD vs VARPRO, max |% diff|',  max(abs(100*(D.concH./D.concV - 1)),[],1), '%9.1f');

say('\n--- F. group CV%% of concentration (n = %d) ---\n', n);
head(nm, 'peak');
row('VARPRO', cv_(D.concV), '%9.1f');
row('HSVD',   cv_(D.concH), '%9.1f');

if ~isempty(D.pairs)
    P = D.pairs;  npair = size(P,1);
    say('\n--- G. test-retest CV%% (%d pairs) ---\n', npair);
    head(nm, 'peak');
    for m = {'concV','concH'}
        T = nan(npair, np);
        for p = 1:npair
            i1 = find(strcmp(D.tags, P{p,1}), 1);
            i2 = find(strcmp(D.tags, P{p,2}), 1);
            if isempty(i1) || isempty(i2), continue; end
            x = D.(m{1})([i1 i2], :);
            T(p,:) = 100*std(x,0,1) ./ abs(mean(x,1));
        end
        lbl = 'VARPRO';  if strcmp(m{1},'concH'), lbl = 'HSVD'; end
        for p = 1:npair
            row(sprintf('  %s %s', lbl, strtok(P{p,1},'_')), T(p,:), '%9.1f');
        end
        row(sprintf('%s MEAN', lbl), mean(T,1,'omitnan'), '%9.1f');
    end
end

iH2 = find(strcmp(nm,'NADH2'),1);  iH6 = find(strcmp(nm,'NADH6'),1);  iH4 = find(strcmp(nm,'NADH4'),1);
if ~isempty(iH2) && ~isempty(iH6) && ~isempty(iH4)
    say('\n--- H. NAD+ self-consistency (target 1.0) ---\n');
    say('%-8s %18s %18s %8s\n', 'fit', 'H6/H2 mean+-SD', 'H4/H2 mean+-SD', 'MAD');
    for m = {'concV','concH'}
        X = D.(m{1});
        r62 = X(:,iH6)./X(:,iH2);  r42 = X(:,iH4)./X(:,iH2);
        ok  = isfinite(r62) & isfinite(r42);
        lbl = 'VARPRO';  if strcmp(m{1},'concH'), lbl = 'HSVD'; end
        say('%-8s %9.3f +- %5.3f %9.3f +- %5.3f %8.3f\n', lbl, ...
            mean(r62,'omitnan'), std(r62,'omitnan'), mean(r42,'omitnan'), std(r42,'omitnan'), ...
            mean([abs(r62(ok)-1); abs(r42(ok)-1)]));
    end
end

say('\n--- I. per session: VARPRO -> HSVD concentration (uM) ---\n');
say('%-42s', 'session');
for k = 1:np, say(' %15s', nm{k}); end
say(' %7s %6s\n', 'wat H/V', 'coher');
for i = 1:n
    say('%-42s', D.tags{i});
    for k = 1:np, say(' %6.0f->%6.0f', D.concV(i,k), D.concH(i,k)); end
    say(' %7.3f %6.2f\n', D.watRatio(i), D.watCoh(i));
end

say('\n--- J. QC flags ---\n');
nf = 0;
for i = 1:n
    msg = {};
    if ~all(D.foundH(i,:) == 1)
        msg{end+1} = ['HSVD missed ' strjoin(nm(D.foundH(i,:) ~= 1), ',')]; %#ok<AGROW>
    end
    if ~isfinite(D.watRatio(i)) || D.watRatio(i) < 0.7 || D.watRatio(i) > 1.4
        msg{end+1} = sprintf('water area ratio %.2f', D.watRatio(i)); %#ok<AGROW>
    end
    if D.watCoh(i) < 0.05
        msg{end+1} = sprintf('water components nearly cancel (coherence %.2f)', D.watCoh(i)); %#ok<AGROW>
    end
    if D.qcFailed(i) == 1
        msg{end+1} = 'coarse-phase QC failed (fell back to the VARPRO estimator)'; %#ok<AGROW>
    elseif D.qcEnabled(i) ~= 1 && D.qcNegFrac(i) > 0.25
        msg{end+1} = sprintf('phased fit window %.0f%% negative, coarse-phase QC not enabled', ...
                             100*D.qcNegFrac(i)); %#ok<AGROW>
    end
    if ~isempty(msg)
        nf = nf + 1;
        say('  %-42s %s\n', D.tags{i}, strjoin(msg, '; '));
    end
end
if nf == 0, say('  none\n'); end
end

function head(nm, lbl)
say('%-34s', lbl);
for k = 1:numel(nm), say(' %9s', nm{k}); end
say('\n');
end

function row(lbl, v, fmt)
say('%-34s', lbl);
for k = 1:numel(v), say([' ' fmt], v(k)); end
say('\n');
end

function txt = say(fmt, varargin)
% Print to the console and keep the text.  say() with no argument returns
% everything said so far and empties the buffer -- the processing pipeline
% closes file identifiers, so the report cannot be held open across a run.
persistent buf
if isempty(buf), buf = {}; end
if nargin == 0
    txt = strjoin(buf, '');
    buf = {};
    return
end
s = sprintf(fmt, varargin{:});
fprintf('%s', s);
buf{end+1} = s;
end

function g = cv_(M)
g = 100*std(M,0,1,'omitnan') ./ abs(mean(M,1,'omitnan'));
end

function d = angdiff_(a, b)
d = mod(a - b + 180, 360) - 180;
end

% =====================================================================
% tables and figures
% =====================================================================
function writePeakTsv(path, R)
f = fopen(path,'w');
fprintf(f, '%s\n', strjoin({'cohort','session','peak','foundHSVD','nCompHSVD', ...
    'areaVARPRO','areaHSVD','areaRatio','ppmVARPRO','ppmHSVD','widthHzVARPRO','widthHzHSVD', ...
    'etaVARPRO','phaseDegVARPRO','phaseDegHSVD','concVARPRO_uM','concHSVD_uM', ...
    'waterAreaRatio'}, '\t'));
for c = 1:numel(R)
    D = R(c);
    for i = 1:numel(D.tags)
        for k = 1:numel(D.names)
            fprintf(f, '%s\t%s\t%s\t%d\t%d', D.label, D.tags{i}, D.names{k}, ...
                D.foundH(i,k) == 1, D.nCompH(i,k));
            fprintf(f, '\t%.6g', D.areaV(i,k), D.areaH(i,k), D.ratio(i,k), ...
                D.ctrV(i,k), D.ctrH(i,k), D.wV(i,k), D.wH(i,k), D.etaV(i,k), ...
                D.phV(i,k), D.phH(i,k), D.concV(i,k), D.concH(i,k), D.watRatio(i));
            fprintf(f, '\n');
        end
    end
end
fclose(f);
end

function writeWaterTsv(path, R)
f = fopen(path,'w');
fprintf(f, '%s\n', strjoin({'cohort','session','areaVARPRO','areaHSVD','areaRatio', ...
    'nCompHSVD','coherence','ppmVARPRO','ppmHSVD','widthHzVARPRO','widthHzHSVD', ...
    'K','lambda'}, '\t'));
for c = 1:numel(R)
    D = R(c);
    for i = 1:numel(D.tags)
        fprintf(f, '%s\t%s', D.label, D.tags{i});
        fprintf(f, '\t%.6g', D.watAreaV(i), D.watAreaH(i), D.watRatio(i), D.watNComp(i), ...
            D.watCoh(i), D.watPpmV(i), D.watPpmH(i), D.watWV(i), D.watWH(i), ...
            D.watK(i), D.watLam(i));
        fprintf(f, '\n');
    end
end
fclose(f);
end

function plotCohort(D, png)
np = numel(D.names);
fig = figure('Visible','off','Position',[100 100 260*(np+1) 300]);
for k = 1:np
    subplot(1, np+1, k);
    x = D.concV(:,k);  y = D.concH(:,k);
    ok = isfinite(x) & isfinite(y);
    plot(x(ok), y(ok), 'o', 'MarkerFaceColor',[0.2 0.4 0.8], 'MarkerEdgeColor','none'); hold on
    lim = [0 max([x(ok); y(ok); eps])*1.1];
    plot(lim, lim, 'k--');
    xlim(lim); ylim(lim); axis square; grid on
    xlabel('VARPRO (uM)'); if k == 1, ylabel('HSVD (uM)'); end
    title(sprintf('%s  (%.2f)', D.names{k}, median(y(ok)./x(ok))));
end
subplot(1, np+1, np+1);
x = D.watAreaV;  y = D.watAreaH;
ok = isfinite(x) & isfinite(y);
plot(x(ok), y(ok), 'o', 'MarkerFaceColor',[0.8 0.3 0.2], 'MarkerEdgeColor','none'); hold on
lim = [0 max([x(ok); y(ok); eps])*1.1];
plot(lim, lim, 'k--');
xlim(lim); ylim(lim); axis square; grid on
xlabel('VARPRO area'); ylabel('HSVD area');
title(sprintf('water  (%.2f)', median(y(ok)./x(ok))));
sgtitle(sprintf('%s: HSVD vs VARPRO', D.label));
print(fig, png, '-dpng', '-r150');
close(fig);
end
