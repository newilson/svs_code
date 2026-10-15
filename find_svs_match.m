function find_svs_match(targetRAS, targetExt)
% FIND_SVS_MATCH  Locate Siemens SVS/ISIS/CSI series whose voxel/slab matches.
%   targetRAS: 1x3 in NIfTI/RAS+ mm. We convert to DICOM/LPS for comparison.
%   targetExt: 1x3 expected slab/voxel size (Sag,Cor,Tra) in mm, optional.

if nargin < 1 || isempty(targetRAS), targetRAS = [-2.52, 9.43, 6.31]; end  % mask centroid
if nargin < 2, targetExt = [200, 180, 75]; end                              % slab extent
targetLPS = [-targetRAS(1), -targetRAS(2), targetRAS(3)];

roots = { ...
    'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_RR\test\2025', ...
    'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_RR\test\2026'};

fprintf('Target center (RAS) = [%+7.3f %+7.3f %+7.3f] mm  -> LPS [%+7.3f %+7.3f %+7.3f]\n', targetRAS, targetLPS);
fprintf('Target slab extent (Sag Cor Tra) = [%g %g %g] mm\n\n', targetExt);

rows = {};   % {dist, sessionRel, lps, ext, source, seriesName}
for r = 1:numel(roots)
    if ~isfolder(roots{r}), continue; end
    sessions = dir(roots{r});
    sessions = sessions([sessions.isdir] & ~startsWith({sessions.name},'.'));
    for s = 1:numel(sessions)
        sdir = fullfile(roots{r}, sessions(s).name);
        seriesList = pickAllLocalizedSeries(sdir);
        for q = 1:numel(seriesList)
            svs1 = seriesList{q};
            dcm  = pickFirstDicom(svs1);
            if isempty(dcm), continue; end
            try
                inf = dicominfo(dcm, 'UseVRHeuristic', false);
            catch
                continue
            end
            [lps, ext, src] = extractVoxel(inf);
            if isempty(lps)
                rows(end+1,:) = {NaN, sessions(s).name, [NaN NaN NaN], [NaN NaN NaN], 'NO POSITION', dirOf(svs1)}; %#ok<AGROW>
                continue
            end
            d = norm(lps - targetLPS);
            rows(end+1,:) = {d, sessions(s).name, lps, ext, src, dirOf(svs1)}; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    fprintf('No SVSSEL* series found under either root.\n'); return
end

% Composite score: center distance + extent mismatch (penalize big extent diffs)
d = cell2mat(rows(:,1));
ext = cell2mat(rows(:,4));
extDiff = sum(abs(ext - targetExt), 2);
score = d + 0.2*extDiff;          % center dominates; extent breaks ties
score(isnan(d)) = Inf;
[~, ord] = sort(score, 'ascend');
rows = rows(ord,:);

fprintf('%-55s  %-8s  %-23s  %-22s  %s\n', 'session', 'd(mm)', 'center LPS', 'slab Sag/Cor/Tra', 'series');
fprintf('%s\n', repmat('-',1,170));
for k = 1:min(size(rows,1), 60)
    if isnan(rows{k,1})
        fprintf('%-55s  %8s  %-23s  %-22s  %s\n', ...
            char(rows{k,2}), '---', '---', '---', char(rows{k,6}));
    else
        e = rows{k,4};
        fprintf('%-55s  %8.3f  [%+6.2f %+6.2f %+6.2f]  [%5.1f %5.1f %5.1f]  %s\n', ...
            char(rows{k,2}), rows{k,1}, rows{k,3}, e(1), e(2), e(3), char(rows{k,6}));
    end
end

best = rows{1,1};
if ~isnan(best) && best < 0.5
    fprintf('\nEXACT MATCH (d < 0.5 mm): %s :: %s\n', char(rows{1,2}), char(rows{1,6}));
elseif ~isnan(best) && best < 5
    fprintf('\nClosest (within 5 mm): %s :: %s (d=%.3f mm)\n', char(rows{1,2}), char(rows{1,6}), best);
else
    fprintf('\nNo series within 5 mm; closest is %s :: %s at %.3f mm.\n', char(rows{1,2}), char(rows{1,6}), best);
end
end

% ---------------- helpers ----------------

function out = pickAllLocalizedSeries(sessionDir)
% Pick one representative folder per *kind* of localized acquisition.
% Within a session, all SVSSEL_* should share the same VoI; ditto ISIS_*,
% CSI_*, etc. So we keep one folder per leading-prefix.
patterns = {'SVSSEL_*','ISIS_*','CSI_*','MRS_*','SVS_*'};
candidates = {};
% direct children
for p = 1:numel(patterns)
    hits = dir(fullfile(sessionDir, patterns{p}));
    hits = hits([hits.isdir]);
    for h = 1:numel(hits), candidates{end+1} = fullfile(sessionDir, hits(h).name); end %#ok<AGROW>
end
% one level deeper
sub = dir(sessionDir);
sub = sub([sub.isdir] & ~startsWith({sub.name},'.'));
for k = 1:numel(sub)
    for p = 1:numel(patterns)
        hits = dir(fullfile(sessionDir, sub(k).name, patterns{p}));
        hits = hits([hits.isdir]);
        for h = 1:numel(hits)
            candidates{end+1} = fullfile(sessionDir, sub(k).name, hits(h).name); %#ok<AGROW>
        end
    end
end
% dedupe by leading prefix (first token before underscore-number suffix)
seen = containers.Map();
out = {};
for k = 1:numel(candidates)
    [~, nm] = fileparts(candidates{k});
    pref = regexprep(nm, '_\d+$','');          % strip trailing _0007
    pref = regexp(pref, '^[^_]+', 'match','once'); % keep prefix (SVSSEL, ISIS, ...)
    if ~seen.isKey(pref)
        seen(pref) = true;
        out{end+1} = candidates{k}; %#ok<AGROW>
    end
end
end

function f = pickFirstDicom(seriesDir)
f = '';
cand = [dir(fullfile(seriesDir,'*.IMA')); dir(fullfile(seriesDir,'*.dcm')); ...
        dir(fullfile(seriesDir,'*.DCM'))];
if ~isempty(cand), f = fullfile(seriesDir, cand(1).name); return; end
% No standard extension: take first non-dot file
cand = dir(seriesDir);
cand = cand(~[cand.isdir] & ~startsWith({cand.name},'.'));
if ~isempty(cand), f = fullfile(seriesDir, cand(1).name); end
end

function n = dirOf(p)
[~, n] = fileparts(p);
end

function [lps, ext, source] = extractVoxel(inf)
% Return voxel/slab center LPS and Sag/Cor/Tra extent (mm). [] if missing.
lps = []; ext = [NaN NaN NaN]; source = '';

for tag = {'Private_0029_1120','Private_0029_1020'}
    if isfield(inf,tag{1})
        ascii = char(inf.(tag{1})(:)');
        m = regexp(ascii, '### ASCCONV BEGIN.*?### ASCCONV END ###','match','once');
        if isempty(m), m = ascii; end
        sag = parseAsccDouble(m,'sSpecPara.sVoI.sPosition.dSag');
        cor = parseAsccDouble(m,'sSpecPara.sVoI.sPosition.dCor');
        tra = parseAsccDouble(m,'sSpecPara.sVoI.sPosition.dTra');
        thk = parseAsccDouble(m,'sSpecPara.sVoI.dThickness');
        ph  = parseAsccDouble(m,'sSpecPara.sVoI.dPhaseFOV');
        ro  = parseAsccDouble(m,'sSpecPara.sVoI.dReadoutFOV');
        if any(~isnan([sag cor tra thk ph ro]))
            if isnan(sag), sag = 0; end
            if isnan(cor), cor = 0; end
            if isnan(tra), tra = 0; end
            lps = [sag, cor, tra];
            % Siemens ASCCONV reports VoI in patient axes only when
            % unrotated; we just report the three numbers in their stored
            % order (sag,cor,tra) regardless of normal orientation.
            ext = [ro, ph, thk];   % readout=sag, phase=cor, thickness=tra (approx)
            source = 'ASCCONV'; return
        end
    end
end

if isfield(inf,'VolumeLocalizationSequence') && ~isempty(inf.VolumeLocalizationSequence)
    s = inf.VolumeLocalizationSequence;
    flds = fieldnames(s);
    for j = 1:numel(flds)
        item = s.(flds{j});
        if isfield(item,'MidSlabPosition')
            lps = double(item.MidSlabPosition(:)'); source = 'VLS.MidSlabPos'; return
        end
    end
end

if isfield(inf,'ImagePositionPatient') && numel(inf.ImagePositionPatient)==3
    lps = double(inf.ImagePositionPatient(:)'); source = 'ImagePosPat'; return
end
end

function v = parseAsccDouble(txt, key)
% Find "<key> = <number>" in an ASCCONV block. Returns NaN if absent.
% Key contains dots that must be escaped for regex.
keyEsc = regexprep(key,'\.','\\.');
tok = regexp(txt, [keyEsc '\s*=\s*([-+0-9.eE]+)'], 'tokens', 'once');
if isempty(tok), v = NaN; else, v = str2double(tok{1}); end
end
