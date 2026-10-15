function repeatabilityU01_2026(inDir)
% REPEATABILITYU01_2026  Test-retest repeatability of the 2026 U01 measurements.
%
%   repeatabilityU01_2026()
%   repeatabilityU01_2026(inDir)
%
% Uses only the sessions marked cohort=='repeat' (those under Repeatability\),
% which are same-subject repeat visits a few days apart.  For every reported
% quantity it computes the standard test-retest statistics from the paired
% differences:
%
%   wsSD  within-subject SD  = sqrt( sum(d_i^2) / (2*n) )   [Bland-Altman]
%   wsCV% = 100 * wsSD / grand mean      <- THE repeatability number
%   RC    repeatability coeff = 2.77 * wsSD
%         (2 measurements on one subject differ by less than RC 95% of the time)
%   bsCV% between-subject CV of the subject means
%   D/W   bsCV / wsCV -- discrimination ratio.  Must exceed 1 for the
%         measurement to separate subjects at all; >2 is genuinely useful.
%
% The group CVs printed by summarizeU01_2026 are NOT repeatability: they mix
% within- and between-subject variance.  These are the ones to quote.
%
% Output: <inDir>\U01_2026_repeatability.tsv, plus phone-readable tables.

if nargin < 1 || isempty(inDir)
    inDir = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026';
end

files = {fullfile(inDir,'U01_2026_brain.tsv'), fullfile(inDir,'U01_2026_calf.tsv')};
parts = {'brain','calf'};

% ---- quantities to assess ---------------------------------------------
qty = {};
qty = [qty; makeQ({'uM_NADH2','uM_NADH6','uM_NADH4','uM_NAD_1H_mean','uM_Trp'}, ...
                  {'NADH2','NADH6','NADH4','NAD mean','Trp'}, '1H [uM]')];
spec  = {'NADplus','NADH','UDP','aATP','PCr'};
for r = {'aATPref','PCrref','extref'}
    lbl = strrep(strrep(r{1},'ref','-ref'),'aATP-ref','aATP-ref');
    qty = [qty; makeQ(strcat(r{1},'_mM_',spec), spec, sprintf('31P %s [mM]', lbl))]; %#ok<AGROW>
end
qty = [qty; makeQ({'f_GM','f_WM','f_CSF'}, {'GM','WM','CSF'}, 'Voxel frac')];

rows = {};
report = struct();

for p = 1:numel(parts)
    part = parts{p};
    if ~isfile(files{p}), report.(part) = []; continue; end
    T = readtable(files{p},'FileType','text','Delimiter','\t','VariableNamingRule','preserve');
    if ~ismember('cohort', T.Properties.VariableNames)
        error(['%s has no cohort column -- re-run collectU01_2026_TSV to add it.'], files{p});
    end
    T = T(strcmp(string(T.cohort),'repeat'), :);
    report.(part).T = T;
    if isempty(T), continue; end

    subs = unique(string(T.subject));
    report.(part).subs = subs;
    for q = 1:size(qty,1)
        col = qty{q,1};
        if ~ismember(col, T.Properties.VariableNames), continue; end
        v = double(T.(col));
        [st, pairs] = pairStats(v, string(T.subject), subs);
        report.(part).stats.(matlab.lang.makeValidName(col)) = st;
        report.(part).pairs.(matlab.lang.makeValidName(col)) = pairs;
        rows{end+1} = {part, qty{q,3}, qty{q,2}, st.nPairs, st.mean, ...
                       st.wsSD, st.wsCV, st.RC, st.bsCV, st.DW}; %#ok<AGROW>
    end
end

if isempty(rows)
    fprintf('No repeatability sessions found (cohort=="repeat").\n'); return
end

S = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'bodypart','block','quantity','n_pairs','mean','wsSD','wsCV_pct', ...
     'RC','bsCV_pct','discrim_ratio'});
outPath = fullfile(inDir,'U01_2026_repeatability.tsv');
writetable(S, outPath, 'FileType','text','Delimiter','\t');

% =====================================================================
% Print
% =====================================================================
W = 46;
for p = 1:numel(parts)
    part = parts{p};
    if ~isfield(report,part) || isempty(report.(part)) || isempty(report.(part).T)
        fprintf('\n%s: no repeatability sessions.\n', upper(part));
        continue
    end
    T = report.(part).T;
    subs = report.(part).subs;
    fprintf('\n%s\n%s\n', ['REPEATABILITY: ' upper(part)], repmat('=',1,20));
    fprintf('%d scans, %d subjects x 2 visits\n', height(T), numel(subs));

    lastBlock = '';
    for q = 1:size(qty,1)
        col = qty{q,1};
        key = matlab.lang.makeValidName(col);
        if ~isfield(report.(part),'stats') || ~isfield(report.(part).stats, key), continue; end
        st = report.(part).stats.(key);
        if ~strcmp(qty{q,3}, lastBlock)
            lastBlock = qty{q,3};
            fprintf('\n%s\n', lastBlock);
            fprintf('%s\n', repmat('-',1,W));
            fprintf('%-10s %8s %7s %8s %6s\n','', 'mean','wsCV%','RC','D/W');
            fprintf('%s\n', repmat('-',1,W));
        end
        fprintf('%-10s %8s %7s %8s %6s\n', qty{q,2}, ...
            fmtNum(st.mean), fmtNum(st.wsCV), fmtNum(st.RC), fmtNum(st.DW));
    end
    fprintf('%s\n', repmat('-',1,W));
end

fprintf('\nwsCV%% = within-subject CV (repeatability).\n');
fprintf('RC = 2.77*wsSD; repeat scans differ by\n');
fprintf('     less than RC 95%% of the time.\n');
fprintf('D/W = between/within CV. >1 to tell\n');
fprintf('      subjects apart, >2 to be useful.\n');
fprintf('\nRepeatability TSV: %s\n', outPath);
end

% =====================================================================
function q = makeQ(cols, labels, block)
q = cell(numel(cols),3);
for i = 1:numel(cols)
    q{i,1} = cols{i}; q{i,2} = labels{i}; q{i,3} = block;
end
end

% =====================================================================
function [st, pairs] = pairStats(v, subj, subs)
% Paired test-retest statistics. Subjects contribute only if BOTH visits
% produced a finite value for this quantity.
d = []; m = []; pairs = [];
for s = 1:numel(subs)
    idx = find(subj == subs(s));
    x = v(idx);
    x = x(isfinite(x));
    if numel(x) < 2, continue; end
    x = x(1:2);                       % first two visits
    d(end+1)  = x(2) - x(1);          %#ok<AGROW>
    m(end+1)  = mean(x);              %#ok<AGROW>
    pairs = [pairs; struct('subject',subs(s),'v1',x(1),'v2',x(2), ...
             'diff',x(2)-x(1),'pctDiff',100*(x(2)-x(1))/mean(x))]; %#ok<AGROW>
end

st = struct('nPairs',numel(d),'mean',NaN,'wsSD',NaN,'wsCV',NaN, ...
            'RC',NaN,'bsCV',NaN,'DW',NaN);
if isempty(d), return; end

st.mean = mean(m);
% Bland-Altman within-subject SD from paired differences
st.wsSD = sqrt(sum(d.^2) / (2*numel(d)));
st.wsCV = 100 * st.wsSD / abs(st.mean);
st.RC   = 2.77 * st.wsSD;
if numel(m) > 1
    st.bsCV = 100 * std(m) / abs(mean(m));
    if st.wsCV > 0, st.DW = st.bsCV / st.wsCV; end
end
end

% =====================================================================
function s = fmtNum(x)
if ~isfinite(x),        s = '-';
elseif abs(x) >= 1000,  s = sprintf('%.0f', x);
elseif abs(x) >= 100,   s = sprintf('%.0f', x);
elseif abs(x) >= 10,    s = sprintf('%.1f', x);
elseif abs(x) >= 1,     s = sprintf('%.2f', x);
else,                   s = sprintf('%.3f', x);
end
end
