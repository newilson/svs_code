function S = statsU01_2026_latest(inDir)
% STATSU01_2026_LATEST  Preliminary NAD+/Trp/NADH statistics, latest protocol only.
%
%   S = statsU01_2026_latest()
%   S = statsU01_2026_latest(inDir)
%
% Reads U01_2026_brain.tsv / U01_2026_calf.tsv (written by collectU01_2026_TSV)
% and reports, for brain and calf separately:
%
%   1) group mean / SD / CV% of
%        1H  : NADH2, NADH6, NADH4, their mean, and Trp        [uM]
%        31P : NAD+ and NADH, INTERNAL reference only          [mM]
%      Group CV is given two ways: over all scans, and over subject means.
%      Several subjects were scanned more than once, so the all-scan CV is
%      inflated by repeated sampling of the same people; the subject-mean CV
%      is the closer thing to a between-subject CV.
%   2) within-subject (test-retest) CV for every subject scanned twice.
%
% PROTOCOL FILTER
%   "Latest protocol" is the one the external-reference sessions were run
%   with, identified from the acquisition headers (prot_*/seq_*/TE_31P_ms
%   columns) rather than folder names:
%       31P isis_svs_20260303 1d, TE 0.1 ms, excitation flip 70 deg,
%           rect NOE pulse (alFree 3/200), TR 5 s, 144 avg
%       1H  svssel_20260101, NAD 9.7 ppm, TE 13 ms, 1024 avg, 20 mm slab
%   The previous protocol used TE 0.3 ms / flip 90 deg / hann NOE (alFree
%   9/500) on the 31P side, so TE_31P_ms separates the two families exactly.
%   Presence of the phantom is NOT part of the filter: sessions run with the
%   protocol but with no usable external reference are included.
%
%   The 1H side barely changed across the revision -- only the sequence binary
%   name (svssel_latest -> svssel_20260101) and, in three July sessions, TR
%   (500 -> 520/600 ms).  So the previous-protocol scans still carry usable 1H
%   data even though their 31P is a different acquisition, which is why the
%   repeatability section reports them separately for 1H.
%
% NO EXTERNAL REFERENCING is used anywhere here: brain is referenced to
% internal aATP (2.7 mM) and calf to internal PCr (33 mM), i.e. the primary
% reference each pipeline already uses.

if nargin < 1 || isempty(inDir)
    inDir = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026';
end

TE31_LATEST = 0.1;              % ms
SEQ1H_LATEST = 'svssel_20260101';
REFSET = struct('brain','aATPref', 'calf','PCrref');

parts = {'brain','calf'};
S = struct();

for p = 1:numel(parts)
    part = parts{p};
    f = fullfile(inDir, sprintf('U01_2026_%s.tsv', part));
    if ~isfile(f), warning('missing %s', f); continue; end
    T = readtable(f,'FileType','text','Delimiter','\t','VariableNamingRule','preserve');

    T.latest31P = abs(double(T.TE_31P_ms) - TE31_LATEST) < 1e-9;
    T.latest1H  = contains(string(T.seq_1H), SEQ1H_LATEST);

    S.(part).all      = T;
    S.(part).latest   = T(T.latest31P & T.latest1H, :);
    S.(part).previous = T(~(T.latest31P & T.latest1H), :);
    S.(part).qty      = quantList(REFSET.(part));
end

% =====================================================================
banner('SCAN INVENTORY');
for p = 1:numel(parts)
    part = parts{p};
    if ~isfield(S,part), continue; end
    fprintf('\n%s\n', upper(part));
    printInventory(S.(part).latest,   'LATEST protocol (31P TE 0.1 ms / flip 70 / rect)');
    printInventory(S.(part).previous, 'previous protocol -- EXCLUDED from group stats');
end

banner('GROUP CONCENTRATIONS -- LATEST PROTOCOL ONLY');
fprintf('internal reference only (brain aATP 2.7 mM, calf PCr 33 mM)\n');
for p = 1:numel(parts)
    part = parts{p};
    if ~isfield(S,part) || isempty(S.(part).latest), continue; end
    S.(part).group = printGroup(S.(part).latest, S.(part).qty, upper(part));
end

banner('TEST-RETEST REPEATABILITY');
fprintf(['wsCV%% = 100*wsSD/mean, wsSD = sqrt(sum(d^2)/(2n)) over paired\n' ...
         'visits (Bland-Altman).  RC = 2.77*wsSD.\n']);

for p = 1:numel(parts)
    part = parts{p};
    if ~isfield(S,part), continue; end
    T = S.(part).all;

    % (a) subjects repeated inside the latest-protocol group
    S.(part).retest_latest = printRetest(S.(part).latest, S.(part).qty, ...
        sprintf('%s -- repeats WITHIN the latest protocol', upper(part)));

    % (b) the dedicated Repeatability\ cohort, on whatever protocol it used
    if ismember('cohort', T.Properties.VariableNames)
        R = T(strcmp(string(T.cohort),'repeat'), :);
    else
        R = T([], :);
    end
    S.(part).retest_cohort = printRetest(R, S.(part).qty, ...
        sprintf('%s -- dedicated Repeatability cohort (all protocols)', upper(part)));

    % (c) same cohort, 1H only, restricted to the current 1H protocol.  Its
    %     31P is the previous acquisition so those rows are dropped here.
    S.(part).retest_cohort_1H = printRetest(R(R.latest1H, :), ...
        quantList1H(), sprintf('%s -- Repeatability cohort, 1H at CURRENT 1H protocol', upper(part)));

    % (d) same cohort, 31P only, restricted to the standard 31P acquisition.
    %     S139's second visit was processed from the TR 10 s / 16-avg variant
    %     rather than the TR 5 s / 144-avg one, so that pair is not a like-for-
    %     like 31P comparison; this block drops it.
    if ismember('TR_31P_ms', R.Properties.VariableNames)
        S.(part).retest_cohort_31P = printRetest(R(double(R.TR_31P_ms) == 5000, :), ...
            quantList31P(REFSET.(part)), ...
            sprintf('%s -- Repeatability cohort, 31P at TR 5 s / 144 avg only', upper(part)));
    end
end

end

% =====================================================================
function q = quantList(refset)
% {column, label, unit, block}
q = [quantList1H(); quantList31P(refset)];
end

function q = quantList31P(refset)
q = {
  [refset '_mM_NADplus'],    'NAD+',     'mM', '31P'
  [refset '_mM_NADH'],       'NADH',     'mM', '31P'
  };
end

function q = quantList1H()
q = {
  'uM_NADH2',                'NAD+ H2',  'uM', '1H'
  'uM_NADH6',                'NAD+ H6',  'uM', '1H'
  'uM_NADH4',                'NAD+ H4',  'uM', '1H'
  'uM_NAD_1H_mean',          'NAD+ mean','uM', '1H'
  'uM_Trp',                  'Trp',      'uM', '1H'
  };
end

% =====================================================================
function printInventory(T, title)
fprintf('  %s: n=%d\n', title, height(T));
for i = 1:height(T)
    fprintf('      %-50s %s  1H:%-42s 31P:%s\n', txt(T.session,i), txt(T.scandate,i), ...
            txt(T.prot_1H,i), txt(T.prot_31P,i));
end
end

% =====================================================================
function G = printGroup(T, qty, title)
subj = string(T.subject);
uS   = unique(subj);
fprintf('\n%s   (%d scans, %d subjects)\n', title, height(T), numel(uS));
line(72);
fprintf('%-11s %-4s %4s %10s %10s %7s %10s %7s\n', ...
        '','unit','n','mean','SD','CV%','subjMean','sCV%');
line(72);
G = struct();
lastBlock = '';
for k = 1:size(qty,1)
    col = qty{k,1};
    if ~ismember(col, T.Properties.VariableNames), continue; end
    if ~strcmp(qty{k,4}, lastBlock)
        lastBlock = qty{k,4};
        fprintf('%s\n', lastBlock);
    end
    v = double(T.(col));
    ok = isfinite(v);
    % per-subject means, so a subject scanned 3x does not count 3x
    sm = arrayfun(@(s) mean(v(ok & subj==s),'omitnan'), unique(subj(ok)));
    sm = sm(isfinite(sm));
    st = struct('n',sum(ok), 'mean',mean(v(ok)), 'sd',std(v(ok)), ...
                'cv',100*std(v(ok))/abs(mean(v(ok))), ...
                'nSubj',numel(sm), 'subjMean',mean(sm), ...
                'subjCV',100*std(sm)/abs(mean(sm)));
    G.(matlab.lang.makeValidName(col)) = st;
    fprintf('  %-9s %-4s %4d %10s %10s %7s %10s %7s\n', qty{k,2}, qty{k,3}, ...
        st.n, num(st.mean), num(st.sd), num(st.cv), num(st.subjMean), num(st.subjCV));
end
line(72);
end

% =====================================================================
function R = printRetest(T, qty, title)
R = struct('nPairs',0);
if isempty(T), fprintf('\n%s: no scans.\n', title); return; end
subj = string(T.subject);
uS   = unique(subj);
rep  = uS(arrayfun(@(s) sum(subj==s) >= 2, uS));
if isempty(rep), fprintf('\n%s: no subject scanned twice.\n', title); return; end

fprintf('\n%s\n', title);
fprintf('  pairs: %s\n', strjoin(cellstr(rep(:)'), ', '));
line(72);
fprintf('%-11s %-4s %6s %10s %10s %7s %10s\n','','unit','nPairs','mean','wsSD','wsCV%','RC');
line(72);
lastBlock = '';
for k = 1:size(qty,1)
    col = qty{k,1};
    if ~ismember(col, T.Properties.VariableNames), continue; end
    if ~strcmp(qty{k,4}, lastBlock)
        lastBlock = qty{k,4};
        fprintf('%s\n', lastBlock);
    end
    v = double(T.(col));
    d = []; m = [];
    for s = 1:numel(rep)
        x = v(subj == rep(s));
        x = x(isfinite(x));
        if numel(x) < 2, continue; end
        x = x(1:2);                       % first two visits, chronological
        d(end+1) = x(2)-x(1);             %#ok<AGROW>
        m(end+1) = mean(x);               %#ok<AGROW>
    end
    st = struct('nPairs',numel(d),'mean',NaN,'wsSD',NaN,'wsCV',NaN,'RC',NaN);
    if ~isempty(d)
        st.mean = mean(m);
        st.wsSD = sqrt(sum(d.^2)/(2*numel(d)));
        st.wsCV = 100*st.wsSD/abs(st.mean);
        st.RC   = 2.77*st.wsSD;
    end
    R.(matlab.lang.makeValidName(col)) = st;
    fprintf('  %-9s %-4s %6d %10s %10s %7s %10s\n', qty{k,2}, qty{k,3}, ...
        st.nPairs, num(st.mean), num(st.wsSD), num(st.wsCV), num(st.RC));
end
R.nPairs = numel(rep);
line(72);
end

% =====================================================================
function line(w), fprintf('%s\n', repmat('-',1,w)); end
function banner(s), fprintf('\n\n%s\n%s\n', s, repmat('=',1,numel(s))); end

function s = txt(col, i)
% readtable types scandate as datetime and session as cellstr; print either.
if iscell(col),          s = col{i};
elseif isdatetime(col),  s = datestr(col(i), 'yyyy-mm-dd'); %#ok<DATST>
else,                    s = char(string(col(i)));
end
end

function s = num(x)
if ~isfinite(x),       s = '-';
elseif abs(x) >= 100,  s = sprintf('%.0f', x);
elseif abs(x) >= 10,   s = sprintf('%.1f', x);
elseif abs(x) >= 1,    s = sprintf('%.2f', x);
else,                  s = sprintf('%.3f', x);
end
end
