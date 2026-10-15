function summarizeU01_2026(inDir)
% SUMMARIZEU01_2026  Group statistics over the 2026 U01 reprocessing TSVs.
%
%   summarizeU01_2026()
%   summarizeU01_2026(inDir)
%
% Reads U01_2026_brain.tsv and U01_2026_calf.tsv (written by
% collectU01_2026_TSV) and produces, for each nucleus x body-part combination:
%
%   1) mean, SD and CV% of every reported concentration
%   2) the per-scan [NAD+]_31P / [NAD+]_1H ratio, averaged within body part,
%      for ALL THREE reference standards (aATP, PCr, external phantom)
%
% The ratio is computed per scan and THEN averaged (not as a ratio of group
% means), so each scan contributes its own internally consistent pairing.
% 1H concentrations are uM and 31P are mM, so the 1H value is divided by 1000
% before the ratio is taken -- a ratio of 1 therefore means the two nuclei
% agree on [NAD+].
%
% The 1H denominator is reported four ways because the three NAD+ resonances
% do not agree: H2, H6, H4 and their mean.  H4 is the least reliable (it sits
% near the fit-window edge and is badly behaved in calf), so the H2 and H6
% columns are usually the ones to trust.
%
% Output: <inDir>\U01_2026_summary.tsv  (long format: one row per statistic)
% and phone-readable tables printed to the command window.  Values are rounded
% to 3 significant figures (CV% to 2) -- full precision lives in the
% per-session brain/calf TSVs.

if nargin < 1 || isempty(inDir)
    inDir = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026';
end

fB = fullfile(inDir,'U01_2026_brain.tsv');
fC = fullfile(inDir,'U01_2026_calf.tsv');
T = struct();
if isfile(fB), T.brain = readTsv(fB); else, T.brain = table(); end
if isfile(fC), T.calf  = readTsv(fC); else, T.calf  = table(); end

parts = {};
if ~isempty(T.brain), parts{end+1} = 'brain'; end
if ~isempty(T.calf),  parts{end+1} = 'calf';  end
assert(~isempty(parts), 'No TSVs found in %s', inDir);

% ---- what to summarize -------------------------------------------------
oneH = {'uM_NADH2','uM_NADH6','uM_NADH4','uM_NAD_1H_mean','uM_Trp'};
spec = {'NADplus','NADH','UDP','aATP','PCr'};
refs = {'aATPref','PCrref','extref'};
refLabel = struct('aATPref','aATP-ref','PCrref','PCr-ref','extref','ext-ref');

rows = {};   % long-format output accumulator

% =====================================================================
% 1) means / CVs
% =====================================================================
for p = 1:numel(parts)
    part = parts{p};
    for i = 1:numel(oneH)
        rows = addStat(rows, part, '1H', 'uM', oneH{i}, getCol(T.(part), oneH{i}));
    end
    for r = 1:numel(refs)
        for i = 1:numel(spec)
            v = getCol(T.(part), [refs{r} '_mM_' spec{i}]);
            rows = addStat(rows, part, ['31P ' refLabel.(refs{r})], 'mM', spec{i}, v);
        end
    end
    if strcmp(part,'brain')
        for f = {'f_GM','f_WM','f_CSF'}
            rows = addStat(rows, part, 'seg', 'frac', f{1}, getCol(T.(part), f{1}));
        end
    end
end

% =====================================================================
% 2) per-scan 31P/1H NAD+ ratio, averaged within body part
% =====================================================================
dens = {'uM_NADH2','uM_NADH6','uM_NADH4','uM_NAD_1H_mean'};
denLbl = {'H2','H6','H4','mean'};
ratio = struct();
for p = 1:numel(parts)
    part = parts{p};
    for r = 1:numel(refs)
        n31 = getCol(T.(part), [refs{r} '_mM_NADplus']);
        for d = 1:numel(dens)
            n1h = getCol(T.(part), dens{d}) / 1000;      % uM -> mM
            v   = n31 ./ n1h;
            v(~isfinite(v)) = NaN;                        % kill /0 and NaN pairs
            ratio.(part).(refs{r}).(denLbl{d}) = v;
            rows = addStat(rows, part, ['ratio 31P/1H ' refLabel.(refs{r})], ...
                           '-', ['vs_' denLbl{d}], v);
        end
    end
end

% =====================================================================
% Write long-format TSV
% =====================================================================
S = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'bodypart','block','quantity','unit','n','mean','sd','cv_pct'});
outPath = fullfile(inDir,'U01_2026_summary.tsv');
writetable(S, outPath, 'FileType','text','Delimiter','\t');

% =====================================================================
% Print phone-readable tables (kept under ~46 chars wide)
% =====================================================================
W = 46;
hr  = @() fprintf('%s\n', repmat('-',1,W));
ttl = @(s) fprintf('\n%s\n%s\n', s, repmat('=',1,numel(s)));

haveB = ismember('brain', parts);

ttl('MEANS AND CVs');
fprintf('n scans:  brain=%d  calf=%d\n', ...
    heightOr(T,'brain'), heightOr(T,'calf'));

printBlock('1H  NAD+/Trp  [uM]', oneH, stripPfx('uM_'));
for r = 1:numel(refs)
    printBlock(sprintf('31P %s  [mM]', refLabel.(refs{r})), ...
               strcat(refs{r},'_mM_',spec), stripPfx([refs{r} '_mM_']));
end
if haveB
    printBlock('Voxel composition  [frac]', {'f_GM','f_WM','f_CSF'}, stripPfx('f_'));
end

ttl('31P/1H NAD+ RATIO');
fprintf('per-scan ratio, then averaged.\n');
fprintf('1 = both nuclei agree on [NAD+].\n');
for r = 1:numel(refs)
    fprintf('\n%s\n', refLabel.(refs{r}));
    hr();
    fprintf('%-8s %16s %16s\n', 'vs 1H', 'BRAIN', 'CALF');
    fprintf('%-8s %7s %8s %7s %8s\n', '', 'mean', 'CV%', 'mean', 'CV%');
    hr();
    for d = 1:numel(denLbl)
        fprintf('%-8s', denLbl{d});
        for part = {'brain','calf'}
            if ismember(part{1}, parts)
                v = ratio.(part{1}).(refs{r}).(denLbl{d});
                fprintf(' %7s %8s', fmtNum(mean(v,'omitnan')), fmtCv(v));
            else
                fprintf(' %7s %8s', '-', '-');
            end
        end
        fprintf('\n');
    end
    hr();
end

fprintf('\nSummary TSV: %s\n', outPath);

% ---------------------------------------------------------------------
    function printBlock(title, names, labelFcn)
        fprintf('\n%s\n', title);
        hr();
        fprintf('%-11s %15s %15s\n', '', 'BRAIN', 'CALF');
        fprintf('%-11s %7s %7s %7s %7s\n', '', 'mean', 'CV%', 'mean', 'CV%');
        hr();
        for q = 1:numel(names)
            lbl = labelFcn(names{q});
            fprintf('%-11s', lbl);
            for pp = {'brain','calf'}
                if ismember(pp{1}, parts)
                    v = getCol(T.(pp{1}), names{q});
                    fprintf(' %7s %7s', fmtNum(mean(v,'omitnan')), fmtCv(v));
                else
                    fprintf(' %7s %7s', '-', '-');
                end
            end
            fprintf('\n');
        end
        hr();
    end
end

% =====================================================================
function f = stripPfx(pfx)
f = @(s) erase(s, pfx);
end

function h = heightOr(T, part)
if isfield(T, part) && ~isempty(T.(part)), h = height(T.(part)); else, h = 0; end
end

function rows = addStat(rows, part, block, unit, quantity, v)
v = v(isfinite(v));
n = numel(v);
if n == 0
    rows{end+1} = {part, block, quantity, unit, 0, NaN, NaN, NaN};
    return
end
mu = mean(v);
sd = std(v);
% The reference metabolite is a fixed constant, so its SD is pure float noise
% (~1e-16). Snap it to zero rather than printing 4.8e-16.
if sd < 1e-10 * abs(mu), sd = 0; end
cv = 100 * sd / abs(mu);
% Round for readability. The underlying scans do not justify more than this,
% and the full-precision values are in the per-session brain/calf TSVs.
rows{end+1} = {part, block, quantity, unit, n, ...
               sigfig(mu,3), sigfig(sd,3), sigfig(cv,2)};
end

function y = sigfig(x, n)
if isfinite(x) && x ~= 0
    y = round(x, n, 'significant');
else
    y = x;
end
end

function v = getCol(T, name)
if isempty(T) || ~ismember(name, T.Properties.VariableNames)
    v = NaN; return
end
c = T.(name);
if iscell(c), c = str2double(c); end
v = double(c);
end

function s = fmtNum(x)
if ~isfinite(x),        s = '-';
elseif abs(x) >= 1000,  s = sprintf('%.0f', x);
elseif abs(x) >= 10,    s = sprintf('%.1f', x);
elseif abs(x) >= 1,     s = sprintf('%.2f', x);
else,                   s = sprintf('%.3f', x);
end
end

function s = fmtCv(v)
v = v(isfinite(v));
if numel(v) < 2, s = '-'; return; end
cv = 100 * std(v) / abs(mean(v));
if ~isfinite(cv), s = '-'; else, s = sprintf('%.0f', cv); end
end

function T = readTsv(f)
T = readtable(f, 'FileType','text', 'Delimiter','\t', ...
              'VariableNamingRule','preserve');
end
