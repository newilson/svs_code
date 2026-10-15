function S = checkU01_2026_relaxation(inDir, dataRoot)
% CHECKU01_2026_RELAXATION  Audit the T1/T2 correction factors behind the plots.
%
%   S = checkU01_2026_relaxation()
%
% For every scan that appears in the four plotU01_2026_scatter figures, reports
% the relaxation assumptions actually baked into that scan's _abs.mat and flags
% any that deviate from the modal value of its group.
%
% WHY THIS IS NOT JUST "COMPARE R":
%   R = (1 - exp(-TR/T1)) * exp(-TE/T2) depends on TR, and the 1H TR was 500,
%   520 or 600 ms across these sessions.  Two scans with identical T1/T2
%   assumptions therefore have different R.  So:
%     31P -- absQ stores T1_ms/T2_ms/nPhos per metabolite; compare those directly.
%     1H  -- absQ stores only the resulting R_M.  T1 and T2 cannot both be
%            recovered from one number, so instead we PREDICT R from the
%            reference table for that scan's own TR/TE and check the stored R_M
%            matches.  A mismatch means that scan used different assumptions;
%            the implied T1 (holding T2 at the table value) is reported so the
%            size of the discrepancy is interpretable.
%
% Reference tables are read from the drivers, so if a driver is edited after
% processing this audit will flag every scan -- which is the intended alarm.

if nargin < 1 || isempty(inDir)
    inDir = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026';
end
if nargin < 2 || isempty(dataRoot)
    dataRoot = fullfile(inDir, '2026');
end

% ---- reference tables, mirroring the drivers ---------------------------
% processU01_Brain_1H.m: met_T1_ms / met_T2_ms
ref1H.brain = struct('NADH2',[205.6 33.6], 'NADH6',[211.6 29.1], ...
                     'NADH4',[237.3 42.3], 'Trp',[75 19]);
% processU01_Calf_1H.m
ref1H.calf  = struct('NADH2',[50 14], 'NADH6',[50 14], 'NADH4',[50 14], ...
                     'NR',[100 11], 'Trp',[100 11]);
% water, per tissue, in driver order
refW.brain = struct('tissues',{{'CSF','GM','WM'}}, 'T1',[4400 1939 1130], 'T2',[300 40 40]);
refW.calf  = struct('tissues',{{'Muscle'}},        'T1',2000,             'T2',25);

% metabolites that actually appear in the plots
want1H = struct('brain',{{'NADH2','NADH6','Trp'}}, 'calf',{{'NADH2','NADH6'}});
want31P = {'NADplus','NADH'};
refName = struct('brain','aATP','calf','PCr');

parts = {'brain','calf'};
S = struct();

for p = 1:numel(parts)
    anat = parts{p};
    T = readtable(fullfile(inDir, sprintf('U01_2026_%s.tsv', anat)), 'FileType','text', ...
                  'Delimiter','\t', 'VariableNamingRule','preserve');
    T = T(abs(double(T.TE_31P_ms)-0.1) < 1e-9 & contains(string(T.seq_1H),'svssel_20260101'), :);
    T = T(~ismember(string(T.session), "S007_07222026_31P-1H-HeadCoil_external_ref"), :);
    T = sortrows(T, {'scandate','session'});
    S.(anat).sessions = string(T.session);

    n = height(T);
    rows31 = {}; rows1H = {};

    for k = 1:n
        sess = T.session{k};
        pd   = fullfile(dataRoot, sess, 'processed');
        if ~isfolder(pd), pd = fullfile(dataRoot,'Repeatability', sess, 'processed'); end

        % ---------------- 31P: stored directly ----------------
        f = fullfile(pd, sprintf('%s_%s_31P_abs.mat', sess, anat));
        if isfile(f)
            m = load(f, 'absQ');
            names = want31P;
            if isfield(m.absQ,'reference') && isfield(m.absQ.reference,'name')
                names = [want31P {m.absQ.reference.name}]; %#ok<AGROW>
            end
            for q = 1:numel(names)
                [t1,t2,np] = get31P(m.absQ, names{q});
                rows31{end+1} = {k, sess, names{q}, t1, t2, np}; %#ok<AGROW>
            end
        end

        % ---------------- 1H: verify R_M against the table ----------------
        f = fullfile(pd, sprintf('%s_%s_1H_abs.mat', sess, anat));
        if isfile(f)
            m = load(f, 'absQ');
            TR = getf(m.absQ.inputs,'TR_ms',NaN);
            TE = getf(m.absQ.inputs,'TE_ms',NaN);
            for q = 1:numel(want1H.(anat))
                nmq = want1H.(anat){q};
                Rm  = get1H(m.absQ, nmq);
                tab = ref1H.(anat).(nmq);
                Rex = (1-exp(-TR/tab(1))) * exp(-TE/tab(2));
                % invert holding T2 at the table value, so a mismatch has units
                arg = 1 - Rm/exp(-TE/tab(2));
                if isfinite(arg) && arg > 0 && arg < 1
                    T1imp = -TR/log(arg);
                else
                    T1imp = NaN;
                end
                rows1H{end+1} = {k, sess, nmq, TR, TE, Rm, Rex, Rm/Rex, tab(1), T1imp}; %#ok<AGROW>
            end
            % water reference relaxation
            rw = getf(m.absQ.water,'Rw_per_tissue',[]);
            TRw = getf(m.absQ.inputs,'TR_ws_ms',NaN);
            TEw = getf(m.absQ.inputs,'TE_ws_ms',NaN);
            w   = refW.(anat);
            rwEx = (1-exp(-TRw./w.T1)) .* exp(-TEw./w.T2);
            S.(anat).water(k).session = sess;
            S.(anat).water(k).Rw      = rw(:).';
            S.(anat).water(k).Rw_exp  = rwEx(:).';
            S.(anat).water(k).ok      = numel(rw)==numel(rwEx) && ...
                                        all(abs(rw(:).'-rwEx(:).') < 1e-9);
        end
    end

    S.(anat).t31 = cell2table(vertcat(rows31{:}), 'VariableNames', ...
        {'scan','session','metab','T1_ms','T2_ms','nPhos'});
    S.(anat).t1H = cell2table(vertcat(rows1H{:}), 'VariableNames', ...
        {'scan','session','metab','TR_ms','TE_ms','R_M','R_expected','R_ratio', ...
         'T1_table_ms','T1_implied_ms'});
end

% =====================================================================
% Report
% =====================================================================
for p = 1:numel(parts)
    anat = parts{p};
    banner(sprintf('%s -- 31P (stored T1/T2)', upper(anat)));
    reportModal(S.(anat).t31, {'T1_ms','T2_ms','nPhos'});

    banner(sprintf('%s -- 1H (R_M vs table prediction)', upper(anat)));
    t = S.(anat).t1H;
    mets = unique(t.metab, 'stable');
    fprintf('%-8s %8s %8s %10s %12s %12s\n','metab','TR set','T2 tbl','T1 tbl','T1 implied','R ratio');
    fprintf('%s\n', repmat('-',1,64));
    for q = 1:numel(mets)
        r = t(strcmp(t.metab, mets{q}), :);
        trs = unique(r.TR_ms).';
        imp = unique(round(r.T1_implied_ms, 6));
        fprintf('%-8s %8s %8s %10.4g %12s %12.6g\n', mets{q}, ...
                strjoin(compose('%g',trs),'/'), '-', r.T1_table_ms(1), ...
                strjoin(compose('%.4g',imp.'),'/'), max(abs(r.R_ratio-1))+1);
        bad = r(abs(r.R_ratio-1) > 1e-9, :);
        for b = 1:height(bad)
            fprintf('    FLAG scan %d %s: R_M %.6g vs expected %.6g (implied T1 %.4g ms)\n', ...
                bad.scan(b), bad.session{b}, bad.R_M(b), bad.R_expected(b), bad.T1_implied_ms(b));
        end
    end

    % water
    w = S.(anat).water;
    okAll = all([w.ok]);
    fprintf('\nwater reference relaxation: %s\n', ternary(okAll, ...
        'all scans match the driver tissue table', 'MISMATCH -- see S.(anat).water'));
    if ~okAll
        for k = 1:numel(w)
            if ~w(k).ok
                fprintf('    FLAG %s: Rw %s vs expected %s\n', w(k).session, ...
                    mat2str(w(k).Rw,5), mat2str(w(k).Rw_exp,5));
            end
        end
    end
end
end

% =====================================================================
function reportModal(t, cols)
mets = unique(t.metab, 'stable');
fprintf('%-10s', 'metab');
for c = 1:numel(cols), fprintf(' %14s', cols{c}); end
fprintf('   %s\n', 'n');
fprintf('%s\n', repmat('-',1,10+15*numel(cols)+6));
for q = 1:numel(mets)
    r = t(strcmp(t.metab, mets{q}), :);
    fprintf('%-10s', mets{q});
    flags = {};
    for c = 1:numel(cols)
        v = r.(cols{c});
        mv = mode(v);
        u  = unique(v);
        if numel(u) == 1
            fprintf(' %14g', mv);
        else
            fprintf(' %13g*', mv);
            odd = find(v ~= mv);
            for o = odd(:).'
                flags{end+1} = sprintf('scan %d %s: %s = %g (modal %g)', ...
                    r.scan(o), r.session{o}, cols{c}, v(o), mv); %#ok<AGROW>
            end
        end
    end
    fprintf('   %d\n', height(r));
    for i = 1:numel(flags), fprintf('    FLAG %s\n', flags{i}); end
end
end

% =====================================================================
function [t1,t2,np] = get31P(absQ, name)
t1 = NaN; t2 = NaN; np = NaN;
if isfield(absQ,'reference') && strcmp(absQ.reference.name, name)
    t1 = absQ.reference.T1_ms; t2 = absQ.reference.T2_ms;
    np = getf(absQ.reference,'nPhos',NaN);
    return
end
M = absQ.metabolites;
for k = 1:numel(M)
    if strcmp(M(k).name, name)
        t1 = getf(M(k),'T1_ms',NaN); t2 = getf(M(k),'T2_ms',NaN);
        np = getf(M(k),'nPhos',NaN);
        return
    end
end
end

function R = get1H(absQ, name)
R = NaN;
M = absQ.metabolites;
for k = 1:numel(M)
    if strcmp(M(k).name, name), R = getf(M(k),'R_M',NaN); return; end
end
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
function banner(s), fprintf('\n\n%s\n%s\n', s, repmat('=',1,numel(s))); end
