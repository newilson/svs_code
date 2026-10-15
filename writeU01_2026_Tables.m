function writeU01_2026_Tables(dataRoot, outDir)
% WRITEU01_2026_TABLES  Regenerate every U01_2026 table and save it to disk.
%
%   writeU01_2026_Tables()
%   writeU01_2026_Tables(dataRoot, outDir)
%
% Runs the whole reporting chain in dependency order and captures the
% formatted (phone-readable) tables into ONE text file, alongside the
% machine-readable TSVs that each step writes anyway:
%
%   U01_2026_brain.tsv          per-session, brain
%   U01_2026_calf.tsv           per-session, calf
%   U01_2026_summary.tsv        group means / CVs / 31P-1H NAD+ ratios
%   U01_2026_repeatability.tsv  test-retest over the Repeatability\ cohort
%   U01_2026_extref.tsv         the external-phantom subset, 3 references
%   U01_2026_tables.txt         <- every printed table, in one file
%   U01_2026_collate.log        per-session collation chatter, kept separate
%
% Re-collation is not optional and always runs first: the TSVs are built from
% the per-session *_abs.mat files, so re-running any processing stage silently
% invalidates them.
%
% Each reporting step is wrapped so that one failure does not cost you the
% others -- the traceback lands in the text file where you will see it.

if nargin < 1 || isempty(dataRoot)
    dataRoot = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026\2026';
end
if nargin < 2 || isempty(outDir)
    outDir = fileparts(dataRoot);
end
assert(isfolder(dataRoot), 'dataRoot not found: %s', dataRoot);

txtPath = fullfile(outDir, 'U01_2026_tables.txt');
logPath = fullfile(outDir, 'U01_2026_collate.log');
stamp   = char(datetime('now','Format','yyyy-MM-dd HH:mm'));

% ---- 1) re-collate the per-session TSVs (noisy -- its own log) ----------
diary off
if isfile(logPath), delete(logPath); end
diary(logPath); diary on
fprintf('collectU01_2026_TSV   %s\n\n', stamp);
try
    collectU01_2026_TSV(dataRoot, outDir);
catch ME
    fprintf('COLLATION FAILED: %s\n%s\n', ME.message, getReport(ME,'basic'));
    diary off
    rethrow(ME)          % nothing downstream is meaningful without the TSVs
end
diary off

% ---- 2) the tables -----------------------------------------------------
if isfile(txtPath), delete(txtPath); end
diary(txtPath); diary on

fprintf('%s\n', repmat('=',1,46));
fprintf('U01_2026 TABLES\n');
fprintf('generated  %s\n', stamp);
fprintf('data root  %s\n', dataRoot);
fprintf('%s\n', repmat('=',1,46));

runStep(@() summarizeU01_2026(outDir),          'GROUP SUMMARY');
runStep(@() repeatabilityU01_2026(outDir),      'REPEATABILITY');
runStep(@() extRefReportU01_2026(dataRoot, outDir), 'EXTERNAL REFERENCE');

fprintf('\n%s\n', repmat('=',1,46));
fprintf('end of tables\n');
diary off

fprintf('Tables:    %s\n', txtPath);
fprintf('Collation: %s\n', logPath);
end

% =====================================================================
function runStep(fcn, name)
fprintf('\n\n%s\n#  %s\n%s\n', repmat('#',1,46), name, repmat('#',1,46));
try
    fcn();
catch ME
    fprintf('\n*** %s FAILED: %s\n', name, ME.message);
    fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
end
end
