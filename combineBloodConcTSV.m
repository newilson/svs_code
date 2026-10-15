function outPath = combineBloodConcTSV(inDir, outName)
% COMBINEBLOODCONCTSV  Combine per-session blood 1H result TSVs into one.
%
%   outPath = combineBloodConcTSV(inDir)
%   outPath = combineBloodConcTSV(inDir, outName)
%
% Scans inDir for every '<stem>_1H_abs.tsv' (stem ends in '_blood') written by
% processU01_Blood_1H.m (one header row + one data row per session:
% dataset, water T2 / consistency columns, <peak>_uM_water, <peak>_uM_TSP)
% and stacks them into one TSV, one row per session.  Columns are the union
% across files, ordered as first encountered; missing values are written NA.
%
% Default outName: 'blood_1H_abs_combined.tsv' (written into inDir).

if nargin < 2 || isempty(outName), outName = 'blood_1H_abs_combined.tsv'; end
assert(isfolder(inDir), 'inDir not found: %s', inDir);

files = dir(fullfile(inDir, '*_blood_1H_abs.tsv'));
files = files(~ismember({files.name}, {outName}));
assert(~isempty(files), 'No *_blood_1H_abs.tsv files found in %s', inDir);

cols = {};
rows = struct('dataset', {}, 'vals', {});
for i = 1:numel(files)
    f = fullfile(inDir, files(i).name);
    T = readtable(f, 'FileType','text', 'Delimiter','\t', 'TextType','char');
    assert(ismember('dataset', T.Properties.VariableNames), ...
        '%s is not a per-session (wide) blood TSV', files(i).name);
    m = containers.Map('KeyType','char','ValueType','double');
    for v = setdiff(T.Properties.VariableNames, {'dataset'}, 'stable')
        m(v{1}) = T.(v{1})(1);
        if ~any(strcmp(cols, v{1})), cols{end+1} = v{1}; end %#ok<AGROW>
    end
    rows(end+1) = struct('dataset', char(string(T.dataset(1))), 'vals', m); %#ok<AGROW>
end

outPath = fullfile(inDir, outName);
fid = fopen(outPath, 'w');
assert(fid ~= -1, 'Could not open %s for writing', outPath);
fprintf(fid, 'dataset');
for j = 1:numel(cols), fprintf(fid, '\t%s', cols{j}); end
fprintf(fid, '\n');
for i = 1:numel(rows)
    fprintf(fid, '%s', rows(i).dataset);
    for j = 1:numel(cols)
        if isKey(rows(i).vals, cols{j})
            fprintf(fid, '\t%.4g', rows(i).vals(cols{j}));
        else
            fprintf(fid, '\tNA');
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('Combined %d sessions -> %s\n', numel(rows), outPath);
end
