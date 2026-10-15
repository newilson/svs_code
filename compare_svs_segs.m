function compare_svs_segs(outDir)
% COMPARE_SVS_SEGS  Stage 3 of testBrainSeg pipeline.
%   Reads brainseg.sh's per-scan SVS-masked stats TSV
%   (<stem>_SPM8bias_hdbet_fast3_first_seg_svsmask.tsv), computes
%   CSF/GM/WM voxel fractions, treats MPRAGE_P2* as ground truth, and
%   ranks every other scan by L1 distance of its [CSF GM WM] vector to
%   the ground truth.
%   Writes outDir/svs_seg_comparison.csv.

M = readtable(fullfile(outDir,'manifest.csv'), ...
    'Delimiter',',', 'ReadVariableNames',true, ...
    'VariableNamingRule','preserve', 'TextType','string');
n = height(M);

frac   = nan(n,3);             % [CSF GM WM]
labels = strings(n,1);

for k = 1:n
    stem = char(M{k,1});       % scanName (column 1)
    tsv  = fullfile(outDir, [stem '_SPM8bias_hdbet_fast3_first_seg_svsmask.tsv']);
    labels(k) = stem;
    if ~isfile(tsv)
        warning('Missing TSV for %s: %s', stem, tsv);
        continue
    end
    % expected columns: labelval, voxels, volume_mm_3_, labelname
    T = readtable(tsv,'FileType','text','Delimiter','\t', ...
                  'VariableNamingRule','preserve');
    voxCol  = find(strcmpi(T.Properties.VariableNames,'voxels'),1);
    nameCol = find(strcmpi(T.Properties.VariableNames,'labelname'),1);
    assert(~isempty(voxCol) && ~isempty(nameCol), ...
        'Unexpected TSV columns in %s: %s', tsv, strjoin(T.Properties.VariableNames,','));
    nm  = string(T{:,nameCol});
    vx  = double(T{:,voxCol});
    csf = sum(vx(nm=="CSF")); gm = sum(vx(nm=="GM")); wm = sum(vx(nm=="WM"));
    tot = csf + gm + wm;
    if tot > 0
        frac(k,:) = [csf gm wm] / tot;
    end
end

% --- ground truth: first scan whose stem starts with MPRAGE_P2 ---
gtIdx = find(startsWith(labels,'MPRAGE_P2'),1);
assert(~isempty(gtIdx), 'No MPRAGE_P2* scan in manifest; cannot compare.');
gtFrac = frac(gtIdx,:);

dist = sum(abs(frac - gtFrac), 2);

% --- print ranked report ---
[~, ord] = sort(dist);
fprintf('\n%-40s  %6s %6s %6s   %s\n','scan','CSF','GM','WM','L1 vs MPRAGE_P2');
fprintf('%s\n', repmat('-',1,80));
for kk = 1:n
    i = ord(kk);
    if i == gtIdx, tag = ' (ground truth)'; else, tag = ''; end
    fprintf('%-40s  %6.3f %6.3f %6.3f   %.4f%s\n', ...
        labels(i), frac(i,1), frac(i,2), frac(i,3), dist(i), tag);
end

% closest non-GT
nonGt = ord(ord~=gtIdx);
if ~isempty(nonGt)
    fprintf('\nClosest non-GT scan: %s (L1=%.4f)\n', labels(nonGt(1)), dist(nonGt(1)));
end

% --- save CSV ---
out = table(labels, frac(:,1), frac(:,2), frac(:,3), dist, ...
    'VariableNames', {'scan','CSF','GM','WM','L1_vs_MPRAGE_P2'});
out = sortrows(out,'L1_vs_MPRAGE_P2');
writetable(out, fullfile(outDir,'svs_seg_comparison.csv'));
fprintf('Wrote %s\n', fullfile(outDir,'svs_seg_comparison.csv'));

end
