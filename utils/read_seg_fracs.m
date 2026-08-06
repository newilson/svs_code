function segFracs = read_seg_fracs(segTsv)
% READ_SEG_FRACS  CSF/GM/WM volume fractions from a brainseg.sh TSV.
%
%   segFracs = read_seg_fracs(segTsv)
%
% Accepts either of the two layouts brainseg.sh has emitted, so sessions
% segmented with the old and new versions of the script can be quantified by
% the same pipeline.
%
% LONG -- one row per FSL label
%
%   labelval  voxels  volume(mm^3)  labelname
%   1         23356   93424.00      CSF
%   2         61358   245432.00     GM
%   3         68398   273592.00     WM
%
% WIDE -- one row per ROI
%
%   ScanID   CSFvol(cc)  GMvol(cc)  WMvol(cc)  CSFfraction ... nonBrain(cc)  ROIindex
%   dummyID  76.9859     217.6992   183.6991   .1609       ... 122.8620      1
%
% Input:
%   segTsv - path to the *_seg_svsmask.tsv written by brainseg.sh
%
% Output:
%   segFracs - [f_CSF f_GM f_WM], summing to 1

T    = readtable(segTsv,'FileType','text','Delimiter','\t', ...
                 'VariableNamingRule','preserve');
cols = T.Properties.VariableNames;

voxCol  = find(strcmpi(cols,'voxels'),    1);
nameCol = find(strcmpi(cols,'labelname'), 1);

if ~isempty(voxCol) && ~isempty(nameCol)
    % ---- LONG ----
    nm = string(T{:,nameCol});
    vx = tonum(T{:,voxCol});
    v  = [sum(vx(nm=="CSF")) sum(vx(nm=="GM")) sum(vx(nm=="WM"))];
else
    % ---- WIDE ----
    volCols  = {'CSFvol(cc)',  'GMvol(cc)',  'WMvol(cc)' };
    fracCols = {'CSFfraction', 'GMfraction', 'WMfraction'};
    iVol  = cellfun(@(c) find(strcmpi(cols,c),1), volCols,  'UniformOutput', false);
    iFrac = cellfun(@(c) find(strcmpi(cols,c),1), fracCols, 'UniformOutput', false);
    if all(~cellfun(@isempty, iVol))
        use = iVol;
    elseif all(~cellfun(@isempty, iFrac))
        use = iFrac;
    else
        error('Unrecognised seg TSV layout in %s: %s', segTsv, strjoin(cols, ','));
    end

    % One row per ROI.  Refuse to guess when several are present rather than
    % silently quantifying against the wrong ROI.
    assert(height(T) == 1, ...
        ['Seg TSV %s has %d ROI rows; this reader expects one. ' ...
         'Pass a single-ROI TSV.'], segTsv, height(T));

    v = [tonum(T{1,use{1}}) tonum(T{1,use{2}}) tonum(T{1,use{3}})];
end

tot = sum(v);
assert(tot > 0, 'No CSF/GM/WM volume in %s', segTsv);
segFracs = v / tot;     % [CSF GM WM]

end

% =====================================================================
function x = tonum(v)
% Numeric value(s) from a table column readtable may have returned as text
% (a single stray character anywhere in the column is enough to do it).

if isnumeric(v)
    x = double(v);
elseif iscell(v)
    x = str2double(v);
else
    x = str2double(string(v));
end
end
