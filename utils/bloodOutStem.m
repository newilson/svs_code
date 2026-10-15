function stem = bloodOutStem(SessionDir)
% BLOODOUTSTEM  Output filename stem for a blood session: the session folder
% name, with '_blood' appended unless it already ends in '_blood'.
%   ...\859875_S01_09012026_blood -> 859875_S01_09012026_blood
%   ...\859875_S01_09012026       -> 859875_S01_09012026_blood

[~, scanID] = fileparts(regexprep(SessionDir, '[\\/]+$', ''));
if endsWith(scanID, '_blood', 'IgnoreCase', true)
    stem = scanID;
else
    stem = [scanID '_blood'];
end
end
