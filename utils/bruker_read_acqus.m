function dic = bruker_read_acqus(dirpath, acqusFiles)
% BRUKER_READ_ACQUS  Read Bruker acquisition parameter files (acqus, acqu2s, …)
%
% Usage:
%   dic = bruker_read_acqus(dirpath)
%   dic = bruker_read_acqus(dirpath, acqusFiles)
%
% Inputs:
%   dirpath    - Directory containing acqus / acqu2s / acqu3s / acqu4s files
%   acqusFiles - (optional) Cell array of explicit filenames or full paths.
%                When omitted the function auto-detects standard file names.
%
% Output:
%   dic  - Struct whose fields are the parameter-file basenames
%          (e.g. dic.acqus, dic.acqu2s …).  Each field is itself a struct
%          whose fields are the parsed JCAMP-DX parameters.
%
% Example:
%   p = bruker_read_acqus('/data/nmr/exp1');
%   td = p.acqus.TD;           % time-domain points
%   sw = p.acqus.SW_h;         % sweep width in Hz
%   sfo1 = p.acqus.SFO1;       % carrier frequency
%
% See also: bruker_read_jcamp, bruker_read

if nargin < 2 || isempty(acqusFiles)
    candidates = {'acqus','acqu2s','acqu3s','acqu4s'};
    acqusFiles = {};
    for k = 1:numel(candidates)
        fp = fullfile(dirpath, candidates{k});
        if isfile(fp)
            acqusFiles{end+1} = fp; %#ok<AGROW>
        end
    end
else
    % resolve relative filenames against dirpath
    for k = 1:numel(acqusFiles)
        if ~isfile(acqusFiles{k})
            acqusFiles{k} = fullfile(dirpath, acqusFiles{k});
        end
    end
end

dic = struct();

for k = 1:numel(acqusFiles)
    [~, fname, ext] = fileparts(acqusFiles{k});
    fieldName = matlab.lang.makeValidName([fname ext]);
    dic.(fieldName) = bruker_read_jcamp(acqusFiles{k});
end

end
