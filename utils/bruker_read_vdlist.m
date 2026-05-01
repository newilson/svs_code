function vdlist = bruker_read_vdlist(dirpath, fname)
% BRUKER_READ_VDLIST  Read a Bruker variable-delay list file (vdlist).
%
% vdlist files store delay times used in NMR relaxation experiments.
% Values may carry unit suffixes: n (ns), u (µs), m (ms), s (s).
% This function converts all delays to seconds.
%
% Usage:
%   vdlist = bruker_read_vdlist(dirpath)
%   vdlist = bruker_read_vdlist(dirpath, fname)
%
% Inputs:
%   dirpath  - Directory containing the vdlist file
%   fname    - (optional) Filename (default: 'vdlist')
%
% Output:
%   vdlist   - Column vector of delay times in seconds (double)
%
% See also: bruker_read

if nargin < 2 || isempty(fname)
    fname = 'vdlist';
end

fullpath = fullfile(dirpath, fname);
if ~isfile(fullpath)
    error('bruker_read_vdlist:fileNotFound', ...
          'vdlist file not found: %s', fullpath);
end

lines = readlines(fullpath, 'EmptyLineRule', 'skip');

vdlist = zeros(numel(lines), 1);
for k = 1:numel(lines)
    tok = strtrim(lines{k});
    if isempty(tok), continue; end

    % replace unit suffix with matching power-of-10 exponent
    tok = regexprep(tok, 'n\s*$', 'e-9');
    tok = regexprep(tok, 'u\s*$', 'e-6');
    tok = regexprep(tok, 'm\s*$', 'e-3');
    tok = regexprep(tok, 's\s*$', 'e0');

    val = str2double(tok);
    if isnan(val)
        error('bruker_read_vdlist:parseFail', ...
              'Cannot parse vdlist value: %s', lines{k});
    end
    vdlist(k) = val;
end

end
