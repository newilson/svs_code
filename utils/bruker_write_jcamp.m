function bruker_write_jcamp(dic, filename, overwrite)
% BRUKER_WRITE_JCAMP  Write a Bruker JCAMP-DX parameter file from a struct.
%
% Usage:
%   bruker_write_jcamp(dic, filename)
%   bruker_write_jcamp(dic, filename, overwrite)
%
% Inputs:
%   dic        - Struct of Bruker parameters (as returned by bruker_read_jcamp)
%   filename   - Full path of the output file
%   overwrite  - (optional) true to overwrite existing file (default: false)
%
% Notes:
%   * The output will differ slightly from Bruker's own files: arrays are
%     always written on a new line after the header, and line breaks may vary.
%   * The special struct fields 'x_coreheader' and 'x_comments' are written
%     first as literal lines; all other fields become ##$KEY= VALUE entries.
%
% See also: bruker_read_jcamp

if nargin < 3, overwrite = false; end

if isfile(filename) && ~overwrite
    warning('bruker_write_jcamp:fileExists', ...
            'File already exists and overwrite=false: %s', filename);
    return;
end

fid = fopen(filename, 'w', 'n', 'UTF-8');
if fid == -1
    error('bruker_write_jcamp:cannotOpen', ...
          'Cannot open file for writing: %s', filename);
end

% ---- write core header lines -------------------------------------------
if isfield(dic, 'x_coreheader')
    for k = 1:numel(dic.x_coreheader)
        fprintf(fid, '%s\n', dic.x_coreheader{k});
    end
end

% ---- write comments -----------------------------------------------------
if isfield(dic, 'x_comments')
    for k = 1:numel(dic.x_comments)
        fprintf(fid, '%s\n', dic.x_comments{k});
    end
end

% ---- write parameter key/value pairs (alphabetical) --------------------
skip = {'x_coreheader', 'x_comments'};
keys = sort(fieldnames(dic));

for ki = 1:numel(keys)
    k = keys{ki};
    if any(strcmp(k, skip)), continue; end
    % Convert valid MATLAB field name back to Bruker parameter name
    % (matlab.lang.makeValidName replaced leading digits / special chars)
    brukerKey = k;   % usually identical for Bruker parameter names
    write_jcamp_pair(fid, brukerKey, dic.(k));
end

fprintf(fid, '##END=\n');
fclose(fid);

end % bruker_write_jcamp

% =========================================================================
function write_jcamp_pair(fid, key, value)
% Write a single ##$KEY= value line (possibly multi-line for arrays).

lineHead = sprintf('##$%s= ', key);

if islogical(value) && isscalar(value)
    if value
        fprintf(fid, '%syes\n', lineHead);
    else
        fprintf(fid, '%sno\n', lineHead);
    end

elseif ischar(value)
    fprintf(fid, '%s<%s>\n', lineHead, value);

elseif isnumeric(value) && isscalar(value)
    fprintf(fid, '%s%s\n', lineHead, num2repr(value));

elseif isnumeric(value) && ~isscalar(value)
    % flatten to 1-D row
    vec = value(:).';
    fprintf(fid, '%s(0..%d)\n', lineHead, numel(vec)-1);
    line = '';
    for vi = 1:numel(vec)
        tok = num2repr(vec(vi));
        if length(line) + length(tok) + 1 > 80
            fprintf(fid, '%s\n', strtrim(line));
            line = '';
        end
        line = [line, tok, ' ']; %#ok<AGROW>
    end
    if ~isempty(strtrim(line))
        fprintf(fid, '%s\n', strtrim(line));
    end

elseif iscell(value)
    fprintf(fid, '%s(0..%d)\n', lineHead, numel(value)-1);
    line = '';
    for vi = 1:numel(value)
        v = value{vi};
        if ischar(v)
            tok = sprintf('<%s>', v);
        elseif isnumeric(v) && isscalar(v)
            tok = num2repr(v);
        else
            tok = mat2str(v);
        end
        if length(line) + length(tok) + 1 > 80
            fprintf(fid, '%s\n', strtrim(line));
            line = '';
        end
        line = [line, tok, ' ']; %#ok<AGROW>
    end
    if ~isempty(strtrim(line))
        fprintf(fid, '%s\n', strtrim(line));
    end

else
    fprintf(fid, '%s%s\n', lineHead, mat2str(value));
end

end % write_jcamp_pair

function s = num2repr(v)
% Return a compact string representation matching Python's repr() style.
if v == floor(v) && abs(v) < 1e15
    s = sprintf('%d', v);
else
    s = sprintf('%.10g', v);
end
end
