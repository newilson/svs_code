function dic = bruker_read_jcamp(filename)
% BRUKER_READ_JCAMP  Parse a Bruker JCAMP-DX parameter file into a struct.
%
% Bruker JCAMP-DX files (acqus, acqu2s, procs, …) contain parameters
% prefixed with '##$'.  This function parses those into a MATLAB struct.
%
% Usage:
%   dic = bruker_read_jcamp(filename)
%
% Input:
%   filename  - Full path to a Bruker JCAMP-DX text file (acqus, procs, …)
%
% Output:
%   dic       - Struct of parsed parameters.  Special fields:
%                 dic.x_coreheader  - cell array of '##' header lines
%                 dic.x_comments    - cell array of '$$' comment lines
%               All Bruker parameters '$FOO' become struct fields 'FOO'.
%               Values are char, double scalar, double array, or logical.
%
% Notes:
%   * Not a fully general JCAMP-DX reader – designed for Bruker acqus files.
%   * Attempts UTF-8 encoding first, falls back to CP-1252 (Windows-1252).
%
% See also: bruker_read_acqus, bruker_read

dic = struct('x_coreheader', {{}}, 'x_comments', {{}});

fid = -1;
encodings = {'UTF-8', 'windows-1252'};
for ei = 1:numel(encodings)
    try
        fid = fopen(filename, 'r', 'n', encodings{ei});
        if fid == -1, continue; end
        dic = parse_jcamp_file(fid, dic);
        fclose(fid);
        return;
    catch ME
        if fid ~= -1
            fclose(fid);
            fid = -1;
        end
        if ei == numel(encodings)
            rethrow(ME);
        end
        % try next encoding
    end
end

end % bruker_read_jcamp

% =========================================================================
function dic = parse_jcamp_file(fid, dic)
% Read every line from an open file handle and populate dic.

while true
    rawLine = fgetl(fid);
    if ~ischar(rawLine)      % end of file
        break;
    end
    line = strtrim(rawLine); % remove trailing whitespace / CR

    if isempty(line)
        continue;
    end

    if strncmp(line, '##END=', 6)
        break;

    elseif strncmp(line, '$$', 2)
        dic.x_comments{end+1} = line;

    elseif strncmp(line, '##', 2) && ~strncmp(line, '##$', 3)
        dic.x_coreheader{end+1} = line;

    elseif strncmp(line, '##$', 3)
        try
            [key, value] = parse_jcamp_line(line, fid);
            safeKey = matlab.lang.makeValidName(key);
            dic.(safeKey) = value;
        catch
            warning('bruker_read_jcamp:parseFail', ...
                    'Unable to parse line: %s', line);
        end
    else
        % extraneous line – silently ignore
    end
end

end % parse_jcamp_file

% =========================================================================
function [key, value] = parse_jcamp_line(line, fid)
% Extract parameter name and value from a single '##$FOO= ...' line.
% May read additional continuation lines from fid for arrays / strings.

eqIdx = strfind(line, '=');
key   = strtrim(line(4 : eqIdx(1)-1));   % strip '##$' prefix
text  = strtrim(line(eqIdx(1)+1 : end));

if ~isempty(text) && text(1) == '<'
    % ----- string value: might span multiple lines until '>' --------
    while ~any(text == '>')
        nextLine = fgetl(fid);
        if ~ischar(nextLine), break; end
        text = [text, newline, strtrim(nextLine)]; %#ok<AGROW>
    end
    % strip surrounding < >
    value = text(2 : find(text == '>', 1) - 1);

elseif ~isempty(text) && text(1) == '('
    % ----- array value ----------------------------------------------
    % header looks like (0..N) with optional data on the same line
    dotdot = strfind(text, '..');
    rp     = strfind(text, ')');
    num    = str2double(text(dotdot(1)+2 : rp(1)-1)) + 1;  % expected length

    % data after the closing paren on the same line
    remainder = strtrim(text(rp(1)+1 : end));
    tokens    = strsplit(remainder);
    value     = {};
    for ti = 1:numel(tokens)
        if ~isempty(tokens{ti})
            value{end+1} = parse_jcamp_value(tokens{ti}); %#ok<AGROW>
        end
    end

    % read continuation lines until we have all elements
    while numel(value) < num
        nextLine = fgetl(fid);
        if ~ischar(nextLine), break; end
        tokens = strsplit(strtrim(nextLine));
        for ti = 1:numel(tokens)
            if ~isempty(tokens{ti})
                value{end+1} = parse_jcamp_value(tokens{ti}); %#ok<AGROW>
            end
        end
    end

    % convert to numeric array if all elements are numeric
    isNum = cellfun(@isnumeric, value);
    if all(isNum)
        value = cell2mat(value);
    end

elseif strcmp(text, 'yes')
    value = true;

elseif strcmp(text, 'no')
    value = false;

else
    % ----- simple scalar value --------------------------------------
    value = parse_jcamp_value(text);
end

end % parse_jcamp_line

% =========================================================================
function value = parse_jcamp_value(text)
% Convert a text token to MATLAB scalar (double or char).

if isempty(text)
    value = [];
    return;
end

if text(1) == '<' && text(end) == '>'
    value = text(2:end-1);   % strip angle brackets -> string
    return;
end

% try integer first (no decimal / exponent markers other than integers)
if ~isempty(regexp(text, '^-?\d+$', 'once'))
    value = str2double(text);  % returns double (exact for integers < 2^53)
    return;
end

% try float
if ~isempty(regexp(text, '^-?[\d.]+([eE][+-]?\d+)?$|inf|nan', 'once', 'ignorecase'))
    value = str2double(text);
    if ~isnan(value) || strcmpi(text, 'nan')
        return;
    end
end

% fallback: return as char
value = text;

end % parse_jcamp_value
