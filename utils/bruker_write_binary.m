function bruker_write_binary(filename, data, big, isfloat, overwrite)
% BRUKER_WRITE_BINARY  Write NMR data to a Bruker binary file (fid / ser).
%
% Usage:
%   bruker_write_binary(filename, data, big, isfloat)
%   bruker_write_binary(filename, data, big, isfloat, overwrite)
%
% Inputs:
%   filename   - Full path to the output file
%   data       - Numeric array (complex or real).
%                Complex → written as interleaved real/imag (int32 or float64).
%                Real    → written directly as int32 or float64.
%   big        - true for big-endian, false for little-endian
%   isfloat    - true for float64 (DTYPA=2), false for int32 (DTYPA=0)
%   overwrite  - (optional) true to allow overwriting (default: false)
%
% See also: bruker_read_binary, bruker_read

if nargin < 5, overwrite = false; end

if isfile(filename) && ~overwrite
    warning('bruker_write_binary:fileExists', ...
            'File exists and overwrite=false: %s', filename);
    return;
end

if big
    machineFormat = 'ieee-be';
else
    machineFormat = 'ieee-le';
end

fid = fopen(filename, 'wb', machineFormat);
if fid == -1
    error('bruker_write_binary:cannotOpen', ...
          'Cannot open file for writing: %s', filename);
end

% ---- uncomplexify if complex -------------------------------------------
if ~isreal(data)
    raw = uncomplexify(data, isfloat);
else
    if isfloat
        raw = double(data);
    else
        raw = int32(data);
    end
end

if isfloat
    fwrite(fid, raw(:), 'float64');
else
    fwrite(fid, raw(:), 'int32');
end

fclose(fid);

end % bruker_write_binary

% =========================================================================
function out = uncomplexify(data, isfloat)
% Pack complex array into interleaved [R0 I0 R1 I1 ...] array.
sz      = size(data);
sz(end) = sz(end) * 2;
if isfloat
    out = zeros(sz, 'double');
else
    out = zeros(sz, 'int32');
end

nd  = ndims(out);
sR  = repmat({':'}, 1, nd);  sI = sR;
sR{nd} = 1:2:sz(end);
sI{nd} = 2:2:sz(end);

if isfloat
    out(sR{:}) = double(real(data));
    out(sI{:}) = double(imag(data));
else
    out(sR{:}) = int32(real(data));
    out(sI{:}) = int32(imag(data));
end

end
