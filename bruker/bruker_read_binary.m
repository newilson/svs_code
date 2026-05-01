function data = bruker_read_binary(filename, shape, cplex, big, isfloat)
% BRUKER_READ_BINARY  Read a Bruker binary FID or SER file into MATLAB.
%
% Usage:
%   data = bruker_read_binary(filename, shape, cplex, big, isfloat)
%
% Inputs:
%   filename  - Full path to the binary file (fid or ser)
%   shape     - Row vector of dimension sizes, e.g. [TD2/2, TD/2] for 2D.
%               The *last* element is the number of complex points in the
%               direct dimension (R+I pairs have already been divided by 2
%               before calling this function if cplex is true).
%   cplex     - true  → direct dimension is complex (interleaved R,I)
%               false → real-only data
%   big       - true  → big-endian byte order
%               false → little-endian byte order
%   isfloat   - true  → float64 (DTYPA=2)
%               false → int32  (DTYPA=0)
%
% Output:
%   data      - Numeric array.
%               complex128 when cplex=true, double/int32 when cplex=false.
%               Shape matches the requested 'shape' argument where possible.
%
% Notes:
%   * Bruker pads each trace to 1024-byte boundaries (256 int32 points or
%     128 float64 points).  This function reads all points and reshapes.
%   * If reshaping fails the full 1-D vector is returned with a warning.
%
% See also: bruker_read, bruker_guess_shape

narginchk(5, 5);

% ---- pick format string -------------------------------------------------
if isfloat
    precision = 'float64';   % 8 bytes per point
else
    precision = 'int32';     % 4 bytes per point
end

if big
    machineFormat = 'ieee-be';
else
    machineFormat = 'ieee-le';
end

% ---- read raw data from file --------------------------------------------
fid = fopen(filename, 'rb', machineFormat);
if fid == -1
    error('bruker_read_binary:cannotOpen', ...
          'Cannot open file: %s', filename);
end
raw = fread(fid, Inf, precision);
fclose(fid);

% convert to double for subsequent arithmetic
raw = double(raw);

% ---- interleave real/imag -> complex ------------------------------------
if cplex
    data = complexify(raw);   % length(raw)/2 complex points
else
    data = raw;
end

% ---- reshape ------------------------------------------------------------
if isempty(shape) || isequal(shape, 0)
    % no shape info – return 1-D
    return;
end

try
    data = reshape(data, fliplr(shape));   % MATLAB is column-major
    % permute so that dimension ordering matches Python/nmrglue convention:
    % first index = slowest, last index = fastest (direct)
    ndim = numel(shape);
    if ndim > 1
        data = permute(data, ndim:-1:1);
    end
catch
    % try reshape as 2D (nTraces x nDirectPoints)
    nDirect = shape(end);
    try
        data = reshape(data, [], nDirect);
        warning('bruker_read_binary:reshapeFallback', ...
            'Could not reshape to %s; returned 2D (%d x %d).', ...
            mat2str(shape), size(data,1), size(data,2));
    catch
        warning('bruker_read_binary:reshape1D', ...
            'Could not reshape data; returning 1D array.');
    end
end

end % bruker_read_binary

% =========================================================================
function c = complexify(real_imag)
% Convert interleaved [R0 I0 R1 I1 ...] vector to complex array.
real_imag = real_imag(:).';        % ensure row vector
c = complex(real_imag(1:2:end), real_imag(2:2:end));
end
