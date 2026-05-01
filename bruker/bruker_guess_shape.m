function [shape, cplex] = bruker_guess_shape(params)
% BRUKER_GUESS_SHAPE  Infer data shape and complexity from acqus parameters.
%
% Usage:
%   [shape, cplex] = bruker_guess_shape(params)
%
% Input:
%   params  - Parameter struct returned by bruker_read_acqus (contains
%             params.acqus, params.acqu2s, etc.) or any struct that has
%             a FILE_SIZE field and acqus sub-struct.
%
% Outputs:
%   shape  - Row vector of dimension sizes (number of *complex* points in
%            the direct dimension, raw points in indirect dimensions).
%            The last element corresponds to the direct (fastest) dimension.
%   cplex  - true when the direct dimension is complex (AQ_mod 1 or 3).
%
% The shape follows the same convention as nmrglue / NMRPipe:
%   1D:   [nDirect]
%   2D:   [nIndirect1, nDirect]
%   3D:   [nIndirect2, nIndirect1, nDirect]
%
% See also: bruker_read, bruker_read_acqus

% ---- direct dimension complexity ----------------------------------------
try
    aq_mod = params.acqus.AQ_mod;
catch
    aq_mod = 0;
end

if ismember(aq_mod, [0, 2])
    cplex = false;
elseif ismember(aq_mod, [1, 3])
    cplex = true;
else
    error('bruker_guess_shape:unknownAQmod', ...
          'Unknown AQ_mod value: %d', aq_mod);
end

% ---- file size ----------------------------------------------------------
if ~isfield(params, 'FILE_SIZE')
    warning('bruker_guess_shape:noFileSize', ...
            'FILE_SIZE not in params; returning default shape.');
    shape = 1;
    cplex = true;
    return;
end
fsize = params.FILE_SIZE;

% ---- data type (determines padding block size) --------------------------
try
    dtypa = params.acqus.DTYPA;
catch
    dtypa = 0;
end

if dtypa == 0
    bytesize = 4;    % int32
elseif dtypa == 2
    bytesize = 8;    % float64
else
    error('bruker_guess_shape:unknownDTYPA', ...
          'Unexpected DTYPA value: %d', dtypa);
end

% ---- time-domain sizes from acqus / acqu2s / acqu3s / acqu4s -----------
try, td0 = double(params.acqus.TD);  catch, td0 = 1024; end
try, td2 = double(params.acqu2s.TD); catch, td2 = 0;    end
try, td1 = double(params.acqu3s.TD); catch, td1 = td2;  end
try, td3 = double(params.acqu4s.TD); catch, td3 = td1;  end

% ---- direct dimension size (padded to 1024-byte blocks) ----------------
pointsize  = 1024.0 / bytesize;
td0_padded = ceil(td0 / pointsize) * pointsize;   % raw (R+I) points

% ---- build shape vector [slowest … fastest] ----------------------------
% Start with placeholders for up to 4D data.
% shape(4) = direct dim (R+I points on disk)
% shape(3) = acqu2s TD
% shape(2) = derived from file size
% shape(1) = derived from file size

shapeVec = zeros(1, 4);
shapeVec(4) = td0_padded;
shapeVec(3) = td2;

if shapeVec(3) ~= 0 && shapeVec(4) ~= 0
    shapeVec(2) = fsize / (shapeVec(4) * shapeVec(3) * bytesize);
    shapeVec(1) = fsize / (shapeVec(4) * shapeVec(3) * max(shapeVec(2),1) * bytesize);
end

% Remove leading dimensions that are <= 1
shapeVec = shapeVec(shapeVec > 1);
if isempty(shapeVec)
    shapeVec = td0_padded;
end

% ---- convert direct dim to complex points if cplex ----------------------
% The caller (bruker_read) divides by 2; here we just return the raw shape.
% Match Python nmrglue behaviour: caller halves the last dim if cplex.
shape = shapeVec;

end
