function data = bruker_remove_digital_filter(params, data, varargin)
% BRUKER_REMOVE_DIGITAL_FILTER  Remove Bruker DSP digital filter from FID.
%
% Bruker spectrometers apply a digital filter during acquisition that
% introduces a group delay. This function corrects for it.
%
% Usage:
%   data = bruker_remove_digital_filter(params, data)
%   data = bruker_remove_digital_filter(params, data, 'Truncate', true)
%
% Inputs:
%   params    - Parameter struct from bruker_read / bruker_read_acqus
%   data      - Complex FID array (last dimension = direct dimension)
%
% Optional Name-Value:
%   'Truncate'  - true (default) to truncate grpdly, false otherwise
%
% Output:
%   data      - FID with digital filter removed
%
% See also: bruker_read

p = inputParser;
addRequired(p, 'params');
addRequired(p, 'data');
addParameter(p, 'Truncate', true, @islogical);
parse(p, params, data, varargin{:});
truncate = p.Results.Truncate;

if ~isfield(params, 'acqus')
    error('bruker_remove_digital_filter:noAcqus', ...
          'params does not contain acqus sub-struct');
end
acq = params.acqus;

if ~isfield(acq, 'DECIM')
    error('bruker_remove_digital_filter:noDECIM', 'DECIM parameter missing');
end
if ~isfield(acq, 'DSPFVS')
    error('bruker_remove_digital_filter:noDSPFVS', 'DSPFVS parameter missing');
end

decim  = acq.DECIM;
dspfvs = acq.DSPFVS;

if isfield(acq, 'GRPDLY')
    grpdly = acq.GRPDLY;
else
    grpdly = 0;
end

data = bruker_rm_dig_filter(data, decim, dspfvs, grpdly, truncate);

end % bruker_remove_digital_filter

% =========================================================================
function data = bruker_rm_dig_filter(data, decim, dspfvs, grpdly, truncate)
% Low-level digital filter removal (mirrors nmrglue rm_dig_filter logic).

% DSP table: bruker_dsp_table{firmware_version}{decimation_rate} = phase_pts
dsp_table = build_dsp_table();

if grpdly > 0
    phase = grpdly;
else
    if isfield(dsp_table, sprintf('v%d', dspfvs))
        fw = dsp_table.(sprintf('v%d', dspfvs));
        key = sprintf('d%d', decim);
        if isfield(fw, key)
            phase = fw.(key);
        else
            error('bruker_rm_dig_filter:unknownDECIM', ...
                  'DECIM=%d not in DSP table for DSPFVS=%d', decim, dspfvs);
        end
    else
        error('bruker_rm_dig_filter:unknownDSPFVS', ...
              'DSPFVS=%d not in DSP table', dspfvs);
    end
end

if truncate
    phase = floor(phase);
end

% number of points to add / skip
add  = max(0, floor(phase));
skip = add + 1;

n = size(data, ndims(data));   % direct dimension length

% frequency shift by 'phase' points using FFT-based approach
data = fft_shift_data(data, phase);

% wrap tail back onto the beginning then truncate
if add > 0 && skip <= n
    idx_add  = 1:add;
    idx_tail = (n - add):(n-1);  % MATLAB 1-based
    if idx_tail(end) <= n
        data_tail = flip_last_dim(data, idx_tail + 1);
        data = add_last_dim(data, idx_add, data_tail);
    end
end

% remove last 'skip' points
if skip <= n
    data = slice_last_dim(data, 1:(n - skip));
end

end % bruker_rm_dig_filter

% ---- helper: FFT-based fractional shift ---------------------------------
function data = fft_shift_data(data, phase)
n    = size(data, ndims(data));
freq = (-n/2 : n/2-1) / n;
freq = ifftshift(freq);
phase_arr = reshape(exp(2j * pi * phase * freq), ones(1, ndims(data)-1)*1, n);
data = ifft(fft(data, [], ndims(data)) .* phase_arr, [], ndims(data));
end

function out = flip_last_dim(data, idx)
nd  = ndims(data);
s   = repmat({':'}, 1, nd);
s{nd} = idx;
chunk = data(s{:});
s{nd} = numel(idx):-1:1;
out   = chunk(s{:});
end

function data = add_last_dim(data, idx_dest, to_add)
nd = ndims(data);
s  = repmat({':'}, 1, nd);
s{nd} = idx_dest;
data(s{:}) = data(s{:}) + to_add;
end

function data = slice_last_dim(data, idx)
nd = ndims(data);
s  = repmat({':'}, 1, nd);
s{nd} = idx;
data = data(s{:});
end

% =========================================================================
function tbl = build_dsp_table()
% Returns a struct mirroring the bruker_dsp_table from nmrglue.
% Keys are firmware version (v10, v11, …) and decimation rate (d2, d4, …).

tbl = struct();

v10.d2    = 44.75;     v10.d3   = 33.5;
v10.d4    = 66.625;    v10.d6   = 59.083333333333333;
v10.d8    = 68.5625;   v10.d12  = 60.375;
v10.d16   = 69.53125;  v10.d24  = 61.020833333333333;
v10.d32   = 70.015625; v10.d48  = 61.34375;
v10.d64   = 70.2578125;v10.d96  = 61.505208333333333;
v10.d128  = 70.37890625;v10.d192= 61.5859375;
v10.d256  = 70.439453125;v10.d384=61.626302083333333;
v10.d512  = 70.4697265625;v10.d768=61.646484375;
v10.d1024 = 70.48486328125;v10.d1536=61.656575520833333;
v10.d2048 = 70.492431640625;
tbl.v10   = v10;

v11.d2    = 46.;       v11.d3   = 36.5;
v11.d4    = 48.;       v11.d6   = 50.166666666666667;
v11.d8    = 53.25;     v11.d12  = 69.5;
v11.d16   = 72.25;     v11.d24  = 70.166666666666667;
v11.d32   = 72.75;     v11.d48  = 70.5;
v11.d64   = 73.;       v11.d96  = 70.666666666666667;
v11.d128  = 72.5;      v11.d192 = 71.333333333333333;
v11.d256  = 72.25;     v11.d384 = 71.666666666666667;
v11.d512  = 72.125;    v11.d768 = 71.833333333333333;
v11.d1024 = 72.0625;   v11.d1536= 71.916666666666667;
v11.d2048 = 72.03125;
tbl.v11   = v11;

v12       = v11;       % same values as v11
tbl.v12   = v12;

v13.d2    = 2.75;      v13.d3   = 2.8333333333333333;
v13.d4    = 2.875;     v13.d6   = 2.9166666666666667;
v13.d8    = 2.9375;    v13.d12  = 2.9583333333333333;
v13.d16   = 2.96875;   v13.d24  = 2.9791666666666667;
v13.d32   = 2.984375;  v13.d48  = 2.9895833333333333;
v13.d64   = 2.9921875; v13.d96  = 2.9947916666666667;
tbl.v13   = v13;

end
