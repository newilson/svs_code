function [params, data] = bruker_read(dirpath, varargin)
% BRUKER_READ  Read Bruker NMR files (fid/ser + acqus) from a directory.
%
% Usage:
%   [params, data] = bruker_read(dirpath)
%   [params, data] = bruker_read(dirpath, 'Name', Value, ...)
%
% Inputs:
%   dirpath       - Path to experiment directory (containing fid/ser and acqus)
%
% Optional Name-Value pairs:
%   'BinFile'     - Binary filename override (default: auto-detect 'fid'/'ser')
%   'AcqusFiles'  - Cell array of acqus filenames (default: auto-detect)
%   'ReadAcqus'   - true/false (default: true)
%   'Big'         - true = big-endian, false = little-endian, [] = auto
%   'IsFloat'     - true = float64, false = int32, [] = auto from DTYPA
%   'Shape'       - Override data shape as a vector (default: auto)
%   'Cplex'       - Override complex flag (default: auto from AQ_mod)
%
% Outputs:
%   params  - Struct of Bruker parameters (acqus fields at top level,
%             plus params.acqus, params.acqu2s, etc.)
%   data    - Complex (or real) NMR data array
%
% Examples:
%   [p, fid] = bruker_read('/data/nmr/exp1');
%   [p, fid] = bruker_read('/data/nmr/exp1', 'ReadAcqus', true);
%
% See also: bruker_read_acqus, bruker_read_binary, bruker_read_jcamp

% ---- parse inputs -------------------------------------------------------
p = inputParser;
addRequired(p,  'dirpath',    @ischar);
addParameter(p, 'BinFile',    '',    @ischar);
addParameter(p, 'AcqusFiles', {},    @iscell);
addParameter(p, 'ReadAcqus',  true,  @islogical);
addParameter(p, 'Big',        [],    @(x) islogical(x) || isempty(x));
addParameter(p, 'IsFloat',    [],    @(x) islogical(x) || isempty(x));
addParameter(p, 'Shape',      [],    @isnumeric);
addParameter(p, 'Cplex',      [],    @(x) islogical(x) || isempty(x));
parse(p, dirpath, varargin{:});
opts = p.Results;

if ~isfolder(opts.dirpath)
    error('bruker_read:badDir', 'Directory does not exist: %s', opts.dirpath);
end

params = struct();

% ---- find binary file ---------------------------------------------------
binFile = opts.BinFile;
searchDir = opts.dirpath;

if isempty(binFile)
    if isfile(fullfile(searchDir, 'fid'))
        binFile = 'fid';
    elseif isfile(fullfile(searchDir, 'ser'))
        binFile = 'ser';
    else
        % look two levels up (e.g. pdata/1 -> experiment root)
        upDir = fileparts(fileparts(searchDir));
        if isfolder(upDir)
            if isfile(fullfile(upDir, 'fid'))
                binFile  = 'fid';
                searchDir = upDir;
            elseif isfile(fullfile(upDir, 'ser'))
                binFile  = 'ser';
                searchDir = upDir;
            end
        end
    end
end

if isempty(binFile)
    error('bruker_read:noFile', ...
          'No Bruker binary file (fid/ser) found in %s', opts.dirpath);
end

% ---- read acqus parameters ----------------------------------------------
if opts.ReadAcqus
    acqusDict = bruker_read_acqus(searchDir, opts.AcqusFiles);
    % merge all acqus sub-structs into params
    fields = fieldnames(acqusDict);
    for k = 1:numel(fields)
        params.(fields{k}) = acqusDict.(fields{k});
    end
end

% ---- file size ----------------------------------------------------------
binPath = fullfile(searchDir, binFile);
finfo = dir(binPath);
params.FILE_SIZE = finfo.bytes;

% ---- determine endianness -----------------------------------------------
big = opts.Big;
if isempty(big)
    big = false;  % default: little-endian
    if isfield(params, 'acqus') && isfield(params.acqus, 'BYTORDA')
        big = (params.acqus.BYTORDA == 1);
    end
end

% ---- determine data type ------------------------------------------------
isfloat = opts.IsFloat;
if isempty(isfloat)
    isfloat = false;  % default: int32
    if isfield(params, 'acqus') && isfield(params.acqus, 'DTYPA')
        isfloat = (params.acqus.DTYPA == 2);
    end
end

% ---- determine shape and complexity -------------------------------------
shape = opts.Shape;
cplex = opts.Cplex;

if isempty(shape) || isempty(cplex)
    [gshape, gcplex] = bruker_guess_shape(params);
    if isempty(shape),  shape = gshape; end
    if isempty(cplex),  cplex = gcplex; end
end

% ---- read binary data ---------------------------------------------------
data = bruker_read_binary(binPath, shape, cplex, big, isfloat);

end
