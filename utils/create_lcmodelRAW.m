function create_lcmodelRAW(filepath, complexfid, varargin)
% create_lcmodelRAW  Write an LCModel-compatible .RAW file.
%
%   create_lcmodelRAW(filepath, complexfid)
%   create_lcmodelRAW(filepath, complexfid, 'Name', Value, ...)
%
%   Inputs:
%     filepath   - Full path for the output .RAW file (e.g., '/data/myscan.RAW')
%     complexfid - Complex time-domain FID vector (column or row)
%
%   Optional Name-Value Parameters:
%     SEQPAR header (all optional; omitted entirely if none are provided):
%       'HZPPPM'  - Field strength in proton MHz (REAL), e.g. 297.2 for 7T
%       'ECHOT'   - Echo time in ms (REAL)
%       'SEQ'     - Localization sequence, 'PRESS' or 'STEAM' (CHARACTER*5)
%       'NUNFIL'  - Number of complex data points (INTEGER)
%       'DELTAT'  - Dwell time in seconds (REAL)

%     NMID header:
%       'ID'      - Identifier string, max 20 chars (default: filename stem)
%       'FMTDAT'  - Fortran format spec (default: '(2E14.5)')
%       'TRAMP'   - Amplitude scaling factor (default: 1.0)
%       'VOLUME'  - Voxel size in mL (default: 1.0)
%
%   Example:
%     create_lcmodelRAW('/data/scan.RAW', fid, 'HZPPPM', 49.9815, ...
%       'NUNFIL', 1024, 'DELTAT', 1.040e-04)

    p = inputParser;
    addRequired(p, 'filepath', @ischar);
    addRequired(p, 'complexfid', @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'HZPPPM', [], @(x) isscalar(x) && isnumeric(x));
    addParameter(p, 'ECHOT', [], @(x) isscalar(x) && isnumeric(x));
    addParameter(p, 'SEQ', '', @ischar);
    addParameter(p, 'NUNFIL', [], @(x) isscalar(x) && isnumeric(x));
    addParameter(p, 'DELTAT', [], @(x) isscalar(x) && isnumeric(x));
    addParameter(p, 'ID', '', @ischar);
    addParameter(p, 'FMTDAT', '(2E14.5)', @ischar);
    addParameter(p, 'TRAMP', 1.0, @(x) isscalar(x) && isnumeric(x));
    addParameter(p, 'VOLUME', 1.0, @(x) isscalar(x) && isnumeric(x));
    parse(p, filepath, complexfid, varargin{:});

    opts = p.Results;

    % Force column vector
    complexfid = complexfid(:);

    % Derive default ID from filename if not provided
    if isempty(opts.ID)
        [~, opts.ID] = fileparts(filepath);
        if length(opts.ID) > 20
            opts.ID = opts.ID(1:20);
        end
    end

    % create output directory if it does not exist
    outdir = fileparts(filepath);
    if ~isempty(outdir) && ~isfolder(outdir)
        mkdir(outdir);
    end

    fid_out = fopen(filepath, 'w');
    if fid_out == -1
        error('create_lcmodelRAW:fileOpen', 'Cannot open file: %s', filepath);
    end

    % --- Write $SEQPAR namelist (only if any SEQPAR parameters are provided) ---
    has_seqpar = ~isempty(opts.HZPPPM) || ~isempty(opts.ECHOT) || ~isempty(opts.SEQ) ...
        || ~isempty(opts.NUNFIL) || ~isempty(opts.DELTAT);
    if has_seqpar
        fprintf(fid_out, ' $SEQPAR\n');
        if ~isempty(opts.HZPPPM)
            fprintf(fid_out, ' HZPPPM=%g\n', opts.HZPPPM);
        end
        if ~isempty(opts.ECHOT)
            fprintf(fid_out, ' ECHOT=%g\n', opts.ECHOT);
        end
        if ~isempty(opts.SEQ)
            fprintf(fid_out, ' SEQ=''%s''\n', opts.SEQ);
        end
        if ~isempty(opts.NUNFIL)
            fprintf(fid_out, ' NUNFIL=%d\n', opts.NUNFIL);
        end
        if ~isempty(opts.DELTAT)
            fprintf(fid_out, ' DELTAT=%g\n', opts.DELTAT);
        end
        fprintf(fid_out, ' $END\n');
    end

    % --- Write $NMID namelist ---
    fprintf(fid_out, ' $NMID ID=''%s'', FMTDAT=''%s''\n', opts.ID, opts.FMTDAT);
    fprintf(fid_out, ' TRAMP=%g, VOLUME=%g $END\n', opts.TRAMP, opts.VOLUME);

    % --- Write time-domain data ---
    data_re = real(complexfid);
    data_im = -imag(complexfid);
    fprintf(fid_out, '%14.5E%14.5E\n', [data_re'; data_im']);

    fclose(fid_out);

    fprintf('Wrote %s (%d complex points)\n', filepath, length(complexfid));

end
