function [complexfid, seqpar, nmid] = read_lcmodelRAW(filepath)
% read_lcmodelRAW  Read an LCModel .RAW file.
%
%   [complexfid, seqpar, nmid] = read_lcmodelRAW(filepath)
%
%   Inputs:
%     filepath   - Full path to the .RAW file
%
%   Outputs:
%     complexfid - Nx1 complex time-domain FID vector
%     seqpar     - struct with $SEQPAR fields (HZPPPM, ECHOT, SEQ, NUNFIL, DELTAT)
%                  Empty struct if no $SEQPAR block is present.
%     nmid       - struct with $NMID fields (ID, FMTDAT, TRAMP, VOLUME)
%
%   Example:
%     [fid, seqpar, nmid] = read_lcmodelRAW('/data/scan.RAW');
%     plot(real(fid));

    if nargin < 1
        [file,path] = uigetfile('*.RAW','Choose RAW file');
        filepath = fullfile(path,file);
    end

    fid_in = fopen(filepath, 'r');
    if fid_in == -1
        error('read_lcmodelRAW:fileOpen', 'Cannot open file: %s', filepath);
    end
    raw = fread(fid_in, '*char')';
    fclose(fid_in);

    % --- Parse $SEQPAR namelist ---
    seqpar = struct();
    seqpar_match = regexp(raw, '\$SEQPAR(.*?)\$END', 'tokens', 'once');
    if ~isempty(seqpar_match)
        block = seqpar_match{1};
        seqpar.HZPPPM = parse_numeric(block, 'HZPPPM');
        seqpar.ECHOT  = parse_numeric(block, 'ECHOT');
        seqpar.SEQ    = parse_string(block, 'SEQ');
        seqpar.NUNFIL = parse_numeric(block, 'NUNFIL');
        seqpar.DELTAT = parse_numeric(block, 'DELTAT');
    end

    % --- Parse $NMID namelist ---
    nmid = struct();
    nmid_match = regexp(raw, '\$NMID(.*?)\$END', 'tokens', 'once');
    if ~isempty(nmid_match)
        block = nmid_match{1};
        nmid.ID     = parse_string(block, 'ID');
        nmid.FMTDAT = parse_string(block, 'FMTDAT');
        nmid.TRAMP  = parse_numeric(block, 'TRAMP');
        nmid.VOLUME = parse_numeric(block, 'VOLUME');
    end

    % --- Find where data begins (after last $END) ---
    end_inds = regexp(raw, '\$END');
    data_str = raw(end_inds(end)+4:end);

    % --- Read time-domain data ---
    vals = sscanf(data_str, '%e');
    if mod(length(vals),2) ~= 0
        warning('read_lcmodelRAW:oddValues', 'Odd number of values in data block');
    end
    data_re = vals(1:2:end);
    data_im = vals(2:2:end);
    complexfid = data_re - 1i*data_im; % undo sign convention from create_lcmodelRAW

    fprintf('Read %s (%d complex points)\n', filepath, length(complexfid));

end


function val = parse_numeric(block, name)
    tok = regexp(block, [name '\s*=\s*([^\s,\n]+)'], 'tokens', 'once');
    if ~isempty(tok)
        val = str2double(tok{1});
    else
        val = [];
    end
end

function val = parse_string(block, name)
    tok = regexp(block, [name '\s*=\s*''([^'']*)'''], 'tokens', 'once');
    if ~isempty(tok)
        val = tok{1};
    else
        val = '';
    end
end
