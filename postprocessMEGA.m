function output = postprocessMEGA(filenames,pars)
% filenames{1} is EDIT OFF
% filenames{2} is EDIT ON

addpath('utils');
if nargin<2 || ~isfield(pars,'lb')
    pars.lb = 3;
end

if ~isfield(pars,'den')
    pars.den = 'none'; % none, hsvd, wavelet
end

if ~isfield(pars,'peak')
    pars.peak = 3.0; % spectral region of interest (eg 3.0 for GABA)
end

if nargin>=1 && ~iscell(filenames)
    warning('input must be a cell array of filenames')
    filenames = [];
end

if nargin<1 || isempty(filenames)
    [fname,pname] = uigetfile({'*.IMA';'*.dcm';'*.rda'},'Choose Edit OFF Dicom/rda');
    filenames{1} = fullfile(pname,fname);
    [fname,pname] = uigetfile({'*.IMA';'*.dcm';'*.rda'},'Choose Edit ON Dicom/rda');
    filenames{2} = fullfile(pname,fname);
end
output.filenames = filenames;

if ~strcmp(filenames{1}(end-3:end),filenames{2}(end-3:end))
    error('files must have the same extension')
end

npts = []; bw = []; f0 = [];
for ii=1:2
    if contains(filenames{ii},'.rda')
        [output.hdr{ii}, output.complex_fid(:,ii)] = Read_rda_file(filenames{ii});
        if ~isempty(npts)
            if ~isequal(npts,output.hdr{ii}.VectorSize) || ~isequal(bw,1e6/output.hdr{ii}.DwellTime) || ~isequal(f0,output.hdr{ii}.MRFrequency)
                error('scans must have the same parameters')
            end
        end
    	npts = output.hdr{ii}.VectorSize;
    	bw = 1e6/output.hdr{ii}.DwellTime;
    	f0 = output.hdr{ii}.MRFrequency;
    elseif contains(filenames{ii},'.IMA') | contains(filenames{ii},'.dcm')
        [output.complex_fid(:,ii), output.hdr{ii}, output.short_hdr{ii}] = readSiemensDicomSpectrum(filenames{ii});
        if ~isempty(npts)
            if ~isequal(npts,output.short_hdr{ii}.npts) || ~isequal(bw,output.short_hdr{ii}.bw) || ~isequal(f0,output.short_hdr{ii}.f0)
                error('scans must have the same parameters')
            end
        end
    	npts = output.short_hdr{ii}.npts;
    	bw = output.short_hdr{ii}.bw;
    	f0 = output.short_hdr{ii}.f0;
    else
        error('unknown file type')
    end
end

% time axis
t = (0:npts-1)*1/bw;

% freq axes
hz = (-1/2:1/npts:1/2-1/npts)*bw;
ppm = 4.72 - hz/f0;

fid = zeros(npts,2);
spec = zeros(npts,4);
for ii=1:2
    fid(:,ii) = output.complex_fid(:,ii);

    % exponential filtering
    if isfield(pars,'lb')
        fid(:,ii) = expFilter(t,pars.lb,fid(:,ii));
        if pars.lb>0
            output.short_hdr{ii}.flag.lb = true;
        end
    end
git gui
    % remove residual water
    if isfield(pars,'wsopts')
        wsopts = pars.wsopts;
        fid(:,ii) = removeResidualWater(fid(:,ii),wsopts);
        output.short_hdr{ii}.flag.ws = true;
    end

    % svd denoising
    if isfield(pars,'den') && ~strcmp(pars.den,'none')
        if strcmp(pars.den,'hsvd')
            f = NWhsvdden(fid(:,ii),t,bw);
            waitfor(f);
            if exist('out','var')
                fid(:,ii) = out; clear out
                output.short_hdr{ii}.flag.hsvdden = true;
            end
        else
            warning('unknown denoising type not performed')
        end
    end

   

    % FT
    spec(:,ii) = fftshift(fft(fid(:,ii),[],1),1);

end
spec(:,3) = spec(:,1) + spec(:,2);
spec(:,4) = spec(:,2) - spec(:,1);

% phase correction on difference spectrum
if ~isfield(pars,'PC') || ~isnan(pars.PC)
    f = NWman_phase(spec(:,4),ppm,'Manual phase correction: save and close when done');
    waitfor(f);
    if exist('out','var')
        spec(:,4) = out; clear out
        output.short_hdr.flag.pc = true;
        pars.PC = parsPC;
    end
end


output.spec = spec;
output.ppm = ppm;
output.time = t;
output.hz = hz;
output.processed_fid = fid;

