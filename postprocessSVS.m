function output = postprocessSVS(filename,pars)

addpath('utils');

if nargin<2 || ~isfield(pars,'fitmode')
    pars.fitmode = 3; % 1: Gaussian, 2: Lorentzian, 3: Complex Lorentzian, 4: Complex Gaussian
end

if nargin<2 || ~isfield(pars,'peaks')
    pars.peaks = [];
end

if nargin<1 || isempty(filename)
    [fname,pname] = uigetfile({'*.IMA';'*.dcm';'*.rda'},'Choose Spectroscopy Dicom/rda');
    filename = fullfile(pname,fname);
end
output.filename = filename;

if contains(filename,'.rda')
    [output.hdr, output.complex_fid] = Read_rda_file(filename);
	npts = output.hdr.VectorSize;
	bw = 1e6/output.hdr.DwellTime;
	f0 = output.hdr.MRFrequency;
elseif contains(filename,'.IMA') | contains(filename,'.dcm')
    [output.complex_fid, output.hdr, output.short_hdr] = readSiemensDicomSpectrum(filename);
	npts = output.hdr.npts;
	bw = output.hdr.bw;
	f0 = output.hdr.f0;
else
    error('unknown file type')
end

fid = output.complex_fid;

% time axis
t = (0:npts-1)*1/bw;

% freq axes
hz = (-1/2:1/npts:1/2-1/npts)*bw;
ppm = 4.72 - hz/f0;

% filtering
if isfield(pars,'lb')
    fid = expFilter(t,pars.lb,fid);
    if pars.lb>0
        output.short_hdr.flag.lb = true;
    end
end

% svd denoising
if isfield(pars,'den') && strcmp(pars.den,'hsvd')
    f = NWhsvdden(fid,t,bw);
    waitfor(f);
    if exist('out','var')
        fid = out; clear out
        output.short_hdr.flag.hsvdden = true;
    end
end

% remove residual water
if isfield(pars,'wsopts')
    wsopts = pars.wsopts;
    fid = removeResidualWater(fid,wsopts);
    output.short_hdr.flag.ws = true;
end

% FT
spec = fftshift(fft(fid,[],1),1);

% denoising
if isfield(pars,'den')
    switch pars.den
        case 'svd'
            % hsvd denoising
            if size(spec,2)>1
                f = NWsvdTS(spec,ppm,'SVD denoising: save and close when done');
                waitfor(f);
                if exist('out','var')
                    spec = out; clear out
                    output.short_hdr.flag.svdden = true;
                end
            end
        case 'wav'
            % wavelet filtering
            if exist('wmaxlev')
                f = NWwavden(spec,ppm,'Wavelet denoising: save and close when done');
                waitfor(f);
                if exist('out','var')
                    spec = out; clear out
                    output.short_hdr.flag.wavden = true;
                end
            end
    end
end

% auto phase/freq correction
if ~isempty(pars.peaks)
    ranges = 0.1*ones(length(pars.peaks),2);
    addlb = 3; % additional line broadening to locate peaks
    [spec,hzshift] = peakSpectrumShift(spec,hz,f0,addlb,pars.peaks,ranges);
    disp(['hz shift = ' num2str(hzshift)])
    spec = peakSpectrumPhase(spec,ppm,t,addlb,pars.peaks,ranges/2); % tighter ranges here
end

% phase correction
if ~isfield(pars,'PC') || ~isnan(pars.PC)
    f = NWman_phase(spec,ppm,'Manual phase correction: save and close when done');
    waitfor(f);
    if exist('out','var')
        spec = out; clear out
        output.short_hdr.flag.pc = true;
        pars.PC = parsPC;
    end
end

% baseline correction
if ~isfield(pars,'base') || ~isnan(pars.base)
    f = NWsemiman_base(spec,ppm,'Semi-manual baseline correction: save and close when done');
    waitfor(f);
    if exist('out','var')
        spec = out; clear out
        output.short_hdr.flag.base = true;
    end
end

% semi manual fitting
if pars.fitmode>0
    [output.fit.yfit,output.fit.n,output.fit.names,output.fit.ampl,output.fit.pos,output.fit.width,output.fit.integral,output.fit.ip,output.fit.fitpars] = curvefitMan(ppm,real(spec),0,pars.fitmode,pars.peaks);
    output.short_hdr.flag.fit = true;
end

output.spec = spec;
output.ppm = ppm;
output.time = t;
output.hz = hz;
