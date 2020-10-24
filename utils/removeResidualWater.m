function fidws = removeResidualWater(fid,opts)
%
% fidws = removeResidualWater(fid,opts)
%
% fid is [npts x whatever]

if nargin<2 || isempty(opts)
    opts.type = 'hsvd'; % 'wavelet'
    opts.plt = false;
end
if strcmp(opts.type,'wavelet')
    addpath(genpath('Wavelab850'));
    if ~isfield(opts,'wavelet') || ~isfield(opts.wavelet,'zf')
        opts.wavelet.zf = 0;
    end
    if ~isfield(opts.wavelet,'scale')
        opts.wavelet.scale = 5;
    end
    if ~isfield(opts.wavelet,'type')
        opts.wavelet.type = 'Daubechies';
    end
    if ~isfield(opts.wavelet,'par')
        opts.wavelet.par = 10;
    end
elseif strcmp(opts.type,'hsvd')
    if ~isfield(opts.hsvd,'bounds')
        opts.hsvd.bounds = [-0.025 0.025];
    end
    if ~isfield(opts.hsvd,'nsin')
        opts.hsvd.nsin = 25;
    end
end

si = size(fid);
fid = reshape(fid,si(1),[]);
fidws = zeros(size(fid));
switch opts.type
    case 'hsvd'
        fidws = transpose(h2osup_big(transpose(fid), opts.hsvd.bounds, opts.hsvd.nsin));
    case 'wavelet'
        for ii=1:size(fid,2)
            fidws(:,ii) = complex(wavewat(real(fid(:,ii)),0,3),wavewat(imag(fid(:,ii)),opts.wavelet.zf,opts.wavelet.scale,opts.wavelet.type,opts.wavelet.par));
        end
end
if opts.plt
    NWplayplot(real(fftshift(fft(fid,[],1),1)),[],'water removal',real(fftshift(fft(fidws,[],1),1)),[]);
end
fidws = reshape(fidws,si);

