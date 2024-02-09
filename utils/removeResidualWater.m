function fidws = removeResidualWater(fid,opts)
%
% fidws = removeResidualWater(fid,opts)
%
% fid is [npts x whatever]

fid = double(fid); % svd requires double

if nargin<2 || isempty(opts)
    opts.type = 'hsvd'; % 'wavelet'
    opts.hsvd = [];
    opts.plt = false;
end
if ~isfield(opts,'plt')
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
elseif strcmp(opts.type,'filt')
    if ~isfield(opts.filt,'N')
        opts.filt.N = 20;
    end
    if ~isfield(opts.filt,'type')
        opts.filt.type = 'average'; % average or gaussian
    end
    if strcmp(opts.filt.type,'gaussian') && ~isfield(opts.filt,'sigma')
        opts.filt.sigma = 5;
    end
elseif strcmp(opts.type,'none')
    fidws = fid;
    return;
end

si = size(fid);
fid = reshape(fid,si(1),[]);
fidws = zeros(size(fid));
switch opts.type
    case 'hsvd'
        fidws = transpose(h2osup_big(transpose(fid), opts.hsvd.bounds, opts.hsvd.nsin));
    case 'wavelet'
        for ii=1:size(fid,2)
            if fid(:,ii)~=0*fid(:,ii)
                temp = complex(wavewat(real(fid(:,ii)),0,3),wavewat(imag(fid(:,ii)),opts.wavelet.zf,opts.wavelet.scale,opts.wavelet.type,opts.wavelet.par));
                fidws(:,ii) = temp(1:si(1),:); % if not a multiple of 2
            end
        end
    case 'filt'
        switch opts.filt.type
            case 'average'
                h = fspecial('average',[opts.filt.N 1]);
            case 'gaussian'
                h = fspecial('gaussian',[opts.filt.N 1],opts.filt.sigma);
            otherwise
                error('unimplemented filter type')
        end
        for ii=1:size(fid,2)
            if fid(:,ii)~=0*fid(:,ii)
                lowfreq = imfilter(real(fid(:,ii)),h,'symmetric','same') + 1i * imfilter(imag(fid(:,ii)),h,'symmetric','same');
                fidws(:,ii) = fid(:,ii) - lowfreq;
            end
        end
end
if opts.plt
    NWplayplot(real(fftshift(fft(fid,[],1),1)),[],'water removal',real(fftshift(fft(fidws,[],1),1)),[]);
end
fidws = reshape(fidws,si);

