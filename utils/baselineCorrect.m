function [spec_cor,opts] = baselineCorrect(spec,ppm,opts)
%
% spec_cor = baselineCorrect(spec,ppm,algo,opts)
% 
if ~isfield(opts,'algo')
    opts.algo = 2;
elseif ( isnumeric(opts.algo) && opts.algo==0 ) || strcmp(opts.algo,'none')
    spec_cor = spec;
    opts.algo = 0;
    return;
end
if ndims(spec)==1
    spec = spec(:);
end
si = size(spec);

spec = reshape(spec,si(1),[]);
spec_cor = zeros(size(spec));

if nargin<2 || isempty(ppm)
    ppm = 1:si(1);
end
if ~isequal(si(1),length(ppm))
    error('inconsistent dimensions')
end

switch opts.algo
    case 1 % Asymmetric least squares
        
        if ~isfield(opts,'lambda')
            opts.lambda = 10^5;
        end
        if ~isfield(opts,'p')
            opts.p = 10^-2;
        end
        if ~isfield(opts,'maxiter')
            opts.maxiter = 6;
        end
        if ~isfield(opts,'ratio')
            opts.ratio = 10^-3;
        end
        m = si(1);
        D = diff(speye(m), 2);
        w = ones(m, 1);
        lDtD = opts.lambda * (D' * D);
        for ii=1:size(spec,2)
            y = real(spec(:,ii));
            if y~=0*y
                for it = 1:opts.maxiter
                    W = spdiags(w, 0, m, m);
                    C = chol(W + lDtD);
                    z = C \ (C' \ (w .* y));
                    wt = opts.p * (y > z) + (1 - opts.p) * (y < z);
                    if norm(w-wt)/norm(w) < opts.ratio, break; end
                    w = wt;
                end
                if isreal(spec)
                    spec_cor(:,ii) = spec(:,ii) - z;
                else
                    spec_cor(:,ii) = spec(:,ii) - z * (1 + 1i);
                end
            end
        end
        
    case 2 % arPLS - asymmetrically reweighted penalized least squares
        
        if ~isfield(opts,'lambda')
            opts.lambda = 10^6;
        end
        if ~isfield(opts,'maxiter')
            opts.maxiter = 6;
        end
        if ~isfield(opts,'ratio')
            opts.ratio = 10^-3;
        end
        N = si(1);
        D = diff(speye(N),2);
        H = opts.lambda * (D' * D);
        w = ones(N,1);
        for ii=1:size(spec,2)
            y = real(spec(:,ii));
            if y~=0*y
                for jj = 1:opts.maxiter
                    W = spdiags(w,0,N,N);
                    C = chol(W + H);
                    z = C \ (C' \ (w .* y));
                    d = y-z;
                    dn = d(d<0);
                    m = mean(dn);
                    s = std(dn);
                    wt = 1 ./ ( 1 + exp( 2 * (d-(2*s-m))/s ) );
                    if norm(w-wt)/norm(w) < opts.ratio, break; end
                    w = wt;
                end
                if isreal(spec)
                    spec_cor(:,ii) = spec(:,ii) - z;
                else
                    spec_cor(:,ii) = spec(:,ii) - z * (1 + 1i);
                end
            end
        end
        
    case 3

        imspec = [];
        if ~isreal(spec)
            spec = real(spec);
            imspec = imag(spec);
            imspec_cor = zeros(size(spec));
        end
        if ~isfield(opts,'method')
            opts.method = 'pchip';
        end
        if ~isfield(opts,'stepsize')
            opts.stepsize = 0.5;
        end
        if ~isfield(opts,'windowsize')
            opts.windowsize = 0.5;
        end
        
        
        for ii=1:size(spec,2)
            if spec(:,ii)~=0*spec(:,ii)
                spec_cor(:,ii) = msbackadj(ppm(:),spec(:,ii),'RegressionMethod',opts.method,...
                    'StepSize',opts.stepsize,'WindowSize',opts.windowsize);
                if ~isempty(imspec)
                    imspec_cor(:,ii) = msbackadj(ppm(:),spec(:,ii),'RegressionMethod',opts.method,...
                        'StepSize',opts.stepsize,'WindowSize',opts.windowsize);
                end
            end
        end
        if ~isempty(imspec)
            spec_cor = spec_cor + 1i*imspec_cor;
        end
end

spec_cor = reshape(spec_cor,si);
