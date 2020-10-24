function spec_cor = baselineCorrect(spec,ppm,opts)
%
% spec_cor = baselineCorrect(spec,ppm,opts)
%

si = size(spec);
imspec = [];
if ~isreal(spec)
    spec = real(spec);
    imspec = imag(spec);
    imspec_cor = zeros(size(spec));
end
if ~isequal(si(1),length(ppm))
    error('inconsistent dimensions')
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

spec = reshape(spec,si(1),[]);
spec_cor = zeros(size(spec));
for ii=1:size(spec,2)
    spec_cor(:,ii) = msbackadj(ppm(:),spec(:,ii),'RegressionMethod',opts.method,...
        'StepSize',opts.stepsize,'WindowSize',opts.windowsize);
    if ~isempty(imspec)
        imspec_cor(:,ii) = msbackadj(ppm(:),spec(:,ii),'RegressionMethod',opts.method,...
            'StepSize',opts.stepsize,'WindowSize',opts.windowsize);
    end
end
if ~isempty(imspec)
    spec_cor = spec_cor + 1i*imspec_cor;
end    
spec_cor = reshape(spec_cor,si);
