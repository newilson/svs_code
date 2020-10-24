function [diff_spec, spec2corr, pars] = NWoptimizeDiffSpec(spec1,fid2,spec_mask,pars0,lb,ub)

oldoptions = optimset('lsqnonlin');
options = optimset(oldoptions, 'TolFun', 1e-12,'TolX', 1e-12,'MaxFunEval',20000*10,'MaxIter', 3000, 'Display','iter' );

spec1 = spec1(:);
fid2 = fid2(:);
t = 1:length(fid2);
t = t(:);

% anonymous function
f = @(pars)mydiff(pars,fid2,spec1,spec_mask);

pars = lsqnonlin(f,pars0,lb,ub,options);

filt = exp(pi*pars(1)*t/length(fid2)).*exp(-pars(2)*t/length(fid2));
fid2corr = filt.*fid2;
spec2corr = fftshift(fft(fid2corr,[],1),1);
spec2corr = fraccircshift(spec2corr,[pars(3),0]);
spec2corr = pars(4)*spec2corr;
phase = exp(-1i*(pars(5)+pars(6)*(t-pars(7))));
spec2corr = spec2corr .* phase;

diff_spec = real(spec1-spec2corr);
end


function diff = mydiff(pars,fid2,spec1,spec_mask)
% Note the mixed input. All processing done on the second spectrum

% allow for lineshape changes and linear and 1st order phase changes
t = 1:length(fid2);
t = t(:);

% Voigt filter fid
filt = exp(-pi*pars(1)*t/length(fid2)).*exp(-pars(2)*t/length(fid2));
fid2 = filt.*fid2;

% fft
spec = fftshift(fft(fid2,[],1),1);

% shift
spec = fraccircshift(spec,[pars(3) 0]);

% scale
spec = pars(4) * spec;

% phase correct spectrum
phase = exp(-1i*(pars(5) + pars(6)*(t-pars(7))));
spec = spec .* phase;

% take real part and mask
diff = spec_mask .* real(spec1-spec);

end