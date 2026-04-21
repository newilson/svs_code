function [cor_fids] = fidfreqshift_fd(time,fids,ppm,ppmRange)
% Frequency-domain analog of fidfreqshift.
% Aligns each FID to the mean spectrum, but only matches points whose
% frequency falls inside ppmRange. This lets the user exclude residual
% water / lipid / other nuisance signal from driving the alignment.
%
% Optimization variables: amp, freq shift (Hz), zero-order phase (rad).
% Correction is applied in the time domain (same expression as
% fidfreqshift), so the returned FIDs are directly compatible.
%
% Inputs:
%   time      [npts x 1]    time vector (s)
%   fids      [npts x navg] complex FIDs
%   ppm       [npts x 1]    ppm axis matching fftshift(fft(fid,[],1),1)
%   ppmRange  [1 x 2]       [ppm_lo ppm_hi] over which to match
%                           (default = full ppm axis)
%
% Output:
%   cor_fids  [npts x navg] freq/phase-corrected FIDs (amp scaling removed)
%
% NW 04/2026

if (nargin < 4) || isempty(ppmRange), ppmRange = [min(ppm) max(ppm)]; end

p = gcp('nocreate');
doPar = ~isempty(p);

fids = double(fids);
time = time(:);
ppm  = ppm(:);

mask  = (ppm >= min(ppmRange)) & (ppm <= max(ppmRange));
nroi  = nnz(mask);
if nroi < 4
    error('fidfreqshift_fd: ppmRange selects fewer than 4 points');
end

% Reference = mean spectrum, restricted to ROI
specref     = fftshift(fft(mean(fids,2),[],1),1);
specref_roi = specref(mask);
ref_cat     = [real(specref_roi); imag(specref_roi)];

% lsqcurvefit needs xdata with length matching ydata; values are unused
xdummy = (1:2*nroi).';

opts = optimoptions(@lsqcurvefit,'Display','off');
par0 = [1 0 0]; % [amp, freq Hz, phase rad]

cor_fids = zeros(size(fids));

if doPar
    parfor ii = 1:size(fids,2)
        fid_i   = fids(:,ii);
        curFunc = @(p,~) fidfreqshift_fd_func(p,time,fid_i,mask);
        pars    = lsqcurvefit(curFunc,par0,xdummy,ref_cat,[],[],opts);
        cor_fids(:,ii) = apply_shift(fid_i,time,pars) / pars(1);
    end
else
    for ii = 1:size(fids,2)
        fid_i   = fids(:,ii);
        curFunc = @(p,~) fidfreqshift_fd_func(p,time,fid_i,mask);
        pars    = lsqcurvefit(curFunc,par0,xdummy,ref_cat,[],[],opts);
        cor_fids(:,ii) = apply_shift(fid_i,time,pars) / pars(1);
    end
end

end

% ------------------------------------------------------------------------
function out = fidfreqshift_fd_func(p,t,fid,mask)
fid_s    = p(1) .* exp(-1i*p(3)) .* exp(-1i*p(2)*2*pi*t) .* fid;
spec_s   = fftshift(fft(fid_s,[],1),1);
spec_roi = spec_s(mask);
out      = [real(spec_roi); imag(spec_roi)];
end

% ------------------------------------------------------------------------
function fid_c = apply_shift(fid,t,p)
fid_c = p(1) .* exp(-1i*p(3)) .* exp(-1i*p(2)*2*pi*t) .* fid;
end
