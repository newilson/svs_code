function [cor_fids] = fidfreqshift(time,fids,Nfit,skippts)
% Find freq shift in Hz and phase to match FIDS to average
% Fits amp, phase, and frequency
%
% NW 03/2026 from ME originally


% check if parallel pools
p = gcp('nocreate');
if isempty(p)
    doPar = false;
else
    doPar = true;
end

if (nargin < 4) || isempty(skippts), skippts = 0; end           % ignore initial FID points
if (nargin < 3) || isempty(Nfit), Nfit    = size(fids,1) - skippts; end % fit to only the first Nfit points if desired

fids = double(fids);

orig_fids = fids; % original, full fids
cor_fids = 0 * orig_fids;
fids       = fids((1+skippts):(Nfit+skippts),:);                                                       % truncated
fidref = mean(fids,2); % fids must be [npts x avgs]
fidref_cat = [real(fidref) ; imag(fidref)]; % vector

time = time(:); % must be a column vector
tr_time       = time((1+skippts):(Nfit+skippts)); % truncated
tr_time = tr_time(:); % column vector

opts = optimoptions(@lsqcurvefit,'Display','off');
par0 = [1 0 0]; % [amp, freq, phase]

if doPar
    parfor ii=1:size(fids,2)

        % fit
        curFunc = @(p,t) fidfreqshift_func(p,t,fids(:,ii)); % anonymous function
        [pars,rnorm] = lsqcurvefit(curFunc,par0,tr_time,fidref_cat,[],[],opts);

        % correct
        corFunc = @(p,t) fidfreqshift_func(p,t,orig_fids(:,ii));
        corfid_cat = corFunc(pars,time) / pars(1); % take off amplitude scaling
        cor_fids(:,ii) = complex(corfid_cat(1:end/2), corfid_cat(end/2+1:end));

    end

else

    for ii=1:size(fids,2)

        % fit
        curFunc = @(p,t) fidfreqshift_func(p,t,fids(:,ii)); % anonymous function
        [pars,rnorm] = lsqcurvefit(curFunc,par0,tr_time,fidref_cat,[],[],opts);

        % correct
        corFunc = @(p,t) fidfreqshift_func(p,t,orig_fids(:,ii));
        corfid_cat = corFunc(pars,time)/ pars(1); % take off amplitude scaling
        cor_fids(:,ii) = complex(corfid_cat(1:end/2), corfid_cat(end/2+1:end));

    end

end


end

% ------------------------------------------------------------------------
function fid_cat = fidfreqshift_func(p,t,fid)
fidguess  = p(1) .* exp(-1i*p(3)) .* exp(-1i*p(2)*2*pi*t) .* fid;
fid_cat   = [real(fidguess)  ; imag(fidguess)];
end

