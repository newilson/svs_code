function c = hsvdComponents(spec, ppmAx, centerPpm, f0, K, npts, lambda)
% One HSVD decomposition of a spectrum, shared by the metabolite coarse fit and
% the water reference.  lambda = [] uses the hsvd_fit ridge heuristic.
    N   = numel(ppmAx);
    bw  = N * f0 * abs(ppmAx(2) - ppmAx(1));   % Hz; ppm step = bw/(N*f0)
    dt  = 1/bw;
    fid = ifft(ifftshift(spec(:)));
    nu  = min(N, npts);
    tv  = (0:nu-1)' * dt;

    [fHz, dHz, amp, ~, ~, lamUsed] = hsvd_fit(tv, fid(1:nu), K, lambda);

    c = struct('fHz',fHz(:), 'dHz',dHz(:), 'amp',amp(:), ...
               'ppm',centerPpm + fHz(:)/f0, 'lwHz',dHz(:)/pi, ...
               'phaseDeg',angle(amp(:))*180/pi, ...
               'dt',dt, 'N',N, 'npts',nu, 'lambda',lamUsed);
