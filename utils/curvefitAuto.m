function output = curvefitAuto(x,y,mode,amp0,center0,width0,ph_range,minw)
% curvefitAuto - nonlinear curve fitting for multiple peaks
%
% Inputs:
%   x        - frequency axis (ppm or Hz)
%   y        - spectral data to fit
%   mode     - lineshape model:
%                1 = Gaussian
%                2 = real Lorentzian
%                3 = complex Lorentzian
%                4 = complex Gaussian
%                5 = complex Voigt (default)
%                6 = magnitude Lorentzian
%                7 = magnitude Voigt
%   amp0     - initial amplitudes [n x 1]
%   center0  - initial peak positions [n x 1]
%   width0   - initial linewidths FWHM [n x 1]
%   ph_range - [min max] phase bounds in degrees (modes 3,4,5 only)
%   minw     - minimum linewidth (optional, default 0)
%
% Output:
%   output.pars     - optimized parameters
%   output.parnames - parameter names (describes order of pars)
%   output.fit      - total fitted spectrum [npts x 1]
%   output.comp     - individual peak components [npts x n]
%   output.areas    - peak areas [1 x n]
%   output.lb       - lower bounds
%   output.ub       - upper bounds
%   output.mode     - fitting mode
%
% Parameter order (output.pars):
%   modes 1,2,6: [amp, pos, width]
%   modes 3,4:   [amp, pos, width, phase]
%   mode 5:      [amp, pos, widthL, phase, widthG, width]
%   mode 7:      [amp, pos, widthL, widthG, width]
%
% Note: For modes 5 and 7, 'width' is the effective FWHM calculated using
%       the Olivero-Longbothum approximation from widthL and widthG.
%
% NW - 01/2026 


x = x(:); y = y(:);
np = length(x);
if (length(y) ~= np)
    error('x and y array size mismatch')
end

% NW
if nargin<8 || isempty(minw)
    minw = 0;
end

if (nargin < 3) || isempty(mode)
    mode = 5; % complex Voigt
end

BW = abs(x(end)-x(1));
dwelltime = 1/BW;
% global time
time = col(dwelltime*(0:length(x)-1));

% Initial guesses
if nargin<6 || isempty(amp0) || isempty(center0) || isempty(width0)
    [~,center,width,amp] = findpeaks(y,x);
    figure, plot(x,y,'-k',center,amp,'*r')
    n = length(center);
else
    amp = amp0(:)';
    center = center0(:)';
    width = width0(:)';
    n = length(amp);
end

centermax = center+width/2;
centermin = center-width/2;

amin = 0.2 * amp; %ag is the y-point chosen by user; search 20% of this value
amax = 1.2 * amp; % NW

% Olivero–Longbothum approximation coefficients for Voigt
Acoef = 0.5346; Bcoef = 0.2166;

pars( 1:n ) = amp;
lb(1:n) = amin;
ub(1:n) = amax;

pars( n + (1:n) ) = center;
lb( n + (1:n) ) = centermin;
ub( n + (1:n) ) = centermax;

if mode==5 % complex Voigt
    widthmin = 1e-6*width;
    widthmax = 2*width;

    widthL = 0.4 * width;

    term = (width - Acoef*widthL).^2 - Bcoef*(widthL.^2);
    widthG = sqrt(max(0, term));

    pars( 2*n  + (1:n) ) = widthL;
    lb( 2*n + (1:n) ) = widthmin;
    ub( 2*n + (1:n) ) = widthmax;

    pars( 4*n  + (1:n) ) = widthG;
    lb( 4*n + (1:n) ) = widthmin;
    ub( 4*n + (1:n) ) = widthmax;
elseif mode==7 % magnitude Voigt
    widthmin = 1e-6*width;
    widthmax = 2*width;

    widthL = 0.4 * width;

    term = (width - Acoef*widthL).^2 - Bcoef*(widthL.^2);
    widthG = sqrt(max(0, term));

    pars( 2*n  + (1:n) ) = widthL;
    lb( 2*n + (1:n) ) = widthmin;
    ub( 2*n + (1:n) ) = widthmax;

    pars( 3*n  + (1:n) ) = widthG;  % no phase, so widthG at 3n
    lb( 3*n + (1:n) ) = widthmin;
    ub( 3*n + (1:n) ) = widthmax;
else
    widthmin = 0.5*width;
    widthmax = 2*width;
    for k = 1:n
        if (widthmin(k) < minw)
            widthmin(k) = minw;
        end
    end

    pars( n + n  + (1:n) ) = width;
    lb( n + n + (1:n) ) = widthmin;
    ub( n + n + (1:n) ) = widthmax;
end


% NW
if mode>2 && mode~=6 && mode~=7 % modes 6,7 are magnitude (no phase)
    if nargin<7 || isempty(ph_range) || length(ph_range)~=2
        pars(3*n + (1:n)) = 0; % phase
        lb(3*n+(1:n)) = -179; % degrees
        ub(3*n+(1:n)) = 179;
    else
        pars(3*n + (1:n)) = mean(ph_range);
        lb(3*n+(1:n)) = min(ph_range);
        ub(3*n+(1:n)) = max(ph_range);
    end
end

%%
if size(x)==size(y'), y = y'; end

spins = ones(n,1); % NW fix this later
oldoptions = optimset('lsqcurvefit');
options = optimset(oldoptions, 'TolFun', 1e-12,'TolX', 1e-12,'MaxFunEval',20000*n,'MaxIter', 12000 );

if (mode == 1) % Gaussian

    parnames={'amp','pos','width'};
    ip = lsqcurvefit(@composite,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip( (1:n) + n);
    width = abs(ip( (1:n) + n + n));
    y_comp = zeros(length(x), n);
    for k = 1:n
        temppars = [ampl(k) pos(k) width(k)];
        y_comp(:,k) = composite(temppars, x);
        areas(k) = ampl(k) * width(k) * 0.6006 * sqrt(pi); % NW - fixed 1/15/26
    end
    yfit = sum(y_comp, 2);

elseif mode==2 % Lorentzian

    parnames={'amp','pos','width'};
    ip = lsqcurvefit(@compositel,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip( (1:n) + n);
    width = abs(ip( (1:n) + n + n));
    y_comp = zeros(length(x), n);
    for k = 1:n
        temppars = [ampl(k) pos(k) width(k)];
        y_comp(:,k) = compositel(temppars, x);
        areas(k) = ampl(k) * width(k) * 1.5708;  % Analytical form of lorentzian integral- a*w*pi/2
    end
    yfit = sum(y_comp, 2);

elseif mode==3 % complex Lorentzian - NW

    parnames={'amp','pos','width','ph'};
    ip = lsqcurvefit(@compositel_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip( (1:n) + n);
    width = abs(ip( (1:n) + n + n));
    phase = ip((1:n)+3*n);
    y_comp = zeros(length(x), n);
    for k = 1:n
        temppars = [ampl(k) pos(k) width(k) phase(k)];
        y_comp(:,k) = compositel_complex(temppars, x);
    end
    yfit = sum(y_comp, 2);
    areas = ampl .* width * pi/2; % same as real Lorentzian

elseif mode==4 % complex Gaussian - NW

    parnames={'amp','pos','width','ph'};
    ip = lsqcurvefit(@compositeG_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip((1:n)+n);
    width = abs(ip((1:n)+2*n));
    phase = ip((1:n)+3*n);
    y_comp = zeros(length(x), n);
    for k = 1:n
        temppars = [ampl(k) pos(k) width(k) phase(k)];
        y_comp(:,k) = compositeG_complex(temppars, x);
    end
    yfit = sum(y_comp, 2);
    areas = ampl .* width * 0.6006 * sqrt(pi); % NW - fixed 1/15/26

elseif mode==5 % complex Voigt - NW

    parnames={'amp','pos','widthL','ph','widthG','width'};
    ip = lsqcurvefit(@compositeV_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip((1:n)+n);
    widthL = abs(ip((1:n)+2*n));
    phase = ip((1:n)+3*n);
    widthG = abs(ip((1:n)+4*n));
    width = Acoef * widthL + sqrt(Bcoef*widthL.^2 + widthG.^2);     % Olivero–Longbothum approximation
    ip = [ip width];  % append effective width
    y_comp = zeros(length(x), n);
    areas = ampl .* width; % relative
    for ii=1:n
        temppars = [ampl(ii) pos(ii) widthL(ii) phase(ii) widthG(ii)];
        y_comp(:,ii) = compositeV_complex(temppars,x);
        % evaluate with 0 phase for area integration
        % temppars_0ph = [ampl(ii) pos(ii) widthL(ii) 0 widthG(ii)];
        % peakfit = compositeV_complex(temppars_0ph,x);
        % areas(ii) = sum(peakfit(:));
    end
    yfit = sum(y_comp, 2);

elseif mode==6 % magnitude Lorentzian - NW

    parnames={'amp','pos','width'};
    ip = lsqcurvefit(@compositel_magnitude,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip((1:n)+n);
    width = abs(ip((1:n)+2*n));
    y_comp = zeros(length(x), n);
    % dx = abs(x(2)-x(1));
    areas = ampl .* width; % relative
    for ii=1:n
        temppars = [ampl(ii) pos(ii) width(ii)];
        y_comp(:,ii) = compositel_magnitude(temppars,x);
        % areas(ii) = sum(y_comp(:,ii)) * dx;
    end
    yfit = sum(y_comp, 2);

elseif mode==7 % magnitude Voigt - NW

    parnames={'amp','pos','widthL','widthG','width'};
    ip = lsqcurvefit(@compositeV_magnitude,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip((1:n)+n);
    widthL = abs(ip((1:n)+2*n));
    widthG = abs(ip((1:n)+3*n));
    width = Acoef * widthL + sqrt(Bcoef*widthL.^2 + widthG.^2);     % Olivero–Longbothum approximation
    ip = [ip width];  % append effective width
    y_comp = zeros(length(x), n);
    % dx = abs(x(2)-x(1));
    areas = ampl .* width; % relative
    for ii=1:n
        temppars = [ampl(ii) pos(ii) widthL(ii) widthG(ii)];
        y_comp(:,ii) = compositeV_magnitude(temppars,x);
        % areas(ii) = sum(y_comp(:,ii)) * dx;
    end
    yfit = sum(y_comp, 2);

end

% output
output.pars = ip;
output.lb = lb;
output.ub = ub;
output.areas = areas;
output.fit = yfit;
output.comp = y_comp;
output.parnames = parnames;
output.mode = mode;

return

function yfit = composite(ip,x)

n = (size(ip,2))/3;
yfit = x*0 ;
for k = 1:n
    a = abs( ip(k) );
    p = ip(k+n);
    w = abs( ip(k+n+n) );
    xpw = (x-p)/w;
    yfit = yfit + a * exp( -( (xpw/0.6006).^2 ) ); % 0.6006 = 0.5/(sqrt(ln(2.0)))
end
return

% function yfit = compositeG_complex(ip,x) % NW
% global time
% n = size(ip,2)/4;
% yfit = x*0;
% for k=1:n
%     a = abs( ip(k));
%     p = ip(k+n);
%     w = abs(ip(k+2*n));
%     ph = ip(k+3*n);
%     alpha = 7*0.6006*w/2; % empiric - NW
%     fid = a * exp(1i*(2*pi*p*time + ph*pi/180) - (alpha*time).^2);
%     fid(1) = fid(1)/2;
%     spec = fftshift(fft(fid,[],1),1)/(2/w * 1.33 * 1000/pi); % empiric scaling
%     yfit = yfit + real(spec);
% end
% return

function yfit = compositeG_complex(ip, x)  % NW (analytic Gaussian with phase)
% Robust version using Dawson function (avoids complex erf)
%
% ip = [a_1..a_n, p_1..p_n, w_1..w_n, ph_1..ph_n]
% phase (degrees), FWHM parameterization same as your code
% x is the frequency axis (Hz or ppm)
%
% yfit returns the real-valued spectrum after applying the phase offset.

n = size(ip,2)/4;
yfit = zeros(size(x));
c = 1/0.6006;            % converts FWHM to Gaussian sigma scaling
deg2rad = pi/180;

for k = 1:n
    a  = abs(ip(k));
    p  = ip(k+n);
    w  = abs(ip(k+2*n));
    ph = ip(k+3*n)*deg2rad;

    z = c * (x - p) / w;

    % Real and imaginary (dispersion) parts for the complex Gaussian
    G_real = exp(-z.^2);
    G_imag = (2/sqrt(pi)) * dawson(z);   % Dawson areas, real for real z

    % Combine with phase rotation
    G = G_real + 1i*G_imag;
    L = a * G * exp(-1i*ph);

    yfit = yfit + real(L);  % final real-valued spectrum (like your Lorentzian case)
end
return



function yfit = compositel(ip,x)

n = (size(ip,2))/3;
yfit = x*0 ;
for k = 1:n
    a = abs( ip(k) );
    p = ip(k+n);
    w = abs( ip(k+n+n) );
    xpw = (x-p)/w;
    yfit = yfit + a ./ ( (4.0*(xpw .*xpw)) + 1.0 );
end
return


function yfit = compositel_complex(ip,x) % NW

n = size(ip,2)/4;
yfit = x*0;
for k=1:n
    a = abs( ip(k));
    p = ip(k+n);
    w = abs(ip(k+2*n));
    ph = ip(k+3*n);
    xpw = (x-p)/w;
    L = a ./ (1 + 1i*2.0*xpw) * exp(-1i*pi/180*ph);
    yfit = yfit + real(L);
end
return

function yfit = compositeV_complex(ip, x)
% compositeV_complex : Voigt profile with phase using fadf(z) from File Exchange (ID 47801)
% ip = [a_1..a_n, p_1..p_n, wL_1..wL_n, ph_1..ph_n, wG_1..wG_n]
% x  : frequency axis (Hz or ppm; units must match p, wL, wG)
% Returns the real-valued spectrum: sum_k real( a_k * Vc_norm(x) * exp(-1i*phi_k) )
%
% Requires fadf.m on path (The Voigt/complex error function, second version).
% Reference: w(z) = exp(-z.^2) * erfc(-1i*z) (Faddeeva function).
%
% Normalization: for each peak we divide by Re{w(z0)} at the center so that,
% with phi=0, the peak height at x=p equals 'a'.

    n = size(ip,2)/5;
    if mod(n,1) ~= 0
        error('ip must be packed as [a.., p.., wL.., wG.., ph..].');
    end

    yfit = zeros(size(x));
    deg2rad = pi/180;

    for k = 1:n
        a   = abs(ip(k));
        p   = ip(k+n);
        wL  = abs(ip(k+2*n));   % Lorentzian FWHM
        ph  = ip(k+3*n)*deg2rad;
        wG  = abs(ip(k+4*n));   % Gaussian FWHM

        % Convert to z for Faddeeva:
        % sigma (Gaussian sigma) and gamma (Lorentzian HWHM)
        sigma = wG / (2*sqrt(2*log(2)));
        gamma = wL / 2;

        if sigma == 0 && gamma == 0
            % degenerate; nothing to add
            continue
        elseif sigma == 0
            % pure Lorentzian with phase (analytic fallback)
            xpw = (x - p) / gamma;  % note gamma is HWHM
            Vc  = 1 ./ (1 + 1i*2*xpw);
            % peak center value (phi=0) is 1, so normalization = 1
            Vc_norm = Vc;
        else
            % general Voigt via Faddeeva
            z = ((x - p) + 1i*gamma) / (sigma*sqrt(2));
            Vc = fadf(z);   % complex Voigt lineshape w(z)

            % Normalize so that Re{Vc_norm(x=p)} == 1 when phi=0
            z0 = (1i*gamma) / (sigma*sqrt(2));     % x=p -> (x-p)=0
            V0 = fadf(z0);                          % complex value at center
            N0 = real(V0);
            if ~isfinite(N0) || N0 == 0
                % extremely narrow/wide corner case; fall back to max of real part near center
                [~, idx0] = min(abs(x - p));
                N0 = real(Vc(idx0));
                if ~isfinite(N0) || N0 == 0
                    N0 = 1; % last resort
                end
            end
            Vc_norm = Vc / N0;
        end

        % Apply amplitude and phase, then take real part
        yfit = yfit + real( a * Vc_norm * exp(-1i*ph) );
    end
return

function yfit = compositel_magnitude(ip,x)
% Magnitude Lorentzian: |L(x)| where L(x) = a/(1 + i*2*(x-p)/w)
% |L(x)| = a / sqrt(1 + 4*(x-p)^2/w^2) = a*w/2 / sqrt((w/2)^2 + (x-p)^2)

n = size(ip,2)/3;
yfit = zeros(size(x));
for k=1:n
    a = abs(ip(k));
    p = ip(k+n);
    w = abs(ip(k+2*n));
    xpw = (x-p)/w;
    yfit = yfit + a ./ sqrt(1 + 4.0*(xpw.^2));
end
return

function yfit = compositeV_magnitude(ip,x)
% Magnitude Voigt: |V(x)| where V(x) is complex Voigt
% ip = [a_1..a_n, p_1..p_n, wL_1..wL_n, wG_1..wG_n]

n = size(ip,2)/4;
if mod(n,1) ~= 0
    error('ip must be packed as [a.., p.., wL.., wG..].');
end

yfit = zeros(size(x));

for k = 1:n
    a   = abs(ip(k));
    p   = ip(k+n);
    wL  = abs(ip(k+2*n));   % Lorentzian FWHM
    wG  = abs(ip(k+3*n));   % Gaussian FWHM

    % Convert to z for Faddeeva
    sigma = wG / (2*sqrt(2*log(2)));
    gamma = wL / 2;

    if sigma == 0 && gamma == 0
        % degenerate case
        continue
    elseif sigma == 0
        % pure Lorentzian magnitude
        gamma_hwhm = wL / 2;
        yfit = yfit + a * gamma_hwhm ./ sqrt(gamma_hwhm^2 + (x-p).^2);
    else
        % general Voigt magnitude via Faddeeva
        z = ((x - p) + 1i*gamma) / (sigma*sqrt(2));
        Vc = fadf(z);   % complex Voigt

        % Normalize by peak center value
        z0 = (1i*gamma) / (sigma*sqrt(2));
        V0 = fadf(z0);
        N0 = abs(V0);
        if ~isfinite(N0) || N0 == 0
            [~, idx0] = min(abs(x - p));
            N0 = abs(Vc(idx0));
            if ~isfinite(N0) || N0 == 0
                N0 = 1;
            end
        end
        Vc_norm = Vc / N0;

        yfit = yfit + a * abs(Vc_norm);
    end
end
return

function out = col(in)
out = in(:);
return
