function [yfit,areas,ip,pars,lb,ub] = curvefitWat(x,y,minw,mode)
% water spectrum should be mostly phase corrected already

%[yfit,ampl,pos,width,integral, pars, n] =
%curvefit(x,y,minw,n,mode,miny,maxy)

% nonlinear curve fitting program for n peaks
% inputs:    x = array of xpoints
%            y = array of spectral points
%            n = number of peaks expected
%            mode 1= gaussian, 2 = lorentzian, 3 = complex lorentzian (NW),
%            4 = complex gaussian (NW)
%            ag0, cg0, wg0 are initial guesses for amplitude, center, and
%            width of peaks
%            ph_range [2x1] sets bounds for the phases in complex fitting 


x = x(:); y = y(:);
np = length(x);
if (length(y) ~= np)
    error('x and y array size mismatch')
end

% NW
if nargin<3 || isempty(minw)
    minw = 0;
end

if (nargin < 4) || isempty(mode)
    mode = 5; % complex Voigt
end

% Initial guesses
[amp0, indC] = max(y);
[~,indW] = min(abs(y - max(y)/2 )); % HWHM
center0 = x(indC(1));
width0 = 2 * abs( x(indW(1))-center0 ); % FWHM

centermax = center0 + width0/2;
centermin = center0 - width0/2;

amin = 0.9 * amp0;
amax = 1.2 * amp0;

pars(1) = amp0;
lb(1) = amin;
ub(1) = amax;

pars(2) = center0;
lb(2) = centermin;
ub(2) = centermax;


if mode==5 % complex Voigt
    widthmin = 1e-6*width0;
    widthmax = 2*width0;

    % Olivero–Longbothum approximation
    Acoef = 0.5346; Bcoef = 0.2166;
    widthL = 0.4 * width0;

    term = (width0 - Acoef*widthL).^2 - Bcoef*(widthL.^2);
    widthG = sqrt(max(0, term));

    pars(3) = widthL;
    lb(3) = widthmin;
    ub(3) = widthmax;

    pars(5) = widthG;
    lb(5) = widthmin;
    ub(5) = widthmax;
else
    widthmin = 0.5*width0;
    widthmax = 2*width0;
    if (widthmin < minw)
        widthmin = minw;
    end

    pars(3) = width0;
    lb(3) = widthmin;
    ub(3) = widthmax;
end

% NW
if mode>2
    pars(4) = 0; % phase
    lb(4) = -15; % degrees
    ub(4) = 15;
end

%%
if size(x)==size(y'), y = y'; end

oldoptions = optimset('lsqcurvefit');
options = optimset(oldoptions, 'TolFun', 1e-12,'TolX', 1e-12,'MaxFunEval',20000,'MaxIter', 12000 );

if (mode == 1) % Gaussian

    ip = lsqcurvefit(@composite,pars,x,y,lb,ub,options);
    ampl = abs(ip(1));
    pos = ip(2);
    width = abs(ip(3));
    xpw = (x-pos)/width;
    yfit = ampl * exp( -( (xpw/0.6006).^2 ) ); % 0.6006 = 0.5/(sqrt(ln(2.0)))
    areas = ampl * width * 1.37362;  % Analytical form of gaussian integral- a*w*sqrt(2pi*0.3003)

elseif mode==2 % Lorentzian

    ip = lsqcurvefit(@compositel,pars,x,y,lb,ub,options);
    ampl = abs(ip(1));
    pos = ip(2);
    width = abs(ip(3));
    xpw = (x-pos)/width;
    yfit = ampl ./ ( (4.0*(xpw .*xpw)) + 1.0 );
    areas = ampl * width * 1.5708;  % Analytical form of lorentzian integral- a*w*pi/2

elseif mode==3 % complex Lorentzian - NW

    ip = lsqcurvefit(@compositel_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1));
    pos = ip(2);
    width = abs(ip(3));
    phase = ip(4);
    yfit = compositel_complex(ip,x);
    areas = ampl .* width * pi/2; % same as real Lorentzian

elseif mode==4 % complex Gaussian - NW

    ip = lsqcurvefit(@compositeG_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1));
    pos = ip(2);
    width = abs(ip(3));
    phase = ip(4);
    yfit = compositeG_complex(ip,x);
    areas = ampl .* width * 1.37362; % same as real Gaussian

elseif mode==5 % complex Voigt - NW

    ip = lsqcurvefit(@compositeV_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1));
    pos = ip(2);
    widthL = abs(ip(3));
    phase = ip(4);
    widthG = abs(ip(5));
    width = Acoef * widthL + sqrt(Bcoef*widthL.^2 + widthG.^2);     % Olivero–Longbothum approximation
    yfit = compositeV_complex(ip,x);
    temppars = [ampl pos widthL 0 widthG];
    peakfit = compositeV_complex(temppars,x);
    areas = sum(peakfit(:));

end

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

function out = col(in)
out = in(:);
return
    