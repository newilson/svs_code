% function [yfit,n,names,ampl,pos,width,areas,ip,pars] = curvefitMan(x,y,minw,mode,peak)
function [yfit,names,areas,ip,pars,lb,ub] = curvefitMan(x,y,minw,mode,peak)


%[yfit,ampl,pos,width,areas, pars, n] =
%curvefit(x,y,minw,n,mode,miny,maxy)

% nonlinear curve fitting program for n peaks
% inputs:    x = array of xpoints
%            y = array of spectral points
%            n = number of peaks expected
%            mode 1= gaussian, 2 = lorentzian, 3 = complex lorentzian (NW),
%            4 = complex gaussian (NW), 5 = complex Voigt (NW)

%input parameters for testing:
% x=xarr(gluindexs:gluindexe);
% y=gluspec(spec_shift,1);
% minw=(lw/cf);
% n=npeaks;
% mode=fitmode;
% nargin = 5

x = x(:); y = y(:);
np = length(x);
if (length(y) ~= np)
    error('x and y array size mismatch')
end

% NW
if nargin<3 || isempty(minw)
    minw = 0;
end

if (nargin < 4)
    mode = 1;
end

if (nargin < 6)
    miny = min(y);
    maxy = max(y);
end

% Use n peaks - work in ppm units
% cg = zeros(n,1);
% wg = cg;
% ag = cg;


cfh = figure;
h1 = plot(x,y,'b-'); axis([min(x) max(x) miny maxy]);
hold on;
yres = y;
h2 = plot(x,yres,'r-');

BW = abs(x(end)-x(1));
dwelltime = 1/BW;
% global time
time = col(dwelltime*(0:length(x)-1));

% First get initial estimates of position, amplitude and widths from user
% interactively.
if nargin<5 || isempty(peak)
    niter = 1000;
else
    niter = length(peak);
end
for k = 1:niter%n NW
    if nargin<5 || isempty(peak)
        str1 = sprintf('Place cursor on top of peak %2i and click left button (right click if done)',k);
    else
        str1 = sprintf('Place cursor on top of peak %s and click left button',peak{k});
    end
    disp(str1);
    but = 0;
    while but~=1
        [xk,yk,but] = ginput(1);
        if but==3
            break;
        elseif but==122 % "z" zoom in
            ax = axis; width = ax(2)-ax(1); height = ax(4)-ax(3);
            axis([xk-width/2 xk+width/2 yk-height/2 yk+height/2])
            zoom(3/2);
        elseif but==120 % "x" zoom out
            ax = axis; width = ax(2)-ax(1); height = ax(4)-ax(3);
            axis([xk-width/2 xk+width/2 yk-height/2 yk+height/2])
            zoom(2/3);
        end
    end
    if but==3, break, end % NW
    tmpind = find(abs(x-xk) == min(abs(x-xk)));
    center(k) = xk; %stores chosen points in an array
    amp(k) = yk; 
    plot(center(k),amp(k),'k*');
    if nargin<5 || isempty(peak)
        str1 = sprintf('Place cursor on FWHM of peak %2i and click left button',k);
    else
        str1 = sprintf('Place cursor on FWHM of peak %s and click left button',peak{k});
    end
    width(k) = 0;
    while (width(k) == 0)
        disp(str1);
        [xk,yk,but] = ginput(1);
        width(k) = abs((xk)-center(k)); % (HWHM - NW)find width by: left point - ampl. point
        if (width(k) < 0.5*minw)
            width(k) = 0.5*minw; %minw input by user when function called, set to lw/cf: 10/399.486
        end
    end

    plot(center(k)-width(k),yk,'k*');
    plot(center(k)+width(k),yk,'k*');
    xpw = (x-center(k))/width(k);
    if mode==1 || mode==4 || mode==5      
        yres = yres - amp(k) * exp( -( (xpw/0.6006).^2 ) );
    elseif mode==2 || mode==3
        yres = yres - amp(k) ./ ( (4.0*(xpw .*xpw)) + 1.0 );
    end
    h2.YData = yres; % update remaining spectrum
end
if nargin<5 || isempty(peak)
    n = k-1; % NW
    for ii=1:n
        str = sprintf('Name of peak %2i : ',ii);
        names{ii} = input(str,'s');
    end
else
    names = peak;
    n = niter;
end


centermax = center+width;
centermin = center-width;
width = 2*width; % FWHM - NW

amin = 0.2 * amp; %ag is the y-point chosen by user; search 20% of this value
amax = 1.1 * amp; % NW

pars( 1:n ) = amp;
lb(1:n) = amin;
ub(1:n) = amax;

pars( n + (1:n) ) = center;
lb( n + (1:n) ) = centermin;
ub( n + (1:n) ) = centermax;

if mode==5 % complex Voigt
    widthmin = 1e-6*width;
    widthmax = 1.4*width;

    % Olivero–Longbothum approximation
    Acoef = 0.5346; Bcoef = 0.2166;
    widthL = 0.4 * width;

    term = (width - Acoef*widthL).^2 - Bcoef*(widthL.^2);
    widthG = sqrt(max(0, term));

    pars( 2*n  + (1:n) ) = widthL;
    lb( 2*n + (1:n) ) = widthmin;
    ub( 2*n + (1:n) ) = widthmax;

    pars( 4*n  + (1:n) ) = widthG;
    lb( 4*n + (1:n) ) = widthmin;
    ub( 4*n + (1:n) ) = widthmax;
else
    widthmin = 0.5*width;
    widthmax = 1.4*width;
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
if mode>2
    pars(3*n + (1:n)) = 0; % phase
    lb(3*n+(1:n)) = -179; % degrees
    ub(3*n+(1:n)) = 179;
end


%%
if size(x)==size(y'), y = y'; end

spins = ones(n,1); % NW fix this later
oldoptions = optimset('lsqcurvefit');
options = optimset(oldoptions, 'TolFun', 1e-12,'TolX', 1e-12,'MaxFunEval',20000*n,'MaxIter', 12000 );

if (mode == 1)
    ip = lsqcurvefit(@composite,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip( (1:n) + n);
    width = abs(ip( (1:n) + n + n));
    yfit = x*0;
    for k = 1:n
        a = abs( ip(k) );
        p = ip(k+n);
        w = abs( ip(k+n+n) );
        xpw = (x-p)/w;
        yfit = yfit + a * exp( -( (xpw/0.6006).^2 ) ); % 0.6006 = 0.5/(sqrt(ln(2.0)))
        areas(k) = a * w * 1.37362;  % Analytical form of gaussian areas- a*w*sqrt(2pi*0.3003)
    end

elseif mode==2
    ip = lsqcurvefit(@compositel,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip( (1:n) + n);
    width = abs(ip( (1:n) + n + n));
    yfit = x*0;
    for k = 1:n
        a = abs( ip(k) );
        p = ip(k+n);
        w = abs( ip(k+n+n) );
        xpw = (x-p)/w;
        yfit = yfit + a ./ ( (4.0*(xpw .*xpw)) + 1.0 );
        areas(k) = a * w * 1.5708;  % Analytical form of lorentzian areas- a*w*pi/2
    end
elseif mode==3 % complex Lorentzian - NW

    ip = lsqcurvefit(@compositel_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip( (1:n) + n);
    width = abs(ip( (1:n) + n + n));
    phase = ip((1:n)+3*n);
    yfit = compositel_complex(ip,x);
    areas = ampl .* width * pi/2; % same as real Lorentzian

elseif mode==4 % complex Gaussian - NW

    ip = lsqcurvefit(@compositeG_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip((1:n)+n);
    width = abs(ip((1:n)+2*n));
    phase = ip((1:n)+3*n);
    yfit = compositeG_complex(ip,x);
    areas = ampl .* width * 1.37362; % same as real Gaussian

elseif mode==5 % complex Voigt - NW

    ip = lsqcurvefit(@compositeV_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip((1:n)+n);
    widthL = abs(ip((1:n)+2*n));
    phase = ip((1:n)+3*n);
    widthG = abs(ip((1:n)+4*n));
    width = Acoef * widthL + sqrt(Bcoef*widthL.^2 + widthG.^2);     % Olivero–Longbothum approximation
    yfit = compositeV_complex(ip,x);
    for ii=1:n % evaluate each peak with 0 phase and integrate
        temppars = [ampl(ii) pos(ii) widthL(ii) 0 widthG(ii)];
        peakfit = compositeV_complex(temppars,x);
        areas(ii) = sum(peakfit(:));
    end

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
    