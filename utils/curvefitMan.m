function [yfit,n,names,ampl,pos,width,integral,ip,pars] = curvefitMan(x,y,minw,mode,peak)

%[yfit,ampl,pos,width,integral, pars, n] =
%curvefit(x,y,minw,n,mode,miny,maxy)

% nonlinear curve fitting program for n peaks
% inputs:    x = array of xpoints
%            y = array of spectral points
%            n = number of peaks expected
%            mode 1= gaussian, 2 = lorentzian, 3 = complex lorentzian (NW),
%            4 = complex gaussian (NW)

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
global time
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
    cg(k) = xk; %stores chosen points in an array
    ag(k) = yk; 
    plot(cg(k),ag(k),'k*');
    if nargin<5 || isempty(peak)
        str1 = sprintf('Place cursor on FWHM of peak %2i and click left button',k);
    else
        str1 = sprintf('Place cursor on FWHM of peak %s and click left button',peak{k});
    end
    wg(k) = 0;
    while (wg(k) == 0)
        disp(str1);
        [xk,yk,but] = ginput(1);
        wg(k) = abs((xk)-cg(k)); % (HWHM - NW)find width by: left point - ampl. point
        if (wg(k) < 0.5*minw)
            wg(k) = 0.5*minw; %minw input by user when function called, set to lw/cf: 10/399.486
        end
    end

    plot(cg(k)-wg(k),yk,'k*');
    plot(cg(k)+wg(k),yk,'k*');
    xpw = (x-cg(k))/wg(k);
    if mode==1        
        yres = yres - ag(k) * exp( -( (xpw/0.6006).^2 ) );
    elseif mode==2 || mode==3
        yres = yres - ag(k) ./ ( (4.0*(xpw .*xpw)) + 1.0 );
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


cgmax = cg+wg;
cgmin = cg-wg;
wg = 2*wg; % FWHM - NW

amin = 0.2 * ag; %ag is the y-point chosen by user; search 20% of this value
amax = 1.1 * ag; % NW

pars( 1:n ) = ag;
lb(1:n) = amin;
ub(1:n) = amax;

pars( n + (1:n) ) = cg;
lb( n + (1:n) ) = cgmin;
ub( n + (1:n) ) = cgmax;


wgmin = 0.5*wg;
wgmax = 2*wg;
for k = 1:n
    if (wgmin(k) < minw)
        wgmin(k) = minw;
    end
end
pars( n + n  + (1:n) ) = wg;
lb( n + n + (1:n) ) = wgmin;
ub( n + n + (1:n) ) = wgmax;

% NW
pars(3*n + (1:n)) = 0; % phase
lb(3*n+(1:n)) = -179; % degrees
ub(3*n+(1:n)) = 179;

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
        integral(k) = a * w * 1.37362;  % Analytical form of gaussian integral- a*w*sqrt(2pi*0.3003)
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
        integral(k) = a * w * 1.5708;  % Analytical form of lorentzian integral- a*w*pi/2
    end
elseif mode==3 % NW
    ip = lsqcurvefit(@compositel_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip( (1:n) + n);
    width = abs(ip( (1:n) + n + n));
    phase = ip((1:n)+3*n);
    yfit = compositel_complex(ip,x);
    integral = 0;
%     yfit = x*0;
%     for k = 1:n
%         a = abs( ip(k) );
%         p = ip(k+n);
%         w = abs( ip(k+n+n) );
%         ph = ip(k+3*n);
%         xpw = (x-p)/w;
%         L = a ./ (1 + 1i*2.0*xpw) * exp(-1i*pi/180*ph);
%         yfit = yfit + real(L);
%         integral(k) = a * w * 1.5708;  % Analytical form of lorentzian integral- a*w*pi/2
%     end
elseif mode==4 % NW
    ip = lsqcurvefit(@compositeG_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip((1:n)+n);
    width = abs(ip((1:n)+2*n));
    phase = ip((1:n)+3*n);
    yfit = compositeG_complex(ip,x);
    integral = 0;
%     yfit = x*0;
%     for k = 1:n
%         a = abs(ip(k));
%         p = ip(k+n);
%         w = abs(ip(k+2*n));
%         ph = ip(k+3*n);
%         alpha = 7*0.6006*w/2; % empiric - NW
%         fid = a * exp(1i*(2*pi*p*time + ph*pi/180) - (alpha*time).^2);
%         fid(1) = fid(1)/2;
%         spec = fftshift(fft(fid,[],1),1)/(2/w * 1.33 * 1000/pi); % empiric scaling
%         yfit = yfit + real(spec);
%         integral(k) = a*w; % fix later
%     end
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

function yfit = compositeG_complex(ip,x) % NW
global time
n = size(ip,2)/4;
yfit = x*0;
for k=1:n
    a = abs( ip(k));
    p = ip(k+n);
    w = abs(ip(k+2*n));
    ph = ip(k+3*n);
    alpha = 7*0.6006*w/2; % empiric - NW
    fid = a * exp(1i*(2*pi*p*time + ph*pi/180) - (alpha*time).^2);
    fid(1) = fid(1)/2;
    spec = fftshift(fft(fid,[],1),1)/(2/w * 1.33 * 1000/pi); % empiric scaling
    yfit = yfit + real(spec);
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

function out = col(in)
out = in(:);
return
    