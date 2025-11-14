function [yfit,n,names,ampl,pos,width,integral,ip,pars] = curvefitAuto(x,y,minw,mode,peak,ag0,cg0,wg0, ph_range)

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
    mode = 3; % complex Lorentzian
end

yres = y;

BW = abs(x(end)-x(1));
dwelltime = 1/BW;
global time
time = col(dwelltime*(0:length(x)-1));

names = [];

% Initial guesses
if nargin<8 || isempty(ag0) || isempty(cg0) || isempty(wg0)
    [~,cg,wg,ag] = findpeaks(y,x);
    figure, plot(x,y,'-k',cg,ag,'*r')
    n = length(cg);
else
    ag = ag0;
    cg = cg0;
    wg = wg0;
end

cgmax = cg+wg/2;
cgmin = cg-wg/2;

amin = 0.2 * ag; %ag is the y-point chosen by user; search 20% of this value
amax = 1.2 * ag; % NW

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
if nargin<9 || isempty(ph_range) || length(ph_range)~=2
    pars(3*n + (1:n)) = 0; % phase
    lb(3*n+(1:n)) = -179; % degrees
    ub(3*n+(1:n)) = 179;
else
    pars(3*n + (1:n)) = mean(ph_range);
    lb(3*n+(1:n)) = min(ph_range);
    ub(3*n+(1:n)) = max(ph_range);
end

%%
if size(x)==size(y'), y = y'; end

spins = ones(n,1); % NW fix this later
oldoptions = optimset('lsqcurvefit');
options = optimset(oldoptions, 'TolFun', 1e-12,'TolX', 1e-12,'MaxFunEval',20000*n,'MaxIter', 12000 );

if (mode == 1) % Gaussian
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

elseif mode==2 % Lorentzian
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
elseif mode==3 % complex Lorentzian - NW
    ip = lsqcurvefit(@compositel_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip( (1:n) + n);
    width = abs(ip( (1:n) + n + n));
    phase = ip((1:n)+3*n);
    yfit = compositel_complex(ip,x);
    integral = 0;
elseif mode==4 % complex Gaussian - NW
    ip = lsqcurvefit(@compositeG_complex,pars,x,y,lb,ub,options);
    ampl = abs(ip(1:n));
    pos = ip((1:n)+n);
    width = abs(ip((1:n)+2*n));
    phase = ip((1:n)+3*n);
    yfit = compositeG_complex(ip,x);
    integral = 0;
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
    