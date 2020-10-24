function [yfit,n,names,ampl,pos,width,integral,spins,pars] = curvefitIMCL(x,y,minw,mode)

%[yfit,ampl,pos,width,integral, pars, n] =
%curvefit(x,y,minw,n,mode,miny,maxy)

% nonlinear curve fitting program for n peaks
% inputs:    x = array of xpoints
%            y = array of spectral points
%            n = number of peaks expected
%            mode 1= gaussian, 2 = lorentzian

%input parameters for testing:
% x=xarr(gluindexs:gluindexe);
% y=gluspec(spec_shift,1);
% minw=(lw/cf);
% n=npeaks;
% mode=fitmode;
% nargin = 5


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
plot(x,y); axis([min(x) max(x) miny maxy]);
hold on;

% First get initial estimates of position, amplitude and widths from user
% interactively.
for k = 1:1000%n NW
    str1 = sprintf('Place cursor on top of peak %2i and click left button (right click if done)',k);
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
    if but==3, break; end
    tmpind = find(abs(x-xk) == min(abs(x-xk)));
    cg(k) = xk; %stores chosen points in an array
    ag(k) = yk; 
    plot(cg(k),ag(k),'k*');
    str1 = sprintf('Place cursor on FWHM of peak %2i and click left button',k);
    wg(k) = 0;
    while (wg(k) == 0)
        disp(str1);
        [xk,yk,but] = ginput(1);
        wg(k) = abs((xk)-cg(k)); %find width by: left point - ampl. point
        if (wg(k) < 0.5*minw)
            wg(k) = 0.5*minw; %minw input by user when function called, set to lw/cf: 10/399.486
        end
    end

    plot(cg(k)-wg(k),yk,'k*');
    plot(cg(k)+wg(k),yk,'k*');
end
n = k-1; % NW
for ii=1:n
    str = sprintf('Name of peak %2i : ',ii); 
    names{ii} = input(str,'s');
end


cgmax = cg+wg;
cgmin = cg-wg;
wg = 2*wg;

amin = 0.2 * ag; %ag is the y-point chosen by user; search 20% of this value
amax = 1.3 * ag; % NW

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

else
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