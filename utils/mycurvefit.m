function [yfit,ampl,pos,hwhm] = mycurvefit(x,yin,ampl0,pos0,hwhm0,mode)

n = length(ampl0);
if ~isequal(n,length(pos0),length(hwhm0))
    error('initial values must have the same size')
end
np = length(x);
if (length(yin) ~= np)
    error('x and y array size mismatch')
end
if ( size(x,1) == size(yin,1) )
    y = double(yin);
else
    y = double(yin)';
end

oldoptions = optimset('lsqcurvefit');
options = optimset(oldoptions, 'TolFun', 1e-12,'TolX', 1e-12,'MaxFunEval',20000,'MaxIter', 500, 'Display', 'off' );

pars0 = zeros(3*n,1);
pars0(1:n) = ampl0;
pars0(n+1:2*n) = pos0;
pars0(2*n+1:3*n) = hwhm0;

if mode==1 % gaussian
    pars = lsqcurvefit(@gaussfit,pars0,x,y,[],[],options);
    ampl = pars(1:n);
    pos = pars(n+1:2*n);
    hwhm = pars(2*n+1:3*n);
    yfit = 0*x;
    integral = 0*ampl;
    for k=1:n
        xpw = (x-pos(k))/hwhm(k);
        integral(k) = ampl(k) * hwhm(k) * sqrt(2*pi*0.3003);
        yfit = yfit + ampl(k) * exp( -( (xpw/0.6006).^2 ) ); % 0.6006 = 0.5/(sqrt(ln(2.0)))
    end
elseif mode==0 % magnitude lorentzian
    pars = lsqcurvefit(@abslorentzfit,pars0,x,y,[],[],options);
    ampl = pars(1:n);
    pos = parsn(n+1:2*n);
    hwhm = pars(2*n+1:3*n);
    yfit = 0*x;
    for k=1:n
        yfit = yfit + ampl(k) ./ sqrt(hwhm(k)^2 + (x-pos(k)).^2);
    end
else % lorentzian
    pars = lsqcurvefit(@lorentzfit,pars0,x,y,[],[],options);
    ampl = pars(1:n);
    pos = pars(n+1:2*n);
    hwhm = pars(2*n+1:3*n);
    yfit = 0*x;
    for k=1:n
        yfit = yfit + ampl(k)*hwhm(k) ./ (hwhm(k)^2 + (x-pos(k)).^2);
    end
end

end

function yfit = gaussfit(ip,x)

n = (size(ip,2))/3;
yfit = x*0 ;
for k = 1:n
    a = abs( ip(k) );
    p = ip(k+n);
    w = abs( ip(k+n+n) );
    xpw = (x-p)/w;
    yfit = yfit + a * exp( -( (xpw/0.6006).^2 ) ); % 0.6006 = 0.5/(sqrt(ln(2.0)))
end
end


function yfit = lorentzfit(ip,x)

n = (size(ip,2))/3;
yfit = x*0 ;
for k = 1:n
    a = abs( ip(k) );
    p = ip(k+n);
    w = abs( ip(k+n+n) );
    yfit = yfit + a * w ./ (w^2 + (x-p).^2);
end

end

function yfit = abslorentzfit(ip,x)

n = size(ip,2)/3;
yfit = x*0;
for k=1:n
    a = abs(ip(k));
    p = ip(k+n);
    w = abs(ip(k+2*n));
    yfit = yfit + a ./ sqrt(w^2 + (x-p).^2);
end

end

