function filt = NWSiemensRawFilter2d(N,W,S)
% N is total length
% W is width of nonzero filter
% S is width of sloped part
% Note that center W-2S is 1s and that outside (N-W)/2 is 0s on each side

if nargin<2 || isempty(W)
    W = N(1);
end
if nargin<3 || isempty(S)
    S = round(W/4);
end

if length(N)>2 || length(W)>2 || length(S)>2
    error('vectors must have 2d only')
end
if min(W-2*S)<0
    error('<0 length flat part')
end
if min(N-W)<0
    error('filter is larger than image')
end

filtr = 0.5 + 0.5*cos(pi*(0:S(1)-1)/S(1));
filtl = 0.5 - 0.5*cos(pi*(1:S(1))/S(1));
filtc = ones(1,W(1)-2*S(1));

zf = round((N(1)-W(1))/2);
filt1 = [zeros(1,zf),filtl,filtc,filtr,zeros(1,zf)];
filt1 = filt1(:);

if isequal(1,length(N),length(W),length(S))
    filt2 = filt1';
else
    filt2r = 0.5 + 0.5*cos(pi*(0:S(end)-1)/S(end));
    filt2l = 0.5 - 0.5*cos(pi*(1:S(end))/S(end));
    filt2c = ones(1,W(end)-2*S(end));
    
    zf = round((N(end)-W(end))/2);
    filt2 = [zeros(1,zf),filt2l,filt2c,filt2r,zeros(1,zf)];
    filt2 = filt2(:)';
end

filt = filt1*filt2;

end