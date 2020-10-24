function [out,filt] = expFilter(t,lb,in)
%
% out = expFilter(t,lb,out)
%
% t is [npts x 1]
% lb is line broadening in Hz

if isequal(lb,0)
    filt = 1;
    out = in;
    return
end

if nargin>2
    npts = size(in,1);
    if ~isequal(length(t),npts)
        warning('number of points do not match')
        in = [];
    end
else
    in = [];
end

out = [];
filt = exp(-2*pi*t(:)*lb);

if ~isempty(in)
    si = size(in);
    in = reshape(in,npts,[]);
    filt2d = repmat(filt,[1 size(in,2)]);
    out = in .* filt2d;
    out = reshape(out,si);
end
