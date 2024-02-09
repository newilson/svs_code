function fidba = blockAverage(fid,block_size)
%
% fidba = blockAverage(fid,block_size)
%
% fid is [npts x navgs x other stuff]

if nargin<2 || isempty(block_size)
    fidba = fid;
    return;
end

si = size(fid);
siba = si;
if si(2)<block_size || mod(si(2),block_size)
    warning('block averaging not performed')
    fidba = fid;
    return;
end
siba(2) = si(2)/block_size;
fidba = zeros(siba);

for ii=1:siba(2)
    fidba(:,ii,:) = mean(fid(:,(ii-1)*block_size+1:ii*block_size,:),2);
end

