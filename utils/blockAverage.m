function fidba = blockAverage(fid,block_size)
%
% fidba = blockAverage(fid,block_size)
%
% fid is [npts x navgs x other stuff]

if nargin<2
    fidba = fid;
    return;
end

si = size(fid);
siba = si;
siba(2) = si(2)/block_size;
fidba = zeros(siba);

for ii=1:siba(2)
    fidba(:,ii,:) = mean(fid(:,(ii-1)*block_size+1:ii*block_size,:),2);
end

