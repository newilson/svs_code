function fid = NWhsvd(fid,nsv,ncad,ncol)
%
% Hankel singular value decomposition with Cadzow enhancement

if ~isvector(fid), warning('not a vector'), return; end
doaug = false;
if mod(length(fid),2)==0, fid(end+1) = 0; doaug = true; end % make fid odd
npts = length(fid);

if nargin<2 || isempty(nsv), nsv = 10; end % number of singular values
if nargin<3 || isempty(ncad), ncad = 3; end % number of Cadzow iterations
if nargin<4 || isempty(ncol), ncol = ceil(npts/2); end % square Hankel matrix

nrow = npts - ncol + 1;

% construct Hankel matrix
H = zeros(nrow,ncol);
for ii=1:nrow
    H(ii,:) = fid(ii:ii+ncol-1);
end

% Cadzow iteration
for zz=1:ncad
    % truncated SVD
    [U,S,V] = svds(H,nsv,'largest','MaxIterations',150);
    
    % get new H
    H = U*S*V';
    
    % average antidiagonals of H to keep Hankel structure    
    H2 = fliplr(H);
    for ii=-(nrow-1):ncol-1
        vec = diag(H2,ii);
        avg = mean(vec);
        H2 = H2 + diag(avg - vec,ii);
    end
    H = fliplr(H2); % flip back
end
        
% get back fid
fid(1:ncol) = H(1,:);
fid(ncol:end) = H(:,end);
if doaug, fid = fid(1:end-1); end % remove extra point 
