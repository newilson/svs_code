function [weights, fidcc, fid_refcc] = coilCombinationNoPC(fid_ref,min_sig_frac,fid)
%
% Coil combination without phase correction
%
% fid_ref is [npts x nc]
% fid is [npts x nc x other stuff]

if nargin<2 || isempty(min_sig_frac)
    min_sig_frac = 0.15;
end

sig = fid_ref(1,:) .* conj(fid_ref(1,:));
weights = sig/sum(sig(:));
weights  = weights/max(weights(:));
weights(weights<min_sig_frac)=0; % zeros out low signal coils

fidcc = []; fid_refcc = [];
if nargout>1 && nargin>2 && ~isempty(fid)
    [npts,nc] = size(fid_ref);
    si = size(fid);
    fid = fid(:,:,:); % make 3d maximum
    if ~isequal(si(2),nc)
        return;
    end
    fid_refcc = zeros(npts,1);
    if length(si)<3
        fidcc = zeros(npts,1);
    else
        fidcc = zeros([si(1) prod(si(3:end))]);
    end
    for ii=1:nc
        fid_refcc = fid_refcc + fid_ref(:,ii)*weights(ii);
        fidcc = fidcc + squeeze(fid(:,ii,:))*weights(ii);
    end
    if ~(length(si)<3)
        fidcc = reshape(fidcc,si([1 3:end]));
    end

end
