function [fid_sup,fid_fit_arr,f_dt, U,S,V] = h2osup_bigPar(fid,bound,p)
% H2OSUP  Water suppression by hsvd
%
%

plt = false;

if ~exist('bound'), bound = [-0.02 0.02]; end
if ~exist('p'), p = 25; end

if plt, figid = figure; end
n = size(fid);
if bound(1)>bound(2),
  fprintf('bound(1)>bound(2): reordering\n');
  bound = bound(end:-1:1);
end
fid_sup = [];
fid_fit_arr = [];

% NW
doParFor = false;
if ~isempty(gcp('nocreate')) && n(1)>10
    doParFor = true;
end

if ~doParFor
    if n(1)>1, h = waitbar(0,'Water suppression'); end

    for ii=1:n(1)
        n_max = 32768;
        n_svd = n(2)*(n(2)<n_max+1)+n_max*(n(2)>n_max);
        %a_arr,b_arr,lb_dt_arr,f_dt_arr,Z_arr,S_arr,fid_fit, U ,S,V
        [alpha,bet,lb_dt,f_dt,Z,S_Arr,fid_fit,U,S,V] = hsvd_big(fid(ii,:),p,n_svd,bound);
        ind = find(f_dt>bound(1)&f_dt<bound(2));
        pp  = length(ind);

        % generate fitted spectrum
        bms = (bet(ind)*ones(1,n(2))).^(ones(pp,1)*[0:n(2)-1]);
        fid_fit = (alpha(ind)*ones(1,n(2))).*bms;
        fid_fit = sum(fid_fit,1);
        fid_fit_arr = [fid_fit_arr ; fid_fit];

        % water suppressed FID
        fid_sup = [fid_sup ; fid(ii,:)-fid_fit];

        if (exist('figid')==1)
            figure(figid);
            plot_hsvd1d(fid(ii,:),fid_fit,f_dt);
            drawnow
        end

        if n(1)>1, waitbar(ii/n(1),h); end
    end
    if n(1)>1, close(h); end

else
    % NW
    n2 = n(2);
    fid_sup = zeros(size(fid));
    parfor ii=1:n(1)
        n_max = 32768;
        n_svd = n2*(n2<n_max+1)+n_max*(n2>n_max);
        %a_arr,b_arr,lb_dt_arr,f_dt_arr,Z_arr,S_arr,fid_fit, U ,S,V
        [alpha,bet,lb_dt,f_dt,Z,S_Arr,fid_fit,U,S,V] = hsvd_big(fid(ii,:),p,n_svd,bound);
        ind = find(f_dt>bound(1)&f_dt<bound(2));
        pp  = length(ind);

        % generate fitted spectrum
        bms = (bet(ind)*ones(1,n(2))).^(ones(pp,1)*[0:n(2)-1]);
        fid_fit = (alpha(ind)*ones(1,n(2))).*bms;
        fid_fit = sum(fid_fit,1);
        fid_fit_arr = [fid_fit_arr ; fid_fit];

        % water suppressed FID
        fid_sup(ii,:) = fid(ii,:) - fid_fit;

    end

end

