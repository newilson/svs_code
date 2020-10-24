function [a_arr,b_arr,lb_dt_arr,f_dt_arr,Z_arr,S_arr,fid_fit, U ,S,V] = ...
    hsvd_big(fid_arr,p,n,bound,figid,returncell)
% HSVD  Hankel singular value decomposition.
%   Decompose a time domain signal into decaying sinusoids, i.e. fit
%   Lorentzian line shapes. The singular matrix yields decay
%   constants and frequencies. A subsequent linear least squares fit
%   yields (complex) amplitudes, hence concentrations and phases.
%   The underlying svd algorithm is the svds function, as the Lanczos
%   algorithm is too imprecise for noisy data.
%
%    [alpha,beta,lb_dt,f_dt,Z,S,fid_fit] = hsvd(fid,p,n,bound,figid,returncell)
%
%    alpha: complex concentrations        
%     beta: indirect decays/frequencies   
%    lb_dt: decay     [normalised units]
%     f_dt: frequency [normalised units]
%        Z: least squares matrix
%        S: singular matrix
%  fid_fit: generated FID of fit (in time domain for truncation)
%   if dim(t1)>1 output variables all cell arrays: cells = rows of fid
%
%      fid: FID (if 2D => hsvd along t2)
%        p: number of decaying sinusoids
%        n: matrix size for SVD
%    bound: boundaries for fitting [normalised units]
%           [-0.5 0.5] (default): fit whole FID
%           [-0.02 0.02]: fit area of water => usable for h2o suppression
%    figid: figure id; if 0 => no plotting
%returncell: always return cell arrays
%
%  Literature:
%    "NMR Data Processing.", JC Hoch and AS Stern
%
% See also SVDS.

if nargin<1, help(mfilename); return; end
error(nargchk(1,6,nargin));


verb = false;

%overwrie plotting - alex
plt = false;



si = size(fid_arr);
if ~exist('p'), p = 10; end		% no. of decaying sinusoids
if ~exist('n'), n = 32768; end		% no. of data points
if ~exist('bound'), bound = [-0.5 0.5]; end;	% boundary for fitting
if ~exist('figid'), 			% plotting
    if plt
        figid = figure;
    else
        figid = false;
    end
end;
if ~exist('returncell'), returncell = false; end



m = floor(2*n/3);		% filter order
l = n-m;			% no. of equations

if l<p,
  error('p too big; n too small: increase n');
end

fid_fit = [];
if ((si(1)>1)&(~figid)), h = waitbar(0,'HSVD'); end

for ll=1:si(1),
  fid = fid_arr(ll,:);

  % data matrix
  D = ones(l,1)*[1:m] + ([1:l].')*ones(1,m)-1;
  D = fid(D);

  t1 = cputime;
  [U,S,V,flag] = svds(D,p);	% svd for only the p strongest signals
  if flag,
    if si(1)==1,
      fprintf('Warning: svds terminated with flag=%g\n',flag);
    else
      fprintf('Warning (iteration %g): svds terminated with flag=%g\n',...
              ll,flag);
    end
  end
  if verb, fprintf('Required time for svd: %g\n',cputime-t1); end
  
  ubr = U(end,:);			% bottom row of U
  Z   = (eye(p)+ubr'*ubr/(1-ubr*ubr'))*...
        U(1:end-1,:)'*U(2:end,:);  
  % least squares fit (left side is pseudo-inverse of right side)

  b  = eig(Z);			% eigenvalues yield the decays/frequencies
  lb_dt = -log(abs(b))/pi;	% (normalised) line widths
  f_dt  = angle(b)/(2*pi);	% (normalised) frequencies
  
  % create index with sensible values and within bound
  ind = find((lb_dt>0)&(f_dt>bound(1))&(f_dt<bound(2)));
  
  if verb,
    fprintf('Excluding: %g/%g\n',...
            p-length(ind),p);
  end
  lb_dt = lb_dt(ind);
  f_dt  = f_dt(ind);
  b  = b(ind);

  pp    = length(lb_dt);
  nn    = length(fid);
  bm    = (b*ones(1,nn)).^(ones(pp,1)*[0:nn-1]);
  switch 2
   case 1,
    a = fid*pinv(bm);		% crops small singular values
   case 2,
    a = bm.'\fid.';
  end

  % generate fitted FID
  if (nargout>6) | (figid),
    fid_fit_tmp = (a*ones(1,nn)).*bm;
    fid_fit_tmp = sum(fid_fit_tmp,1);
    fid_fit = [fid_fit ; fid_fit_tmp];
  end

  % write into array form
  a_arr{ll}     = a;
  b_arr{ll}     = b;
  lb_dt_arr{ll} = lb_dt;
  f_dt_arr{ll}  = f_dt;
  Z_arr{ll}     = Z;
  S_arr{ll}     = diag(S);
  
  % plotting
  if figid, 
    figure(figid); 
    plot_hsvd1d(fid,fid_fit(ll,:),f_dt);
    title(['HSVD: ' num2str(ll) '/' num2str(si(1))]);
    drawnow;
  else
    % progress waitbar
    if si(1)>1, waitbar(ll/si(1),h); drawnow; end
  end
end

% if input fid 1d => overwrite hash arrays
if si(1)==1,
  if ~returncell,
    a_arr     = a;
    b_arr     = b;
    lb_dt_arr = lb_dt;
    f_dt_arr  = f_dt;
    Z_arr     = Z;
    S_arr     = diag(S);
  end
else
  if ~figid, close(h); end
end
