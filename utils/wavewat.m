function mat_out = wavewat(mat_in,ZF,M,QMF_type,QMF_par,MI,nndp)
%
%WAVEWAT - Suppression of on-resonance water signals as described by
%          Guenther and Ludwig, JMR 2001.
%
%          call:   mat_out = wavewat(mat,ZF,M,QMF_type,QMF_par,MI,nndp)
%                  
%
%
%
  
% ULG - 6.3.2001
% Copyright (c) U.L. Günther, C. Ludwig 2001

global NMRDAT
global NMRPAR

if ~exist('ZF'), ZF = 2; end,
if ~exist('M'),   M = 7; end,
if ~exist('QMF_type'),  QMF_type = 'Daubechies'; end,
if ~exist('QMF_par'),   QMF_par  =  10; end,
if ~exist('MI'),        MI=1; end,
if ~exist('nndp'), nndp=0; end,

L    = 1;    % ?????????????????????????????????????????
QMF = MakeONFilter(QMF_type,QMF_par);

if 0
  [mat_in, scaleFactor] = NormNoise1(mat_in',QMF);
  mat_in=mat_in';
end

[d1size d2size] = size(mat_in);
TD = d1size;

% Make it a column if it is a 1D vector
if ((d2size==1) | (d1size==1))
  mat_in = mat_in(:);
  [d1size d2size] = size(mat_in);
  TD = d1size;
  WB_ON = logical(0);
else
  WB_ON=NMRPAR.WB_ON;
end

if exist('nndp') & nndp~=0
  nndp=abs(round(nndp));
  filt_mat = mat_in(1:nndp,:);
  mat      = mat_in(nndp+1:d1size,:); 
else
  nndp=0;
  mat = mat_in;
end  

% WAVELAB wavelet routines MUST be powers of two!!! 
[msiz1 msiz2] = size(mat);
origlength  = msiz1;
pow2length  = 2^nextpow2(length(mat));
if origlength<pow2length
  mat(pow2length,1)=0;
end

if WB_ON, 
  h=waitbar(0,['WAVEWAT 1 ... ' num2str(d2size) '.']); 
end,

for nn=1:d2size
  if WB_ON, waitbar(nn/d2size); end
  X = conj(permute(mat(:,nn),[2 1]));
  % X = (mat(:,nn))';
  orignpts = length(X);

  if exist('ZF') & ZF>0
	last_point     = X(length(X));
	X              = X - last_point;   % remove DC offset
	X(length(X)*ZF) = 0;
  end
  if MI
    X = [X(length(X):-1:1), X];          % add fid mirror image
  end

  wc    = FWT_PO(X,L,QMF);

  wavecoef = wc(:)';
  [n,J] = dyadlength(wavecoef);
  w = zeros(1,n);
  w(2^M+1:n) = wavecoef(2^M+1:n);
  X_new = IWT_PO(w,M,QMF);

  npts = length(X_new);
  if MI
    fidnew = X_new(npts/2+1:npts/2+orignpts);
  else
    fidnew = X_new(1:orignpts);
  end
  
  nmat(:,nn) = fidnew(:);
end

if exist('h') & ishandle(h)
  close(h);
  
  nmat(:,nn) = fidnew(:);
end

if exist('h') & ishandle(h)
  close(h);
end

if nndp==0
  mat_out = nmat;
else
  mat_out = filt_mat;
  if origlength<pow2length
	nmat = nmat(1:origlength,:);
  end
  mat_out(nndp+1:TD,:) = nmat;
end  

return

