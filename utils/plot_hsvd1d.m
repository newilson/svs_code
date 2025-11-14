function plot_hsvd1d(fid,fid_fit,f_dt)
% PLOT_HSVD1D  Plotting interface for hsvd
%
% See also HSVD.

if nargin<1, help(mfilename); return; end
error(nargchk(2,3,nargin));


n = length(fid);
ws = [-n/2:(n/2-1)]/n;


spec = fftshift(fft(fid,n));
spec_fit = fftshift(fft(fid_fit,n));
spec_sup = spec-spec_fit;

spec = abs(spec);
spec_fit = abs(spec_fit);
spec_sup = abs(spec_sup);


plot(ws,spec,'b',...
     ws,spec_fit,'g',...
     ws,spec_sup,'r');
xlabel('b: original/ g: fit/ r: residual');
if exist('f_dt'),
  hold on
  plot(f_dt,zeros(1,length(f_dt)),'gx');
  hold off
end
