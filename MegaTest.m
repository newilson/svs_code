ws_dat = 'C:\Users\z0043etb\Desktop\MSLASER\TestData\meas_MID00017_FID09301_WIP_1033_svs_slaser_mega_ws.dat';
% [ws, ~, ppm] = dat2mat_svsMEGA(ws_dat,ws_dat);
[ws, ~, ppm] = dat2mat_svsMEGA(ws_dat,[]);


spe1 = 'C:\Users\z0043etb\Desktop\MSLASER\TestData\WriteToFile_00001.spe';
spe2 = 'C:\Users\z0043etb\Desktop\MSLASER\TestData\WriteToFile_00002.spe';
spe3 = 'C:\Users\z0043etb\Desktop\MSLASER\TestData\WriteToFile_00003.spe';
fid1 = read_spe(spe1);
fid2 = read_spe(spe2);
fid3 = read_spe(spe3);
fid = cat(2,fid1,fid2,fid3);
fid = conj(fid);

% spec1 = fftshift(fft(fid1));
% spec2 = fftshift(fft(fid2));
% spec3 = fftshift(fft(fid3));
% spec = cat(2,spec1,spec2,spec3);

spec = fftshift(fft(fid,[],1),1);

sc1 = max(abs(ws(:)));
sc2 = max(abs(spec(:)));

ws = ws/sc1;
spec = spec/sc2;

% remove OS
npts = length(ppm);
ppm = ppm(npts/4+1:3*npts/4);
ws = ws(npts/4+1:3*npts/4,:);

% figure, plot(ppm,flip(real(ws(:,1))),'r',ppm,real(spec(:,1)),'b')
NWplayplot(real(ws),ppm,[],real(spec),ppm);