obj = mapVBVD;
if iscell(obj)
    data = obj{end}.image('');
else
    data = obj.image('');
end
size(data)

data = squeeze(data);
data = permute(data,[1 4 3 2]);
data = fftshift(fft(ifftshift(data,1),[],1),1);
data = fftshift(fft(ifftshift(data,2),[],2),2);
data = fftshift(fft(data,[],3),3);
data = squeeze(sum(abs(data),4));
bw = 1000; % Hz
hzperppm = 123.23; % 1H at 3T
hz = linspace(-bw/2,bw/2,size(data,3));
ppm = hz / hzperppm;
% NWinteractiveCSI(squeeze(spec(:,:,:,1)),ppm);



% unsorted - with reps
obj = mapVBVD;
if iscell(obj)
    data = obj{end}.image.unsorted();
else
    data = obj.image.unsorted();
end
size(data)

nt = obj.image.NLin;
ny = obj.image.NSeg;
nx = obj.image.NCol;
nch = obj.image.NCha;
acqs = obj.image.NAcq;

reps = acqs / (nt*ny);

Yinds = obj.image.Seg(1:nt:nt*ny);
[~,indsNew] = sort(Yinds);

data = reshape(data,nx,nch,nt,ny,reps);
data = permute(data,[1 4 3 2 5]);
data = data(:,indsNew,:,:,:,:);
data = fftshift(fft(ifftshift(data,1),[],1),1);
data = fftshift(fft(ifftshift(data,2),[],2),2);
data = fftshift(fft(data,[],3),3);
data = squeeze(sum(abs(data),4));
bw = 1000; % Hz
hzperppm = 123.23; % 1H at 3T
hz = linspace(-bw/2,bw/2,size(data,3));
ppm = hz / hzperppm;
NWinteractiveCSI3d(permute(data,[1 2 4 3]),ppm);