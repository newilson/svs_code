function [sig,offsets] = dat2mat_svsselRFsweep(datFile,offsetFile)

addpath('mapVBVD')
addpath('utils')

% read in offsets first
fid = fopen(offsetFile);
if fid~=-1
    offsets = [];
    tline = fgetl(fid);
    while ischar(tline)
        if ~strcmp(tline(1),'#') % test for comment
            offsets = [offsets str2double(tline)];
        end
        tline = fgetl(fid);
    end
end
fclose(fid);

% read in data
obj = mapVBVD(datFile);
if iscell(obj)
    obj = obj{end};
end
obj.image.flagDoAverage = true;
bw = 1e9/obj.hdr.Config.DwellTimeSig;
f0 = 1e-6*obj.hdr.Config.Frequency;
nch = obj.image.NCha;
npts = obj.image.NCol;

dat = double(obj.image{''});

% time axis
t = (0:npts-1)*1/bw;

% freq axis
hz = (-1/2:1/npts:1/2-1/npts)*bw;
ppm = 4.72 + hz/f0;

% apodization
lb = 3;
dat = expFilter(t,lb,dat);

% water signal only
wsopts.type = 'hsvd';
wsopts.hsvd.bounds = [-0.025 0.025]; % normalized bounds for water
wsopts.hsvd.nsin = 25; % number of decaying sinusoids
ws = removeResidualWater(dat,wsopts);
wat = dat - ws;
clear dat ws

% derive coil weights from 0 ppm scan
ind = find(offsets-4.72==min(abs(offsets-4.72)));
ph = angle(wat(1,:,ind(1)));
weights = exp(-1i*ph);

% combine coils
for ii=1:length(weights)
    wat(:,ii,:) = wat(:,ii,:) * weights(ii);
end
wat = squeeze(sum(wat,2));

% signal - first time point of water-only fid
sig = wat(1,:);



