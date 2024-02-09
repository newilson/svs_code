addpath('utils')

mainDir = 'C:\Users\CMROI\Documents\MATLAB\neil\svs_code\NAD T1 meas\calf\S000_AS_21_09_10-15_59_21-DST-1_3_12_2_1107_5_2_0_79106\RESEARCH_RAVIN_20210910_161623_932000';
fold{1} = 'SVSSEL_NAD_TR1300MS_SR5MS_0020';
fold{2} = 'SVSSEL_NAD_TR1350MS_SR50MS_0021';
fold{3} = 'SVSSEL_NAD_TR1400MS_SR100MS_0022';
fold{4} = 'SVSSEL_NAD_TR1550MS_SR250MS_0023';
fold{5} = 'SVSSEL_NAD_TR1800MS_SR500MS_0024';
fold{6} = 'SVSSEL_NAD_TR2300MS_SR1000MS_0025';

recovDur = [5 50 100 250 500 1000];

for ii=1:length(fold)
    file = dir([mainDir filesep fold{ii} filesep '*.IMA']);
    if length(file)>1
        error('too many files')
    end
    fname = fullfile(file.folder,file.name);
    [fid,fullhdr,h] = readSiemensDicomSpectrum(fname);
        
    if ii==1
        np = size(fid,1);
        fids = complex(zeros(np,length(fold)), zeros(np,length(fold)));
        bw = h.bw;
        t = ( 0:np-1 ) / bw;
        hz = 0:np-1;
        hz = (hz'/np-0.5)*bw;
        ppm = -hz/h.f0 + 4.72;
    end
    
    fids(:,ii) = fid;
    
end

% apodization
opts.lb = 2;
fids = expFilter(t,opts.lb,fids);

% FT
spec = fftshift(fft(fids,[],1),1);

% display
NWplayplot(real(spec),ppm);

% % phase correction
% f = NWman_phase(spec,ppm,'Manual phase correction: save and close when done');
% waitfor(f);
% if exist('out','var')
%     spec = out; clear out
%     opts.flag.pc = true;
%     opts.PC = parsPC;
% end

% semi manual curvefit last point
fitmode = 3;
peaks = {'car','atp1','atp2','nad91','nad93'};
ind = (ppm<10 & ppm >7.8);
[yfit,n,names,ampl,pos,width,integral,ip,~,lb,ub] = curvefitMan(ppm(ind),real(spec(ind,end)),0,fitmode,peaks);
figure, plot(ppm(ind),real(spec(ind,end)),'k--',ppm(ind),yfit,'r-')
A(:,size(spec,2)) = ampl(:);
P(:,size(spec,2)) = pos(:);
fit(:,size(spec,2)) = yfit(:);

for ii=1:size(spec,2)-1
    [yfit,ampl,pos,width,integral,ip] = curvefitAutoFromMan(ppm(ind),real(spec(ind,ii)),fitmode,n,ip,ub,lb);
    A(:,ii) = ampl(:);
    P(:,ii) = pos(:);
    fit(:,ii) = yfit(:);
end
