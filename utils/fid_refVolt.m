function [refVolt,sig,FA,fit,FAinterp] = fid_refVolt(folderName,opt)
%
% Calibrates reference voltage based on FA calibration option of fid_test1D
% pulse sequence. Signal is fit to a sine curve.
% 08/2025 - NW

DEBUG = false;

if nargin<1 || isempty(folderName)
    folderName = uigetdir('Choose directory with FA calibration data.');
    disp(folderName)
end

if nargin<2 || (~strcmp(opt,'t') && ~strcmp(opt,'f')) % time domain or frequency domain
    opt = 'f'; % default
end

files = dir(folderName);
nFiles = length(files);

% read data
count = 0;
for ii=1:length(files)
    fullFile = fullfile(files(ii).folder,files(ii).name);
    if ~isfolder(fullFile) && isdicom(fullFile)
        count = count+1;
        [curFid, ~, shortHdr,txtname] = readSiemensDicomSpectrum(fullFile);
        if opt=='f'
            curSpec = fftshift(fft(curFid,[],1),1);
            sig(count) = max(abs(curSpec)); % signal is maximum of absolute spectrum
        else
            sig(count) = abs(curFid(4)); % take 4th time point of fid to avoid filter ringing
        end
        if DEBUG
            figure(count), subplot(1,2,1), plot(abs(fftshift(fft(curFid,[],1),1)));
            subplot(1,2,2), plot(abs(curFid))
        end

    end
end

% read in pulse voltage for excitation pulse specific fid_test1d
ftxt = fopen(txtname,'r');
str1 = fgetl(ftxt);
RFvolt = zeros(1,count);
while ( ~isempty(str1) && ~contains(str1,'ASCCONV END') )
    str1 = fgetl(ftxt);

    tnamestr = sprintf('sTXSPEC.aRFPULSE');
    if contains(str1,tnamestr)
        for aa=1:count % arbitrary high number > 16 which is the max in the sequence
            tnamestr = sprintf('sTXSPEC.aRFPULSE[%i].flAmplitude',aa-1);
            if contains(str1,tnamestr) && ~contains(str1,[tnamestr 'NL'])
                tindex1 = strfind(str1,'=')+1;
                tname = str1(tindex1:end);
                value = sscanf(tname,'%f',1);
                RFvolt(aa) = value;
            end
        end
    end
end

sig = sig/max(sig); % scale to 1

% relavant header values
minFA = shortHdr.WIPlong(4);
maxFA = shortHdr.flipAngle;
stepFA = (maxFA - minFA) / (count-1);
TxVolt = shortHdr.refVolt;
FA = maxFA:-1*stepFA:minFA;

% check for pulse clipping
[p,S] = polyfit(FA,RFvolt,1);
if DEBUG || S.rsquared < 1 - 1e-3
    lnfit = polyval(p,0:maxFA);
    figure, plot(0:maxFA,lnfit,'-b',FA,RFvolt,'o'), xlabel('Flip Angle'), ylabel('Voltage (V)')
    text1 = sprintf('Rsq = %3.3f',S.rsquared);
    legend(text1,'Location','NorthWest')
end

% do fitting
fitFunc = @(par,FAdeg) par(1) * sind(par(2) * FAdeg);
par0 = [1 1];
par = lsqcurvefit(fitFunc,par0,FA,sig);

true90 = 90/par(2);

FAinterp = 0:1:max([FA true90]);
fit = fitFunc(par,FAinterp);

refVolt = TxVolt / par(2); 

% plot
figure,
plot(FA,sig,'bo',FAinterp,fit,'k--')
hold on
plot([true90 true90],[0,1],'-r')
xlabel('Flip Angle'), ylabel('Signal')

text0 = sprintf('Prescribed Ref Voltage = %6.1f V', TxVolt);
text1 = sprintf('Signal max at FA = %6.1f deg',true90);
text2 = sprintf('Calibrated Ref Voltage = %6.1f V',refVolt);

text(0,0.8,{text0,text1});
title(text2,'fontweight','bold','fontsize',16);

