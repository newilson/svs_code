%% demo_SVS_CrCEST
clear,clc

%% scan 1 - AM 08/15/2025
pars.fitmode = 1;
pars.peak = 0; % look at water peak - note that peak is in z-spectrum
pars.dummyShots = 1; % removes 1st time point pre and post exercise

mainDir = 'Y:\users\nwilson\S077\RESEARCH_RARMBRUSTER_20250815_152315_224000';

preScan = 'PREP_SVS_CRCEST_BASE_0012';
postScan = 'PREP_SVS_CRCEST_POST_0013';

preData = postprocessSVS_CrCEST([mainDir filesep preScan],pars);
postData = postprocessSVS_CrCEST([mainDir filesep postScan],pars);

CrCEST = cat(2,preData.CrCEST,postData.CrCEST);
pos = cat(2,preData.pos,postData.pos);
neg = cat(2,preData.neg,postData.neg);

figure, plot(1:length(CrCEST),CrCEST,':o');

figure, plot(1:length(pos),pos,':ob',1:length(neg),neg,':or');

% individual spectra
ppm = preData.ppm;
preSpec = fftshift(fft(preData.fid,[],1),1);
postSpec = fftshift(fft(postData.fid,[],1),1);

preSpecP = preSpec(:,1:2:end);
preSpecN = preSpec(:,2:2:end);

postSpecP = postSpec(:,1:2:end);
postSpecN = postSpec(:,2:2:end);

specP = cat(2,preSpecP,0*preSpecP(:,1),postSpecP);
specN = cat(2,preSpecN,0*preSpecN(:,1),postSpecN);

NWplayplot(abs(specP),ppm,[],abs(specN),ppm);


%% scan 2 - AS 08/28/2025
pars.fitmode = 1.5;
pars.peak = 0;
pars.dummyShots = 1;
pars.denoiseSeries = 3;
% pars.denoiseSeries = false;

mainDir = 'Y:\users\rarmbruster\data\S092_20250828\dcm';

isRandomExercise = false;
switch isRandomExercise
    case false
        preScan = 'PREP_SVS_CRCEST_BASE_MG_250V_0011';
        postScan = 'PREP_SVS_CRCEST_POST_MG_250V_0012';

    case true
        preScan = 'PREP_SVS_CRCEST_BASE_MG_250V_2_0016';
        postScan = 'PREP_SVS_CRCEST_POST_MG_250V_RANDOM_0017';

end

preData = postprocessSVS_CrCEST([mainDir filesep preScan],pars);
postData = postprocessSVS_CrCEST([mainDir filesep postScan],pars);

CrCEST = cat(2,preData.CrCEST,postData.CrCEST);
pos = cat(2,preData.pos,postData.pos);
neg = cat(2,preData.neg,postData.neg);

figure, plot(1:length(CrCEST),CrCEST,':o');

figure, plot(1:length(pos),pos,':ob',1:length(neg),neg,':or');

% individual spectra
ppm = preData.ppm;
preSpec = fftshift(fft(preData.fid,[],1),1);
postSpec = fftshift(fft(postData.fid,[],1),1);

preSpecP = preSpec(:,1:2:end);
preSpecN = preSpec(:,2:2:end);

postSpecP = postSpec(:,1:2:end);
postSpecN = postSpec(:,2:2:end);

specP = cat(2,preSpecP,0*preSpecP(:,1),postSpecP);
specN = cat(2,preSpecN,0*preSpecN(:,1),postSpecN);

% NWplayplot(abs(specP),ppm,[],abs(specN),ppm);

%% scan 3 - BB 08/29/2025
pars.fitmode = 1;
pars.peak = 0;
pars.dummyShots = 0;

mainDir = 'Y:\users\rarmbruster\data\S068_20250829\dcm';

isRandomExercise = true;
isLargeVoxel = false;
switch isRandomExercise
    case false
        switch isLargeVoxel
            case true
                preScan = 'PREP_SVS_CRCEST_PRE_MG_LARGESHIMVOL_0013';
                postScan = 'PREP_SVS_CRCEST_POST_MG_LARGEVOXSHIM_0015';
            case false
                preScan = 'PREP_SVS_CRCEST_PRE_MG_LARGESHIMVOL_0013';
                postScan = 'PREP_SVS_CRCEST_PRE_MG_SMALLSHIMVOL_0014';
        end
    case true
        preScan = 'PREP_SVS_CRCEST_PRE_MG_LARGESHIMVOL_RANDOM_0019';
        postScan = 'PREP_SVS_CRCEST_POST_MG_LARGEVOXSHIM_RANDOM_0020';
end

preData = postprocessSVS_CrCEST([mainDir filesep preScan],pars);
postData = postprocessSVS_CrCEST([mainDir filesep postScan],pars);

CrCEST = cat(2,preData.CrCEST,postData.CrCEST);
pos = cat(2,preData.pos,postData.pos);
neg = cat(2,preData.neg,postData.neg);

figure, plot(1:length(CrCEST),CrCEST,':o');

figure, plot(1:length(pos),pos,':ob',1:length(neg),neg,':or');

% individual spectra
ppm = preData.ppm;
preSpec = fftshift(fft(preData.fid,[],1),1);
postSpec = fftshift(fft(postData.fid,[],1),1);

preSpecP = preSpec(:,1:2:end);
preSpecN = preSpec(:,2:2:end);

postSpecP = postSpec(:,1:2:end);
postSpecN = postSpec(:,2:2:end);

specP = cat(2,preSpecP,0*preSpecP(:,1),postSpecP);
specN = cat(2,preSpecN,0*preSpecN(:,1),postSpecN);

NWplayplot(abs(specP),ppm,[],abs(specN),ppm);
