%% demo_DFS_test
clear,clc

%% 091625 using slice selective svssel

mainDir = 'Y:\users\rnanga\U01_RR\test\S095_09162025_28Ch-KneeCoil\datfiles';

nwsAll = 'S095_09162025_meas_MID00440_FID23669_svssel_latest_M0_Water_1avg_tra.dat';
ws30_1 = 'S095_09162025_meas_MID00457_FID23686_svssel_latest_30ppm_1024avg_tra.dat';
ws30_2 = 'S095_09162025_meas_MID00460_FID23689_svssel_latest_30ppm_8kBW_2048avg_tra.dat';
ws78_1 = 'S095_09162025_meas_MID00458_FID23687_svssel_latest_78ppm_2048avg_tra.dat';
ws78_2 = 'S095_09162025_meas_MID00459_FID23688_svssel_latest_78ppm_8kBW_2048avg_tra.dat';

pars.lb = 20;
pars.eccopt = -1;

scan = 3.1;
switch scan
    case 1 % original ws30_1
        ws_dat = fullfile(mainDir,ws30_1);
        nws_dat = [];
    case 1.1 % with water scan
        ws_dat = fullfile(mainDir,ws30_1);
        nws_dat = fullfile(mainDir,nwsAll);
    case 2 % original ws30_2
        ws_dat = fullfile(mainDir,ws30_2);
        nws_dat = [];
    case 2.1 % with water scan
        ws_dat = fullfile(mainDir,ws30_2);
        nws_dat = fullfile(mainDir,nwsAll);
    case 3 % original ws78_1
        ws_dat = fullfile(mainDir,ws78_1);
        nws_dat = [];
    case 3.1
        ws_dat = fullfile(mainDir,ws78_1);
        nws_dat = fullfile(mainDir,nwsAll);
    case 4
        ws_dat = fullfile(mainDir,ws78_2);
        nws_dat = [];
    case 4.1
        ws_dat = fullfile(mainDir,ws78_2);
        nws_dat = fullfile(mainDir,nwsAll);
        
end

recon  = 1;
switch recon 
    case 0
        if isempty(nws_dat)
            [ws,nws,ppm] = dat2mat_svs(ws_dat,ws_dat,pars);
        else
            [ws,nws,ppm] = dat2mat_svs(ws_dat,nws_dat,pars);
        end
        figure, plot(ppm,real(ws),'k'), set(gca,'Xdir','reverse'), title(num2str(scan))

    case 1
        [output,pars] = dat2mat_svsNAD(pars,ws_dat,nws_dat);
        
        figure, plot(output.ppm,real(output.met.spec),'k'), set(gca,'Xdir','reverse'), title(num2str(scan))
end

if scan<3
    xlim([10 50])
else
    xlim([60 100])
end

