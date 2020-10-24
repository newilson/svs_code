scan = 1;
switch scan
    case 1
        mainDir = 'C:\Users\CMROI\Documents\MATLAB\test_data\Spectro4Ravi\08292018_DT\datfiles';
        nws_dat = 'meas_MID03544_FID25202_svsglu_wref_172V_8avg_dacc.dat';
        ws_dat = 'meas_MID03545_FID25203_svsglu_ws_172V_64avg_dacc.dat';
    case 2
        mainDir = 'C:\Users\CMROI\Documents\MATLAB\test_data\Spectro4Ravi\08302018_DT\datfiles';
        nws_dat = 'meas_MID03711_FID25369_svsglu_wref_260V_8avg_dacc.dat';
        ws_dat = 'meas_MID03712_FID25370_svsglu_ws_180V_512avg.dat';
end
ws_path = fullfile(mainDir,ws_dat);
nws_path = fullfile(mainDir,nws_dat);

[ws,nws,ppm] = dat2mat_svs(ws_path,nws_path);

figure(scan), plot(ppm,real(ws),'k'), xlim([1 7]), set(gca,'Xdir','reverse')
