%% 31-P / 1-H birdcage coil calibration

mainDir = 'Y:\users\nwilson\RESEARCH_PAUL_20250819_154147_436000';

scan  = 2;
switch scan
    case 0
        % 1-H normal flips
        subDir = 'FID_1HFACALIB_30_140_8_200US_0003';

    case 1
        % 1-H low flips
        subDir = 'FID_1HFACALIB_0_70_16_200US_0004';

    case 2
        % 31-P normal flips
        subDir = 'FID_31P_FACALIB_30_140_16_500US_0006';

    case 3
        % 31-P low flips
        subDir = 'FID_31P_FACALIB_0_70_16_500US_0005';

    case 4
        mainDir = 'Y:\users\nwilson';
        subDir = 'FID_TEST_1H_FACALIB_0079';
end

opt = {'t','f'};
for ii=1:length(opt)
    [refVolt, sig, FA, fit, FAinterp] = fid_refVolt([mainDir filesep subDir],opt{ii});
end