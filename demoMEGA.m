% demo of postprocessMEGA.m

scan = 6;

pars.lb = 2; % 2 Hz line broadening
pars.den = 'none'; % no additional denoising
pars.PC = nan; % no phase correction

wsopts.type = 'none'; % options: 'none' or 'hsvd'
wsopts.hsvd.bounds = [-0.025 0.025]; % normalized bounds for water
wsopts.hsvd.nsin = 25; % number of decaying sinusoids
wsopts.plt = false;
pars.wsopts = wsopts;

switch scan 
    case 0
        mainDir = 'Y:\users\rnanga\853707_CEW\Testrun\S000_12212023\rdafiles';
        filenames{1} = [mainDir filesep 'S000_gaba_edit-off.rda'];
        filenames{2} = [mainDir filesep 'S000_gaba_edit-on.rda'];
    case 1
        mainDir = 'Y:\users\rnanga\853707_CEW\Testrun\S00000_01112024\rdafiles';
        filenames{1} = [mainDir filesep 'S00000_01112024_GABA_editoff.rda'];
        filenames{2} = [mainDir filesep 'S00000_01112024_GABA_editon.rda'];
    case 2
        mainDir = 'Y:\users\rnanga\853707_CEW\Testrun\S108_12202023\rdafiles';
        filenames{1} = [mainDir filesep 'S108_12202023_gaba_edit-off.rda'];
        filenames{2} = [mainDir filesep 'S108_12202023_gaba_edit-on.rda'];
        filenames{3} = [mainDir filesep 'S108_12202023_gabamm_edit-off.rda'];
        filenames{4} = [mainDir filesep 'S108_12202023_gabamm_edit-on.rda'];
    case 3
        mainDir = 'Y:\users\rnanga\853707_CEW\Testrun\S000_12212023\E100\EJA_SVS_MSLASER_DACC_0018';
        filenames{1} = [mainDir filesep 'S000.MR.CAMRIS_BACKEDUP_USETHIS_WIERS.0018.0001.2023.12.21.12.45.31.774424.119190910.IMA'];
        filenames{2} = [mainDir filesep 'S000.MR.CAMRIS_BACKEDUP_USETHIS_WIERS.0018.0002.2023.12.21.12.45.31.774424.119190921.IMA'];
    case 4
        mainDir = 'Y:\users\rnanga\853707_CEW\Testrun\S092_02222024_GABAonly\rdafiles';
        filenames{1} = [mainDir filesep 'S092_02222024_dacc_gaba-editoff.rda'];
        filenames{2} = [mainDir filesep 'S092_02222024_dacc_gaba-editon.rda'];
        filenames{3} = [mainDir filesep 'S092_02222024_dacc_symMM_gaba-editoff.rda'];
        filenames{4} = [mainDir filesep 'S092_02222024_dacc_symMM_gaba-editon.rda'];
    case 5
        mainDir = 'Y:\users\rnanga\853707_CEW\Testrun\S106_02212024_GABAonly\rdafiles';
        filenames{1} = [mainDir filesep 'S106_02212024_dacc_gaba-editoff.rda'];
        filenames{2} = [mainDir filesep 'S106_02212024_dacc_gaba-editon.rda'];
        filenames{3} = [mainDir filesep 'S106_02212024_dacc_symMM_gaba-editoff.rda'];
        filenames{4} = [mainDir filesep 'S106_02212024_dacc_symMM_gaba-editon.rda'];
    case 6 % 1st subject scan1
        mainDir = '/Volumes/camipm/users/rnanga/853707_CEW/853707_S01_03212024/rdafiles';
        filenames{1} = [mainDir filesep '853707_S01_03212024_gabaMM_editoff.rda'];
        filenames{2} = [mainDir filesep '853707_S01_03212024_gabaMM_editon.rda'];
    case 7 % 1st subject scan2
        mainDir = '/Volumes/camipm/users/rnanga/853707_CEW/853707_S01_04112024/rdafiles';
        filenames{1} = [mainDir filesep '853707_S01_04112024_gabaMM_editoff.rda'];
        filenames{2} = [mainDir filesep '853707_S01_04112024_gabaMM_editon.rda'];
    case 8 % 2nd subject scan1
        mainDir = '/Volumes/camipm/users/rnanga/853707_CEW/853707_S02_04082024/rdafiles';
        filenames{1} = [mainDir filesep '853707_S02_04082024_gabaMM_editoff.rda'];
        filenames{2} = [mainDir filesep '853707_S02_04082024_gabaMM_editon.rda'];
    case 9 % 3rd subject scan1
        mainDir = '/Volumes/camipm/users/rnanga/853707_CEW/853707_S03_04152024/rdafiles';
        filenames{1} = [mainDir filesep '853707_S03_04152024_gabaMM_editoff.rda'];
        filenames{2} = [mainDir filesep '853707_S03_04152024_gabaMM_editon.rda'];
end

output = postprocessMEGA(filenames(1:2),pars);
if length(filenames)>2
    output2 = postprocessMEGA(filenames(3:4),pars);
    figure, plot(output2.ppm, real(output2.spec(:,4))), xlim([2.5 3.5]), set(gca,'XDir','reverse'), title('symMM')
end

figure, plot(output.ppm, real(output.spec(:,4))), xlim([2.5 3.5]), set(gca,'XDir','reverse')
% figure, plot(output.ppm, real(output.spec(:,1))), set(gca,'Xdir','reverse')
