% demo of postprocessMEGA.m

scan = 3;

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
        mainDir = 'Y:\users\rnanga\853707\Testrun\S000_12212023\rdafiles';
        filenames{2} = [mainDir filesep 'S000_gaba_edit-on.rda'];
        filenames{1} = [mainDir filesep 'S000_gaba_edit-off.rda'];
    case 1
        mainDir = 'Y:\users\rnanga\853707\Testrun\S108_12202023\rdafiles';
        filenames{2} = [mainDir filesep 'S108_gaba_edit-on.rda'];
        filenames{1} = [mainDir filesep 'S108_gaba_edit-off.rda'];
    case 2
        mainDir = 'Y:\users\rnanga\853707\Testrun\S108_12202023\rdafiles';
        filenames{2} = [mainDir filesep 'S108_gabamm_edit-on.rda'];
        filenames{1} = [mainDir filesep 'S108_gabamm_edit-off.rda'];
    case 3
        mainDir = 'Y:\users\rnanga\853707\Testrun\S000_12212023\E100\EJA_SVS_MSLASER_DACC_0018';
        filenames{1} = [mainDir filesep 'S000.MR.CAMRIS_BACKEDUP_USETHIS_WIERS.0018.0001.2023.12.21.12.45.31.774424.119190910.IMA'];
        filenames{2} = [mainDir filesep 'S000.MR.CAMRIS_BACKEDUP_USETHIS_WIERS.0018.0002.2023.12.21.12.45.31.774424.119190921.IMA'];
end

output = postprocessMEGA(filenames,pars);

figure, plot(output.ppm, real(output.spec(:,4))), xlim([2.5 3.5])