function [hdr,cest_fid,twix_obj] = mapraw_cestsvs(fullFileName)
%
% Generalization of readraw_cestsvs using mapVBVD
% NW

if (nargin == 0)
    [FileName,PathName] = uigetfile(fullfile('./','*.dat'),'Select meas.dat file from cestzsvs32_ec to open');
    fullFileName = [PathName FileName];
    if ~exist(fullFileName,'file')
        disp(['Error, file: ' fullFileName ' not found']);
        return
    end
    hdr.impath = PathName;
    hdr.FileName = FileName; %Name of the file
end

twix_obj = mapVBVD(fullFileName);
if length(twix_obj)>1
    for ii=1:length(twix_obj)
        if isfield(twix_obj{ii},'image')
            twix_obj = twix_obj{ii};
            twix_obj.image.flagDoAverage = 1;
            twix_obj.image.flagremoveOS = 1;
            cest_fid = twix_obj.image{''};
            break
        end
    end
else
    twix_obj.image.flagDoAverage = 1;
    twix_obj.image.flagRemoveOS = 1;
    cest_fid = twix_obj.image{''};
end
cest_fid = conj(cest_fid); % to match output from readraw_cestsvs

hdr.Seqname = twix_obj.hdr.MeasYaps.tSequenceFileName;
hdr.ProtocolName = twix_obj.hdr.MeasYaps.tProtocolName;
hdr.BaselineString = twix_obj.hdr.Dicom.SoftwareVersions;
% hdr.patname = twix_obj.hdr.Meas.tPatientsName;
hdr.MRS_CF = twix_obj.hdr.Config.Frequency;
hdr.MRS_DWus = twix_obj.hdr.Config.DwellTime/1000;
hdr.MRS_N = twix_obj.image.NCol;
hdr.VOIPositionLR = twix_obj.hdr.Config.VoI_Position_Sag;
hdr.VOIPositionAP = twix_obj.hdr.Config.VoI_Position_Cor;
hdr.VOIPositionHF = twix_obj.hdr.Config.VoI_Position_Tra;
hdr.VOINormalTra = twix_obj.hdr.Config.VoI_Normal_Tra;
hdr.VOINormalSag = twix_obj.hdr.Config.VoI_Normal_Sag;
hdr.VOINormalCor = twix_obj.hdr.Config.VoI_Normal_Cor;
hdr.cestzsvs_averages = twix_obj.image.NAve;
hdr.reps = twix_obj.image.NRep;
hdr.VOIThickness = twix_obj.hdr.Config.VoI_SliceThickness;
hdr.VOIPhaseFOV = twix_obj.hdr.Config.VoI_PeFOV;
hdr.VOIReadoutFOV = twix_obj.hdr.Config.VoI_RoFOV;
% hdr.nrcvrs = twix_obj.hdr.Meas.MaxNoOfRxChannels; % double check this
hdr.reqreps = hdr.reps;
hdr.actreps = hdr.reps;

hdr.MRS_BW  = 1e6/hdr.MRS_DWus;

if strfind(hdr.BaselineString,'VB17') % add option for VE
    software_ver = 'VB17';
elseif strfind(hdr.BaselineString,'VD13')
    software_ver = 'VD13';
elseif strfind(hdr.BaselineString,'E12')
    software_ver = 'VE12';
else
    error('unknown software version')
end

hdr.software_ver = software_ver;
switch software_ver %add option for VE
    case 'VB17'
        hdr.WIPlong = zeros(64,1); % have to initialize or else 0 values are not included
        hdr.WIPlong(1:length(twix_obj.hdr.Phoenix.sWiPMemBlock.alFree)) = cell2mat(twix_obj.hdr.Phoenix.sWiPMemBlock.alFree);
        hdr.WIPdbl = zeros(16,1);
        hdr.WIPdbl(1:length(twix_obj.hdr.Phoenix.sWiPMemBlock.adFree)) = cell2mat(twix_obj.hdr.Phoenix.sWiPMemBlock.adFree);

        hdr.acqtype = hdr.WIPlong(1);
        hdr.cestpultype = hdr.WIPlong(4);
        hdr.cestpwms = hdr.WIPlong(5);
        
        hdr.cestppmrange = hdr.WIPdbl(4);
        hdr.cestppmstep = hdr.WIPdbl(5);
        hdr.cestb1 = hdr.WIPdbl(6);
        hdr.cestpw1ms = hdr.WIPdbl(7);
        hdr.cestdc = hdr.WIPdbl(8);
    case 'VD13'
        hdr.WIPlong = zeros(64,1); % have to initialize or else 0 values are not included
        for ii=1:length(twix_obj.hdr.Phoenix.sWipMemBlock.alFree)
            if isempty(twix_obj.hdr.Phoenix.sWipMemBlock.alFree{ii})
                twix_obj.hdr.Phoenix.sWipMemBlock.alFree{ii} = 0;
            end
        end
        hdr.WIPlong(1:length(twix_obj.hdr.Phoenix.sWipMemBlock.alFree)) = cell2mat(twix_obj.hdr.Phoenix.sWipMemBlock.alFree);
        hdr.WIPdbl = zeros(16,1);
        for ii=1:length(twix_obj.hdr.Phoenix.sWipMemBlock.adFree)
            if isempty(twix_obj.hdr.Phoenix.sWipMemBlock.adFree{ii})
                twix_obj.hdr.Phoenix.sWipMemBlock.adFree{ii} = 0;
            end
        end
        hdr.WIPdbl(1:length(twix_obj.hdr.Phoenix.sWipMemBlock.adFree)) = cell2mat(twix_obj.hdr.Phoenix.sWipMemBlock.adFree);
        
        hdr.acqtype = hdr.WIPlong(1);
        hdr.cestpwms = hdr.WIPlong(2);
        hdr.cestpw1ms = hdr.WIPlong(3);
        
        hdr.cestppmrange = hdr.WIPdbl(1);
        hdr.cestppmstep = hdr.WIPdbl(2);
        hdr.cestb1 = hdr.WIPdbl(3);
        hdr.cestdc = hdr.WIPdbl(4);
    otherwise
        hdr.WIPlong = zeros(64,1); % have to initialize or else 0 values are not included
        for ii=1:length(twix_obj.hdr.Phoenix.sWipMemBlock.alFree)
            if isempty(twix_obj.hdr.Phoenix.sWipMemBlock.alFree{ii})
                twix_obj.hdr.Phoenix.sWipMemBlock.alFree{ii} = 0;
            end
        end
        hdr.WIPlong(1:length(twix_obj.hdr.Phoenix.sWipMemBlock.alFree)) = cell2mat(twix_obj.hdr.Phoenix.sWipMemBlock.alFree);
        hdr.WIPdbl = zeros(16,1);
        for ii=1:length(twix_obj.hdr.Phoenix.sWipMemBlock.adFree)
            if isempty(twix_obj.hdr.Phoenix.sWipMemBlock.adFree{ii})
                twix_obj.hdr.Phoenix.sWipMemBlock.adFree{ii} = 0;
            end
        end
        hdr.WIPdbl(1:length(twix_obj.hdr.Phoenix.sWipMemBlock.adFree)) = cell2mat(twix_obj.hdr.Phoenix.sWipMemBlock.adFree);
        
end

