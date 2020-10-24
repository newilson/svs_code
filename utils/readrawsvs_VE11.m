function [hdr,outfid] = readrawsvs_VE11(fullFileName)
%
%
% Hari Hariharan, 2016_04_09
% Reads measure.dat files VD13D and returns header and raw data
% 

if (nargin == 0)
    [FileName,PathName] = uigetfile(fullfile('.','*.dat'),'Select meas.dat file to open');
    fullFileName = fullfile(PathName,FileName);
    if ~exist(fullFileName,'file')
        disp(['Error, file: ' fullFileName ' not found']);
        return
    end
    hdr.impath = PathName;
    hdr.FileName = FileName; %Name of the file
else
    a = strfind(fullFileName,filesep);
    hdr.impath = fullFileName(1:(a(length(a)-1)));
    hdr.FileName = fullFileName((a(length(a))+1):length(fullFileName));
end

a = mapVBVD(fullFileName);
b = squeeze(a.image());


hdr.ProtocolName = a.hdr.Phoenix.ProtocolName;
hdr.SeqName      = a.hdr.Phoenix.tSequenceFileName;
hdr.swversion    = a.hdr.Phoenix.sProtConsistencyInfo.tMeasuredBaselineString;
hdr.Nucleus      = a.hdr.Phoenix.sTXSPEC.asNucleusInfo{1}.tNucleus;
hdr.contrasts    = a.hdr.Phoenix.lContrasts;
hdr.avg          = a.hdr.Phoenix.lAverages;
hdr.echoes       = hdr.contrasts;
for i = 1:hdr.echoes
    hdr.TE(i)    = a.hdr.Phoenix.alTE{i};
end
hdr.TR           = a.hdr.Phoenix.alTR{1};
hdr.WIPlong      = cell2mat(a.hdr.Phoenix.sWipMemBlock.alFree);
hdr.WIPdbl       = cell2mat(a.hdr.Phoenix.sWipMemBlock.adFree);
hdr.sf           = a.hdr.Phoenix.sTXSPEC.asNucleusInfo{1}.lFrequency * 1.0e-6;
hdr.dwus         = a.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1.0e-3;
hdr.sw           = 1.0e6/hdr.dwus;
hdr.MRS_N        = a.hdr.Phoenix.sSpecPara.lVectorSize * 2;
hdr.VOIthick     = a.hdr.Phoenix.sSpecPara.sVoI.dThickness;
hdr.VOIread      = a.hdr.Phoenix.sSpecPara.sVoI.dReadoutFOV;
hdr.VOIphase     = a.hdr.Phoenix.sSpecPara.sVoI.dPhaseFOV;

tname = hdr.Nucleus;
if (~isempty(strfind(tname,'1H')))
    hdr.gamma = 42.5756;
elseif (~isempty(strfind(tname,'23N')))
    hdr.gamma = 11.2620;
elseif (~isempty(strfind(tname,'17O')))
    hdr.gamma = 5.7716;
elseif (~isempty(strfind(tname,'13C')))
    hdr.gamma = 10.7063;
else
    hdr.gamma = 42.5756;
end

hdr.MRS_BW  = 1e6/hdr.dwus;

[npts,nchan,navg,nreps] = size(b);
hdr.nrcvrs = nchan;
hdr.nreps = nreps;

if ( (npts ~= hdr.MRS_N) )
    fprintf(' Inconsistent data file \n');
    outfid = [];
    return;
end
outfid = double(conj(b) * 32767.0);


hdrfilename = regexprep(fullFileName,'\.dat','\.txt');
ftxt = fopen(hdrfilename,'w');
fprintf(ftxt,'hdr.FileName = %s\n', hdr.FileName);
fprintf(ftxt,'hdr.ProtocolName = %s\n', hdr.ProtocolName);
fprintf(ftxt,'hdr.SeqName = %s\n', hdr.SeqName);
fprintf(ftxt,'hdr.SWVesion = %s\n', hdr.swversion);
fprintf(ftxt,'hdr.Nucleus = %s\n', hdr.Nucleus);
fprintf(ftxt,'hdr.echoes = %i\n', hdr.echoes);
fprintf(ftxt,'hdr.TEms = '); fprintf(ftxt,'%.1f ', hdr.TE * 0.001); fprintf(ftxt, '\n');
fprintf(ftxt,'hdr.TRms = %.1f\n', hdr.TR * 0.001);
fprintf(ftxt,'hdr.WIPlong = '); fprintf(ftxt,'%i ', hdr.WIPlong); fprintf(ftxt, '\n');
fprintf(ftxt,'hdr.WIPdouble = '); fprintf(ftxt,'%.4f ', hdr.WIPdbl); fprintf(ftxt, '\n');
fprintf(ftxt,'hdr.sf = %.6f\n', hdr.sf);
fprintf(ftxt,'hdr.dwus = %i\n', hdr.dwus);
fprintf(ftxt,'hdr.sw = %.1f\n', hdr.sw);
fprintf(ftxt,'hdr.MRS_N = %i\n', hdr.MRS_N);
fprintf(ftxt,'hdr.VOIthick = %i\n', hdr.VOIthick);
fprintf(ftxt,'hdr.VOIphase = %i\n', hdr.VOIphase);
fprintf(ftxt,'hdr.VOIread = %i\n', hdr.VOIread);
fprintf(ftxt,'hdr.gamma = %.6f\n', hdr.gamma);
fprintf(ftxt,'hdr.avg = %i\n', hdr.avg);
fprintf(ftxt,'hdr.reps = %i\n', hdr.nreps);
fprintf(ftxt,'hdr.nrcvrs = %i\n', hdr.nrcvrs);

fclose(ftxt);

return
end
