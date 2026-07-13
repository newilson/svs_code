function [hdrout,images,dicomhdr,pathname] = readdicomfiles2d(pathname)


if (nargin < 1)
    str1 = 'Choose directory with 2D dicom images';
    pathname = uigetdir(str1); % file chooserend
end
files = dir(pathname);
nfiles = length(files);
files.name;
i1 = 0;
for ifile = 1:nfiles
    ascstr = char(files(ifile).name);
    if ( ascstr(1) ~= '.' && ascstr(1) ~= '$')
        fullfilen = fullfile(pathname,ascstr);
        if ( ~isfolder(fullfilen) && isdicom( fullfilen) )
            i1 = i1+1;
            filenames{i1}= files(ifile).name;
        end
    end
    clear ascstr
end

nfiles = i1;

% Read header from first file
fullfilen = fullfile( pathname ,char(filenames{1}) );


tmphdr = dicominfo(fullfilen);
A = dicomread(fullfilen);
hdr.n1 = size(A,1);
hdr.n2 = size(A,2);
hdr.trajectory = 'cartesian';

% NW
if isfield(tmphdr,'SoftwareVersions')
    hdr.swversion = tmphdr.SoftwareVersions;
elseif isfield(tmphdr,'ImplementationVersionName') && contains(tmphdr.ImplementationVersionName,'XA')
    hdr.swversion = tmphdr.ImplementationVersionName;
end

% % NW check for enhanced dicom
% isEnhanced = false;
% if isfield(tmphdr,'ImplementationVersionName') && contains(tmphdr.ImplementationVersionName, 'XA')
%     isEnhanced = true;
%     hdr.swversion = tmphdr.ImplementationVersionName;
% end

if ( strfind(fullfilen,'.dcm') )
    txtfilen = regexprep(fullfilen,'\.dcm','\.txt');
elseif ( strfind(fullfilen,'.IMA') )
    txtfilen = regexprep(fullfilen,'\.IMA','\.txt');
else
    txtfilen = [pathname filesep 'PVTHDR.txt'];
end
% txtfilen

if( fopen(txtfilen,'r') == -1 )
    ftxt = fopen(txtfilen,'w');
    
    % NW
    if isfield(tmphdr,'SharedFunctionalGroupsSequence') && isfield(tmphdr.SharedFunctionalGroupsSequence,'Item_1') && isfield(tmphdr.SharedFunctionalGroupsSequence.Item_1,'Private_0021_10fe') && isa(tmphdr.SharedFunctionalGroupsSequence.Item_1.Private_0021_10fe.Item_1.Private_0021_1019,'uint8')
        pvthdr = char(tmphdr.SharedFunctionalGroupsSequence.Item_1.Private_0021_10fe.Item_1.Private_0021_1019)';
    elseif isfield(tmphdr,'Private_0021_1019') && isa(tmphdr.Private_0021_1019,'uint8')
        pvthdr = char(tmphdr.Private_0021_1019)';
    elseif isfield(tmphdr,'Private_0029_1020') && isa(tmphdr.Private_0029_1020,'uint8')
        pvthdr = char(tmphdr.Private_0029_1020)';
    else
        error('did not find private header')
    end
    indasc1 = strfind(pvthdr,'### ASCCONV BEGIN');
    indasc2 = strfind( pvthdr(indasc1(1):length(pvthdr)), '### ASCCONV END' );
    if ~isequal(ftxt,-1)
        fwrite( ftxt, pvthdr(indasc1(1):indasc1(1)+indasc2(1)+19) );
        fclose(ftxt);
    end
    clear pvthdr;
end
% Initialize WIP values
for i = 1:64
    hdr.WIPlong(i) = 0;
    if (i < 17)
        hdr.WIPdbl(i)=0;
    end
end

% fprintf('Reading headers from %s ...\n',txtfilen);

hdr.ndim = 2; % initialization

fid2 = fopen(txtfilen,'r');
str1 = fgetl(fid2);
while( ~isempty(str1) &&  ~contains(str1,'ASCCONV END')  )
    str1 = fgetl(fid2);
    
    tnamestr = 'tSequenceFileName';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'""')+2;
        tindex2 = length(str1)-2;
        tname = str1(tindex1(1):tindex2);
        hdr.SeqName = tname;
    end
    
    tnamestr = 'tProtocolName';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'""')+2;
        tindex2 = length(str1)-2;
        tname = str1(tindex1(1):tindex2);
        tname = regexprep(tname,'+AF8-','-');
        hdr.ProtocolName = tname;
    end
    
    tnamestr = 'sProtConsistencyInfo.tBaselineString';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'""')+2;
        tindex2 = length(str1)-2;
        tname = str1(tindex1(1):tindex2);
        hdr.swversion = tname;
    end
    
    tnamestr = 'sProtConsistencyInfo.tMeasuredBaselineString'; % NW added 11/07/16
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'""')+2;
        tindex2 = length(str1)-2;
        tname = str1(tindex1(1):tindex2);
        if isfield(hdr,'swversion') % already exists
            warning('multiple definitions of software version')
            disp(hdr.swversion)
            disp(tname)
        end
        hdr.swversion = tname;
    end
    
    tnamestr = 'sTXSPEC.asNucleusInfo[0].tNucleus';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'""')+2;
        tindex2 = length(str1)-2;
        tname = str1(tindex1(1):tindex2);
        hdr.Nucleus = tname;
        if (contains(tname,'1H'))
            hdr.gamma = 42.5756;
        elseif (contains(tname,'23N'))
            hdr.gamma = 11.2620;
        elseif (contains(tname,'17O'))
            hdr.gamma = 5.7716;
        elseif (contains(tname,'13C'))
            hdr.gamma = 10.7063;
        else
            hdr.gamma = 42.5756;
        end
    end
    
    tnamestr = 'sTXSPEC.asNucleusInfo[0].lFrequency';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.sf = (sscanf(tname,'%i',1))* 1.0e-6;
    end

    % NW 
    tnamestr = 'sProtConsistencyInfo.flNominalB0';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.B0 = (sscanf(tname,'%f',1));
    end
    
    tnamestr = 'sRXSPEC.alDwellTime[0]';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.dwus = 0.001*(sscanf(tname,'%i',1));
        hdr.sw = 1.0e6 / hdr.dwus;
    end
    
    tnamestr = 'lContrasts';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.contrasts = sscanf(tname,'%i',1);
    end
    
    tnamestr = 'alTR[0]';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.TR = sscanf(tname,'%i',1);
    end
    
    for iTE = 1:10
        tnamestr = sprintf('alTE[%i]',iTE-1);
        if (contains(str1,tnamestr))
            tindex1 = strfind(str1,'=')+1;
            tname = str1(tindex1(1):numel(str1));
            hdr.TEms(iTE) = double(sscanf(tname,'%i',1))*0.001;
            hdr.echoes = iTE;
        end
    end
    
    tnamestr = 'sSliceArray.asSlice[0].dThickness';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.slthk = (sscanf(tname,'%f',1));
    end
    
    tnamestr = 'sSliceArray.asSlice[0].dPhaseFOV';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.phfov = (sscanf(tname,'%f',1));
    end
    
    tnamestr = 'sSliceArray.asSlice[0].dReadoutFOV';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.rdfov = (sscanf(tname,'%f',1));
    end

    % NW
    tnamestr = 'sSliceArray.asSlice[0].sPosition.dSag';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1:numel(str1));
        hdr.posSag = (sscanf(tname,'%f',1));
    end

    tnamestr = 'sSliceArray.asSlice[0].sPosition.dCor';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1:numel(str1));
        hdr.posCor = (sscanf(tname,'%f',1));
    end

    tnamestr = 'sSliceArray.asSlice[0].sPosition.dTra';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1:numel(str1));
        hdr.posTra = (sscanf(tname,'%f',1));
    end
    
    tnamestr = 'sSliceArray.asSlice[0].dInPlaneRot';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1:numel(str1));
        hdr.rot = (sscanf(tname,'%f',1));
    end

    tnamestr = 'sSliceArray.lSize';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.nslices = sscanf(tname,'%i',1);
    end
    
    tnamestr = 'sKSpace.lBaseResolution';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.nx = sscanf(tname,'%i',1);
    end
    
    tnamestr = 'sKSpace.lPhaseEncodingLines';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.ny = sscanf(tname,'%i',1);
    end
    
    tnamestr = 'sKSpace.lImagesPerSlab';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.nz = sscanf(tname,'%i',1);
    end
    
    tnamestr = 'sKSpace.lRadialViews';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.nrad = sscanf(tname,'%i',1);
    end
    
%     hdr.ndim = 2;
    tnamestr = 'sKSpace.ucDimension';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        if (contains(tname,'0x'))
            tname1 = tname(4:numel(tname));
            ndimflag = hex2dec(tname1);
            hdr.ndim = 1;
            while(ndimflag > 1)
                hdr.ndim = hdr.ndim+1;
                ndimflag = ndimflag/2;
            end
        else
            hdr.ndim = sscanf(tname,'%i',1);
        end
    end
    
    tnamestr = 'sKSpace.ucTrajectory';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        if (contains(tname,'0x'))
            tname1 = tname(4:numel(tname));
            traj = hex2dec(tname1);
        else
            traj = sscanf(tname,'%i',1);
        end
        if (traj < 2)
            hdr.trajectory = 'cartesian';
        else
            hdr.trajectory = 'radial';
        end
    end
    
    tnamestr = 'lAverages';
    if (contains(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1(1):numel(str1));
        hdr.averages = (sscanf(tname,'%i',1));
    end
    
    tnamestr = 'lRepetitions';
    if (contains(str1,tnamestr))
        if ( ~contains(str1,'alRepetitions') )
            tindex1 = strfind(str1,'=')+1;
            tname = str1(tindex1(1):numel(str1));
            hdr.reps = (sscanf(tname,'%i',1)) + 1;
        end
    end
    
    tnamestr = sprintf('sWiPMemBlock.alFree');
    if ( contains(str1,tnamestr) )
        for ilong = 1:64
            tnamestr = sprintf('sWiPMemBlock.alFree[%i]',ilong-1);
            if (contains(str1,tnamestr))
                tindex1 = strfind(str1,'=')+1;
                tname = str1(tindex1(1):numel(str1));
                value = sscanf(tname,'%f',1);
                hdr.WIPlong(ilong) = value;
            end
        end
    end
    tnamestr = sprintf('sWiPMemBlock.adFree');
    if ( contains(str1,tnamestr) )
        for idbl = 1:16
            tnamestr = sprintf('sWiPMemBlock.adFree[%i]',idbl-1);
            if (contains(str1,tnamestr))
                tindex1 = strfind(str1,'=')+1;
                tname = str1(tindex1(1):numel(str1));
                value = sscanf(tname,'%f',1);
                hdr.WIPdbl(idbl) = value;
            end
        end
    end

    tnamestr = sprintf('sWipMemBlock.alFree');
    if ( contains(str1,tnamestr) )
        for ilong = 1:64
            tnamestr = sprintf('sWipMemBlock.alFree[%i]',ilong-1);
            if (contains(str1,tnamestr))
                tindex1 = strfind(str1,'=')+1;
                tname = str1(tindex1(1):numel(str1));
                value = sscanf(tname,'%f',1);
                hdr.WIPlong(ilong) = value;
            end
        end
    end
    tnamestr = sprintf('sWipMemBlock.adFree');
    if ( contains(str1,tnamestr) )
        for idbl = 1:16
            tnamestr = sprintf('sWipMemBlock.adFree[%i]',idbl-1);
            if (contains(str1,tnamestr))
                tindex1 = strfind(str1,'=')+1;
                tname = str1(tindex1(1):numel(str1));
                value = sscanf(tname,'%f',1);
                hdr.WIPdbl(idbl) = value;
            end
        end
    end
    
end  % While ASCCONV
fclose(fid2);

if( ~isfield(hdr,'echoes') )
    hdr.echoes = 1;
end

hdr.nTE = hdr.echoes;

if ( isempty( strfind(hdr.trajectory, 'radi') ) )
    hdr.nrad = 0;
end
if (hdr.ndim < 3)
    hdr.nz = 1;
end
if ( ~isfield(hdr, 'reps'))
    hdr.reps = 1;
end

hdrout = hdr;

% if (nfiles ~= hdr.reps * hdr.nTE )
%     warndlg(sprintf('Number of files(%i) is not consistent with reps(%i)*echoes(%i). No results found.',nfiles,hdr.reps,hdr.nTE));
%     return
% end
% 
% images = zeros(hdr.n1,hdr.n2,hdr.reps);
images = zeros(hdr.n1,hdr.n2,nfiles); % NW since B0MC has 3 files but only 1 rep

hwaitbar = waitbar( 0,sprintf('Reading %i image files ...',nfiles) );
% for irep = 1: hdr.reps
for irep = 1:nfiles % NW
        ind1 = irep;
        waitbar(ind1/nfiles);
        fullfilen = fullfile( pathname,char(filenames{ind1}) );
        if ( isdicom(fullfilen) )
            if (ind1 == 1)
                dicomhdr = dicominfo(fullfilen);
            end
            A = double(dicomread(fullfilen));
            [~,fname,~] = fileparts(fullfilen);
            if isnan(str2double(fname)) % XA names files by rep number (NW)
                images(:,:,irep) = A;
            else
                images(:,:,str2double(fname)) = A;
            end
        end
end % for irep
close(hwaitbar);

% disp('readdicomimage2d done');
end

