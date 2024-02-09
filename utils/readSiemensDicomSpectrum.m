function [complex_fid,info,hdr] = readSiemensDicomSpectrum(filename)
% NW

if nargin<1 || isempty(filename)
%     [filename,pathname] = uigetfile({'*.IMA','*.dcm'},'Choose Spectrum Dicom file');
    [filename,pathname] = uigetfile;
    filename = fullfile(pathname,filename);
end

%% data
fd = dicom_open(filename);
complex_fid = dicom_get_spectrum_siemens(fd);

%% header
info = dicominfo(filename);

if nargout>2
    if ( strfind(filename,'.dcm') )
        txtname = regexprep(filename,'.dcm','.txt');
    elseif ( strfind(filename,'.IMA') )
        txtname = regexprep(filename,'.IMA','.txt');
    else
        txtname = [pathname filesep 'PVTHDR.txt'];
    end
    
    % write private header to text file
    if( fopen(txtname,'r') == -1 )
        ftxt = fopen(txtname,'w');
        if isfield(info,'Private_0029_1120')
            pvthdr = char(info.Private_0029_1120)'; % NW changed from 1020 to 1120
        elseif isfield(info,'Private_0029_1020')
            pvthdr = char(info.Private_0029_1020)';
        else
            error('unknown private header')
        end
        indasc1 = strfind(pvthdr,'### ASCCONV BEGIN');
        indasc2 = strfind( pvthdr(indasc1(1):length(pvthdr)), '### ASCCONV END' );
        fwrite( ftxt, pvthdr(indasc1(1):indasc1(1)+indasc2(1)+19) );
        fclose(ftxt);
        clear pvthdr;
    end
    
    % Initialize WIP values
    for i = 1:64
        hdr.WIPlong(i) = 0;
        if (i < 17)
            hdr.WIPdbl(i)=0;
        end
    end
    
    fid2 = fopen(txtname,'r');
    str1 = fgetl(fid2);
    nloop = 0;
    while( ~isempty(str1) && isempty( strfind(str1,'ASCCONV END') ) )
        str1 = fgetl(fid2);
        nloop = nloop+1;
        
        tnamestr = 'sRXSPEC.alDwellTime[0]';
        if (~isempty(strfind(str1,tnamestr)))
            tindex1 = strfind(str1,'=')+1;
            tname = str1(tindex1:numel(str1));
            hdr.dwus = 0.001*(sscanf(tname,'%i',1));
            hdr.bw = 1.0e6 / hdr.dwus;
        end
        
        tnamestr = 'sTXSPEC.asNucleusInfo[0].lFrequency';
        if (~isempty(strfind(str1,tnamestr)))
            tindex1 = strfind(str1,'=')+1;
            tname = str1(tindex1:numel(str1));
            hdr.f0 = (sscanf(tname,'%i',1))* 1.0e-6;
        end
        
        tnamestr = 'sSpecPara.lVectorSize';
        if (~isempty(strfind(str1,tnamestr)))
            tindex1 = strfind(str1,'=')+1;
            tname = str1(tindex1:numel(str1));
            hdr.npts = (sscanf(tname,'%i',1));
        end
        
        tnamestr = 'sSpecPara.lFinalMatrixSizePhase';
        if (~isempty(strfind(str1,tnamestr)))
            tindex1 = strfind(str1,'=')+1;
            tname = str1(tindex1:numel(str1));
            hdr.nx = (sscanf(tname,'%i',1));
        end
        
        tnamestr = 'sSpecPara.lFinalMatrixSizeRead';
        if (~isempty(strfind(str1,tnamestr)))
            tindex1 = strfind(str1,'=')+1;
            tname = str1(tindex1:numel(str1));
            hdr.ny = (sscanf(tname,'%i',1));
        end
        
        tnamestr = 'sSpecPara.lFinalMatrixSizeSlice';
        if (~isempty(strfind(str1,tnamestr)))
            tindex1 = strfind(str1,'=')+1;
            tname = str1(tindex1:numel(str1));
            hdr.nz = (sscanf(tname,'%i',1));
        end
        
%         sSpecPara.ucRemoveOversampling	 = 	0x0
        tnamestr = 'sSpecPara.ucRemoveOversampling';
        if (~isempty(strfind(str1,tnamestr)))
            tindex1 = strfind(str1,'=')+1;
            tname = str1(tindex1:numel(str1));
            tindex2 = strfind(tname,'0x');
            if (~isempty(tindex2))
                tname1 = tname(tindex2+2:end);
                hdr.removeOS = hex2dec(tname1);
            end
            hdr.removeOS = (sscanf(tname,'%i',1));
        end
        
        tnamestr = sprintf('sWiPMemBlock.alFree');
        if ( ~isempty(strfind(str1,tnamestr)) )
            for ilong = 1:64
                tnamestr = sprintf('sWiPMemBlock.alFree[%i]',ilong-1);
                if (~isempty(strfind(str1,tnamestr)))
                    tindex1 = strfind(str1,'=')+1;
                    tname = str1(tindex1:numel(str1));
                    value = sscanf(tname,'%f',1);
                    hdr.WIPlong(ilong) = value;
                end
            end
        end
        tnamestr = sprintf('sWiPMemBlock.adFree');
        if ( ~isempty(strfind(str1,tnamestr)) )
            for idbl = 1:16
                tnamestr = sprintf('sWiPMemBlock.adFree[%i]',idbl-1);
                if (~isempty(strfind(str1,tnamestr)))
                    tindex1 = strfind(str1,'=')+1;
                    tname = str1(tindex1:numel(str1));
                    value = sscanf(tname,'%f',1);
                    hdr.WIPdbl(idbl) = value;
                end
            end
        end
        
        tnamestr = sprintf('sWipMemBlock.alFree');
        if ( ~isempty(strfind(str1,tnamestr)) )
            for ilong = 1:64
                tnamestr = sprintf('sWipMemBlock.alFree[%i]',ilong-1);
                if (~isempty(strfind(str1,tnamestr)))
                    tindex1 = strfind(str1,'=')+1;
                    tname = str1(tindex1:numel(str1));
                    value = sscanf(tname,'%f',1);
                    hdr.WIPlong(ilong) = value;
                end
            end
        end
        tnamestr = sprintf('sWipMemBlock.adFree');
        if ( ~isempty(strfind(str1,tnamestr)) )
            for idbl = 1:16
                tnamestr = sprintf('sWipMemBlock.adFree[%i]',idbl-1);
                if (~isempty(strfind(str1,tnamestr)))
                    tindex1 = strfind(str1,'=')+1;
                    tname = str1(tindex1:numel(str1));
                    value = sscanf(tname,'%f',1);
                    hdr.WIPdbl(idbl) = value;
                end
            end
        end
        
    end  % While ASCCONV
    fclose(fid2);
    if isfield(hdr,'removeOS') % NW
        if hdr.removeOS
            hdr.bw = hdr.bw/2;
            hdr.dwus = hdr.dwus*2;
		else
			hdr.npts = hdr.npts*2;
        end
    end
        
end