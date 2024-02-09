function hdr = dicom_hdr(fname)
% NW

if contains(fname,'.dcm')
    txtname = regexprep(fname,'.dcm','.txt');
elseif contains(fname,'.IMA')
    txtname = regexprep(fname,'.IMA','.txt');
end

if fopen(txtname,'r')==-1
    ftxt = fopen(txtname,'w');
    info = dicominfo(fname);
    if isfield(info,'Private_0029_1120')
        pvthdr = char(info.Private_0029_1120)';
    elseif isfield(info,'Private_0029_1020')
        pvthdr = char(info.Private_0029_1020)';
    end
    indasc1 = strfind(pvthdr,'### ASCCONV BEGIN');
    indasc2 = strfind( pvthdr(indasc1(1):length(pvthdr)), '### ASCCONV END' );
    fwrite( ftxt, pvthdr(indasc1(1):indasc1(1)+indasc2(1)+19) );
    fclose(ftxt);
    clear pvthdr
end
fid2 = fopen(txtname,'r');

str1 = fgetl(fid2);
while( ~isempty(str1) && isempty( strfind(str1,'ASCCONV END')) )
    str1 = fgetl(fid2);
%     disp(str1);
    
    tnamestr = 'tSequenceFileName';
    if ~isempty(strfind(str1,tnamestr))
        tindex1 = strfind(str1,'""')+2;
        tindex2 = length(str1)-2;
        tname = str1(tindex1:tindex2);
        hdr.SeqName = tname;
    end
    
    tnamestr = 'tProtocolName';
    if ~isempty(strfind(str1,tnamestr))
        tindex1 = strfind(str1,'""')+2;
        tindex2 = length(str1)-2;
        tname = str1(tindex1:tindex2);
        tname = regexprep(tname,'+AF8-','-');
        hdr.ProtocolName = tname;
    end
    
    tnamestr = 'sProtConsistencyInfo.tBaselineString';
    if ~isempty(strfind(str1,tnamestr))
        tindex1 = strfind(str1,'""')+2;
        tindex2 = length(str1)-2;
        tname = str1(tindex1:tindex2);
        hdr.swversion = tname;
    end
    
    tnamestr = 'sProtConsistencyInfo.tMeasuredBaselineString'; % NW added 11/07/16
    if ~isempty(strfind(str1,tnamestr))
        tindex1 = strfind(str1,'""')+2;
        tindex2 = length(str1)-2;
        tname = str1(tindex1:tindex2);
        hdr.swversion = tname;
    end
    
    tnamestr = 'sRXSPEC.alDwellTime[0]';
    if ~isempty(strfind(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1:numel(str1));
        hdr.dwus = 0.001*(sscanf(tname,'%i',1));
        hdr.sw = 1.0e6 / hdr.dwus;
    end
    
    tnamestr = 'sTXSPEC.asNucleusInfo[0].lFrequency';
    if ~isempty(strfind(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1:numel(str1));
        hdr.sf = (sscanf(tname,'%i',1))* 1.0e-6;
    end
    
    for iTE = 1:10
        tnamestr = sprintf('alTE[%i]',iTE-1);
        if ~isempty(strfind(str1,tnamestr))
            tindex1 = strfind(str1,'=')+1;
            tname = str1(tindex1:numel(str1));
            hdr.TEms(iTE) = double(sscanf(tname,'%i',1))*0.001;
            hdr.echoes = iTE;
        end
    end
    
    tnamestr = 'sSpecPara.sVoI.dThickness';
    if ~isempty(strfind(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1:numel(str1));
        hdr.voisl = (sscanf(tname,'%i',1));
    end
    
    tnamestr = 'sSpecPara.sVoI.dPhaseFOV';
    if ~isempty(strfind(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1:numel(str1));
        hdr.voiph = (sscanf(tname,'%i',1));
    end
    
    tnamestr = 'sSpecPara.sVoI.dReadoutFOV';
    if ~isempty(strfind(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1:numel(str1));
        hdr.voiro = (sscanf(tname,'%i',1));
    end
    
    tnamestr = 'lAverages';
    if ~isempty(strfind(str1,tnamestr))
        tindex1 = strfind(str1,'=')+1;
        tname = str1(tindex1:numel(str1));
        hdr.averages = (sscanf(tname,'%i',1));
    end
    
    tnamestr = 'lRepetitions';
    if ~isempty(strfind(str1,tnamestr))
        if isempty(strfind(str1,'alRepetitions') )
            tindex1 = strfind(str1,'=')+1;
            tname = str1(tindex1:numel(str1));
            hdr.reps = (sscanf(tname,'%i',1)) + 1;
        end
    end
    
    tnamestr = sprintf('sWiPMemBlock.alFree');
    if ~isempty(strfind(str1,tnamestr) )
        for ilong = 1:64
            tnamestr = sprintf('sWiPMemBlock.alFree[%i]',ilong-1);
            if ~isempty(strfind(str1,tnamestr))
                tindex1 = strfind(str1,'=')+1;
                tname = str1(tindex1:numel(str1));
                value = sscanf(tname,'%f',1);
                hdr.WIPlong(ilong) = value;
            end
        end
    end
    tnamestr = sprintf('sWiPMemBlock.adFree');
    if ~isempty(strfind(str1,tnamestr) )
        for idbl = 1:16
            tnamestr = sprintf('sWiPMemBlock.adFree[%i]',idbl-1);
            if ~isempty(strfind(str1,tnamestr))
                tindex1 = strfind(str1,'=')+1;
                tname = str1(tindex1:numel(str1));
                value = sscanf(tname,'%f',1);
                hdr.WIPdbl(idbl) = value;
            end
        end
    end

    tnamestr = sprintf('sWipMemBlock.alFree');
    if ~isempty(strfind(str1,tnamestr) )
        for ilong = 1:64
            tnamestr = sprintf('sWipMemBlock.alFree[%i]',ilong-1);
            if ~isempty(strfind(str1,tnamestr))
                tindex1 = strfind(str1,'=')+1;
                tname = str1(tindex1:numel(str1));
                value = sscanf(tname,'%f',1);
                hdr.WIPlong(ilong) = value;
            end
        end
    end
    tnamestr = sprintf('sWipMemBlock.adFree');
    if ~isempty(strfind(str1,tnamestr) )
        for idbl = 1:16
            tnamestr = sprintf('sWipMemBlock.adFree[%i]',idbl-1);
            if ~isempty(strfind(str1,tnamestr))
                tindex1 = strfind(str1,'=')+1;
                tname = str1(tindex1:numel(str1));
                value = sscanf(tname,'%f',1);
                hdr.WIPdbl(idbl) = value;
            end
        end
    end
    
end  % While ASCCONV
fclose(fid2);