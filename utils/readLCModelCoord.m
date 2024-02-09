function [out, filename] = readLCModelCoord(filename)
% NW

if nargin<1 || isempty(filename)
    [filename,pathname] = uigetfile({'.coord'},'Choose coord file');
    filename = fullfile(pathname,filename);
end

strtext1 = 'lines in following concentration table';
strtext2 = 'lines in following misc. output table';
strtext3 = 'points on ppm-axis';
strtext4 = 'Conc. =';

npts = [];

fid = fopen(filename);
if fid==-1
    error('unable to open file')
end

tline = fgetl(fid);
metCount = 0;

while 1

    if ischar(tline) && contains(tline,strtext1)
        ind = strfind(tline,strtext1);
        Lines = str2double(tline(1:ind(1)-1));

        for ii=1:Lines
            tline = fgetl(fid);
            if ii>1
                tempstr = split(strtrim(tline),' ');
                tempstr = tempstr(~cellfun('isempty',tempstr));
                if length(tempstr)~=4
                    warning('unexpected number of table elements')
                end
                out.conc(ii-1) = str2double(tempstr{1});
                out.crlb(ii-1) = str2double(tempstr{2}(1:end-1));
                out.ratio(ii-1) = str2double(tempstr{3});
                out.met{ii-1} = tempstr{4};
            end
        end
    end

    if ischar(tline) && contains(tline,strtext2)
        ind = strfind(tline,strtext2);
        Lines = str2double(tline(1:ind(1)-1));
        for ii=1:Lines
            tline = fgetl(fid);
        end
    end

    if ischar(tline) && contains(tline,strtext3)
        ind = strfind(tline,strtext3);
        npts = str2double(tline(1:ind(1)-1));
        data = [];
        tline = fgetl(fid);
        while 1
            temp = split(strtrim(tline),' ');
            temp = str2double(temp(~cellfun('isempty',temp)));
            data = cat(1,data,temp);
            if length(data)==npts
                break;
            end
            tline = fgetl(fid);
        end
        out.fit.ppm = data;

        data = [];
        fgetl(fid);
        tline = fgetl(fid);
        while 1
            temp = split(strtrim(tline),' ');
            temp = str2double(temp(~cellfun('isempty',temp)));
            data = cat(1,data,temp);
            if length(data)==npts
                break;
            end
            tline = fgetl(fid);
        end
        out.fit.specRaw = data;

        data = [];
        fgetl(fid);
        tline = fgetl(fid);
        while 1
            temp = split(strtrim(tline),' ');
            temp = str2double(temp(~cellfun('isempty',temp)));
            data = cat(1,data,temp);
            if length(data)==npts
                break;
            end
            tline = fgetl(fid);
        end
        out.fit.specFit = data;

        data = [];
        fgetl(fid);
        tline = fgetl(fid);
        while 1
            temp = split(strtrim(tline),' ');
            temp = str2double(temp(~cellfun('isempty',temp)));
            data = cat(1,data,temp);
            if length(data)==npts
                break;
            end
            tline = fgetl(fid);
        end
        out.fit.specBack = data;
    end

    if ischar(tline) && contains(tline,strtext4)
        metCount = metCount + 1;

        ind = strfind(tline, ' ');
        met = tline(ind(1)+1:ind(2)-1);

        data = [];
        tline = fgetl(fid);
        while 1
            temp = split(strtrim(tline),' ');
            temp = str2double(temp(~cellfun('isempty',temp)));
            data = cat(1,data,temp);
            if length(data)==npts
                break;
            end
            tline = fgetl(fid);
        end
        out.fit.met{metCount} = met;
        out.fit.specMet{metCount} = data;
    end


    tline = fgetl(fid);

    if feof(fid)
        break
    end

end