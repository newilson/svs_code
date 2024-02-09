function writeRDAfromTemplate(data,rdaFname,rdaTemplateFname)
%
% writes RDA file copying header from existing rda file
% data should be [time x X x Y x Z] complex array of time/spatial data

if nargin<2
    rdaFname = fullfile(pwd,['Spectro' date '.rda']);
end
    
if nargin<3
    [filename,pathname] = uigetfile('*.rda','Select an RDA file');
    rdaTemplateFname = fullfile(pathname,filename);
end

head_start_text = '>>> Begin of header <<<';
head_end_text = '>>> End of header <<<';

fid1 = fopen(rdaTemplateFname);
if fid1==-1
    error('unable to open new file')
end

fid2 = fopen(rdaFname,'w');
if fid2==-1
    fclose(fid1);
    error('unable to open template file')
end


isVec = isvector(data);
if isVec
    data = data(:); % ensure column vector
    nx = 1;
    ny = 1;
    nz = 1;
elseif ndims(data)==3
    nx = size(data,2);
    ny = size(data,3);
    nz = 1;
elseif ndims(data)==4
    nx = size(data,2);
    ny = size(data,3);
    nz = size(data,4);
else
    error('unexpected data size')
end
npts = size(data,1);


% write header

tline = fgets(fid1);
fwrite(fid2,tline);

if ~contains(tline,head_start_text)
    warning('unexpected header opening')
end

while ~contains(tline,head_end_text)
        
    tline = fgets(fid1);
    
    if contains(tline,'VectorSize')
        ind = strfind(tline,':');
        ind = ind(1); % first instance only
        VectorSize = str2double(tline(ind+1:end));
        if ~isequal(VectorSize,npts)
            warning('VectorSize not equal to template')
            tline = [tline(1:ind) num2str(npts)];
        end
    end
    
    if contains(tline,'CSIMatrixSize[0]')
        ind = strfind(tline,':');
        ind = ind(1); % first instance only
        CSIMatrixSize0 = str2double(tline(ind+1:end));
        if ~isequal(CSIMatrixSize0,nx)
            warning('CSIMatrixSize0 not equal to template')
            tline = [tline(1:ind) num2str(nx)];
        end
    end
    
    if contains(tline,'CSIMatrixSize[1]')
        ind = strfind(tline,':');
        ind = ind(1); % first instance only
        CSIMatrixSize1 = str2double(tline(ind+1:end));
        if ~isequal(CSIMatrixSize1,ny)
            warning('CSIMatrixSize1 not equal to template')
            tline = [tline(1:ind) num2str(ny)];
        end
    end
    
    if contains(tline,'CSIMatrixSize[2]')
        ind = strfind(tline,':');
        ind = ind(1); % first instance only
        CSIMatrixSize2 = str2double(tline(ind+1:end));
        if ~isequal(CSIMatrixSize2,nz)
            warning('CSIMatrixSize2 not equal to template')
            tline = [tline(1:ind) num2str(nz)];
        end
    end    
    
    fwrite(fid2,tline,'char');
%     fprintf(fid2,'%s',tline);
    
end
fclose(fid1);
        
% write data
dataRI = zeros(2,npts,nx,ny,nz);
dataRI(1,:,:,:,:) = real(data);
dataRI(2,:,:,:,:) = imag(data);
fwrite(fid2,dataRI,'double');
fclose(fid2);




