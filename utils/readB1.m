function [B1map,hdr,alpha,images,pathname1] = readB1(pathname1,dofilt,roimask)

if nargin<1 || isempty(pathname1)
    [hdr,images,dicomhdr,pathname1] = readdicomfiles2d;
elseif isstruct(pathname1)
    hdr = pathname1.hdr;
    images = pathname1.images;
    if isfield(pathname1,'dicomhdr')
        dicomhdr = pathname1.dicomhdr;
    else
        dicomhdr = [];
    end
else
    [hdr,images,dicomhdr] = readdicomfiles2d(pathname1);
end
if nargin<2, dofilt = false; end

nx = size(images,1); ny = size(images,2);
nz = size(images,3)/hdr.reps;

if nargin<3 || isempty(roimask)
    if nz>1
        roimask = ones(nx,ny,nz);
    else
        roimask = ones(nx,ny);
    end
end

if isstr(pathname1) && contains(lower(pathname1),'flip30') && isequal(size(images,3),1)
    alpha = pi/6.0;
    pathparts = strsplit(pathname1,'_');
    for jj=1:length(pathparts)
        if strcmp(pathparts{jj},'flip30')
            pathparts{jj} = 'flip60';
        elseif strcmp(pathparts{jj},'FLIP30')
            pathparts{jj} = 'FLIP60';
        end
    end
    pathparts{end} = num2str(str2double(pathparts{end})+1,'%04d');
    pathname2 = strjoin(pathparts,'_');
    [~,image2] = readdicomfiles2d(pathname2);
    images = cat(3,images,image2);
    if dofilt
        si = size(images);
        if length(dofilt)==3
            if dofilt(1)~=si(1)
                dofilt(1) = si(1);
            end
            filt = repmat(NWSiemensRawFilter2d(dofilt(1),dofilt(2),dofilt(3)),[1 1 si(3)]);
        elseif length(dofilt)==6
            if dofilt(1)~=si(1)
                dofilt(1) = si(1);
            end
            if dofilt(2)~=si(1)
                dofilt(2) = si(1);
            end
            filt = repmat(NWSiemensRawFilter2d(dofilt(1:2),dofilt(3:4),dofilt(5:6)),[1 1 si(3)]);
        else
            filt = repmat(NWSiemensRawFilter2d(si(1:2),si(1:2),si(1:2)/2),[1 1 si(3)]);
        end
        images = abs(fft2c(ifft2c(images).*filt));
    end
    if nz>1
        images = reshape(images,nx,ny,nz,[]);
        B1map = zeros(nx,ny,nz);
        for jj=1:nz
            B1map(:,:,jj) = calc_B1map(images(:,:,jj,1),images(:,:,jj,2),ones(size(squeeze(images(:,:,jj,1)))),alpha);
        end
    else
        B1map = calc_B1map(images(:,:,1),images(:,:,2),ones(size(squeeze(images(:,:,1)))),alpha);
    end
elseif isstr(pathname1) && contains(lower(pathname1),'flip60') && isequal(size(images,3),1)
    alpha = pi/6.0;
    pathparts = strsplit(pathname1,'_');
    for jj=1:length(pathparts)
        if strcmp(pathparts{jj},'flip60')
            pathparts{jj} = 'flip30';
        elseif strcmp(pathparts{jj},'FLIP60')
            pathparts{jj} = 'FLIP30';
        end
    end
    pathparts{end} = num2str(str2double(pathparts{end})-1,'%04d');
    pathname2 = strjoin(pathparts,'_');
    [~,image2] = readdicomfiles2d(pathname2);
    images = cat(3,image2,images);
    if dofilt
        si = size(images);
        if length(dofilt)==3
            if dofilt(1)~=si(1)
                dofilt(1) = si(1);
            end
            filt = repmat(NWSiemensRawFilter2d(dofilt(1),dofilt(2),dofilt(3)),[1 1 si(3)]);
        elseif length(dofilt)==6
            if dofilt(1)~=si(1)
                dofilt(1) = si(1);
            end
            if dofilt(2)~=si(1)
                dofilt(2) = si(1);
            end
            filt = repmat(NWSiemensRawFilter2d(dofilt(1:2),dofilt(3:4),dofilt(5:6)),[1 1 si(3)]);
        else
            filt = repmat(NWSiemensRawFilter2d(si(1:2),si(1:2),si(1:2)/2),[1 1 si(3)]);
        end
        images = abs(fft2c(ifft2c(images).*filt));
    end
    if nz>1
        images = reshape(images,nx,ny,nz,[]);
        B1map = zeros(nx,ny,nz);
        for jj=1:nz
            B1map(:,:,jj) = calc_B1map(images(:,:,jj,1),images(:,:,jj,2),ones(size(squeeze(images(:,:,jj,1)))),alpha);
        end
    else
        B1map = calc_B1map(images(:,:,1),images(:,:,2),ones(size(squeeze(images(:,:,1)))),alpha);
    end
else
    if dofilt
        si = size(images);
        if length(dofilt)==3
            if dofilt(1)~=si(1)
                dofilt(1) = si(1);
            end
            filt = repmat(NWSiemensRawFilter2d(dofilt(1),dofilt(2),dofilt(3)),[1 1 si(3)]);
        elseif length(dofilt)==6
            if dofilt(1)~=si(1)
                dofilt(1) = si(1);
            end
            if dofilt(2)~=si(1)
                dofilt(2) = si(1);
            end
            filt = repmat(NWSiemensRawFilter2d(dofilt(1:2),dofilt(3:4),dofilt(5:6)),[1 1 si(3)]);
        else
            filt = repmat(NWSiemensRawFilter2d(si(1:2),si(1:2),si(1:2)/2),[1 1 si(3)]);
        end
        images = abs(fft2c(ifft2c(images).*filt));
    end
    alpha = hdr.WIPdbl(2)*pi/180;
    if nz>1
        images = reshape(images,nx,ny,nz,[]);
        B1map = zeros(nx,ny,nz);
        for jj=1:nz
            if isequal(size(images,4),2)
                B1map(:,:,jj) = calc_B1map(images(:,:,jj,1),images(:,:,jj,2),roimask(:,:,jj),alpha);
            elseif isequal(size(images,4),3)
                B1map(:,:,jj) = calc_B1map_new(images(:,:,jj,1),images(:,:,jj,2),images(:,:,jj,3),roimask(:,:,jj),alpha);
            else
                B1map = 0;
            end
        end
    else        
        if isequal(size(images,3),2)
            B1map = calc_B1map(images(:,:,1),images(:,:,2),ones(size(squeeze(images(:,:,1)))),alpha);
        elseif isequal(size(images,3),3)
            B1map = calc_B1map_new(images(:,:,1),images(:,:,2),images(:,:,3),ones(size(squeeze(images(:,:,1)))),alpha);
        else
            B1map = 0;
        end        
    end
end

if nz>1
    for jj=1:nz
        B1map(:,:,jj) = anisodiff(B1map(:,:,jj),20,50,0.03,1);
    end
else
    B1map = anisodiff(B1map,20,50,0.03,1);
end


end