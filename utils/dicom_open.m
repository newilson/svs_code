function fd = dicom_open(strFilename)

%
% Gets a file descriptor to a DICOM file and check
% it vaguely for DICOM like properties...
%

fd = fopen(strFilename, 'rb');

if -1 == fd 
    fprintf('\nFailed to open DICOM file.\n');
    return;
end

% lose the header
hdr = fread(fd, 128, 'uchar');

if strcmp(sprintf('%c', fread(fd, 4, 'schar')),'DICM')
    fprintf('\nFile appears to be valid DICOM.\n');
else
    fprintf('\nFile does NOT appear to be valid DICOM.\n');
    fclose(fd);
    fd = -1;
    return;
end