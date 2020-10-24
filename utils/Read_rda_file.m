function [rda,Time_domain_data] = Read_rda_file(rda_filename)
%
% Read spectroscopy data from Siemens machine
%
% As input take the rda file
% return the rda structure, which includes all the information about the rda file
%
% Read a .rda file
%
%
% MDR 11 June 2002
%

if(nargin < 1)
[filename , pathname ] = uigetfile('*.rda' , 'Select an RDA file');
rda_filename = [pathname , filename]; %'c:/data/spectroscopy/spec raw data/MrSpec.20020531.160701.rda'
end

fid = fopen(rda_filename);

head_start_text = '>>> Begin of header <<<';
head_end_text = '>>> End of header <<<';

tline = fgets(fid);

while (isempty(strfind(tline , head_end_text)))

tline = fgets(fid);

if ( isempty(strfind (tline , head_start_text)) + isempty(strfind (tline , head_end_text )) == 2)


% Store this data in the appropriate format

occurence_of_colon = findstr(':',tline);
variable = tline(1:occurence_of_colon-1) ;
value = tline(occurence_of_colon+1 : length(tline)) ;

switch variable
case { 'PatientName' }
% This is included for patient confidentiality reasons.
% It will always be possible to get the name from cross referencing the patient ID
rda.PatientName = 'Anonymous';
case { 'PatientID' , 'StudyDescription' , 'PatientBirthDate' , 'StudyDate' , 'StudyTime' , 'PatientAge' , 'SeriesDate' , ...
'SeriesTime' , 'SeriesDescription' , 'ProtocolName' , 'PatientPosition' , 'ModelName' , 'StationName' , 'InstitutionName' , ...
'DeviceSerialNumber', 'InstanceDate' , 'InstanceTime' , 'InstanceComments' , 'SequenceName' , 'SequenceDescription' , 'Nucleus' ,...
'TransmitCoil' }
eval(['rda.' , variable , ' = value ;']);
case { 'PatientSex' }
% Sex converter! (int to M,F,U)
switch value
case 0
rda.sex = 'Unknown';
case 1
rda.sex = 'Male';
case 2

rda.sex = 'Female';
end

case { 'SeriesNumber' , 'InstanceNumber' , 'AcquisitionNumber' , 'NumOfPhaseEncodingSteps' , 'NumberOfRows' , 'NumberOfColumns' , 'VectorSize' }
%Integers
eval(['rda.' , variable , ' = str2num(value) ;']);
case { 'PatientWeight' , 'TR' , 'TE' , 'TM' , 'DwellTime' , 'NumberOfAverages' , 'MRFrequency' , 'MagneticFieldStrength' , 'FlipAngle' , ...
'SliceThickness' , 'FoVHeight' , 'FoVWidth' , 'PercentOfRectFoV' , 'PixelSpacingRow' , 'PixelSpacingCol'}
%Floats
eval(['rda.' , variable , ' = str2num(value) ;']);
case {'SoftwareVersion[0]' }
rda.software_version = value;
case {'CSIMatrixSize[0]' }
rda.CSIMatrix_Size(1) = str2num(value);
case {'CSIMatrixSize[1]' }
rda.CSIMatrix_Size(2) = str2num(value);
case {'CSIMatrixSize[2]' }
rda.CSIMatrix_Size(3) = str2num(value);
case {'PositionVector[0]' }
rda.PositionVector(1) = str2num(value);
case {'PositionVector[1]' }
rda.PositionVector(2) = str2num(value);
case {'PositionVector[2]' }
rda.PositionVector(3) = str2num(value);
case {'RowVector[0]' }
rda.RowVector(1) = str2num(value);
case {'RowVector[1]' }
rda.RowVector(2) = str2num(value);
case {'RowVector[2]' }
rda.RowVector(3) = str2num(value);
case {'ColumnVector[0]' }
rda.ColumnVector(1) = str2num(value);
case {'ColumnVector[1]' }
rda.ColumnVector(2) = str2num(value);
case {'ColumnVector[2]' }
rda.ColumnVector(3) = str2num(value);

otherwise
% We don't know what this variable is. Report this just to keep things clear
disp(['Unrecognised variable ' , variable ]);
end

else
% Don't bother storing this bit of the output
end

end

%
% So now we should have got to the point after the header text
%
% Siemens documentation suggests that the data should be in a double complex format (8bytes for real, and 8 for imaginary?)
%

bytes_per_point = 16;
complex_data = fread(fid , rda.CSIMatrix_Size(1) * rda.CSIMatrix_Size(1) *rda.CSIMatrix_Size(1) *rda.VectorSize * 2 , 'double');
% Now convert this data into something meaningful

%Reshape so that we can get the real and imaginary separated
hmm = reshape(complex_data, 2 , rda.VectorSize , rda.CSIMatrix_Size(1) , rda.CSIMatrix_Size(2) , rda.CSIMatrix_Size(3) );

%Combine the real and imaginary into the complex matrix
hmm_complex = complex(hmm(1,:,:,:,:),hmm(2,:,:,:,:));

%Remove the redundant first element in the array
Time_domain_data = reshape(hmm_complex, rda.VectorSize , rda.CSIMatrix_Size(1) , rda.CSIMatrix_Size(2) , rda.CSIMatrix_Size(3));

%fread(fid , 1, 'double'); %This was a check to confirm that we had read all the data (it passed!)
fclose(fid);