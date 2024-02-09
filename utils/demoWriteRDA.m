%% demo write rda files
clear,clc

%% SVS - same file
newName = 'MyTempRDA.rda';
info = dir('F:\lactate_bigvoxel\3T_lactate\s077\exercise\20220808\post\*.rda');

scan = 1;
fname = fullfile(info(scan).folder,info(scan).name);

[~,data1,hdr1] = Read_rda_file(fname);

writeRDAfromTemplate(data1,newName,fname);

[~,data2,hdr2] = Read_rda_file(newName);

isequal(data1,data2)
isequal(hdr1,hdr2)

%% CSI - same file

%% Modified data dimensions