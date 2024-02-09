% Run dat2matsvs first to get ws data

si = size(ws);
templateName = fullfile(infoRDA(2).folder,infoRDA(2).name);

for ii=1:si(2)
    newName = ['TimePoint' num2str(ii) 'of' num2str(si(2)) '.rda'];
    tempdata = conj(ifft(ifftshift(squeeze(ws(:,ii)),1),[],1));
    writeRDAfromTemplate(tempdata,newName,templateName);
end
disp('DONE')