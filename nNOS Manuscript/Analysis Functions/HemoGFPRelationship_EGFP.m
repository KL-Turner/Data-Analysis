function [Results_IntSig_GCaMP] = HemoGFPRelationship_EGFP(animalID,group,set,rootFolder,delim,Results_IntSig_GCaMP)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
LH_blueChannel = []; RH_blueChannel = [];
LH_greenChannel = []; RH_greenChannel = [];
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procCataFileIDs(aa,:);
    load(procDataFileID)
    LH_blueChannel = cat(1,LH_blueChannel,ProcData.data.GCaMP7s.LH);
    RH_blueChannel = cat(1,RH_blueChannel,ProcData.data.GCaMP7s.RH);
    LH_greenChannel = cat(1,LH_greenChannel,ProcData.data.CBV.LH);
    RH_greenChannel = cat(1,RH_greenChannel,ProcData.data.CBV.RH);
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_GFP.mat','Results_GFP')
cd([rootFolder delim 'Data'])