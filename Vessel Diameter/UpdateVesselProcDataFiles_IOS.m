function [] = UpdateVesselProcDataFiles_IOS(procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

for aa = 1:size(procDataFileIDs,1)
     procDataFileID = procDataFileIDs(aa,:);
    [animalID,~,fileID] = GetFileInfo_IOS(procDataFileID);
    load(procDataFileID)
    vesselDataFileID = [animalID '_' fileID '_VesselDiam.mat'];
    load(vesselDataFileID)
    ProcData.data.CBV.Veinous = VesselData.V1';
    disp(['Updating ProcData file ' num2str(aa) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    save(procDataFileID,'ProcData')
end

end