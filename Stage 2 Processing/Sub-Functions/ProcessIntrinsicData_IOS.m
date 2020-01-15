function [] = ProcessIntrinsicData_IOS(animalID,imagingType,rawDataFileIDs,procDataFileIDs)
% CBV from ROIs.
if strcmp(imagingType,'bilateral') == true
    ROInames = {'LH','RH','LH_Cement','RH_Cement','Cement'};
elseif strcmp(imagingType,'single') == true
    ROInames = {'Barrels','Cement'};
end
ROIFileDir = dir('*_ROIs.mat');
if isempty(ROIFileDir) == true
    ROIs = [];
else
    ROIFileName = {ROIFileDir.name}';
    ROIFileID = char(ROIFileName);
    load(ROIFileID);
end
[ROIs] = CheckROIDates_IOS(animalID,ROIs,ROInames);

%% BLOCK PURPOSE: [2] Extract CBV data from each ROI for each RawData file in the directory that hasn't been processed yet.
disp('Analyzing Block [2] Extracting mean reflectance data from each ROI.'); disp(' ')
ExtractCBVData_IOS(ROIs,ROInames,rawDataFileIDs)

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    load(procDataFileID)
    CBVfields = fieldnames(ProcData.data.CBV);
    for b = 1:length(CBVfields)
        ProcData.data.CBV.(CBVfields{b}(1:end - 6)) = RawData.data.CBV.(CBVfields{b})(1:end - 1);
    end
    CheckForNaNs_IOS(ProcData);
end

end