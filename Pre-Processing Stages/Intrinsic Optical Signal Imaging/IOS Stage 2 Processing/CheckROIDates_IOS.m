function [ROIs] = CheckROIDates_IOS(animalID,ROIs,ROInames,lensMag,imagingType,imagingColors)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Create/Update ROIs.mat structure to verify all ROIs are drawn
%________________________________________________________________________________________________________________________

% character list of all WindowCam files
windowCamFilesDir = dir('*_WindowCam.bin');
windowCamDataFiles = {windowCamFilesDir.name}';
windowCamDataFileIDs = char(windowCamDataFiles);
% establish the number of unique days based on file IDs
[~,fileDates,~] = GetFileInfo_IOS(windowCamDataFileIDs);
[uniqueDays,~,DayID] = GetUniqueDays_IOS(fileDates);
firstsFileOfDay = cell(1,length(uniqueDays));
for a = 1:length(uniqueDays)
    FileInd = DayID == a;
    dayFilenames = windowCamDataFileIDs(FileInd,:);
    firstsFileOfDay(a) = {dayFilenames(1,:)};
end
% load existing ROI structure if it exists
ROIFileDir = dir('*_ROIs.mat');
ROIFileName = {ROIFileDir.name}';
ROIFileID = char(ROIFileName);
if exist(ROIFileID,'file')
    load(ROIFileID);
else
    ROIs = [];
end
% create the desired window ROI for each day if it doesn't yet exist
for b = 1:length(firstsFileOfDay)
    fileID = firstsFileOfDay{1,b};
    strDay = ConvertDate_IOS(fileID);
    for c = 1:length(ROInames)
        ROIname = [ROInames{1,c} '_' strDay];
        if ~isfield(ROIs,(ROIname))
            if any(strcmp(ROInames{1,c},{'LH','RH','frontalLH','frontalRH','Barrels'})) == true
                if strcmp(imagingColors,'RGB') == true
                    [ROIs] = PlaceTriWavelengthROIs_IOS(animalID,fileID,ROIs,ROInames,lensMag,imagingType);
                elseif strcmpi(imagingColors,'GB') == true
                    [ROIs] = PlaceDualWavelengthROIs_IOS(animalID,fileID,ROIs,ROInames,lensMag,imagingType);
                else
                    [ROIs] = CalculateROICorrelationMatrix_IOS(animalID,strDay,fileID,ROIs,lensMag,imagingType);
                end
            elseif strcmp(ROInames{1,c},{'Cement'}) == true || strcmp(ROInames{1,c},'SSS') == true
                % ROIS drawn free-hand for cement and SSS
                [frames] = ReadDalsaBinary_IOS(animalID,fileID);
                [ROIs] = CreateFreeHandROIs_IOS(frames{3},ROIname,animalID,ROIs);            
            end
            save([animalID '_ROIs.mat'],'ROIs');
        end
    end
end

end
