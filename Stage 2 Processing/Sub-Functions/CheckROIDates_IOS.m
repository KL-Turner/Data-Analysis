function [ROIs] = CheckROIDates_IOS(animalID, ROIs, ROInames)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs: 
%
%   Last Revised: June 26th, 2019
%________________________________________________________________________________________________________________________

windowCamFilesDir = dir('*_WindowCam.bin');
windowCamDataFiles = {windowCamFilesDir.name}';
windowCamDataFileIDs = char(windowCamDataFiles);

[~, fileDates, ~] = GetFileInfo_IOS(windowCamDataFileIDs);
[uniqueDays, ~, DayID] = GetUniqueDays_IOS(fileDates);
firstsFileOfDay = cell(1,length(uniqueDays));
for a = 1:length(uniqueDays)
    FileInd = DayID == a;
    dayFilenames = windowCamDataFileIDs(FileInd,:);
    firstsFileOfDay(a) = {dayFilenames(1,:)};
end

ROIFileDir = dir('*_ROIs.mat');
ROIFileName = {ROIFileDir.name}';
ROIFileID = char(ROIFileName);
if exist(ROIFileID)
    load(ROIFileID);
else
    ROIs = [];
end

for b = 1:length(firstsFileOfDay)
    fileID = firstsFileOfDay{1, b};
    strDay = ConvertDate_IOS(fileID);
    checkROI = [ROInames{1,1} '_' strDay];
    if ~isfield(ROIs, (checkROI))
        [frames] = ReadDalsaBinary_IOS(fileID, 256, 256);
        for c = 1:length(ROInames)
            ROIname = ROInames{1, c};
            [ROIs] = CreateBilateralROIs_IOS(frames{2},[ROIname '_' strDay], animalID, ROIs);
        end
        save([animalID '_ROIs.mat'], 'ROIs');
    else
        disp(['ROIs for ' strDay ' already exist. Continuing...']); disp(' ')
    end
end

for b = 1:length(firstsFileOfDay)
    fileID = firstsFileOfDay{1, b};
    strDay = ConvertDate_IOS(fileID);
    [frames] = ReadDalsaBinary_IOS(fileID, 256, 256);
    ROIname = 'Cement';
    [ROIs] = CreateBilateralROIs_IOS(frames{2},[ROIname '_' strDay], animalID, ROIs);
    save([animalID '_ROIs.mat'], 'ROIs');
end

end