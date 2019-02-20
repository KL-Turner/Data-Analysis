function DrawROIs(animal, hem, ROInames)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs: 
%
%   Last Revised: August 15th, 2018
%________________________________________________________________________________________________________________________

windowCamFiles = ls('*_WindowCam.bin');
[~, ~, fileDates, ~] = GetFileInfo(windowCamFiles);
[uniqueDays, dayIndex, DayID] = GetUniqueDays(fileDates);
firstFileOfDay = cell(1,length(uniqueDays));
for uD = 1:length(uniqueDays)
    FileInd = DayID == uD;
    dayFilenames = windowCamFiles(FileInd,:);
    firstFileOfDay(uD) = {dayFilenames(1,:)};
end

for ii = 1:length(firstFileOfDay)
    fileID = firstFileOfDay{1, ii};
    strDay = ConvertDate(fileID);
    [frames] = ReadDalsaBinary(fileID, 256, 256);
    for iii = 1:length(ROInames)
        ROIname = ROInames{1, iii};
        mask = CreateROI(frames{2},[ROIname '_' strDay], animal, hem);
        close all;
    end
end

end