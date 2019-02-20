function [] = AddROIIntensityToRawdataFile(ROIname, rawDataFiles)
%___________________________________________________________________________________________________
% Edited by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%
%   Last Revised: August 8th, 2018
%___________________________________________________________________________________________________
%
% Author: Aaron Winder
% Affiliation: Engineering Science and Mechanics, Penn State University
% https://github.com/awinde
%
% DESCRIPTION: Create an ROI from camera frames and average the
% reflectance from the ROI to get a timeseries. Add the result into the
% *rawdata.mat file.
%   
%_______________________________________________________________
% PARAMETERS:             
%             ROIname - [string] Description of ROI                        
%_______________________________________________________________
% RETURN:                     
%                               
%_______________________________________________________________

for fileNumber = 1:size(rawDataFiles, 1)
    fileName = rawDataFiles(fileNumber, :);
    disp(['Analyzing file ' num2str(fileNumber) ' of ' num2str(size(rawDataFiles, 1)) '...']); disp(' ')
    [animal, hem, fileDate, fileID] = GetFileInfo(fileName);
    ROIFile = ls('*ROIs.mat');
    if not(isempty(ROIFile))
        load(ROIFile)
    else
        ROIs = [];
    end
    
    strDay = ConvertDate(fileDate);
    load(fileName)
    
    [isok] = CheckROI([ROIname '_' strDay]);
    [frames] = ReadDalsaBinary([fileID '_WindowCam.bin'], RawData.Notes.CBVCamPixelHeight, RawData.Notes.CBVCamPixelWidth);
    
    if not(isok)
        mask = CreateROI(frames{2},[ROIname '_' strDay], animal, hem);
    elseif isok
        mask = GetROI(frames{1},[ROIname '_' strDay], animal, hem);
    end
    
    meanIntensity = BinToIntensity([fileID '_WindowCam.bin'], mask, frames);
    RawData.Data.CBV.(ROIname) = meanIntensity;
    save(fileName, 'RawData')
end

end

