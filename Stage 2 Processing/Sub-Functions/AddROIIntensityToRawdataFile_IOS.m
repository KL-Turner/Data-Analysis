function AddROIIntensityToRawdataFile_IOS(ROIname, rawDataFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Create an ROI from camera frames and average the reflectance from the ROI to get a timeseries. 
%            Add the result into the RawData.mat file.
%________________________________________________________________________________________________________________________
%
%   Inputs: ROIname - typically 'LH', 'RH', or an electrode ROI, and a list of the rawDataFiles in the current folder.
%
%   Outputs: None - but saves the ROI pixel intensity over time as a (1 x n) array in all listed RawData.mat files.
%
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

for fileNumber = 1:size(rawDataFiles, 1)
    fileName = rawDataFiles(fileNumber, :);
    disp(['Analyzing file ' num2str(fileNumber) ' of ' num2str(size(rawDataFiles, 1)) '...']); disp(' ')
    [animal, hem, fileDate, fileID] = GetFileInfo_IOS(fileName);
    ROIFile = ls('*ROIs.mat');
    if not(isempty(ROIFile))
        load(ROIFile)
    else
        ROIs = [];
    end
    
    strDay = ConvertDate_IOS(fileDate);
    load(fileName)
    
    [isok] = CheckROI_IOS([ROIname '_' strDay]);
    [frames] = ReadDalsaBinary_IOS([fileID '_WindowCam.bin'], RawData.Notes.CBVCamPixelHeight, RawData.Notes.CBVCamPixelWidth);
    
    if not(isok)
        mask = CreateROI_IOS(frames{2},[ROIname '_' strDay], animal, hem);
    elseif isok
        mask = GetROI_IOS(frames{1},[ROIname '_' strDay], animal, hem);
    end
    
    meanIntensity = BinToIntensity_IOS([fileID '_WindowCam.bin'], mask, frames);
    RawData.Data.CBV.(ROIname) = meanIntensity;
    save(fileName, 'RawData')
end

end

