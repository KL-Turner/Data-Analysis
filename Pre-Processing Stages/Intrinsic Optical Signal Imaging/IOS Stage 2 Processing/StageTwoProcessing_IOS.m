%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: - Set resting and whisking thresholds using whisker acceleration and force sensor
%          - Downsample and filter many of the data types in 'RawData' creating new 'ProcData' structures
%          - Place ROIs for IOS reflectance and/or fluorescence changes depending on wavelength(s)
%          - Extract the animal's heart rate if green or lime reflectance data is at least 30 Hz
%________________________________________________________________________________________________________________________

% clear the workspace, variables, and command window.
zap;
% character list of all RawData files in the directory from StageOneProcessing_IOS.m
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
% animal ID
[animalID,~,~] = GetFileInfo_IOS(rawDataFileIDs(1,:));
% select imaging type
imagingOptions = {'Single ROI (SI)','Single ROI (SS)','Bilateral ROI (SI)','Bilateral ROI (SI,FC)'};
imagingType = SelectImagingType_IOS(imagingOptions);
% select imaging type
wavelengthOptions = {'1 wavelength (G/L)','1 wavelength (B)','2 wavelengths (G/L & B)','3 wavelengths (R, G/L, & B)'};
imagingColors = SelectWavelengthType_IOS(wavelengthOptions);
% process the RawData structure -> Create Threshold data structure and ProcData structure.
ProcessRawDataFiles_IOS(rawDataFileIDs)
% character list of all ProcData files in the directory from ProcessRawDataFiles_IOS.m
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% process pixel data from each region of interest
ProcessIntrinsicData_IOS(animalID,imagingType,imagingColors,rawDataFileIDs,procDataFileIDs)
% add Heart Rate to the ProcData structures.
ExtractHeartRate_IOS(procDataFileIDs,imagingType)
