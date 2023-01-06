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

zap; 
% character list of all RawData files in the directory from StageOneProcessing_IOS.m
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
% process the RawData structure -> Create Threshold data structure and ProcData structure
ProcessRawDataFiles_IOS(rawDataFileIDs)
% character list of all ProcData files in the directory from ProcessRawDataFiles_IOS.m
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% process pixel data from each region of interest
ProcessIntrinsicData_IOS(rawDataFileIDs,procDataFileIDs)
% add Heart Rate to the ProcData structures
ExtractHeartRate_IOS(procDataFileIDs)