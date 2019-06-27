%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 1) Generate bilateral ROIs for CBV analysis.
%            2) Create ProcData structure using threshholds for the observed data.
%________________________________________________________________________________________________________________________
%
%   Inputs: 1) Selects all _RawData files from all days. Draw all left ROIs first, then the right.
%           2) Selects all _RawData files from all days. Follow the command window prompts for each threshold value.
%
%   Outputs: 1) An animal_hemisphere_ROIs.mat file with the xi, yi coordinates for each day.
%            2) A ProcData.mat structure for each inputed rawdata.mat file, as well as a Thresholds.mat file that serves
%               as a record for each variable/day's set threshold.        
%
%   Last Revised: June 26th, 2019    
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
clc;
clear;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')

% Character list of all RawData files
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFiles = char(rawDataFiles);
[animalID, ~, ~] = GetFileInfo_IOS(rawDataFiles(1,:));

%% BLOCK PURPOSE: [1] Create bilateral regions of interest for the windows
disp('Analyzing Block [1] Creating bilateral regions of interest.'); disp(' ')
ROInames = {'LH', 'RH', 'LH_Electrode', 'RH_Electrode'};
ROIFileDir = dir('*_ROIs.mat');
if isempty(ROIFileDir) == true
    ROIs = [];
else
    ROIFileName = {ROIFileDir.name}';
    ROIFileID = char(ROIFileName);
    load(ROIFileID);
end
[ROIs] = CheckROIDates_IOS(animalID, ROIs, ROInames);

%% BLOCK PURPOSE: [2] Extract CBV data from each ROI for each RawData file in the directory that hasn't been processed yet.
disp('Analyzing Block [2] Extracting cerebral blood volume data from each ROI.'); disp(' ')
ExtractCBVData_IOS(ROIs, ROInames, rawDataFiles)

%% BLOCK PURPOSE: [2] Process the RawData structure -> Create Threshold data structure and ProcData structure.
disp('Analyzing Block [2] Create ProcData files and analyze neural data.'); disp(' ')
for fileNumber = 1:size(rawDataFiles, 1)
    fileName = rawDataFiles(fileNumber, :);
    disp(['Analyzing file ' num2str(fileNumber) ' of ' num2str(size(rawDataFiles, 1)) '...']); disp(' ')
    ProcessRawDataFile_IOS(fileName)
    close all;
end

%% BLOCK PURPOSE: [3] Add Heart Rate to the ProcData structures.
disp('Analyzing Block [3] Add heart rate to ProcData files.'); disp(' ')
for fileNumber = 1:size(procDataFiles, 1)
    fileName = procDataFiles(fileNumber, :);
    disp(['Adding the heart rate to ProcData file ' num2str(fileNumber) ' of ' num2str(size(procDataFiles, 1)) '...']); disp(' ')
    AddHeartRate(fileName)
end

disp('Stage Two Processing - Complete.'); disp(' ')
