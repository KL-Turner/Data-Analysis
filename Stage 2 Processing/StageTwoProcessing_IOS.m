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
rawDataFileIDs = char(rawDataFiles);
[animalID, ~, ~] = GetFileInfo_IOS(rawDataFileIDs(1,:));

%% BLOCK PURPOSE: [1] Create bilateral regions of interest for the windows
disp('Analyzing Block [1] Creating bilateral regions of interest.'); disp(' ')
ROInames = {'LH', 'RH', 'LH_Electrode', 'RH_Electrode', 'Cement'};
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
ExtractCBVData_IOS(ROIs, ROInames, rawDataFileIDs)

%% BLOCK PURPOSE: [3] Process the RawData structure -> Create Threshold data structure and ProcData structure.
disp('Analyzing Block [2] Create ProcData files and process analog data.'); disp(' ')
ProcessRawDataFiles_IOS(rawDataFileIDs)

%% BLOCK PURPOSE: [4] Add Heart Rate to the ProcData structures.
disp('Analyzing Block [4] Add heart rate to ProcData files.'); disp(' ')
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
ExtractHeartRate_IOS(procDataFileIDs)

%% BLOCK PURPOSE: [5] Check/Correct pixel drift 
CheckPixelDrift(procDataFileIDs)
CorrectPixelDrift(procDataFileIDs)

disp('Stage Two Processing - Complete.'); disp(' ')
