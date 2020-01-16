%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 1) Draw ROIs for reflectance analysis
%            2) Extract the average pixel reflectance changes within those ROIs and save to RawData
%            3) Create ProcData structure using threshholds for the observed data
%            4) Use spectral analysis of the reflectance data to pull out the animal's heart rate
%            5) Regress out the pixel drift over the cement from the reflectance data
%            6) Analyze single vessel diameter changes from the IOS movies
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
[animalID,~,~] = GetFileInfo_IOS(rawDataFileIDs(1,:));

curDir = cd;
dirBreaks = strfind(curDir,'\');
curFolder = curDir(dirBreaks(end) + 1:end);
if strcmp(curFolder,'Combined Imaging') == true
    imagingType = 'bilateral';
elseif strcmp(curFolder,'Single Hemisphere') == true
    imagingType = 'single';
end

%% BLOCK PURPOSE: [3] Process the RawData structure -> Create Threshold data structure and ProcData structure.
disp('Analyzing Block [2] Create ProcData files and process analog data.'); disp(' ')
ProcessRawDataFiles_IOS(rawDataFileIDs,imagingType)

%% BLOCK PURPOSE: [4]
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
ProcessIntrinsicData_IOS(animalID,imagingType,rawDataFileIDs,procDataFileIDs)

%% BLOCK PURPOSE: [4] Add Heart Rate to the ProcData structures.
disp('Analyzing Block [4] Add heart rate to ProcData files.'); disp(' ')
ExtractHeartRate_IOS(procDataFileIDs,imagingType)

%% BLOCK PURPOSE: [5] Check/Correct pixel drift
if strcmp(imagingType,'bilateral') == true
    CorrectBilateralPixelDrift_IOS(procDataFileIDs)
elseif strcmp(imagingType,'single') == true
    CorrectPixelDrift_IOS(procDataFileIDs)
end

disp('Stage Two Processing - Complete.'); disp(' ')
