%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 1) Create ProcData structure using threshholds for the observed data
%            2) Extract the average pixel reflectance changes within those ROIs and save to RawData/ProcData files
%            3) Use spectral analysis of the reflectance data to pull out the animal's heart rate
%            4) Remove pixel-drift trends from each file, if necessary
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
% Identify whether this analysis is for bilateral or single hemisphere imaging
curDir = cd;
dirBreaks = strfind(curDir,'\');
curFolder = curDir(dirBreaks(end) + 1:end);
if strcmp(curFolder,'Bilateral Imaging') == true
    imagingType = 'bilateral';
elseif strcmp(curFolder,'Isoflurane Trials') == true
    imagingType = 'bilateral';
elseif strcmp(curFolder,'Single Hemisphere') == true
    imagingType = 'single';
end

%% BLOCK PURPOSE: [1] Process the RawData structure -> Create Threshold data structure and ProcData structure.
disp('Analyzing Block [2] Creating ProcData files and processing analog data.'); disp(' ')
ProcessRawDataFiles_IOS(rawDataFileIDs)

%% BLOCK PURPOSE: [2] Process IOS pixel data from each ROI.
disp('Analyzing Block [2] Proccesing IOS pixel data and ROI analysis.'); disp(' ')
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
ProcessIntrinsicData_IOS(animalID,imagingType,rawDataFileIDs,procDataFileIDs)

%% BLOCK PURPOSE: [3] Add Heart Rate to the ProcData structures.
disp('Analyzing Block [3] Adding heart rate to ProcData files.'); disp(' ')
ExtractHeartRate_IOS(procDataFileIDs,imagingType)

%% BLOCK PURPOSE: [4] Check/Correct IOS pixel drift.
disp('Analyzing Block [4] Correcting pixel drift.'); disp(' ')
if strcmp(imagingType,'bilateral') == true
    CorrectBilateralPixelDrift_IOS(procDataFileIDs)
elseif strcmp(imagingType,'single') == true
    CorrectPixelDrift_IOS(procDataFileIDs)
end

%% fin.
disp('IOS Stage Two Processing - Complete.'); disp(' ')
