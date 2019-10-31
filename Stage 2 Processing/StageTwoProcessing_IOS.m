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

%% BLOCK PURPOSE: [1] Create regions of interest for the windows
disp('Analyzing Block [1] Creating regions of interest for reflectance data.'); disp(' ')
imagingType = input('Imaging Type: (bilateral/single): ','s'); disp(' ')
if strcmp(imagingType,'bilateral') == true
    ROInames = {'LH','RH','LH_Cement','RH_Cement','Cement'};
elseif strcmp(imagingType,'single') == true
    ROInames = {'Barrels','Cement'};
end
ROIFileDir = dir('*_ROIs.mat');
if isempty(ROIFileDir) == true
    ROIs = [];
else
    ROIFileName = {ROIFileDir.name}';
    ROIFileID = char(ROIFileName);
    load(ROIFileID);
end
[ROIs] = CheckROIDates_IOS(animalID,ROIs,ROInames);

%% BLOCK PURPOSE: [2] Extract CBV data from each ROI for each RawData file in the directory that hasn't been processed yet.
disp('Analyzing Block [2] Extracting mean reflectance data from each ROI.'); disp(' ')
ExtractCBVData_IOS(ROIs,ROInames,rawDataFileIDs)

%% BLOCK PURPOSE: [3] Process the RawData structure -> Create Threshold data structure and ProcData structure.
disp('Analyzing Block [2] Create ProcData files and process analog data.'); disp(' ')
ProcessRawDataFiles_IOS(rawDataFileIDs)

%% BLOCK PURPOSE: [4] Add Heart Rate to the ProcData structures.
disp('Analyzing Block [4] Add heart rate to ProcData files.'); disp(' ')
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
ExtractHeartRate_IOS(procDataFileIDs,imagingType)

%% BLOCK PURPOSE: [5] Check/Correct pixel drift
if strcmp(imagingType,'bilateral') == true
    CorrectBilateralPixelDrift_IOS(procDataFileIDs)
elseif strcmp(imagingType,'single') == true
    CorrectPixelDrift_IOS(procDataFileIDs)
end

%% BLOCK PURPOSE: [6] IOS vessel diameter analysis
% rawDataFileIDs = rawDataFileIDs(1,:);
% if strcmp(imagingType,'single') == true
%     CreateIntrinsicTiffStacks_IOS(rawDataFileIDs)
%     DiamCalcSurfaceVessel_IOS(rawDataFileIDs)
%     ExtractTiffData_IOS(rawDataFileIDs)
% end

disp('Stage Two Processing - Complete.'); disp(' ')
