%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised: July 26th, 2019
%________________________________________________________________________________________________________________________

%% Clear workspace/Load in file names for various analysis
clear; clc; close all
disp('Loading necessary file names...'); disp(' ')

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110'};
driveLetters = {'E','E','F','F','F','D','D'};

%% BLOCK PURPOSE [1] Create training data set for a specific animal
% Character list of all '*_ProcData.mat' files
baselineType = 'manualSelection';
curDir = cd;
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

baselineDir = uigetdir;
cd(baselineDir)
% Character name for the '*_RestingBaselines.mat' structure
baselinesFileStruct = dir('*_RestingBaselines.mat');
baselinesFile = {baselinesFileStruct.name}';
baselinesFileID = char(baselinesFile);
load(baselinesFileID)
cd(curDir)
AddSleepParameters_SVM(procDataFileIDs,RestingBaselines,baselineType)
CreateModelDataSet_SVM(procDataFileIDs)
CreateTrainingDataSet_SVM(procDataFileIDs,RestingBaselines,baselineType)
UpdateTrainingDataSets_SVM(procDataFileIDs)

%% BLOCK PURPOSE [2] Train SVM Model - cycle through each data set and update any necessary parameters
baselineType = 'manualSelection';
TrainModel_SVM(animalIDs,driveLetters,baselineType);

%% BLOCK PURPOSE [3] Validate SVM Model - cycle through each data set and check model accuracy against second training set
saveFigs = 'n';
baselineType = 'manualSelection';
animalIDs = {'T111'};
driveLetters = {'D'};
VerifyModelPredictions_SVM(animalIDs,driveLetters,saveFigs,baselineType)

%% BLOCK PURPOSE [4] Sleep score an animal's data set and create a SleepData.mat structure for SVM classification 
% Load SVM model, Use SVM model to sleep score new data
curDir = cd;
modelLocation = 'C:\Users\klt8\Documents\';
cd(modelLocation)
load('SVM_SleepScoringModel.mat')
cd(curDir)

baselinesFileStruct = dir('*_RestingBaselines.mat');
baselinesFile = {baselinesFileStruct.name}';
baselinesFileID = char(baselinesFile);
load(baselinesFileID)

procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

baselineType = 'manualSelection';
AddSleepParameters_SVM(procDataFileIDs,RestingBaselines,baselineType)
CreateModelDataSet_SVM(procDataFileIDs)

modelDataFileStruct = dir('*_ModelData.mat');
modelDataFiles = {modelDataFileStruct.name}';
modelDataFileIDs = char(modelDataFiles);

[SVMResults] = PredictBehaviorEvents_SVM(modelDataFileIDs,SVMModel);
ApplySleepLogical_SVM(procDataFileIDs, SVMResults)
sleepTime = 30;   % seconds
[SleepData] = CreateSleepData_SVM(procDataFileIDs,sleepTime);
