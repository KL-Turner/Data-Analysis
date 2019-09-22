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

animalIDs = {'T101', 'T102', 'T103', 'T105', 'T108', 'T109'};
driveLetters = {'E', 'E', 'E', 'F', 'F', 'F'};

%% BLOCK PURPOSE [1] Create training data set for a specific animal
% Character list of all '*_ProcData.mat' files
curDir = cd;
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

cd ..
% Character name for the '*_RestingBaselines.mat' structure
baselinesFileStruct = dir('*_RestingBaselines.mat');
baselinesFile = {baselinesFileStruct.name}';
baselinesFileID = char(baselinesFile);
load(baselinesFileID)
cd(curDir)
AddSleepParameters_SVM(procDataFileIDs, RestingBaselines)
CreateModelDataSet_SVM(procDataFileIDs)
CreateTrainingDataSet_SVM(procDataFileIDs, RestingBaselines)
UpdateTrainingDataSets(procDataFileIDs)

%% BLOCK PURPOSE [2] Train SVM Model - cycle through each data set and update any necessary parameters
TrainModel_SVM(animalIDs, driveLetters);

%% BLOCK PURPOSE [3] Validate SVM Model - cycle through each data set and check model accuracy against second training set
saveFigs = 'y';
VerifyModelPredictions(animalIDs, driveLetters, saveFigs)

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
AddSleepParameters_SVM(procDataFileIDs, RestingBaselines, baselineType)
CreateModelDataSet_SVM(procDataFileIDs)

modelDataFileStruct = dir('*_ModelData.mat');
modelDataFiles = {modelDataFileStruct.name}';
modelDataFileIDs = char(modelDataFiles);

[SVMResults] = PredictBehaviorEvents_SVM(modelDataFileIDs, SVMModel);
ApplySleepLogical_SVM(procDataFileIDs, SVMResults)
sleepTime = 60;   % seconds
[SleepData] = CreateSleepData_SVM(procDataFileIDs, sleepTime);
