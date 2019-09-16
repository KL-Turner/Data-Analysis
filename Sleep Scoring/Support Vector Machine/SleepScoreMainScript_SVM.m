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
% Character list of all '*_ProcData.mat' files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

% Character list of all '*_TrainingData.mat' files
trainingDataFileStruct = dir('*_TrainingData.mat');
trainingDataFiles = {trainingDataFileStruct.name}';
trainingDataFileIDs = char(trainingDataFiles);

% Character list of all '*_ModelData.mat' files
modelDataFileStruct = dir('*_ModelData.mat');
modelDataFiles = {modelDataFileStruct.name}';
modelDataFileIDs = char(modelDataFiles);

% Character name for the '*_RestingBaselines.mat' structure
baselinesFileStruct = dir('*_RestingBaselines.mat');
baselinesFile = {baselinesFileStruct.name}';
baselinesFileID = char(baselinesFile);
load(baselinesFileID)

modelFileStruct = dir('*_SleepScoringModel.mat');
modelFileName = {modelFileStruct.name}';
modelFile = char(modelFileName);
load(modelFile)

%% Add a 'sleep' folder with a 'parameters' field to each '*_ProcData.mat' file in the directory. 
% This needs to be run first for all data related to sleep scoring
AddSleepParameters_SVM(procDataFileIDs, RestingBaselines)

%% Create a '*_ModelData.mat' file for each '*_ProcData.mat' file in the directory.
% This needs to be run for all testing/unseen data run through the SVM model
CreateModelDataSet_SVM(procDataFileIDs)

%% Update all '*_TrainingData.mat' files with potentially updated ModelData sets
UpdateTrainingDataSets(procDataFileIDs)

%% Create '*_TrainingData.mat' files for each selected '*_ProcData.mat' file (uigetfile).
% Run this to create manually scored behavioral states for SVM model training
CreateTrainingDataSet_SVM(procDataFileIDs,RestingBaselines)

%% Create SVM Model using manually-scored training data from '*_TrainingData.mat' files.
% Saves model to directory containing the training data
TrainModel_SVM(trainingDataFileIDs);
PredictBehaviorEvents_SVM(modelDataFileIDs, SVMModel)
CompareTrainingSet(modelDataFileIDs, RestingBaselines)

%% Load SVM model, Use SVM model to sleep score new data
disp('Select file location of the support vector machine classifier'); disp(' ')
modelLocation = uigetdir;
curDir = cd;
cd(modelLocation)
modelFileStruct = dir('*_SleepScoringModel.mat');
modelFileName = {modelFileStruct.name}';
modelFile = char(modelFileName);
load(modelFile)
cd(curDir)

PredictBehaviorEvents_SVM(modelDataFileIDs, SVMModel)

%% Check model predictions
VerifyModelPredictions(modelDataFileIDs, RestingBaselines)

