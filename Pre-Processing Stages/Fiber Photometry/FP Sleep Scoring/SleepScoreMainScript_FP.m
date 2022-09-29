%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose:
%________________________________________________________________________________________________________________________

% Clear workspace/Load in file names for various analysis
clear; clc; close all
disp('Loading necessary file names...'); disp(' ')
baselineType = 'manualSelection';
startingDirectory = cd;
%% Create training data set for each animal
% cd to the animal's bilateral imaging folder to load the baseline structure
baselineDirectory = [startingDirectory '\Bilateral Imaging\'];
cd(baselineDirectory)
% load the baseline structure
baselinesFileStruct = dir('*_RestingBaselines.mat');
baselinesFile = {baselinesFileStruct.name}';
baselinesFileID = char(baselinesFile);
load(baselinesFileID)
cd(startingDirectory)
% cd to the animal's training set folde
trainingDirectory = [startingDirectory '\Training Data\'];
cd(trainingDirectory)
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% add sleep parameters (each behavior we care about during sleep)
AddSleepParameters_FP(procDataFileIDs,RestingBaselines,baselineType)
% create a table of values for sleep scoring model
CreateModelDataSet_FP(procDataFileIDs)
% create manual decisions for each 5 second bin
CreateTrainingDataSet_FP(procDataFileIDs,RestingBaselines,baselineType)
% combine the existing training set decisions with any sleep parameter changes
UpdateTrainingDataSets_FP(procDataFileIDs)
cd(startingDirectory)
% Train Models - cycle through each data set and update any necessary parameters
[animalID] = TrainSleepModels_FP;
%% Sleep score an animal's data set and create a SleepData.mat structure for classification
modelNames = {'SVM','Ensemble','Forest','Manual'};
SleepData = [];
% cd to the animal's bilateral imaging folder
cd(baselineDirectory)
% load the baseline structure
baselinesFileStruct = dir('*_RestingBaselines.mat');
baselinesFile = {baselinesFileStruct.name}';
baselinesFileID = char(baselinesFile);
load(baselinesFileID)
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% add sleep parameters (each behavior we care about during sleep)
AddSleepParameters_FP(procDataFileIDs,RestingBaselines,baselineType)
% create a table of values for sleep scoring model
CreateModelDataSet_FP(procDataFileIDs)
% character list of all ModelData files
modelDataFileStruct = dir('*_ModelData.mat');
modelDataFiles = {modelDataFileStruct.name}';
modelDataFileIDs = char(modelDataFiles);
for c = 1:length(modelNames)
    modelName = modelNames{1,c};
    [ScoringResults] = PredictBehaviorEvents_FP(animalID,startingDirectory,baselineDirectory,modelDataFileIDs,modelName);
    ApplySleepLogical_FP(startingDirectory,trainingDirectory,baselineDirectory,modelName,ScoringResults)
    NREMsleepTime = 30;   % seconds
    REMsleepTime = 60;   % seconds
    [SleepData] = CreateSleepData_FP(startingDirectory,trainingDirectory,baselineDirectory,NREMsleepTime,REMsleepTime,modelName,SleepData);
end
cd(baselineDirectory)
save([animalID '_SleepData.mat'],'SleepData')
cd(startingDirectory)

disp('Sleep Scoring analysis complete'); disp(' ')
