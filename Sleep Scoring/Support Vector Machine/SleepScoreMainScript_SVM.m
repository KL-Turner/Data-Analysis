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
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111'};
baselineType = 'manualSelection';
startingDirectory = cd;

%% BLOCK PURPOSE [1] Create training data set for each animal
for a = 1:size(animalIDs,2)
    % cd to the animal's bilateral imaging folder to load the baseline structure
    baselineDirectory = [animalIDs{1,a} '\Bilateral Imaging\'];
    cd(baselineDirectory)
    baselinesFileStruct = dir('*_RestingBaselines.mat');
    baselinesFile = {baselinesFileStruct.name}';
    baselinesFileID = char(baselinesFile);
    load(baselinesFileID)
    cd(startingDirectory)
    % cd to the animal's SVM training set folder
    trainingDirectory = [animalIDs{1,a} '\SVM Training Set\'];
    cd(trainingDirectory)
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    AddSleepParameters_SVM(procDataFileIDs,RestingBaselines,baselineType)
    CreateModelDataSet_SVM(procDataFileIDs)
    CreateTrainingDataSet_SVM(procDataFileIDs,RestingBaselines,baselineType)
    UpdateTrainingDataSets_SVM(procDataFileIDs)
    cd(startingDirectory)
end

%% BLOCK PURPOSE [2] Train SVM Model - cycle through each data set and update any necessary parameters
TrainModel_SVM(animalIDs);

%% BLOCK PURPOSE [3] Validate SVM Model - cycle through each data set and check model accuracy against second training set
% saveFigs = 'n';
% VerifyModelPredictions_SVM(animalIDs,driveLetters,saveFigs,baselineType)

%% BLOCK PURPOSE [4] Sleep score an animal's data set and create a SleepData.mat structure for SVM classification 
% Load SVM model, Use SVM model to sleep score new data
for a = 1:size(animalIDs,2)
    % cd to the animal's bilateral imaging folder to load the baseline structure
    load('IOS_SVM_SleepScoringModel.mat')
    % cd to the animal's SVM training set folder
    animalDirectory = [animalIDs{1,a} '\Bilateral Imaging\'];
    cd(animalDirectory)
    baselinesFileStruct = dir('*_RestingBaselines.mat');
    baselinesFile = {baselinesFileStruct.name}';
    baselinesFileID = char(baselinesFile);
    load(baselinesFileID)
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    AddSleepParameters_SVM(procDataFileIDs,RestingBaselines,baselineType)
    CreateModelDataSet_SVM(procDataFileIDs)
    modelDataFileStruct = dir('*_ModelData.mat');
    modelDataFiles = {modelDataFileStruct.name}';
    modelDataFileIDs = char(modelDataFiles);
    [SVMResults] = PredictBehaviorEvents_SVM(modelDataFileIDs,SVMModel);
    ApplySleepLogical_SVM(procDataFileIDs,SVMResults)
    sleepTime = 30;   % seconds
    [SleepData] = CreateSleepData_SVM(procDataFileIDs,sleepTime);
    cd(startingDirectory)
end
