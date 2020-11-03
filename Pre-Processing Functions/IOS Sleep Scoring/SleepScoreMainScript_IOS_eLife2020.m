%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose:
%________________________________________________________________________________________________________________________

%% Clear workspace/Load in file names for various analysis
clear; clc; close all
disp('Loading necessary file names...'); disp(' ')
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
baselineType = 'manualSelection';
startingDirectory = cd;

%% BLOCK PURPOSE [1] Create training data set for each animal
for aa = 1:size(animalIDs,2)
    % cd to the animal's bilateral imaging folder to load the baseline structure
    baselineDirectory = [animalIDs{1,aa} '\Bilateral Imaging\'];
    cd(baselineDirectory)
    % load the baseline structure
    baselinesFileStruct = dir('*_RestingBaselines.mat');
    baselinesFile = {baselinesFileStruct.name}';
    baselinesFileID = char(baselinesFile);
    load(baselinesFileID)
    cd(startingDirectory)
    % cd to the animal's training set folde
    trainingDirectory = [animalIDs{1,aa} '\Training Data\'];
    cd(trainingDirectory)
    % character list of all ProcData files
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    % add sleep parameters (each behavior we care about during sleep)
    AddSleepParameters_IOS_eLife2020(procDataFileIDs,RestingBaselines,baselineType)
    % create a table of values for sleep scoring model
    CreateModelDataSet_IOS_eLife2020(procDataFileIDs)
    % create manual decisions for each 5 second bin
    CreateTrainingDataSet_IOS_eLife2020(procDataFileIDs,RestingBaselines,baselineType)
    % combine the existing training set decisions with any sleep parameter changes
    UpdateTrainingDataSets_IOS_eLife2020(procDataFileIDs)
    cd(startingDirectory)
end

%% BLOCK PURPOSE [2] Train Models - cycle through each data set and update any necessary parameters
TrainSleepModels_IOS_eLife2020(animalIDs);
ComparePredictionAccuracy_IOS_eLife2020

%% BLOCK PURPOSE [3] Sleep score an animal's data set and create a SleepData.mat structure for classification 
modelNames = {'SVM','Ensemble','Forest','Manual'};
for bb = 1:size(animalIDs,2)
    SleepData = [];
    % cd to the animal's bilateral imaging folder
    animalDirectory = [animalIDs{1,bb} '\Bilateral Imaging\'];
    trainingDirectory = [animalIDs{1,bb} '\Training Data\'];
    cd(animalDirectory)
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
    AddSleepParameters_IOS_eLife2020(procDataFileIDs,RestingBaselines,baselineType)
    % create a table of values for sleep scoring model
    CreateModelDataSet_IOS_eLife2020(procDataFileIDs)
    % character list of all ModelData files
    modelDataFileStruct = dir('*_ModelData.mat');
    modelDataFiles = {modelDataFileStruct.name}';
    modelDataFileIDs = char(modelDataFiles);
    for c = 1:length(modelNames)
        modelName = modelNames{1,c};
        [ScoringResults] = PredictBehaviorEvents_IOS_eLife2020(animalIDs{1,bb},startingDirectory,animalDirectory,modelDataFileIDs,modelName);
        ApplySleepLogical_IOS_eLife2020(startingDirectory,trainingDirectory,animalDirectory,modelName,ScoringResults)
        NREMsleepTime = 30;   % seconds
        REMsleepTime = 60;   % seconds
        [SleepData] = CreateSleepData_IOS_eLife2020(startingDirectory,trainingDirectory,animalDirectory,NREMsleepTime,REMsleepTime,modelName,SleepData);
    end
    save([animalIDs{1,bb} '_SleepData.mat'],'SleepData')
    cd(startingDirectory)
end
