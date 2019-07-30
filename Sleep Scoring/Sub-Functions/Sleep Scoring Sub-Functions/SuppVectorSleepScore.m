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

%% Character list of all ProcData files
clear; clc; close all
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

% Character name for the RestingBaselines structure
baselinesFileStruct = dir('*_RestingBaselines.mat');
baselinesFile = {baselinesFileStruct.name}';
baselinesFileID = char(baselinesFile);
load(baselinesFileID)

trainModel = input('Train SVM model? (y/n): ','s'); disp(' ')
createTrainingSet = input('Add additional data to SVM training set? (y/n): ','s'); disp(' ')
if strcmp(trainModel,'n') == true
    disp('Select trained model location:'); disp(' ')
    modelLocation = uigetdir;
    disp(['Loading model from: ' modelLocation]); disp(' ')
    AddSleepParameters_SVM(procDataFileIDs, RestingBaselines)
    CreateModelDataSet_SVM(procDataFileIDs)
else
    AddSleepParameters_SVM(procDataFileIDs, RestingBaselines)
    CreateModelDataSet_SVM(procDataFileIDs)
    if createTrainingSet == true
        CreateTrainingDataSet_SVM(procDataFileIDs, RestingBaselines)
    end
    trainingTableFileStruct = dir('*_TrainingTable.mat');
    trainingTableFiles = {trainingTableFileStruct.name}';
    trainingTableFileIDs = char(trainingTableFiles);
    [SVMModel] = TrainModel_SVM(trainingTableFileIDs);
end

%%
if strcmp(trainModel,'y') == true
    PredictBehaviorEvents_SVM(procDataFileIDs, SVMModel)
else
    curDir = cd;
    cd(modelLocation)
    modelFileStruct = dir('*_TrainingTable.mat');
    modelFileName = {modelFileStruct.name}';
    modelFile = char(modelFileName);
    load(modelFile)
    cd(curDir)
    modelDataFileStruct = dir('*_ModelData.mat');
    modelDataFiles = {modelDataFileStruct.name}';
    modelDataFileIDs = char(modelDataFiles);
    PredictBehaviorEvents_SVM(modelDataFileIDs, SVMModel)
end

