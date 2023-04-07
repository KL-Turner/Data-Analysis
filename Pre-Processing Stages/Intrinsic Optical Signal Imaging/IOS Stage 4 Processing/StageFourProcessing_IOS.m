function [] = StageFourProcessing_IOS()
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Sleep score manual training files, training classification models, and create sleep data structure
%----------------------------------------------------------------------------------------------------------
zap;
disp('Loading necessary file names...'); disp(' ')
% load the baseline structure
baselinesFileStruct = dir('*_RestingBaselines.mat');
baselinesFile = {baselinesFileStruct.name}';
baselinesFileID = char(baselinesFile);
load(baselinesFileID)
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% check and load TrainingFileDates
trainingDatesFileStruct = dir('*_TrainingFileDates.mat');
trainingDatesFile = {trainingDatesFileStruct.name}';
trainingDatesFileID = char(trainingDatesFile);
if isempty(trainingDatesFileID) == true
    [TrainingFiles] = SelectTrainingDates_IOS(procDataFileIDs);
else
    load(trainingDatesFileID,'-mat')
end
% animal ID
[animalID,~,~] = GetFileInfo_IOS(procDataFileIDs(1,:));
% add sleep parameters (each behavior we care about during sleep)
AddSleepParameters_IOS(procDataFileIDs,RestingBaselines,'manualSelection')
% % create a table of values for sleep scoring model
CreateModelDataSet_IOS(procDataFileIDs,imagingType)
% % create manual decisions for each 5 second bin
CreateTrainingDataSet_IOS(procDataFileIDs,RestingBaselines,'manualSelection',TrainingFiles)
% % combine the existing training set decisions with any sleep parameter changes
UpdateTrainingDataSets_IOS(procDataFileIDs,TrainingFiles)
% % train Models - cycle through each data set and update any necessary parameters
TrainSleepModels_IOS();
% sleep score an animal's data set and create a SleepData.mat structure for classification
modelNames = {'SVM','Ensemble','Forest','Manual'};
SleepData = [];
% character list of all ModelData files
modelDataFileStruct = dir('*_ModelData.mat');
modelDataFiles = {modelDataFileStruct.name}';
modelDataFileIDs = char(modelDataFiles);
for c = 1:length(modelNames)
    modelName = modelNames{1,c};
    AddSleepParameters_IOS(procDataFileIDs,RestingBaselines,'manualSelection')
    [ScoringResults] = PredictBehaviorEvents_IOS(animalID,modelDataFileIDs,modelName);
    ApplySleepLogical_IOS(modelName,TrainingFiles,ScoringResults)
    NREMsleepTime = 30; % seconds
    REMsleepTime = 60; % seconds
    if strcmpi(imagingType,'GCaMP') == true
        [SleepData] = CreateSleepData_GCaMP_IOS(NREMsleepTime,REMsleepTime,modelName,TrainingFiles,SleepData);
    else
        [SleepData] = CreateSleepData_IOS(NREMsleepTime,REMsleepTime,modelName,TrainingFiles,SleepData);
    end
end
save([animalID '_SleepData.mat'],'SleepData')
disp('Sleep Scoring analysis complete'); disp(' ')
