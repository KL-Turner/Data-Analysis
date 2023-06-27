function [] = SleepScoreMainScript_IOS()
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
zap;
[TrainingFiles,procDataFileIDs] = SelectTrainingDates_IOS();
% add sleep parameters (each behavior we care about during sleep)
AddSleepParameters_IOS(procDataFileIDs)
% create a table of values for sleep scoring model
CreateModelDataSet_IOS(procDataFileIDs,imagingType)
% create manual decisions for each 5 second bin
CreateTrainingDataSet_IOS(procDataFileIDs,RestingBaselines,baselineType,imagingType,TrainingFiles)
% train Models - cycle through each data set and update any necessary parameters
TrainSleepModels_IOS();
% sleep score data set and create structure for classification
modelNames = {'SVM','Ensemble','Forest','Manual'}; SleepData = [];
for aa = 1:length(modelNames)
    modelName = modelNames{1,aa};
    [ScoringResults] = PredictBehaviorEvents_IOS(modelName);
    ApplySleepLogical_IOS(modelName,TrainingFiles,ScoringResults)
    [SleepData] = CreateSleepData_IOS(modelName,TrainingFiles,SleepData);
end