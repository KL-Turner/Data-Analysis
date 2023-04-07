zap
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
imagingOptions = {'Single ROI (SI)','Single ROI (SSS)','Bilateral ROI (SI)','Bilateral ROI (SI,FC)'};
imagingType = SelectImagingType_IOS(imagingOptions);
% select imaging type
wavelengthOptions = {'Green','Lime','Blue','Green & Blue','Lime & Blue','Red, Green, & Blue','Red, Lime, & Blue'};
imagingWavelengths = SelectImagingType_IOS(wavelengthOptions);
% edit fieldname
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID);
    disp(num2str(aa))
    ProcData.notes.imagingType = imagingType;
    ProcData.notes.imagingWavelengths = imagingWavelengths;
    try
        ProcData.data.HbT = ProcData.data.CBV_HbT;
        ProcData.data = rmfield(ProcData.data,'CBV_HbT');
        ProcData.sleep.parameters = rmfield(ProcData.sleep.parameters,'CBV');
    catch
    end
    save(procDataFileID,'ProcData')
end
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
% check and load TrainingFileDates
trainingDatesFileStruct = dir('*_TrainingFileDates.mat');
trainingDatesFile = {trainingDatesFileStruct.name}';
trainingDatesFileID = char(trainingDatesFile);
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
    [SleepData] = CreateSleepData_IOS(NREMsleepTime,REMsleepTime,modelName,TrainingFiles,SleepData);
end
save([animalID '_SleepData.mat'],'SleepData')
disp('Sleep Scoring analysis complete'); disp(' ')