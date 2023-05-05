zap;
animalIDs = {'T224','T228','T233','T263','T264','T265','T267'};
curDir = cd;
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    dataDir = [curDir '/' animalID '/Imaging/'];
    cd(dataDir)
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    [animalID,~,~] = GetFileInfo_IOS(procDataFileIDs(1,:));

    restFile = ls('*RestData.mat');
    load(restFile)
    baselineFile = ls('*RestingBaselines.mat');
    load(baselineFile)
    % pixel-wise resting baselines
    [RestingBaselines] = CalculatePixelWiselRestingBaselines_IOS(procDataFileIDs,RestingBaselines,'y');
    % correct GCaMP attenuation
    CalculateFullSpectroscopy_IOS(procDataFileIDs,RestingBaselines)
    % re-create the RestData structure now that HbT (and/or corrected GCaMP) is available
    [RestData] = ExtractRestingData_IOS(procDataFileIDs,2);
    % create the EventData structure for CBV and neural data
    [EventData] = ExtractEventTriggeredData_IOS(procDataFileIDs);
    [RestingBaselines,ManualDecisions] = CalculateManualRestingBaselinesTimeIndeces_IOS(animalID,procDataFileIDs,RestData,RestingBaselines,'reflectance');
    % normalize RestData structures by the resting baseline
    [RestData] = NormRestDataStruct_IOS(animalID,RestData,RestingBaselines,'manualSelection');
    % normalize EventData structures by the resting baseline
    [EventData] = NormEventDataStruct_IOS(animalID,EventData,RestingBaselines,'manualSelection');

    % check and load TrainingFileDates
    trainingDatesFileStruct = dir('*_TrainingFileDates.mat');
    trainingDatesFile = {trainingDatesFileStruct.name}';
    trainingDatesFileID = char(trainingDatesFile);
    if isempty(trainingDatesFileID) == true
        [TrainingFiles] = SelectTrainingDates_IOS(procDataFileIDs);
    else
        load(trainingDatesFileID,'-mat')
    end
    % check and load TrainingFileDates
    trainingDatesFileStruct = dir('*_TrainingFileDates.mat');
    trainingDatesFile = {trainingDatesFileStruct.name}';
    trainingDatesFileID = char(trainingDatesFile);
    % sleep score an animal's data set and create a SleepData.mat structure for classification
    modelNames = {'SVM','Ensemble','Forest','Manual'};
    SleepData = [];
    % character list of all ModelData files
    modelDataFileStruct = dir('*_ModelData.mat');
    modelDataFiles = {modelDataFileStruct.name}';
    modelDataFileIDs = char(modelDataFiles);
    AddSleepParameters_IOS(procDataFileIDs,RestingBaselines,'manualSelection')
    for c = 1:length(modelNames)
        modelName = modelNames{1,c};
        [ScoringResults] = PredictBehaviorEvents_IOS(animalID,modelDataFileIDs,modelName);
        ApplySleepLogical_IOS(modelName,TrainingFiles,ScoringResults)
        NREMsleepTime = 30; % seconds
        REMsleepTime = 60; % seconds
        [SleepData] = CreateSleepData_GCaMP_IOS(NREMsleepTime,REMsleepTime,modelName,TrainingFiles,SleepData);
    end
    save([animalID '_SleepData.mat'],'SleepData')

    cd(curDir)
end