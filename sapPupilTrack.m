zap;
animalIDs = {'T240','T241','T242','T243','T259','T260','T261','T262'};
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
    % normalize RestData structures by the resting baseline
    [RestData] = NormRestDataStruct_IOS(RestData,RestingBaselines,'manualSelection');
    % normalize EventData structures by the resting baseline
    [EventData] = NormEventDataStruct_IOS(EventData,RestingBaselines,'manualSelection');

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

    cd(curDir)
end