function [Results_BilatCoher_Ephys] = AnalyzeBilateralCoherence_Ephys(animalID,group,set,rootFolder,delim,Results_BilatCoher_Ephys)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% only run analysis for valid animal IDs
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Bilateral Imaging'];
cd(dataLocation)
% character list of all ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID,'-mat')
% find and load manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID,'-mat')
% find and load RestingBaselines.mat struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load SleepData.mat struct
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
% find and load Forest_ScoringResults.mat struct
forestScoringResultsFileID = [animalID '_Forest_ScoringResults.mat'];
load(forestScoringResultsFileID,'-mat')
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
% go through each valid data type for arousal based coherence analysis
dataTypes = {'HbT','gammaBandPower'};
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    %% Rest
    clear LH_finalRestData RH_finalRestData LH_ProcRestData RH_ProcRestData LH_restData RH_restData
    if strcmp(dataType,'HbT') == true
        samplingRate = RestData.(dataType).LH.CBVCamSamplingRate;
        [restLogical] = FilterEvents_IOS(RestData.(dataType).LH,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.(dataType).LH,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.(dataType).LH.fileIDs(combRestLogical,:);
        restEventTimes = RestData.(dataType).LH.eventTimes(combRestLogical,:);
        restDurations = RestData.(dataType).LH.durations(combRestLogical,:);
        LH_RestingData = RestData.(dataType).LH.data(combRestLogical,:);
        RH_RestingData = RestData.(dataType).RH.data(combRestLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [LH_finalRestData,~,~,~] = RemoveInvalidData_IOS(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    else
        samplingRate = RestData.cortical_LH.(dataType).CBVCamSamplingRate;
        [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
        restEventTimes = RestData.cortical_LH.(dataType).eventTimes(combRestLogical,:);
        restDurations = RestData.cortical_LH.(dataType).durations(combRestLogical,:);
        LH_RestingData = RestData.cortical_LH.(dataType).NormData(combRestLogical,:);
        RH_RestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [LH_finalRestData,~,~,~] = RemoveInvalidData_IOS(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    end
    % detrend and truncate data to minimum length to match events
    for bb = 1:length(LH_finalRestData)
        if length(LH_finalRestData{bb,1}) < params.minTime.Rest*samplingRate
            restChunkSampleDiff = params.minTime.Rest*samplingRate - length(LH_finalRestData{bb,1});
            LH_restPad = (ones(1,restChunkSampleDiff))*LH_finalRestData{bb,1}(end);
            RH_restPad = (ones(1,restChunkSampleDiff))*RH_finalRestData{bb,1}(end);
            LH_ProcRestData{bb,1} = horzcat(LH_finalRestData{bb,1},LH_restPad);
            RH_ProcRestData{bb,1} = horzcat(RH_finalRestData{bb,1},RH_restPad);
            LH_ProcRestData{bb,1} = detrend(LH_ProcRestData{bb,1},'constant');
            RH_ProcRestData{bb,1} = detrend(RH_ProcRestData{bb,1},'constant');
        else
            LH_ProcRestData{bb,1} = detrend(LH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            RH_ProcRestData{bb,1} = detrend(RH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
        end
    end
    % pre-allocate coherence matrix
    LH_restData = zeros(length(LH_ProcRestData{1,1}),length(LH_ProcRestData));
    RH_restData = zeros(length(RH_ProcRestData{1,1}),length(RH_ProcRestData));
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for cc = 1:length(LH_ProcRestData)
        LH_restData(:,cc) = LH_ProcRestData{cc,1};
        RH_restData(:,cc) = RH_ProcRestData{cc,1};
    end
    % parameters for coherencyc - information available in function
    params.tapers = [1,1]; % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;
    params.fpass = [0,1]; % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2,0.05];
    % calculate the coherence between desired signals
    [C_RestData,~,~,~,~,f_RestData,confC_RestData,~,cErr_RestData] = coherencyc(LH_restData,RH_restData,params);
    % save results
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).Rest.C = C_RestData;
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).Rest.f = f_RestData;
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).Rest.confC = confC_RestData;
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).Rest.cErr = cErr_RestData;
    %% Alert
    clear LH_AlertData LH_ProcAlertData LH_alertData RH_AlertData RH_ProcAlertData RH_alertData scoringLabels
    zz = 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144 % 36 bins (180 total) or 3 minutes of asleep
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                if strcmp(dataType,'HbT') == true
                    LH_AlertData{zz,1} = ProcData.data.(dataType).LH;
                    RH_AlertData{zz,1} = ProcData.data.(dataType).RH;
                    zz = zz + 1;
                else
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        LH_AlertData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                        RH_AlertData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                        zz = zz + 1;
                    end
                end
            end
        end
    end
    % detrend data
    if isempty(LH_AlertData) == false
        for bb = 1:length(LH_AlertData)
            LH_ProcAlertData{bb,1} = detrend(LH_AlertData{bb,1},'constant');
            RH_ProcAlertData{bb,1} = detrend(RH_AlertData{bb,1},'constant');
        end
        % preallocate coherence matrix
        LH_alertData = zeros(length(LH_ProcAlertData{1,1}),length(LH_ProcAlertData));
        RH_alertData = zeros(length(RH_ProcAlertData{1,1}),length(RH_ProcAlertData));
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        for cc = 1:length(LH_ProcAlertData)
            LH_alertData(:,cc) = LH_ProcAlertData{cc,1};
            RH_alertData(:,cc) = RH_ProcAlertData{cc,1};
        end
        % calculate the coherence between desired signals
        params.tapers = [10,19]; % Tapers [n, 2n - 1]
        [C_AlertData,~,~,~,~,f_AlertData,confC_AlertData,~,cErr_AlertData] = coherencyc(LH_alertData,RH_alertData,params);
        % save results
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Alert.C = C_AlertData;
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Alert.f = f_AlertData;
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Alert.confC = confC_AlertData;
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Alert.cErr = cErr_AlertData;
    else
        % save results
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Alert.C = [];
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Alert.f = [];
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Alert.confC = [];
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Alert.cErr = [];
    end
    %% analyze bilateral coherence during periods of aasleep
    clear LH_AsleepData LH_ProcAsleepData LH_asleepData RH_AsleepData RH_ProcAsleepData RH_asleepData scoringLabels
    zz = 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) < 36 % 36 bins (180 total) or 3 minutes of alert
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                if strcmp(dataType,'HbT') == true || strcmp(dataType,'Deoxy') == true
                    LH_AsleepData{zz,1} = ProcData.data.(dataType).LH;
                    RH_AsleepData{zz,1} = ProcData.data.(dataType).RH;
                    zz = zz + 1;
                else
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        LH_AsleepData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                        RH_AsleepData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                        zz = zz + 1;
                    end
                end
            end
        end
    end
    % detrend data
    if isempty(LH_AsleepData) == false
        for bb = 1:length(LH_AsleepData)
            LH_ProcAsleepData{bb,1} = detrend(LH_AsleepData{bb,1},'constant');
            RH_ProcAsleepData{bb,1} = detrend(RH_AsleepData{bb,1},'constant');
        end
        % preallocate coherence matrix
        LH_asleepData = zeros(length(LH_ProcAsleepData{1,1}),length(LH_ProcAsleepData));
        RH_asleepData = zeros(length(RH_ProcAsleepData{1,1}),length(RH_ProcAsleepData));
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        for cc = 1:length(LH_ProcAsleepData)
            LH_asleepData(:,cc) = LH_ProcAsleepData{cc,1};
            RH_asleepData(:,cc) = RH_ProcAsleepData{cc,1};
        end
        % calculate the coherence between desired signals
        params.tapers = [10,19]; % Tapers [n, 2n - 1]
        [C_AsleepData,~,~,~,~,f_AsleepData,confC_AsleepData,~,cErr_AsleepData] = coherencyc(LH_asleepData,RH_asleepData,params);
        % save results
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Asleep.C = C_AsleepData;
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Asleep.f = f_AsleepData;
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Asleep.confC = confC_AsleepData;
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Asleep.cErr = cErr_AsleepData;
    else
        % save results
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Asleep.C = [];
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Asleep.f = [];
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Asleep.confC = [];
        Results_BilatCoher_Ephys.(group).(animalID).(dataType).Asleep.cErr = [];
    end
    %% All
    clear LH_AllData LH_ProcAllData LH_allData RH_AllData RH_ProcAllData RH_allData
    zz = 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,~] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        load(procDataFileID,'-mat')
        puffs = ProcData.data.stimulations.LPadSol;
        % don't include trials with stimulation
        if isempty(puffs) == true
            if strcmp(dataType,'HbT') == true
                LH_AllData{zz,1} = ProcData.data.(dataType).LH;
                RH_AllData{zz,1} = ProcData.data.(dataType).RH;
                zz = zz + 1;
            else
                motionArtifact = ProcData.notes.motionArtifact;
                if motionArtifact == false
                    LH_AllData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                    RH_AllData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                    zz = zz + 1;
                end
            end
        end
    end
    % detrend data
    for bb = 1:length(LH_AllData)
        LH_ProcAllData{bb,1} = detrend(LH_AllData{bb,1},'constant');
        RH_ProcAllData{bb,1} = detrend(RH_AllData{bb,1},'constant');
    end
    % preallocate coherence matrix
    LH_allData = zeros(length(LH_ProcAllData{1,1}),length(LH_ProcAllData));
    RH_allData = zeros(length(RH_ProcAllData{1,1}),length(RH_ProcAllData));
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for cc = 1:length(LH_ProcAllData)
        LH_allData(:,cc) = LH_ProcAllData{cc,1};
        RH_allData(:,cc) = RH_ProcAllData{cc,1};
    end
    % calculate the coherence between desired signals
    params.tapers = [10,19]; % Tapers [n, 2n - 1]
    [C_AllData,~,~,~,~,f_AllData,confC_AllData,~,cErr_AllData] = coherencyc(LH_allData,RH_allData,params);
    % save results
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).All.C = C_AllData;
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).All.f = f_AllData;
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).All.confC = confC_AllData;
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).All.cErr = cErr_AllData;
    %% NREM
    clear LH_nremData LH_nrem RH_nremData RH_nrem
    if strcmp(dataType,'HbT') == true
        [LH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    else
        [LH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_LH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    end
    % detrend and truncate data to minimum length to match events
    for ee = 1:length(LH_nremData)
        LH_nremData{ee,1} = detrend(LH_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        RH_nremData{ee,1} = detrend(RH_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
    end
    % pre-allocate coherence matrix
    LH_nrem = zeros(length(LH_nremData{1,1}),length(LH_nremData));
    RH_nrem = zeros(length(RH_nremData{1,1}),length(RH_nremData));
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for ff = 1:length(LH_nremData)
        LH_nrem(:,ff) = LH_nremData{ff,1};
        RH_nrem(:,ff) = RH_nremData{ff,1};
    end
    % calculate the coherence between desired signals
    params.tapers = [3,5]; % Tapers [n, 2n - 1]
    [C_nrem,~,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc(LH_nrem,RH_nrem,params);
    % save results
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).NREM.C = C_nrem;
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).NREM.f = f_nrem;
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).NREM.confC = confC_nrem;
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).NREM.cErr = cErr_nrem;
    %% REM
    clear LH_remData LH_rem RH_remData RH_rem
    if strcmp(dataType,'HbT') == true
        [LH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    else
        [LH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_LH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    end
    % detrend and truncate data to minimum length to match events
    for ee = 1:length(LH_remData)
        LH_remData{ee,1} = detrend(LH_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
        RH_remData{ee,1} = detrend(RH_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
    end
    % pre-allocate coherence matrix
    LH_rem = zeros(length(LH_remData{1,1}),length(LH_remData));
    RH_rem = zeros(length(RH_remData{1,1}),length(RH_remData));
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for ff = 1:length(LH_remData)
        LH_rem(:,ff) = LH_remData{ff,1};
        RH_rem(:,ff) = RH_remData{ff,1};
    end
    % calculate the coherence between desired signals
    params.tapers = [5,9]; % Tapers [n, 2n - 1]
    [C_rem,~,~,~,~,f_rem,confC_rem,~,cErr_rem] = coherencyc(LH_rem,RH_rem,params);
    % save results
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).REM.C = C_rem;
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).REM.f = f_rem;
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).REM.confC = confC_rem;
    Results_BilatCoher_Ephys.(group).(animalID).(dataType).REM.cErr = cErr_rem;
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_BilatCoher_Ephys.mat','Results_BilatCoher_Ephys')
cd([rootFolder delim 'Data'])