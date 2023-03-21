function [Results_PowerSpec_Ephys] = AnalyzePowerSpectrum_Ephys(animalID,group,set,rootFolder,delim,Results_PowerSpec_Ephys)
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
% determine data types based on experimental group
dataTypes = {'CBV_HbT','gammaBandPower'};
% go through each valid data type for arousal based coherence analysis
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    % analyze bilateral coherence during periods of rest
    if strcmp(dataType,'CBV_HbT') == true
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
    clear LH_ProcRestData RH_ProcRestData
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
    % calculate the power spectra of the desired signals
    [LH_rest_S,LH_rest_f,LH_rest_sErr] = mtspectrumc(LH_restData,params);
    [RH_rest_S,RH_rest_f,RH_rest_sErr] = mtspectrumc(RH_restData,params);
    % save results
    Results_PowerSpec_Ephys.(group).(animalID).Rest.(dataType).LH.S = LH_rest_S;
    Results_PowerSpec_Ephys.(group).(animalID).Rest.(dataType).LH.f = LH_rest_f;
    Results_PowerSpec_Ephys.(group).(animalID).Rest.(dataType).LH.sErr = LH_rest_sErr;
    Results_PowerSpec_Ephys.(group).(animalID).Rest.(dataType).RH.S = RH_rest_S;
    Results_PowerSpec_Ephys.(group).(animalID).Rest.(dataType).RH.f = RH_rest_f;
    Results_PowerSpec_Ephys.(group).(animalID).Rest.(dataType).RH.sErr = RH_rest_sErr;
    % analyze bilateral coherence during periods of alert
    zz = 1;
    clear LH_AlertData RH_AlertData LH_ProcAlertData RH_ProcAlertData
    LH_AlertData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        scoringLabels = [];
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
                if strcmp(dataType,'CBV_HbT') == true
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
        % calculate the power spectra of the desired signals
        [LH_alert_S,LH_alert_f,LH_alert_sErr] = mtspectrumc(LH_alertData,params);
        [RH_alert_S,RH_alert_f,RH_alert_sErr] = mtspectrumc(RH_alertData,params);
        % save results
        Results_PowerSpec_Ephys.(group).(animalID).Alert.(dataType).LH.S = LH_alert_S;
        Results_PowerSpec_Ephys.(group).(animalID).Alert.(dataType).LH.f = LH_alert_f;
        Results_PowerSpec_Ephys.(group).(animalID).Alert.(dataType).LH.sErr = LH_alert_sErr;
        Results_PowerSpec_Ephys.(group).(animalID).Alert.(dataType).RH.S = RH_alert_S;
        Results_PowerSpec_Ephys.(group).(animalID).Alert.(dataType).RH.f = RH_alert_f;
        Results_PowerSpec_Ephys.(group).(animalID).Alert.(dataType).RH.sErr = RH_alert_sErr;
    else
        % save results
        Results_PowerSpec_Ephys.(group).(animalID).Alert.(dataType).LH.S = [];
        Results_PowerSpec_Ephys.(group).(animalID).Alert.(dataType).LH.f = [];
        Results_PowerSpec_Ephys.(group).(animalID).Alert.(dataType).LH.sErr = [];
        Results_PowerSpec_Ephys.(group).(animalID).Alert.(dataType).RH.S = [];
        Results_PowerSpec_Ephys.(group).(animalID).Alert.(dataType).RH.f = [];
        Results_PowerSpec_Ephys.(group).(animalID).Alert.(dataType).RH.sErr = [];
    end
    % analyze bilateral coherence during periods of aasleep
    zz = 1;
    clear LH_AsleepData RH_AsleepData LH_ProcAsleepData RH_ProcAsleepData
    LH_AsleepData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        scoringLabels = [];
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of alert
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                if strcmp(dataType,'CBV_HbT') == true
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
        % calculate the power spectra of the desired signals
        params.tapers = [10,19]; % Tapers [n, 2n - 1]
        [LH_asleep_S,LH_asleep_f,LH_asleep_sErr] = mtspectrumc(LH_asleepData,params);
        [RH_asleep_S,RH_asleep_f,RH_asleep_sErr] = mtspectrumc(RH_asleepData,params);
        % save results
        Results_PowerSpec_Ephys.(group).(animalID).Asleep.(dataType).LH.S = LH_asleep_S;
        Results_PowerSpec_Ephys.(group).(animalID).Asleep.(dataType).LH.f = LH_asleep_f;
        Results_PowerSpec_Ephys.(group).(animalID).Asleep.(dataType).LH.sErr = LH_asleep_sErr;
        Results_PowerSpec_Ephys.(group).(animalID).Asleep.(dataType).RH.S = RH_asleep_S;
        Results_PowerSpec_Ephys.(group).(animalID).Asleep.(dataType).RH.f = RH_asleep_f;
        Results_PowerSpec_Ephys.(group).(animalID).Asleep.(dataType).RH.sErr = RH_asleep_sErr;
    else
        % save results
        Results_PowerSpec_Ephys.(group).(animalID).Asleep.(dataType).LH.S = [];
        Results_PowerSpec_Ephys.(group).(animalID).Asleep.(dataType).LH.f = [];
        Results_PowerSpec_Ephys.(group).(animalID).Asleep.(dataType).LH.sErr = [];
        Results_PowerSpec_Ephys.(group).(animalID).Asleep.(dataType).RH.S = [];
        Results_PowerSpec_Ephys.(group).(animalID).Asleep.(dataType).RH.f = [];
        Results_PowerSpec_Ephys.(group).(animalID).Asleep.(dataType).RH.sErr = [];
    end
    % analyze bilateral coherence during periods of all data
    zz = 1;
    clear LH_AllData RH_AllData LH_ProcAllData RH_ProcAllData
    LH_AllData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,~] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        load(procDataFileID,'-mat')
        puffs = ProcData.data.stimulations.LPadSol;
        % don't include trials with stimulation
        if isempty(puffs) == true
            if strcmp(dataType,'CBV_HbT') == true
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
    % calculate the power spectra of the desired signals
    params.tapers = [10,19]; % Tapers [n, 2n - 1]
    [LH_all_S,LH_all_f,LH_all_sErr] = mtspectrumc(LH_allData,params);
    [RH_all_S,RH_all_f,RH_all_sErr] = mtspectrumc(RH_allData,params);
    % save results
    Results_PowerSpec_Ephys.(group).(animalID).All.(dataType).LH.S = LH_all_S;
    Results_PowerSpec_Ephys.(group).(animalID).All.(dataType).LH.f = LH_all_f;
    Results_PowerSpec_Ephys.(group).(animalID).All.(dataType).LH.sErr = LH_all_sErr;
    Results_PowerSpec_Ephys.(group).(animalID).All.(dataType).RH.S = RH_all_S;
    Results_PowerSpec_Ephys.(group).(animalID).All.(dataType).RH.f = RH_all_f;
    Results_PowerSpec_Ephys.(group).(animalID).All.(dataType).RH.sErr = RH_all_sErr;
    % analyze bilateral coherence during periods of NREM asleep
    % pull data from SleepData.mat structure
    if strcmp(dataType,'CBV_HbT') == true
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
    % calculate the power spectra of the desired signals
    params.tapers = [3,5]; % Tapers [n, 2n - 1]
    [LH_nrem_S,LH_nrem_f,LH_nrem_sErr] = mtspectrumc(LH_nrem,params);
    [RH_nrem_S,RH_nrem_f,RH_nrem_sErr] = mtspectrumc(RH_nrem,params);
    % save results
    Results_PowerSpec_Ephys.(group).(animalID).NREM.(dataType).LH.S = LH_nrem_S;
    Results_PowerSpec_Ephys.(group).(animalID).NREM.(dataType).LH.f = LH_nrem_f;
    Results_PowerSpec_Ephys.(group).(animalID).NREM.(dataType).LH.sErr = LH_nrem_sErr;
    Results_PowerSpec_Ephys.(group).(animalID).NREM.(dataType).RH.S = RH_nrem_S;
    Results_PowerSpec_Ephys.(group).(animalID).NREM.(dataType).RH.f = RH_nrem_f;
    Results_PowerSpec_Ephys.(group).(animalID).NREM.(dataType).RH.sErr = RH_nrem_sErr;
    % analyze bilateral coherence during periods of REM asleep
    % pull data from SleepData.mat structure
    if strcmp(dataType,'CBV_HbT') == true
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
    % calculate the power spectra of the desired signals
    params.tapers = [5,9]; % Tapers [n, 2n - 1]
    [LH_rem_S,LH_rem_f,LH_rem_sErr] = mtspectrumc(LH_rem,params);
    [RH_rem_S,RH_rem_f,RH_rem_sErr] = mtspectrumc(RH_rem,params);
    % save results
    Results_PowerSpec_Ephys.(group).(animalID).REM.(dataType).LH.S = LH_rem_S;
    Results_PowerSpec_Ephys.(group).(animalID).REM.(dataType).LH.f = LH_rem_f;
    Results_PowerSpec_Ephys.(group).(animalID).REM.(dataType).LH.sErr = LH_rem_sErr;
    Results_PowerSpec_Ephys.(group).(animalID).REM.(dataType).RH.S = RH_rem_S;
    Results_PowerSpec_Ephys.(group).(animalID).REM.(dataType).RH.f = RH_rem_f;
    Results_PowerSpec_Ephys.(group).(animalID).REM.(dataType).RH.sErr = RH_rem_sErr;
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_PowerSpec_Ephys.mat','Results_PowerSpec_Ephys')
cd([rootFolder delim 'Data'])
