function [Results_PowerSpec] = AnalyzePowerSpectrum(animalID,group,rootFolder,delim,Results_PowerSpec)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs
dataLocation = [rootFolder delim group delim animalID delim 'Bilateral Imaging'];
cd(dataLocation)
% character list of all RawData file IDs
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
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
% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load SleepData.mat strut
SleepDataFileStruct = dir('*_SleepData.mat');
SleepDataFile = {SleepDataFileStruct.name}';
SleepDataFileID = char(SleepDataFile);
load(SleepDataFileID,'-mat')
% find and load Forest_ScoringResults.mat struct
forestScoringResultsFileID = [animalID '_Forest_ScoringResults.mat'];
load(forestScoringResultsFileID,'-mat')
samplingRate = RestData.CBV_HbT.adjLH.CBVCamSamplingRate;
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
% go through each valid data type for behavior-based power spectrum analysis
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    %% analyze power spectra during periods of rest
    % pull data from RestData.mat structure
    if strcmp(dataType,'CBV_HbT') == true
        [restLogical] = FilterEvents_IOS(RestData.(dataType).adjLH,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.(dataType).adjLH,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.(dataType).adjLH.fileIDs(combRestLogical,:);
        restEventTimes = RestData.(dataType).adjLH.eventTimes(combRestLogical,:);
        restDurations = RestData.(dataType).adjLH.durations(combRestLogical,:);
        LH_unstimRestingData = RestData.(dataType).adjLH.data(combRestLogical,:);
        RH_unstimRestingData = RestData.(dataType).adjRH.data(combRestLogical,:);
    else
        [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
        restEventTimes = RestData.cortical_LH.(dataType).eventTimes(combRestLogical,:);
        restDurations = RestData.cortical_LH.(dataType).durations(combRestLogical,:);
        LH_unstimRestingData =RestData.cortical_LH.(dataType).NormData(combRestLogical,:);
        RH_unstimRestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical,:);
        Hip_unstimRestingData = RestData.hippocampus.(dataType).NormData(combRestLogical,:);
    end
    % keep only the data that occurs within the manually-approved awake regions
    [LH_finalRestData,~,~,~] = RemoveInvalidData_IOS(LH_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS(RH_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    if strcmp(dataType,'CBV_HbT') == false
        [Hip_finalRestData,~,~,~] = RemoveInvalidData_IOS(Hip_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    end
    clear LH_ProcRestData RH_ProcRestData Hip_ProcRestData
    % detrend and truncate data to minimum length to match events
    for bb = 1:length(LH_finalRestData)
        if length(LH_finalRestData{bb,1}) < params.minTime.Rest*samplingRate
            restChunkSampleDiff = params.minTime.Rest*samplingRate - length(LH_finalRestData{bb,1});
            LH_restPad = (ones(1,restChunkSampleDiff))*LH_finalRestData{bb,1}(end);
            RH_restPad = (ones(1,restChunkSampleDiff))*RH_finalRestData{bb,1}(end);
            LH_ProcRestData{bb,1} = horzcat(LH_finalRestData{bb,1},LH_restPad); %#ok<*AGROW>
            RH_ProcRestData{bb,1} = horzcat(RH_finalRestData{bb,1},RH_restPad);
            LH_ProcRestData{bb,1} = detrend(LH_ProcRestData{bb,1},'constant');
            RH_ProcRestData{bb,1} = detrend(RH_ProcRestData{bb,1},'constant');
            if strcmp(dataType,'CBV_HbT') == false
                Hip_restPad = (ones(1,restChunkSampleDiff))*Hip_finalRestData{bb,1}(end);
                Hip_ProcRestData{bb,1} = horzcat(Hip_finalRestData{bb,1},Hip_restPad);
                Hip_ProcRestData{bb,1} = detrend(Hip_ProcRestData{bb,1},'constant');
            end
        else
            LH_ProcRestData{bb,1} = detrend(LH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            RH_ProcRestData{bb,1} = detrend(RH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            if strcmp(dataType,'CBV_HbT') == false
                Hip_ProcRestData{bb,1} = detrend(Hip_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            end
        end
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
    LH_restData = zeros(length(LH_ProcRestData{1,1}),length(LH_ProcRestData));
    RH_restData = zeros(length(RH_ProcRestData{1,1}),length(RH_ProcRestData));
    if strcmp(dataType,'CBV_HbT') == false
        Hip_restData = zeros(length(Hip_ProcRestData{1,1}),length(Hip_ProcRestData));
    end
    for cc = 1:length(LH_ProcRestData)
        LH_restData(:,cc) = LH_ProcRestData{cc,1};
        RH_restData(:,cc) = RH_ProcRestData{cc,1};
        if strcmp(dataType,'CBV_HbT') == false
            Hip_restData(:,cc) = Hip_ProcRestData{cc,1};
        end
    end
    % parameters for mtspectrumc - information available in function
    params.tapers = [5,9];   % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;
    params.fpass = [0,1];   % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2,0.05];
    % calculate the power spectra of the desired signals
    [LH_rest_S,LH_rest_f,LH_rest_sErr] = mtspectrumc(LH_restData,params);
    [RH_rest_S,RH_rest_f,RH_rest_sErr] = mtspectrumc(RH_restData,params);
    if strcmp(dataType,'CBV_HbT') == false
        [Hip_rest_S,Hip_rest_f,Hip_rest_sErr] = mtspectrumc(Hip_restData,params);
    end
    % save results
    Results_PowerSpec.(animalID).Rest.(dataType).adjLH.S = LH_rest_S;
    Results_PowerSpec.(animalID).Rest.(dataType).adjLH.f = LH_rest_f;
    Results_PowerSpec.(animalID).Rest.(dataType).adjLH.sErr = LH_rest_sErr;
    Results_PowerSpec.(animalID).Rest.(dataType).adjRH.S = RH_rest_S;
    Results_PowerSpec.(animalID).Rest.(dataType).adjRH.f = RH_rest_f;
    Results_PowerSpec.(animalID).Rest.(dataType).adjRH.sErr = RH_rest_sErr;
    if strcmp(dataType,'CBV_HbT') == false
        Results_PowerSpec.(animalID).Rest.(dataType).Hip.S = Hip_rest_S;
        Results_PowerSpec.(animalID).Rest.(dataType).Hip.f = Hip_rest_f;
        Results_PowerSpec.(animalID).Rest.(dataType).Hip.sErr = Hip_rest_sErr;
    end
    %% analyze power spectra during periods of alert
    zz = 1;
    clear LH_AwakeData RH_AwakeData Hip_AwakeData LH_ProcAwakeData RH_ProcAwakeData Hip_ProcAwakeData
    LH_AwakeData = [];
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
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144   % 36 bins (180 total) or 3 minutes of Asleep
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                if strcmp(dataType,'CBV_HbT') == true
                    LH_AwakeData{zz,1} = ProcData.data.(dataType).adjLH;
                    RH_AwakeData{zz,1} = ProcData.data.(dataType).adjRH;
                    zz = zz + 1;
                else
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        LH_AwakeData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                        RH_AwakeData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                        Hip_AwakeData{zz,1} = (ProcData.data.hippocampus.(dataType) - RestingBaselines.manualSelection.hippocampus.(dataType).(strDay))./RestingBaselines.manualSelection.hippocampus.(dataType).(strDay);
                        zz = zz + 1;
                    end
                end
            end
        end
    end
    if isempty(LH_AwakeData) == false
        % detrend data
        for bb = 1:length(LH_AwakeData)
            LH_ProcAwakeData{bb,1} = detrend(LH_AwakeData{bb,1},'constant');
            RH_ProcAwakeData{bb,1} = detrend(RH_AwakeData{bb,1},'constant');
            if strcmp(dataType,'CBV_HbT') == false
                Hip_ProcAwakeData{bb,1} = detrend(Hip_AwakeData{bb,1},'constant');
            end
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        LH_awakeData = zeros(length(LH_ProcAwakeData{1,1}),length(LH_ProcAwakeData));
        RH_awakeData = zeros(length(RH_ProcAwakeData{1,1}),length(RH_ProcAwakeData));
        if strcmp(dataType,'CBV_HbT') == false
            Hip_awakeData = zeros(length(Hip_ProcAwakeData{1,1}),length(Hip_ProcAwakeData));
        end
        for cc = 1:length(LH_ProcAwakeData)
            LH_awakeData(:,cc) = LH_ProcAwakeData{cc,1};
            RH_awakeData(:,cc) = RH_ProcAwakeData{cc,1};
            if strcmp(dataType,'CBV_HbT') == false
                Hip_awakeData(:,cc) = Hip_ProcAwakeData{cc,1};
            end
        end
        % calculate the power spectra of the desired signals
        [LH_awake_S,LH_awake_f,LH_awake_sErr] = mtspectrumc(LH_awakeData,params);
        [RH_awake_S,RH_awake_f,RH_awake_sErr] = mtspectrumc(RH_awakeData,params);
        if strcmp(dataType,'CBV_HbT') == false
            [Hip_awake_S,Hip_awake_f,Hip_awake_sErr] = mtspectrumc(Hip_awakeData,params);
        end
        % save results
        Results_PowerSpec.(animalID).Awake.(dataType).adjLH.S = LH_awake_S;
        Results_PowerSpec.(animalID).Awake.(dataType).adjLH.f = LH_awake_f;
        Results_PowerSpec.(animalID).Awake.(dataType).adjLH.sErr = LH_awake_sErr;
        Results_PowerSpec.(animalID).Awake.(dataType).adjRH.S = RH_awake_S;
        Results_PowerSpec.(animalID).Awake.(dataType).adjRH.f = RH_awake_f;
        Results_PowerSpec.(animalID).Awake.(dataType).adjRH.sErr = RH_awake_sErr;
        if strcmp(dataType,'CBV_HbT') == false
            Results_PowerSpec.(animalID).Awake.(dataType).Hip.S = Hip_awake_S;
            Results_PowerSpec.(animalID).Awake.(dataType).Hip.f = Hip_awake_f;
            Results_PowerSpec.(animalID).Awake.(dataType).Hip.sErr = Hip_awake_sErr;
        end
    else
        % save results
        Results_PowerSpec.(animalID).Awake.(dataType).adjLH.S = [];
        Results_PowerSpec.(animalID).Awake.(dataType).adjLH.f = [];
        Results_PowerSpec.(animalID).Awake.(dataType).adjLH.sErr = [];
        Results_PowerSpec.(animalID).Awake.(dataType).adjRH.S = [];
        Results_PowerSpec.(animalID).Awake.(dataType).adjRH.f = [];
        Results_PowerSpec.(animalID).Awake.(dataType).adjRH.sErr = [];
        if strcmp(dataType,'CBV_HbT') == false
            Results_PowerSpec.(animalID).Awake.(dataType).Hip.S = [];
            Results_PowerSpec.(animalID).Awake.(dataType).Hip.f = [];
            Results_PowerSpec.(animalID).Awake.(dataType).Hip.sErr = [];
        end
    end
    %% analyze power spectra during periods of Asleep
    zz = 1;
    clear LH_AsleepData RH_AsleepData Hip_AsleepData LH_ProcAsleepData RH_ProcAsleepData Hip_ProcAsleepData
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
        if sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of awake
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                if strcmp(dataType,'CBV_HbT') == true
                    LH_AsleepData{zz,1} = ProcData.data.(dataType).adjLH;
                    RH_AsleepData{zz,1} = ProcData.data.(dataType).adjRH;
                    zz = zz + 1;
                else
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        LH_AsleepData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                        RH_AsleepData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                        Hip_AsleepData{zz,1} = (ProcData.data.hippocampus.(dataType) - RestingBaselines.manualSelection.hippocampus.(dataType).(strDay))./RestingBaselines.manualSelection.hippocampus.(dataType).(strDay);
                        zz = zz + 1;
                    end
                end
            end
        end
    end
    if isempty(LH_AsleepData) == false
        % detrend data
        for bb = 1:length(LH_AsleepData)
            LH_ProcAsleepData{bb,1} = detrend(LH_AsleepData{bb,1},'constant');
            RH_ProcAsleepData{bb,1} = detrend(RH_AsleepData{bb,1},'constant');
            if strcmp(dataType,'CBV_HbT') == false
                Hip_ProcAsleepData{bb,1} = detrend(Hip_AsleepData{bb,1},'constant');
            end
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        LH_AsleepData = zeros(length(LH_ProcAsleepData{1,1}),length(LH_ProcAsleepData));
        RH_AsleepData = zeros(length(RH_ProcAsleepData{1,1}),length(RH_ProcAsleepData));
        if strcmp(dataType,'CBV_HbT') == false
            Hip_AsleepData = zeros(length(Hip_ProcAsleepData{1,1}),length(Hip_ProcAsleepData));
        end
        for cc = 1:length(LH_ProcAsleepData)
            LH_AsleepData(:,cc) = LH_ProcAsleepData{cc,1};
            RH_AsleepData(:,cc) = RH_ProcAsleepData{cc,1};
            if strcmp(dataType,'CBV_HbT') == false
                Hip_AsleepData(:,cc) = Hip_ProcAsleepData{cc,1};
            end
        end
        % calculate the power spectra of the desired signals
        [LH_Asleep_S,LH_Asleep_f,LH_Asleep_sErr] = mtspectrumc(LH_AsleepData,params);
        [RH_Asleep_S,RH_Asleep_f,RH_Asleep_sErr] = mtspectrumc(RH_AsleepData,params);
        if strcmp(dataType,'CBV_HbT') == false
            [Hip_Asleep_S,Hip_Asleep_f,Hip_Asleep_sErr] = mtspectrumc(Hip_AsleepData,params);
        end
        % save results
        Results_PowerSpec.(animalID).Asleep.(dataType).adjLH.S = LH_Asleep_S;
        Results_PowerSpec.(animalID).Asleep.(dataType).adjLH.f = LH_Asleep_f;
        Results_PowerSpec.(animalID).Asleep.(dataType).adjLH.sErr = LH_Asleep_sErr;
        Results_PowerSpec.(animalID).Asleep.(dataType).adjRH.S = RH_Asleep_S;
        Results_PowerSpec.(animalID).Asleep.(dataType).adjRH.f = RH_Asleep_f;
        Results_PowerSpec.(animalID).Asleep.(dataType).adjRH.sErr = RH_Asleep_sErr;
        if strcmp(dataType,'CBV_HbT') == false
            Results_PowerSpec.(animalID).Asleep.(dataType).Hip.S = Hip_Asleep_S;
            Results_PowerSpec.(animalID).Asleep.(dataType).Hip.f = Hip_Asleep_f;
            Results_PowerSpec.(animalID).Asleep.(dataType).Hip.sErr = Hip_Asleep_sErr;
        end
    else
        % save results
        Results_PowerSpec.(animalID).Asleep.(dataType).adjLH.S = [];
        Results_PowerSpec.(animalID).Asleep.(dataType).adjLH.f = [];
        Results_PowerSpec.(animalID).Asleep.(dataType).adjLH.sErr = [];
        Results_PowerSpec.(animalID).Asleep.(dataType).adjRH.S = [];
        Results_PowerSpec.(animalID).Asleep.(dataType).adjRH.f = [];
        Results_PowerSpec.(animalID).Asleep.(dataType).adjRH.sErr = [];
        if strcmp(dataType,'CBV_HbT') == false
            Results_PowerSpec.(animalID).Asleep.(dataType).Hip.S = [];
            Results_PowerSpec.(animalID).Asleep.(dataType).Hip.f = [];
            Results_PowerSpec.(animalID).Asleep.(dataType).Hip.sErr = [];
        end
    end
    %% analyze power spectra during periods of all data
    zz = 1;
    clear LH_AllUnstimData RH_AllUnstimData Hip_AllUnstimData LH_ProcAllUnstimData RH_ProcAllUnstimData Hip_ProcAllUnstimData
    LH_AllUnstimData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allUnstimDataFileDate,~] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allUnstimDataFileDate);
        load(procDataFileID,'-mat')
        puffs = ProcData.data.stimulations.LPadSol;
        % don't include trials with stimulation
        if isempty(puffs) == true
            if strcmp(dataType,'CBV_HbT') == true
                LH_AllUnstimData{zz,1} = ProcData.data.(dataType).adjLH;
                RH_AllUnstimData{zz,1} = ProcData.data.(dataType).adjRH;
                zz = zz + 1;
            else
                motionArtifact = ProcData.notes.motionArtifact;
                if motionArtifact == false
                    LH_AllUnstimData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                    RH_AllUnstimData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                    Hip_AllUnstimData{zz,1} = (ProcData.data.hippocampus.(dataType) - RestingBaselines.manualSelection.hippocampus.(dataType).(strDay))./RestingBaselines.manualSelection.hippocampus.(dataType).(strDay);
                    zz = zz + 1;
                end
            end
        end
    end
    if isempty(LH_AllUnstimData) == false
        % detrend data
        for bb = 1:length(LH_AllUnstimData)
            LH_ProcAllUnstimData{bb,1} = detrend(LH_AllUnstimData{bb,1},'constant');
            RH_ProcAllUnstimData{bb,1} = detrend(RH_AllUnstimData{bb,1},'constant');
            if strcmp(dataType,'CBV_HbT') == false
                Hip_ProcAllUnstimData{bb,1} = detrend(Hip_AllUnstimData{bb,1},'constant');
            end
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        LH_allUnstimData = zeros(length(LH_ProcAllUnstimData{1,1}),length(LH_ProcAllUnstimData));
        RH_allUnstimData = zeros(length(RH_ProcAllUnstimData{1,1}),length(RH_ProcAllUnstimData));
        if strcmp(dataType,'CBV_HbT') == false
            Hip_allUnstimData = zeros(length(Hip_ProcAllUnstimData{1,1}),length(Hip_ProcAllUnstimData));
        end
        for cc = 1:length(LH_ProcAllUnstimData)
            LH_allUnstimData(:,cc) = LH_ProcAllUnstimData{cc,1};
            RH_allUnstimData(:,cc) = RH_ProcAllUnstimData{cc,1};
            if strcmp(dataType,'CBV_HbT') == false
                Hip_allUnstimData(:,cc) = Hip_ProcAllUnstimData{cc,1};
            end
        end
        % calculate the power spectra of the desired signals
        [LH_allUnstim_S,LH_allUnstim_f,LH_allUnstim_sErr] = mtspectrumc(LH_allUnstimData,params);
        [RH_allUnstim_S,RH_allUnstim_f,RH_allUnstim_sErr] = mtspectrumc(RH_allUnstimData,params);
        if strcmp(dataType,'CBV_HbT') == false
            [Hip_allUnstim_S,Hip_allUnstim_f,Hip_allUnstim_sErr] = mtspectrumc(Hip_allUnstimData,params);
        end
        % save results
        Results_PowerSpec.(animalID).All.(dataType).adjLH.S = LH_allUnstim_S;
        Results_PowerSpec.(animalID).All.(dataType).adjLH.f = LH_allUnstim_f;
        Results_PowerSpec.(animalID).All.(dataType).adjLH.sErr = LH_allUnstim_sErr;
        Results_PowerSpec.(animalID).All.(dataType).adjRH.S = RH_allUnstim_S;
        Results_PowerSpec.(animalID).All.(dataType).adjRH.f = RH_allUnstim_f;
        Results_PowerSpec.(animalID).All.(dataType).adjRH.sErr = RH_allUnstim_sErr;
        if strcmp(dataType,'CBV_HbT') == false
            Results_PowerSpec.(animalID).All.(dataType).Hip.S = Hip_allUnstim_S;
            Results_PowerSpec.(animalID).All.(dataType).Hip.f = Hip_allUnstim_f;
            Results_PowerSpec.(animalID).All.(dataType).Hip.sErr = Hip_allUnstim_sErr;
        end
    end
    %% analyze power spectra during periods of NREM
    % pull data from SleepData.mat structure
    if strcmp(dataType,'CBV_HbT') == true
        [LH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    else
        [LH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_LH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [Hip_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.hippocampus.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    end
    % detrend and truncate data to minimum length to match events
    for dd = 1:length(LH_nremData)
        LH_nremData{dd,1} = detrend(LH_nremData{dd,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        RH_nremData{dd,1} = detrend(RH_nremData{dd,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        if strcmp(dataType,'CBV_HbT') == false
            Hip_nremData{dd,1} = detrend(Hip_nremData{dd,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        end
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
    LH_nrem = zeros(length(LH_nremData{1,1}),length(LH_nremData));
    RH_nrem = zeros(length(RH_nremData{1,1}),length(RH_nremData));
    if strcmp(dataType,'CBV_HbT') == false
        Hip_nrem = zeros(length(Hip_nremData{1,1}),length(Hip_nremData));
    end
    for ee = 1:length(LH_nremData)
        LH_nrem(:,ee) = LH_nremData{ee,1};
        RH_nrem(:,ee) = RH_nremData{ee,1};
        if strcmp(dataType,'CBV_HbT') == false
            Hip_nrem(:,ee) = Hip_nremData{ee,1};
        end
    end
    % calculate the power spectra of the desired signals
    [LH_nrem_S,LH_nrem_f,LH_nrem_sErr] = mtspectrumc(LH_nrem,params);
    [RH_nrem_S,RH_nrem_f,RH_nrem_sErr] = mtspectrumc(RH_nrem,params);
    if strcmp(dataType,'CBV_HbT') == false
        [Hip_nrem_S,Hip_nrem_f,Hip_nrem_sErr] = mtspectrumc(Hip_nrem,params);
    end
    % save results
    Results_PowerSpec.(animalID).NREM.(dataType).adjLH.S = LH_nrem_S;
    Results_PowerSpec.(animalID).NREM.(dataType).adjLH.f = LH_nrem_f;
    Results_PowerSpec.(animalID).NREM.(dataType).adjLH.sErr = LH_nrem_sErr;
    Results_PowerSpec.(animalID).NREM.(dataType).adjRH.S = RH_nrem_S;
    Results_PowerSpec.(animalID).NREM.(dataType).adjRH.f = RH_nrem_f;
    Results_PowerSpec.(animalID).NREM.(dataType).adjRH.sErr = RH_nrem_sErr;
    if strcmp(dataType,'CBV_HbT') == false
        Results_PowerSpec.(animalID).NREM.(dataType).Hip.S = Hip_nrem_S;
        Results_PowerSpec.(animalID).NREM.(dataType).Hip.f = Hip_nrem_f;
        Results_PowerSpec.(animalID).NREM.(dataType).Hip.sErr = Hip_nrem_sErr;
    end
    %% analyze power spectra during periods of REM
    % pull data from SleepData.mat structure
    if strcmp(dataType,'CBV_HbT') == true
        [LH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    elseif strcmp(dataType,'LFP') == false
        [LH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_LH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [Hip_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.hippocampus.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    end
    % detrend and truncate data to minimum length to match events
    for ff = 1:length(LH_remData)
        LH_remData{ff,1} = detrend(LH_remData{ff,1}(1:(params.minTime.REM*samplingRate)),'constant');
        RH_remData{ff,1} = detrend(RH_remData{ff,1}(1:(params.minTime.REM*samplingRate)),'constant');
        if strcmp(dataType,'CBV_HbT') == false
            Hip_remData{ff,1} = detrend(Hip_remData{ff,1}(1:(params.minTime.REM*samplingRate)),'constant');
        end
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
    LH_rem = zeros(length(LH_remData{1,1}),length(LH_remData));
    RH_rem = zeros(length(RH_remData{1,1}),length(RH_remData));
    if strcmp(dataType,'CBV_HbT') == false
        Hip_rem = zeros(length(Hip_remData{1,1}),length(Hip_remData));
    end
    for gg = 1:length(LH_remData)
        LH_rem(:,gg) = LH_remData{gg,1};
        RH_rem(:,gg) = RH_remData{gg,1};
        if strcmp(dataType,'CBV_HbT') == false
            Hip_rem(:,gg) = Hip_remData{gg,1};
        end
    end
    % calculate the power spectra of the desired signals
    [LH_rem_S,LH_rem_f,LH_rem_sErr] = mtspectrumc(LH_rem,params);
    [RH_rem_S,RH_rem_f,RH_rem_sErr] = mtspectrumc(RH_rem,params);
    if strcmp(dataType,'CBV_HbT') == false
        [Hip_rem_S,Hip_rem_f,Hip_rem_sErr] = mtspectrumc(Hip_rem,params);
    end
    %save results
    Results_PowerSpec.(animalID).REM.(dataType).adjLH.S = LH_rem_S;
    Results_PowerSpec.(animalID).REM.(dataType).adjLH.f = LH_rem_f;
    Results_PowerSpec.(animalID).REM.(dataType).adjLH.sErr = LH_rem_sErr;
    Results_PowerSpec.(animalID).REM.(dataType).adjRH.S = RH_rem_S;
    Results_PowerSpec.(animalID).REM.(dataType).adjRH.f = RH_rem_f;
    Results_PowerSpec.(animalID).REM.(dataType).adjRH.sErr = RH_rem_sErr;
    if strcmp(dataType,'CBV_HbT') == false
        Results_PowerSpec.(animalID).REM.(dataType).Hip.S = Hip_rem_S;
        Results_PowerSpec.(animalID).REM.(dataType).Hip.f = Hip_rem_f;
        Results_PowerSpec.(animalID).REM.(dataType).Hip.sErr = Hip_rem_sErr;
    end
end
%% analyze LFP power spectra during alert/Asleep/all
behavFields = {'Alert','Asleep','All'};
Data.Alert = []; Data.Asleep = []; Data.All = [];
dataTypes = {'LH','RH','Hip'};
xx = 1; yy = 1; zz = 1;
analogFs = 20000;
dsFs = 1000;
params.Fs = dsFs;
params.fpass = [1,100];
for bb = 1:size(rawDataFileIDs,1)
    rawDataFileID = rawDataFileIDs(bb,:);
    procDataFileID = procDataFileIDs(bb,:);
    [~,~,allDataFileID] = GetFileInfo_IOS(procDataFileID);
    scoringLabels = [];
    for cc = 1:length(ScoringResults.fileIDs)
        if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
            scoringLabels = ScoringResults.labels{cc,1};
        end
    end
    load(procDataFileID,'-mat')
    puffs = ProcData.data.stimulations.LPadSol;
    % don't include trials with stimulation
    if isempty(puffs) == true
        load(rawDataFileID,'-mat')
        Data.All.LH{xx,1} = resample(RawData.data.cortical_LH,dsFs,analogFs);
        Data.All.RH{xx,1} = resample(RawData.data.cortical_RH,dsFs,analogFs);
        Data.All.Hip{xx,1} = resample(RawData.data.hippocampus,dsFs,analogFs);
        xx = xx + 1;
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144   % 36 bins (180 total) or 3 minutes of Asleep
            Data.Alert.LH{yy,1} = resample(RawData.data.cortical_LH,dsFs,analogFs);
            Data.Alert.RH{yy,1} = resample(RawData.data.cortical_RH,dsFs,analogFs);
            Data.Alert.Hip{yy,1} = resample(RawData.data.hippocampus,dsFs,analogFs);
            yy = yy + 1;
        elseif sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of awake
            Data.Asleep.LH{zz,1} = resample(RawData.data.cortical_LH,dsFs,analogFs);
            Data.Asleep.RH{zz,1} = resample(RawData.data.cortical_RH,dsFs,analogFs);
            Data.Asleep.Hip{zz,1} = resample(RawData.data.hippocampus,dsFs,analogFs);
            zz = zz + 1;
        end
    end
end
%% Calculate LFP power spectrum
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        if isempty(Data.(behavField)) == false
            % detrend data
            procData = {};
            for cc = 1:length(Data.(behavField).(dataType))
                procData{cc,1} = detrend(Data.(behavField).(dataType){cc,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            data = zeros(length(procData{1,1}),length(procData));
            for dd = 1:length(procData)
                data(:,dd) = procData{dd,1};
            end
            % calculate the power spectra of the desired signals
            [S,f,sErr] = mtspectrumc(data,params);
            % save results
            AnalysisResults.(animalID).(behavField).LFP.(dataType).S = S;
            AnalysisResults.(animalID).(behavField).LFP.(dataType).f = f;
            AnalysisResults.(animalID).(behavField).LFP.(dataType).sErr = sErr;
        else
            % save results
            AnalysisResults.(animalID).(behavField).LFP.(dataType).S = [];
            AnalysisResults.(animalID).(behavField).LFP.(dataType).f = [];
            AnalysisResults.(animalID).(behavField).LFP.(dataType).sErr = [];
        end
    end
end
% save data
cd(rootFolder)
save('Results_PowerSpec.mat','Results_PowerSpec')

end
