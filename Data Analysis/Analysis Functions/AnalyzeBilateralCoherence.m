function [Results_BilatCoher] = AnalyzeBilateralCoherence(animalID,group,rootFolder,delim,Results_BilatCoher)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the spectral coherence between bilateral hemodynamic [HbT] and neural signals (IOS)
%________________________________________________________________________________________________________________________

modelType = 'Manual';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% only run analysis for valid animal IDs
dataLocation = [rootFolder delim 'Data' delim group delim animalID delim 'Bilateral Imaging'];
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
strBreaks = strfind(group,delim);
groupName = group(strBreaks + 1:end);
if strcmp(groupName,'IOS_Ephys') == true
    dataTypes = {'CBV_HbT','gammaBandPower'};
elseif strcmp(groupName,'IOS_GCaMP7s') == true
    dataTypes = {'CBV_HbT','GCaMP7s','Deoxy'};
end
% go through each valid data type for arousal based coherence analysis
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    %% analyze bilateral coherence during periods of rest
    if strcmp(dataType,'CBV_HbT') == true || strcmp(dataType,'Deoxy') == true
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
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_RestingData = RestData.(dataType).frontalLH.data(combRestLogical,:);
            fRH_RestingData = RestData.(dataType).frontalRH.data(combRestLogical,:);
            % keep only the data that occurs within the manually-approved alert regions
            [fLH_finalRestData,~,~,~] = RemoveInvalidData_IOS(fLH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
            [fRH_finalRestData,~,~,~] = RemoveInvalidData_IOS(fRH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        end
    elseif strcmp(dataType,'GCaMP7s') == true
        samplingRate = RestData.(dataType).corLH.CBVCamSamplingRate;
        [restLogical] = FilterEvents_IOS(RestData.(dataType).corLH,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.(dataType).corLH,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.(dataType).corLH.fileIDs(combRestLogical,:);
        restEventTimes = RestData.(dataType).corLH.eventTimes(combRestLogical,:);
        restDurations = RestData.(dataType).corLH.durations(combRestLogical,:);
        LH_RestingData = RestData.(dataType).corLH.data(combRestLogical,:);
        RH_RestingData = RestData.(dataType).corRH.data(combRestLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [LH_finalRestData,~,~,~] = RemoveInvalidData_IOS(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        fLH_RestingData = RestData.(dataType).corFrontalLH.data(combRestLogical,:);
        fRH_RestingData = RestData.(dataType).corFrontalRH.data(combRestLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [fLH_finalRestData,~,~,~] = RemoveInvalidData_IOS(fLH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [fRH_finalRestData,~,~,~] = RemoveInvalidData_IOS(fRH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
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
    clear LH_ProcRestData RH_ProcRestData fLH_ProcRestData fRH_ProcRestData
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
            if strcmp(groupName,'IOS_GCaMP7s') == true
                fLH_ProcRestData{bb,1} = horzcat(fLH_finalRestData{bb,1},LH_restPad);
                fRH_ProcRestData{bb,1} = horzcat(fRH_finalRestData{bb,1},RH_restPad);
                fLH_ProcRestData{bb,1} = detrend(fLH_ProcRestData{bb,1},'constant');
                fRH_ProcRestData{bb,1} = detrend(fRH_ProcRestData{bb,1},'constant');
            end
        else
            LH_ProcRestData{bb,1} = detrend(LH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            RH_ProcRestData{bb,1} = detrend(RH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            if strcmp(groupName,'IOS_GCaMP7s') == true
                fLH_ProcRestData{bb,1} = detrend(fLH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
                fRH_ProcRestData{bb,1} = detrend(fRH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            end
        end
    end
    % pre-allocate coherence matrix
    LH_restData = zeros(length(LH_ProcRestData{1,1}),length(LH_ProcRestData));
    RH_restData = zeros(length(RH_ProcRestData{1,1}),length(RH_ProcRestData));
    if strcmp(groupName,'IOS_GCaMP7s') == true
        fLH_restData = zeros(length(fLH_ProcRestData{1,1}),length(fLH_ProcRestData));
        fRH_restData = zeros(length(fRH_ProcRestData{1,1}),length(fRH_ProcRestData));
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for cc = 1:length(LH_ProcRestData)
        LH_restData(:,cc) = LH_ProcRestData{cc,1};
        RH_restData(:,cc) = RH_ProcRestData{cc,1};
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_restData(:,cc) = fLH_ProcRestData{cc,1};
            fRH_restData(:,cc) = fRH_ProcRestData{cc,1};
        end
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
    Results_BilatCoher.(animalID).Rest.(dataType).C = C_RestData;
    Results_BilatCoher.(animalID).Rest.(dataType).f = f_RestData;
    Results_BilatCoher.(animalID).Rest.(dataType).confC = confC_RestData;
    Results_BilatCoher.(animalID).Rest.(dataType).cErr = cErr_RestData;
    if strcmp(groupName,'IOS_GCaMP7s') == true
        % calculate the coherence between desired signals
        [fC_RestData,~,~,~,~,ff_RestData,fconfC_RestData,~,fcErr_RestData] = coherencyc(fLH_restData,fRH_restData,params);
        % save results
        Results_BilatCoher.(animalID).Rest.(dataType).fC = fC_RestData;
        Results_BilatCoher.(animalID).Rest.(dataType).ff = ff_RestData;
        Results_BilatCoher.(animalID).Rest.(dataType).fconfC = fconfC_RestData;
        Results_BilatCoher.(animalID).Rest.(dataType).fcErr = fcErr_RestData;
    end
    %% analyze bilateral coherence during periods of alert
    zz = 1;
    clear LH_AlertData RH_AlertData LH_ProcAlertData RH_ProcAlertData
    clear fLH_AlertData fRH_AlertData fLH_ProcAlertData fRH_ProcAlertData
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
                if strcmp(dataType,'CBV_HbT') == true || strcmp(dataType,'Deoxy') == true
                    LH_AlertData{zz,1} = ProcData.data.(dataType).LH;
                    RH_AlertData{zz,1} = ProcData.data.(dataType).RH;
                    if strcmp(groupName,'IOS_GCaMP7s') == true
                        fLH_AlertData{zz,1} = ProcData.data.(dataType).frontalLH;
                        fRH_AlertData{zz,1} = ProcData.data.(dataType).frontalRH;
                    end
                    zz = zz + 1;
                elseif strcmp(dataType,'GCaMP7s') == true
                    LH_AlertData{zz,1} = ProcData.data.(dataType).corLH;
                    RH_AlertData{zz,1} = ProcData.data.(dataType).corRH;
                    fLH_AlertData{zz,1} = ProcData.data.(dataType).corFrontalLH;
                    fRH_AlertData{zz,1} = ProcData.data.(dataType).corFrontalRH;
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
            if strcmp(groupName,'IOS_GCaMP7s') == true
                fLH_ProcAlertData{bb,1} = detrend(fLH_AlertData{bb,1},'constant');
                fRH_ProcAlertData{bb,1} = detrend(fRH_AlertData{bb,1},'constant');
            end
        end
        % preallocate coherence matrix
        LH_alertData = zeros(length(LH_ProcAlertData{1,1}),length(LH_ProcAlertData));
        RH_alertData = zeros(length(RH_ProcAlertData{1,1}),length(RH_ProcAlertData));
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_alertData = zeros(length(fLH_ProcAlertData{1,1}),length(fLH_ProcAlertData));
            fRH_alertData = zeros(length(fRH_ProcAlertData{1,1}),length(fRH_ProcAlertData));
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        for cc = 1:length(LH_ProcAlertData)
            LH_alertData(:,cc) = LH_ProcAlertData{cc,1};
            RH_alertData(:,cc) = RH_ProcAlertData{cc,1};
            if strcmp(groupName,'IOS_GCaMP7s') == true
                fLH_alertData(:,cc) = fLH_ProcAlertData{cc,1};
                fRH_alertData(:,cc) = fRH_ProcAlertData{cc,1};
            end
        end
        % calculate the coherence between desired signals
        params.tapers = [10,19]; % Tapers [n, 2n - 1]
        [C_AlertData,~,~,~,~,f_AlertData,confC_AlertData,~,cErr_AlertData] = coherencyc(LH_alertData,RH_alertData,params);
        % save results
        Results_BilatCoher.(animalID).Alert.(dataType).C = C_AlertData;
        Results_BilatCoher.(animalID).Alert.(dataType).f = f_AlertData;
        Results_BilatCoher.(animalID).Alert.(dataType).confC = confC_AlertData;
        Results_BilatCoher.(animalID).Alert.(dataType).cErr = cErr_AlertData;
        if strcmp(groupName,'IOS_GCaMP7s') == true
            % calculate the coherence between desired signals
            [fC_AlertData,~,~,~,~,ff_AlertData,fconfC_AlertData,~,fcErr_AlertData] = coherencyc(fLH_alertData,fRH_alertData,params);
            % save results
            Results_BilatCoher.(animalID).Alert.(dataType).fC = fC_AlertData;
            Results_BilatCoher.(animalID).Alert.(dataType).ff = ff_AlertData;
            Results_BilatCoher.(animalID).Alert.(dataType).fconfC = fconfC_AlertData;
            Results_BilatCoher.(animalID).Alert.(dataType).fcErr = fcErr_AlertData;
        end
    else
        % save results
        Results_BilatCoher.(animalID).Alert.(dataType).C = [];
        Results_BilatCoher.(animalID).Alert.(dataType).f = [];
        Results_BilatCoher.(animalID).Alert.(dataType).confC = [];
        Results_BilatCoher.(animalID).Alert.(dataType).cErr = [];
        if strcmp(groupName,'IOS_GCaMP7s') == true
            % save results
            Results_BilatCoher.(animalID).Alert.(dataType).fC = [];
            Results_BilatCoher.(animalID).Alert.(dataType).ff = [];
            Results_BilatCoher.(animalID).Alert.(dataType).fconfC = [];
            Results_BilatCoher.(animalID).Alert.(dataType).fcErr = [];
        end
    end
    %% analyze bilateral coherence during periods of aasleep
    zz = 1;
    clear LH_AsleepData RH_AsleepData LH_ProcAsleepData RH_ProcAsleepData
    clear fLH_AsleepData fRH_AsleepData fLH_ProcAsleepData fRH_ProcAsleepData
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
        if sum(strcmp(scoringLabels,'Not Sleep')) < 36 % 36 bins (180 total) or 3 minutes of alert
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                if strcmp(dataType,'CBV_HbT') == true || strcmp(dataType,'Deoxy') == true
                    LH_AsleepData{zz,1} = ProcData.data.(dataType).LH;
                    RH_AsleepData{zz,1} = ProcData.data.(dataType).RH;
                    if strcmp(groupName,'IOS_GCaMP7s') == true
                        fLH_AsleepData{zz,1} = ProcData.data.(dataType).frontalLH;
                        fRH_AsleepData{zz,1} = ProcData.data.(dataType).frontalRH;
                    end
                    zz = zz + 1;
                elseif strcmp(dataType,'GCaMP7s') == true
                    LH_AsleepData{zz,1} = ProcData.data.(dataType).corLH;
                    RH_AsleepData{zz,1} = ProcData.data.(dataType).corRH;
                    fLH_AsleepData{zz,1} = ProcData.data.(dataType).corFrontalLH;
                    fRH_AsleepData{zz,1} = ProcData.data.(dataType).corFrontalRH;
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
            if strcmp(groupName,'IOS_GCaMP7s') == true
                fLH_ProcAsleepData{bb,1} = detrend(fLH_AsleepData{bb,1},'constant');
                fRH_ProcAsleepData{bb,1} = detrend(fRH_AsleepData{bb,1},'constant');
            end
        end
        % preallocate coherence matrix
        LH_asleepData = zeros(length(LH_ProcAsleepData{1,1}),length(LH_ProcAsleepData));
        RH_asleepData = zeros(length(RH_ProcAsleepData{1,1}),length(RH_ProcAsleepData));
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_asleepData = zeros(length(fLH_ProcAsleepData{1,1}),length(fLH_ProcAsleepData));
            fRH_asleepData = zeros(length(fRH_ProcAsleepData{1,1}),length(fRH_ProcAsleepData));
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        for cc = 1:length(LH_ProcAsleepData)
            LH_asleepData(:,cc) = LH_ProcAsleepData{cc,1};
            RH_asleepData(:,cc) = RH_ProcAsleepData{cc,1};
            if strcmp(groupName,'IOS_GCaMP7s') == true
                fLH_asleepData(:,cc) = fLH_ProcAsleepData{cc,1};
                fRH_asleepData(:,cc) = fRH_ProcAsleepData{cc,1};
            end
        end
        % calculate the coherence between desired signals
        params.tapers = [10,19]; % Tapers [n, 2n - 1]
        [C_AsleepData,~,~,~,~,f_AsleepData,confC_AsleepData,~,cErr_AsleepData] = coherencyc(LH_asleepData,RH_asleepData,params);
        % save results
        Results_BilatCoher.(animalID).Asleep.(dataType).C = C_AsleepData;
        Results_BilatCoher.(animalID).Asleep.(dataType).f = f_AsleepData;
        Results_BilatCoher.(animalID).Asleep.(dataType).confC = confC_AsleepData;
        Results_BilatCoher.(animalID).Asleep.(dataType).cErr = cErr_AsleepData;
        if strcmp(groupName,'IOS_GCaMP7s') == true
            % calculate the coherence between desired signals
            [fC_AsleepData,~,~,~,~,ff_AsleepData,fconfC_AsleepData,~,fcErr_AsleepData] = coherencyc(fLH_asleepData,fRH_asleepData,params);
            % save results
            Results_BilatCoher.(animalID).Asleep.(dataType).fC = fC_AsleepData;
            Results_BilatCoher.(animalID).Asleep.(dataType).ff = ff_AsleepData;
            Results_BilatCoher.(animalID).Asleep.(dataType).fconfC = fconfC_AsleepData;
            Results_BilatCoher.(animalID).Asleep.(dataType).fcErr = fcErr_AsleepData;
        end
    else
        % save results
        Results_BilatCoher.(animalID).Asleep.(dataType).C = [];
        Results_BilatCoher.(animalID).Asleep.(dataType).f = [];
        Results_BilatCoher.(animalID).Asleep.(dataType).confC = [];
        Results_BilatCoher.(animalID).Asleep.(dataType).cErr = [];
        if strcmp(groupName,'IOS_GCaMP7s') == true
            % save results
            Results_BilatCoher.(animalID).Asleep.(dataType).fC = [];
            Results_BilatCoher.(animalID).Asleep.(dataType).ff = [];
            Results_BilatCoher.(animalID).Asleep.(dataType).fconfC = [];
            Results_BilatCoher.(animalID).Asleep.(dataType).fcErr = [];
        end
    end
    %% analyze bilateral coherence during periods of all data
    zz = 1;
    clear LH_AllData RH_AllData LH_ProcAllData RH_ProcAllData
    clear fLH_AllData fRH_AllData fLH_ProcAllData fRH_ProcAllData
    LH_AllData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,~] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        load(procDataFileID,'-mat')
        puffs = ProcData.data.stimulations.LPadSol;
        % don't include trials with stimulation
        if isempty(puffs) == true
            if strcmp(dataType,'CBV_HbT') == true || strcmp(dataType,'Deoxy') == true
                LH_AllData{zz,1} = ProcData.data.(dataType).LH;
                RH_AllData{zz,1} = ProcData.data.(dataType).RH;
                if strcmp(groupName,'IOS_GCaMP7s') == true
                    fLH_AllData{zz,1} = ProcData.data.(dataType).frontalLH;
                    fRH_AllData{zz,1} = ProcData.data.(dataType).frontalRH;
                end
                zz = zz + 1;
            elseif strcmp(dataType,'GCaMP7s') == true
                LH_AllData{zz,1} = ProcData.data.(dataType).corLH;
                RH_AllData{zz,1} = ProcData.data.(dataType).corRH;
                fLH_AllData{zz,1} = ProcData.data.(dataType).corFrontalLH;
                fRH_AllData{zz,1} = ProcData.data.(dataType).corFrontalRH;
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
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_ProcAllData{bb,1} = detrend(fLH_AllData{bb,1},'constant');
            fRH_ProcAllData{bb,1} = detrend(fRH_AllData{bb,1},'constant');
        end
    end
    % preallocate coherence matrix
    LH_allData = zeros(length(LH_ProcAllData{1,1}),length(LH_ProcAllData));
    RH_allData = zeros(length(RH_ProcAllData{1,1}),length(RH_ProcAllData));
    if strcmp(groupName,'IOS_GCaMP7s') == true
        fLH_allData = zeros(length(fLH_ProcAllData{1,1}),length(fLH_ProcAllData));
        fRH_allData = zeros(length(fRH_ProcAllData{1,1}),length(fRH_ProcAllData));
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for cc = 1:length(LH_ProcAllData)
        LH_allData(:,cc) = LH_ProcAllData{cc,1};
        RH_allData(:,cc) = RH_ProcAllData{cc,1};
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_allData(:,cc) = fLH_ProcAllData{cc,1};
            fRH_allData(:,cc) = fRH_ProcAllData{cc,1};
        end
    end
    % calculate the coherence between desired signals
    params.tapers = [10,19]; % Tapers [n, 2n - 1]
    [C_AllData,~,~,~,~,f_AllData,confC_AllData,~,cErr_AllData] = coherencyc(LH_allData,RH_allData,params);
    % save results
    Results_BilatCoher.(animalID).All.(dataType).C = C_AllData;
    Results_BilatCoher.(animalID).All.(dataType).f = f_AllData;
    Results_BilatCoher.(animalID).All.(dataType).confC = confC_AllData;
    Results_BilatCoher.(animalID).All.(dataType).cErr = cErr_AllData;
    if strcmp(groupName,'IOS_GCaMP7s') == true
        % calculate the coherence between desired signals
        [fC_AllData,~,~,~,~,ff_AllData,fconfC_AllData,~,fcErr_AllData] = coherencyc(fLH_allData,fRH_allData,params);
        % save results
        Results_BilatCoher.(animalID).All.(dataType).fC = fC_AllData;
        Results_BilatCoher.(animalID).All.(dataType).ff = ff_AllData;
        Results_BilatCoher.(animalID).All.(dataType).fconfC = fconfC_AllData;
        Results_BilatCoher.(animalID).All.(dataType).fcErr = fcErr_AllData;
    end
    %% analyze bilateral coherence during periods of NREM asleep
    % pull data from SleepData.mat structure
    if strcmp(dataType,'CBV_HbT') == true || strcmp(dataType,'GCaMP7s') == true || strcmp(dataType,'Deoxy') == true
        [LH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        if strcmp(groupName,'IOS_GCaMP7s') == true
            [fLH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).frontalLH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            [fRH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).frontalRH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        end
    else
        [LH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_LH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    end
    % detrend and truncate data to minimum length to match events
    for ee = 1:length(LH_nremData)
        LH_nremData{ee,1} = detrend(LH_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        RH_nremData{ee,1} = detrend(RH_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_nremData{ee,1} = detrend(fLH_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
            fRH_nremData{ee,1} = detrend(fRH_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        end
    end
    % pre-allocate coherence matrix
    LH_nrem = zeros(length(LH_nremData{1,1}),length(LH_nremData));
    RH_nrem = zeros(length(RH_nremData{1,1}),length(RH_nremData));
    if strcmp(groupName,'IOS_GCaMP7s') == true
        fLH_nrem = zeros(length(fLH_nremData{1,1}),length(fLH_nremData));
        fRH_nrem = zeros(length(fRH_nremData{1,1}),length(fRH_nremData));
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for ff = 1:length(LH_nremData)
        LH_nrem(:,ff) = LH_nremData{ff,1};
        RH_nrem(:,ff) = RH_nremData{ff,1};
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_nrem(:,ff) = fLH_nremData{ff,1};
            fRH_nrem(:,ff) = fRH_nremData{ff,1};
        end
    end
    % calculate the coherence between desired signals
    params.tapers = [3,5]; % Tapers [n, 2n - 1]
    [C_nrem,~,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc(LH_nrem,RH_nrem,params);
    % save results
    Results_BilatCoher.(animalID).NREM.(dataType).C = C_nrem;
    Results_BilatCoher.(animalID).NREM.(dataType).f = f_nrem;
    Results_BilatCoher.(animalID).NREM.(dataType).confC = confC_nrem;
    Results_BilatCoher.(animalID).NREM.(dataType).cErr = cErr_nrem;
    if strcmp(groupName,'IOS_GCaMP7s') == true
        % calculate the coherence between desired signals
        [fC_nrem,~,~,~,~,ff_nrem,fconfC_nrem,~,fcErr_nrem] = coherencyc(fLH_nrem,fRH_nrem,params);
        % save results
        Results_BilatCoher.(animalID).NREM.(dataType).fC = fC_nrem;
        Results_BilatCoher.(animalID).NREM.(dataType).ff = ff_nrem;
        Results_BilatCoher.(animalID).NREM.(dataType).fconfC = fconfC_nrem;
        Results_BilatCoher.(animalID).NREM.(dataType).fcErr = fcErr_nrem;
    end
    %% analyze bilateral coherence during periods of REM asleep
    % pull data from SleepData.mat structure
    if strcmp(dataType,'CBV_HbT') == true || strcmp(dataType,'GCaMP7s') == true || strcmp(dataType,'Deoxy') == true
        [LH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        if strcmp(groupName,'IOS_GCaMP7s') == true
            [fLH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).frontalLH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            [fRH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).frontalRH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        end
    else
        [LH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_LH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    end
    % detrend and truncate data to minimum length to match events
    for ee = 1:length(LH_remData)
        LH_remData{ee,1} = detrend(LH_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
        RH_remData{ee,1} = detrend(RH_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_remData{ee,1} = detrend(fLH_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
            fRH_remData{ee,1} = detrend(fRH_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
        end
    end
    % pre-allocate coherence matrix
    LH_rem = zeros(length(LH_remData{1,1}),length(LH_remData));
    RH_rem = zeros(length(RH_remData{1,1}),length(RH_remData));
    if strcmp(groupName,'IOS_GCaMP7s') == true
        fLH_rem = zeros(length(fLH_remData{1,1}),length(fLH_remData));
        fRH_rem = zeros(length(fRH_remData{1,1}),length(fRH_remData));
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for ff = 1:length(LH_remData)
        LH_rem(:,ff) = LH_remData{ff,1};
        RH_rem(:,ff) = RH_remData{ff,1};
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_rem(:,ff) = fLH_remData{ff,1};
            fRH_rem(:,ff) = fRH_remData{ff,1};
        end
    end
    % calculate the coherence between desired signals
    params.tapers = [5,9]; % Tapers [n, 2n - 1]
    [C_rem,~,~,~,~,f_rem,confC_rem,~,cErr_rem] = coherencyc(LH_rem,RH_rem,params);
    % save results
    Results_BilatCoher.(animalID).REM.(dataType).C = C_rem;
    Results_BilatCoher.(animalID).REM.(dataType).f = f_rem;
    Results_BilatCoher.(animalID).REM.(dataType).confC = confC_rem;
    Results_BilatCoher.(animalID).REM.(dataType).cErr = cErr_rem;
    if strcmp(groupName,'IOS_GCaMP7s') == true
        % calculate the coherence between desired signals
        [fC_rem,~,~,~,~,ff_rem,fconfC_rem,~,fcErr_rem] = coherencyc(fLH_rem,fRH_rem,params);
        % save results
        Results_BilatCoher.(animalID).REM.(dataType).fC = fC_rem;
        Results_BilatCoher.(animalID).REM.(dataType).ff = ff_rem;
        Results_BilatCoher.(animalID).REM.(dataType).fconfC = fconfC_rem;
        Results_BilatCoher.(animalID).REM.(dataType).fcErr = fcErr_rem;
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_BilatCoher.mat','Results_BilatCoher')

end
