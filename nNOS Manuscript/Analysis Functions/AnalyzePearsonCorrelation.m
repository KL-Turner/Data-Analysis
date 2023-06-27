function [Results_PearsonCorr] = AnalyzePearsonCorrelation(animalID,group,rootFolder,delim,Results_PearsonCorr)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze Pearson's correlation coefficient between bilateral hemodynamic [HbT] and neural signals (IOS)
%________________________________________________________________________________________________________________________

% function parameters
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.Whisk = 7;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% only run analysis for valid animal IDs
dataLocation = [rootFolder delim group delim animalID   delim 'Bilateral Imaging'];
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
% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID,'-mat')
% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load AsleepData.mat strut
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
% find and load Forest_ScoringResults.mat struct
forestScoringResultsFileID = [animalID '_Forest_ScoringResults.mat'];
load(forestScoringResultsFileID,'-mat')
% criteria for whisking
WhiskCriteria.Fieldname = {'duration','puffDistance'};
WhiskCriteria.Comparison = {'gt','gt'};
WhiskCriteria.Value = {5,5};
WhiskPuffCriteria.Fieldname = {'puffDistance'};
WhiskPuffCriteria.Comparison = {'gt'};
WhiskPuffCriteria.Value = {5};
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
% go through each valid data type for arousal-based correlation analysis
for a = 1:length(dataTypes)
    dataType = dataTypes{1,a};
    %% analyze Pearson's correlation coefficient during periods of rest
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
    clear LH_ProcRestData RH_ProcRestData
    clear fLH_ProcRestData fRH_ProcRestData
    % lowpass filter
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    % filter, detrend, and truncate data to minimum length to match events
    for gg = 1:length(LH_finalRestData)
        LH_ProcRestData{gg,1} = detrend(filtfilt(sos,g,LH_finalRestData{gg,1}(1:params.minTime.Rest*samplingRate)),'constant');
        RH_ProcRestData{gg,1} = detrend(filtfilt(sos,g,RH_finalRestData{gg,1}(1:params.minTime.Rest*samplingRate)),'constant');
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_ProcRestData{gg,1} = detrend(filtfilt(sos,g,fLH_finalRestData{gg,1}(1:params.minTime.Rest*samplingRate)),'constant');
            fRH_ProcRestData{gg,1} = detrend(filtfilt(sos,g,fRH_finalRestData{gg,1}(1:params.minTime.Rest*samplingRate)),'constant');
        end
    end
    % analyze correlation coefficient of resting epochs
    for n = 1:length(LH_ProcRestData)
        rest_CC = corrcoef(LH_ProcRestData{n,1},RH_ProcRestData{n,1});
        rest_R(n,1) = rest_CC(2,1);
        if strcmp(groupName,'IOS_GCaMP7s') == true
            rest_fCC = corrcoef(fLH_ProcRestData{n,1},fRH_ProcRestData{n,1});
            rest_fR(n,1) = rest_fCC(2,1);
        end
    end
    meanRest_R = mean(rest_R);
    stdRest_R = std(rest_R,0,1);
    % save results
    Results_PearsonCorr.(animalID).Rest.(dataType).R = rest_R;
    Results_PearsonCorr.(animalID).Rest.(dataType).meanR = meanRest_R;
    Results_PearsonCorr.(animalID).Rest.(dataType).stdR = stdRest_R;
    if strcmp(groupName,'IOS_GCaMP7s') == true
        meanRest_fR = mean(rest_fR);
        stdRest_fR = std(rest_fR,0,1);
        % save results
        Results_PearsonCorr.(animalID).Rest.(dataType).fR = rest_fR;
        Results_PearsonCorr.(animalID).Rest.(dataType).meanfR = meanRest_fR;
        Results_PearsonCorr.(animalID).Rest.(dataType).stdfR = stdRest_fR;
    end
    %% analyze Pearson's correlation coefficient during periods of moderate whisking (2-5 seconds)
    if strcmp(dataType,'CBV_HbT') == true || strcmp(dataType,'Deoxy') == true
        [whiskLogical] = FilterEvents_IOS(EventData.(dataType).LH.whisk,WhiskCriteria);
        [puffLogical] = FilterEvents_IOS(EventData.(dataType).LH.whisk,WhiskPuffCriteria);
        combWhiskLogical = logical(whiskLogical.*puffLogical);
        whiskFileIDs = EventData.(dataType).LH.whisk.fileIDs(combWhiskLogical,:);
        whiskEventTimes = EventData.(dataType).LH.whisk.eventTime(combWhiskLogical,:);
        whiskDurations = EventData.(dataType).LH.whisk.duration(combWhiskLogical,:);
        LH_whiskData = EventData.(dataType).LH.whisk.data(combWhiskLogical,:);
        RH_whiskData = EventData.(dataType).RH.whisk.data(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [LH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(LH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        [RH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_whiskData = EventData.(dataType).frontalLH.whisk.data(combWhiskLogical,:);
            fRH_whiskData = EventData.(dataType).frontalRH.whisk.data(combWhiskLogical,:);
            % keep only the data that occurs within the manually-approved alert regions
            [fLH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(fLH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
            [fRH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(fRH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        end
    elseif strcmp(dataType,'GCaMP7s') == true
        [whiskLogical] = FilterEvents_IOS(EventData.(dataType).corLH.whisk,WhiskCriteria);
        [puffLogical] = FilterEvents_IOS(EventData.(dataType).corLH.whisk,WhiskPuffCriteria);
        combWhiskLogical = logical(whiskLogical.*puffLogical);
        whiskFileIDs = EventData.(dataType).corLH.whisk.fileIDs(combWhiskLogical,:);
        whiskEventTimes = EventData.(dataType).corLH.whisk.eventTime(combWhiskLogical,:);
        whiskDurations = EventData.(dataType).corLH.whisk.duration(combWhiskLogical,:);
        LH_WhiskData = EventData.(dataType).corLH.whisk.data(combWhiskLogical,:);
        RH_WhiskData = EventData.(dataType).corRH.whisk.data(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [LH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(LH_WhiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        [RH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(RH_WhiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        fLH_whiskData = EventData.(dataType).corFrontalLH.whisk.data(combWhiskLogical,:);
        fRH_whiskData = EventData.(dataType).corFrontalRH.whisk.data(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [fLH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(fLH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        [fRH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(fRH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    else
        [whiskLogical] = FilterEvents_IOS(EventData.cortical_LH.(dataType).whisk,WhiskCriteria);
        [puffLogical] = FilterEvents_IOS(EventData.cortical_LH.(dataType).whisk,WhiskPuffCriteria);
        combWhiskLogical = logical(whiskLogical.*puffLogical);
        whiskFileIDs = EventData.cortical_LH.(dataType).whisk.fileIDs(combWhiskLogical,:);
        whiskEventTimes = EventData.cortical_LH.(dataType).whisk.eventTime(combWhiskLogical,:);
        whiskDurations = EventData.cortical_LH.(dataType).whisk.duration(combWhiskLogical,:);
        LH_whiskData = EventData.cortical_LH.(dataType).whisk.NormData(combWhiskLogical,:);
        RH_whiskData = EventData.cortical_RH.(dataType).whisk.NormData(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [LH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(LH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        [RH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    end
    clear LH_ProcWhiskData RH_ProcWhiskData
    clear fLH_ProcWhiskData fRH_ProcWhiskData
    % filter, detrend, and take data from whisk onset through 5 seconds
    for gg = 1:size(LH_finalWhiskData,1)
        LH_ProcWhiskData(gg,:) = detrend(filtfilt(sos,g,LH_finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
        RH_ProcWhiskData(gg,:) = detrend(filtfilt(sos,g,RH_finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_ProcWhiskData(gg,:) = detrend(filtfilt(sos,g,fLH_finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
            fRH_ProcWhiskData(gg,:) = detrend(filtfilt(sos,g,fRH_finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
        end
    end
    % analyze correlation coefficient between epochs
    for n = 1:size(LH_ProcWhiskData,1)
        whisk_CC = corrcoef(LH_ProcWhiskData(n,:),RH_ProcWhiskData(n,:));
        whisk_R(n,1) = whisk_CC(2,1);
        if strcmp(groupName,'IOS_GCaMP7s') == true
            whisk_fCC = corrcoef(fLH_ProcWhiskData(n,:),fRH_ProcWhiskData(n,:));
            whisk_fR(n,1) = whisk_fCC(2,1);
        end
    end
    meanWhisk_R = mean(whisk_R);
    stdWhisk_R = std(whisk_R,0,1);
    % save results
    Results_PearsonCorr.(animalID).Whisk.(dataType).R = whisk_R;
    Results_PearsonCorr.(animalID).Whisk.(dataType).meanR = meanWhisk_R;
    Results_PearsonCorr.(animalID).Whisk.(dataType).stdR = stdWhisk_R;
    if strcmp(groupName,'IOS_GCaMP7s') == true
        % analyze correlation coefficient between epochs
        meanWhisk_fR = mean(whisk_fR);
        stdWhisk_fR = std(whisk_fR,0,1);
        % save results
        Results_PearsonCorr.(animalID).Whisk.(dataType).fR = whisk_fR;
        Results_PearsonCorr.(animalID).Whisk.(dataType).meanfR = meanWhisk_fR;
        Results_PearsonCorr.(animalID).Whisk.(dataType).stdfR = stdWhisk_fR;
    end
    %% analyze Pearson's correlation coefficient during periods of alert
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
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144 % 36 bins (180 total) or 3 minutes of aAasleep
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
    if isempty(LH_AlertData) == false
        % filter and detrend data
        for gg = 1:length(LH_AlertData)
            LH_ProcAlertData{gg,1} = detrend(filtfilt(sos,g,LH_AlertData{gg,1}),'constant');
            RH_ProcAlertData{gg,1} = detrend(filtfilt(sos,g,RH_AlertData{gg,1}),'constant');
            if strcmp(groupName,'IOS_GCaMP7s') == true
                fLH_ProcAlertData{gg,1} = detrend(filtfilt(sos,g,fLH_AlertData{gg,1}),'constant');
                fRH_ProcAlertData{gg,1} = detrend(filtfilt(sos,g,fRH_AlertData{gg,1}),'constant');
            end
        end
        % analyze correlation coefficient between epochs
        for n = 1:length(LH_ProcAlertData)
            alert_CC = corrcoef(LH_ProcAlertData{n,1},RH_ProcAlertData{n,1});
            alert_R(n,1) = alert_CC(2,1);
            if strcmp(groupName,'IOS_GCaMP7s') == true
                alert_fCC = corrcoef(fLH_ProcAlertData{n,1},fRH_ProcAlertData{n,1});
                alert_fR(n,1) = alert_fCC(2,1);
            end
        end
        meanAlert_R = mean(alert_R);
        stdAlert_R = std(alert_R,0,1);
        % save results
        Results_PearsonCorr.(animalID).Alert.(dataType).R = alert_R;
        Results_PearsonCorr.(animalID).Alert.(dataType).meanR = meanAlert_R;
        Results_PearsonCorr.(animalID).Alert.(dataType).stdR = stdAlert_R;
        if strcmp(groupName,'IOS_GCaMP7s') == true
            % analyze correlation coefficient between epochs
            meanAlert_fR = mean(alert_fR);
            stdAlert_fR = std(alert_fR,0,1);
            % save results
            Results_PearsonCorr.(animalID).Alert.(dataType).fR = alert_fR;
            Results_PearsonCorr.(animalID).Alert.(dataType).meanfR = meanAlert_fR;
            Results_PearsonCorr.(animalID).Alert.(dataType).stdfR = stdAlert_fR;
        end
    else
        % save results
        Results_PearsonCorr.(animalID).Alert.(dataType).R = [];
        Results_PearsonCorr.(animalID).Alert.(dataType).meanR = [];
        Results_PearsonCorr.(animalID).Alert.(dataType).stdR = [];
        if strcmp(groupName,'IOS_GCaMP7s') == true
            % save results
            Results_PearsonCorr.(animalID).Alert.(dataType).fR = [];
            Results_PearsonCorr.(animalID).Alert.(dataType).meanfR = [];
            Results_PearsonCorr.(animalID).Alert.(dataType).stdfR = [];
        end
    end
    %% analyze Pearson's correlation coefficient during periods of Aasleep
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
    if isempty(LH_AsleepData) == false
        % filter and detrend data
        for gg = 1:length(LH_AsleepData)
            LH_ProcAsleepData{gg,1} = detrend(filtfilt(sos,g,LH_AsleepData{gg,1}),'constant');
            RH_ProcAsleepData{gg,1} = detrend(filtfilt(sos,g,RH_AsleepData{gg,1}),'constant');
            if strcmp(groupName,'IOS_GCaMP7s') == true
                fLH_ProcAsleepData{gg,1} = detrend(filtfilt(sos,g,fLH_AsleepData{gg,1}),'constant');
                fRH_ProcAsleepData{gg,1} = detrend(filtfilt(sos,g,fRH_AsleepData{gg,1}),'constant');
            end
        end
        % analyze correlation coefficient between epochs
        for n = 1:length(LH_ProcAsleepData)
            Asleep_CC = corrcoef(LH_ProcAsleepData{n,1},RH_ProcAsleepData{n,1});
            Asleep_R(n,1) = Asleep_CC(2,1);
            if strcmp(groupName,'IOS_GCaMP7s') == true
                Asleep_fCC = corrcoef(fLH_ProcAsleepData{n,1},fRH_ProcAsleepData{n,1});
                Asleep_fR(n,1) = Asleep_fCC(2,1);
            end
        end
        meanAsleep_R = mean(Asleep_R);
        stdAsleep_R = std(Asleep_R,0,1);
        % save results
        Results_PearsonCorr.(animalID).Asleep.(dataType).R = Asleep_R;
        Results_PearsonCorr.(animalID).Asleep.(dataType).meanR = meanAsleep_R;
        Results_PearsonCorr.(animalID).Asleep.(dataType).stdR = stdAsleep_R;
        if strcmp(groupName,'IOS_GCaMP7s') == true
            % analyze correlation coefficient between epochs
            meanAsleep_fR = mean(Asleep_fR);
            stdAsleep_fR = std(Asleep_fR,0,1);
            % save results
            Results_PearsonCorr.(animalID).Asleep.(dataType).fR = Asleep_fR;
            Results_PearsonCorr.(animalID).Asleep.(dataType).meanfR = meanAsleep_fR;
            Results_PearsonCorr.(animalID).Asleep.(dataType).stdfR = stdAsleep_fR;
        end
    else
        % save results
        Results_PearsonCorr.(animalID).Asleep.(dataType).R = [];
        Results_PearsonCorr.(animalID).Asleep.(dataType).meanR = [];
        Results_PearsonCorr.(animalID).Asleep.(dataType).stdR = [];
        if strcmp(groupName,'IOS_GCaMP7s') == true
            % save results
            Results_PearsonCorr.(animalID).Asleep.(dataType).fR = [];
            Results_PearsonCorr.(animalID).Asleep.(dataType).meanfR = [];
            Results_PearsonCorr.(animalID).Asleep.(dataType).stdfR = [];
        end
    end
    %% analyze Pearson's correlation coefficient during periods of all data
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
    % filter and detrend data
    for gg = 1:length(LH_AllData)
        LH_ProcAllData{gg,1} = detrend(filtfilt(sos,g,LH_AllData{gg,1}),'constant');
        RH_ProcAllData{gg,1} = detrend(filtfilt(sos,g,RH_AllData{gg,1}),'constant');
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_ProcAllData{gg,1} = detrend(filtfilt(sos,g,fLH_AllData{gg,1}),'constant');
            fRH_ProcAllData{gg,1} = detrend(filtfilt(sos,g,fRH_AllData{gg,1}),'constant');
        end
    end
    % analyze correlation coefficient between epochs
    for n = 1:length(LH_ProcAllData)
        all_CC = corrcoef(LH_ProcAllData{n,1},RH_ProcAllData{n,1});
        all_R(n,1) = all_CC(2,1);
        if strcmp(groupName,'IOS_GCaMP7s') == true
            all_fCC = corrcoef(fLH_ProcAllData{n,1},fRH_ProcAllData{n,1});
            all_fR(n,1) = all_fCC(2,1);
        end
    end
    meanAll_R = mean(all_R);
    stdAll_R = std(all_R,0,1);
    % save results
    Results_PearsonCorr.(animalID).All.(dataType).R = all_R;
    Results_PearsonCorr.(animalID).All.(dataType).meanR = meanAll_R;
    Results_PearsonCorr.(animalID).All.(dataType).stdR = stdAll_R;
    if strcmp(groupName,'IOS_GCaMP7s') == true
        % analyze correlation coefficient between epochs
        meanAll_fR = mean(all_fR);
        stdAll_fR = std(all_fR,0,1);
        % save results
        Results_PearsonCorr.(animalID).All.(dataType).fR = all_fR;
        Results_PearsonCorr.(animalID).All.(dataType).meanfR = meanAll_fR;
        Results_PearsonCorr.(animalID).All.(dataType).stdfR = stdAll_fR;
    end
    %% analyze Pearson's correlation coefficient during periods of NREM
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
    % filter, detrend, and truncate data to data to minimum length to match events
    for j = 1:length(LH_nremData)
        LH_nremData{j,1} = detrend(filtfilt(sos,g,LH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
        RH_nremData{j,1} = detrend(filtfilt(sos,g,RH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_nremData{j,1} = detrend(filtfilt(sos,g,fLH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
            fRH_nremData{j,1} = detrend(filtfilt(sos,g,fRH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
        end
    end
    % analyze correlation coefficient between epochs
    for n = 1:length(LH_nremData)
        nrem_CC = corrcoef(LH_nremData{n,1},RH_nremData{n,1});
        nrem_R(n,1) = nrem_CC(2,1);
        if strcmp(groupName,'IOS_GCaMP7s') == true
            nrem_fCC = corrcoef(fLH_nremData{n,1},fRH_nremData{n,1});
            nrem_fR(n,1) = nrem_fCC(2,1);
        end
    end
    meanNREM_R = mean(nrem_R);
    stdNREM_R = std(nrem_R,0,1);
    % save results
    Results_PearsonCorr.(animalID).NREM.(dataType).R = nrem_R;
    Results_PearsonCorr.(animalID).NREM.(dataType).meanR = meanNREM_R;
    Results_PearsonCorr.(animalID).NREM.(dataType).stdR = stdNREM_R;
    if strcmp(groupName,'IOS_GCaMP7s') == true
        % analyze correlation coefficient between epochs
        meanNREM_fR = mean(nrem_fR);
        stdNREM_fR = std(nrem_fR,0,1);
        % save results
        Results_PearsonCorr.(animalID).NREM.(dataType).fR = nrem_fR;
        Results_PearsonCorr.(animalID).NREM.(dataType).meanfR = meanNREM_fR;
        Results_PearsonCorr.(animalID).NREM.(dataType).stdfR = stdNREM_fR;
    end
    %% analyze Pearson's correlation coefficient during periods of REM
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
    % filter, detrend, and truncate data to data to minimum length to match events
    for m = 1:length(LH_remData)
        LH_remData{m,1} = detrend(filtfilt(sos,g,LH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
        RH_remData{m,1} = detrend(filtfilt(sos,g,RH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
        if strcmp(groupName,'IOS_GCaMP7s') == true
            fLH_remData{m,1} = detrend(filtfilt(sos,g,fLH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
            fRH_remData{m,1} = detrend(filtfilt(sos,g,fRH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
        end
    end
    % analyze correlation coefficient between epochs
    for n = 1:length(LH_remData)
        rem_CC = corrcoef(LH_remData{n,1},RH_remData{n,1});
        rem_R(n,1) = rem_CC(2,1);
        if strcmp(groupName,'IOS_GCaMP7s') == true
            rem_fCC = corrcoef(fLH_remData{n,1},fRH_remData{n,1});
            rem_fR(n,1) = rem_fCC(2,1);
        end
    end
    meanREM_R = mean(rem_R);
    stdREM_R = std(rem_R,0,1);
    % save results
    Results_PearsonCorr.(animalID).REM.(dataType).R = rem_R;
    Results_PearsonCorr.(animalID).REM.(dataType).meanR = meanREM_R;
    Results_PearsonCorr.(animalID).REM.(dataType).stdR = stdREM_R;
    if strcmp(groupName,'IOS_GCaMP7s') == true
        % analyze correlation coefficient between epochs
        meanREM_fR = mean(rem_fR);
        stdREM_fR = std(rem_fR,0,1);
        % save results
        Results_PearsonCorr.(animalID).REM.(dataType).fR = rem_fR;
        Results_PearsonCorr.(animalID).REM.(dataType).meanfR = meanREM_fR;
        Results_PearsonCorr.(animalID).REM.(dataType).stdfR = stdREM_fR;
    end
end
% save data
cd(rootFolder)
save('Results_PearsonCorr.mat','Results_PearsonCorr')

end
