function [Results_PearsonCorrEphys] = AnalyzePearsonCorrelation_Ephys(animalID,group,rootFolder,delim,Results_PearsonCorrEphys)
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
dataTypes = {'CBV_HbT','gammaBandPower'};
% go through each valid data type for arousal-based correlation analysis
for a = 1:length(dataTypes)
    dataType = dataTypes{1,a};
    %% analyze Pearson's correlation coefficient during periods of rest
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
    % lowpass filter
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    % filter, detrend, and truncate data to minimum length to match events
    for gg = 1:length(LH_finalRestData)
        LH_ProcRestData{gg,1} = detrend(filtfilt(sos,g,LH_finalRestData{gg,1}(1:params.minTime.Rest*samplingRate)),'constant');
        RH_ProcRestData{gg,1} = detrend(filtfilt(sos,g,RH_finalRestData{gg,1}(1:params.minTime.Rest*samplingRate)),'constant');
    end
    % analyze correlation coefficient of resting epochs
    for n = 1:length(LH_ProcRestData)
        rest_CC = corrcoef(LH_ProcRestData{n,1},RH_ProcRestData{n,1});
        rest_R(n,1) = rest_CC(2,1);
    end
    meanRest_R = mean(rest_R);
    stdRest_R = std(rest_R,0,1);
    % save results
    Results_PearsonCorrEphys.(animalID).Rest.(dataType).R = rest_R;
    Results_PearsonCorrEphys.(animalID).Rest.(dataType).meanR = meanRest_R;
    Results_PearsonCorrEphys.(animalID).Rest.(dataType).stdR = stdRest_R;
    %% analyze Pearson's correlation coefficient during periods of moderate whisking (2-5 seconds)
    if strcmp(dataType,'CBV_HbT') == true
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
    % filter, detrend, and take data from whisk onset through 5 seconds
    for gg = 1:size(LH_finalWhiskData,1)
        LH_ProcWhiskData(gg,:) = detrend(filtfilt(sos,g,LH_finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
        RH_ProcWhiskData(gg,:) = detrend(filtfilt(sos,g,RH_finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
    end
    % analyze correlation coefficient between epochs
    for n = 1:size(LH_ProcWhiskData,1)
        whisk_CC = corrcoef(LH_ProcWhiskData(n,:),RH_ProcWhiskData(n,:));
        whisk_R(n,1) = whisk_CC(2,1);
    end
    meanWhisk_R = mean(whisk_R);
    stdWhisk_R = std(whisk_R,0,1);
    % save results
    Results_PearsonCorrEphys.(animalID).Whisk.(dataType).R = whisk_R;
    Results_PearsonCorrEphys.(animalID).Whisk.(dataType).meanR = meanWhisk_R;
    Results_PearsonCorrEphys.(animalID).Whisk.(dataType).stdR = stdWhisk_R;
    %% analyze Pearson's correlation coefficient during periods of alert
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
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144 % 36 bins (180 total) or 3 minutes of aAasleep
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
    if isempty(LH_AlertData) == false
        % filter and detrend data
        for gg = 1:length(LH_AlertData)
            LH_ProcAlertData{gg,1} = detrend(filtfilt(sos,g,LH_AlertData{gg,1}),'constant');
            RH_ProcAlertData{gg,1} = detrend(filtfilt(sos,g,RH_AlertData{gg,1}),'constant');
        end
        % analyze correlation coefficient between epochs
        for n = 1:length(LH_ProcAlertData)
            alert_CC = corrcoef(LH_ProcAlertData{n,1},RH_ProcAlertData{n,1});
            alert_R(n,1) = alert_CC(2,1);
        end
        meanAlert_R = mean(alert_R);
        stdAlert_R = std(alert_R,0,1);
        % save results
        Results_PearsonCorrEphys.(animalID).Alert.(dataType).R = alert_R;
        Results_PearsonCorrEphys.(animalID).Alert.(dataType).meanR = meanAlert_R;
        Results_PearsonCorrEphys.(animalID).Alert.(dataType).stdR = stdAlert_R;
    else
        % save results
        Results_PearsonCorrEphys.(animalID).Alert.(dataType).R = [];
        Results_PearsonCorrEphys.(animalID).Alert.(dataType).meanR = [];
        Results_PearsonCorrEphys.(animalID).Alert.(dataType).stdR = [];
    end
    %% analyze Pearson's correlation coefficient during periods of Aasleep
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
        if sum(strcmp(scoringLabels,'Not Sleep')) < 36 % 36 bins (180 total) or 3 minutes of alert
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
    if isempty(LH_AsleepData) == false
        % filter and detrend data
        for gg = 1:length(LH_AsleepData)
            LH_ProcAsleepData{gg,1} = detrend(filtfilt(sos,g,LH_AsleepData{gg,1}),'constant');
            RH_ProcAsleepData{gg,1} = detrend(filtfilt(sos,g,RH_AsleepData{gg,1}),'constant');
        end
        % analyze correlation coefficient between epochs
        for n = 1:length(LH_ProcAsleepData)
            Asleep_CC = corrcoef(LH_ProcAsleepData{n,1},RH_ProcAsleepData{n,1});
            Asleep_R(n,1) = Asleep_CC(2,1);
        end
        meanAsleep_R = mean(Asleep_R);
        stdAsleep_R = std(Asleep_R,0,1);
        % save results
        Results_PearsonCorrEphys.(animalID).Asleep.(dataType).R = Asleep_R;
        Results_PearsonCorrEphys.(animalID).Asleep.(dataType).meanR = meanAsleep_R;
        Results_PearsonCorrEphys.(animalID).Asleep.(dataType).stdR = stdAsleep_R;
    else
        % save results
        Results_PearsonCorrEphys.(animalID).Asleep.(dataType).R = [];
        Results_PearsonCorrEphys.(animalID).Asleep.(dataType).meanR = [];
        Results_PearsonCorrEphys.(animalID).Asleep.(dataType).stdR = [];
    end
    %% analyze Pearson's correlation coefficient during periods of all data
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
    % filter and detrend data
    for gg = 1:length(LH_AllData)
        LH_ProcAllData{gg,1} = detrend(filtfilt(sos,g,LH_AllData{gg,1}),'constant');
        RH_ProcAllData{gg,1} = detrend(filtfilt(sos,g,RH_AllData{gg,1}),'constant');
    end
    % analyze correlation coefficient between epochs
    for n = 1:length(LH_ProcAllData)
        all_CC = corrcoef(LH_ProcAllData{n,1},RH_ProcAllData{n,1});
        all_R(n,1) = all_CC(2,1);
    end
    meanAll_R = mean(all_R);
    stdAll_R = std(all_R,0,1);
    % save results
    Results_PearsonCorrEphys.(animalID).All.(dataType).R = all_R;
    Results_PearsonCorrEphys.(animalID).All.(dataType).meanR = meanAll_R;
    Results_PearsonCorrEphys.(animalID).All.(dataType).stdR = stdAll_R;
    %% analyze Pearson's correlation coefficient during periods of NREM
    if strcmp(dataType,'CBV_HbT') == true
        [LH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    else
        [LH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_LH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    end
    % filter, detrend, and truncate data to data to minimum length to match events
    for j = 1:length(LH_nremData)
        LH_nremData{j,1} = detrend(filtfilt(sos,g,LH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
        RH_nremData{j,1} = detrend(filtfilt(sos,g,RH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
    end
    % analyze correlation coefficient between epochs
    for n = 1:length(LH_nremData)
        nrem_CC = corrcoef(LH_nremData{n,1},RH_nremData{n,1});
        nrem_R(n,1) = nrem_CC(2,1);
    end
    meanNREM_R = mean(nrem_R);
    stdNREM_R = std(nrem_R,0,1);
    % save results
    Results_PearsonCorrEphys.(animalID).NREM.(dataType).R = nrem_R;
    Results_PearsonCorrEphys.(animalID).NREM.(dataType).meanR = meanNREM_R;
    Results_PearsonCorrEphys.(animalID).NREM.(dataType).stdR = stdNREM_R;
    %% analyze Pearson's correlation coefficient during periods of REM
    if strcmp(dataType,'CBV_HbT') == true
        [LH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    else
        [LH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_LH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    end
    % filter, detrend, and truncate data to data to minimum length to match events
    for m = 1:length(LH_remData)
        LH_remData{m,1} = detrend(filtfilt(sos,g,LH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
        RH_remData{m,1} = detrend(filtfilt(sos,g,RH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
    end
    % analyze correlation coefficient between epochs
    for n = 1:length(LH_remData)
        rem_CC = corrcoef(LH_remData{n,1},RH_remData{n,1});
        rem_R(n,1) = rem_CC(2,1);
    end
    meanREM_R = mean(rem_R);
    stdREM_R = std(rem_R,0,1);
    % save results
    Results_PearsonCorrEphys.(animalID).REM.(dataType).R = rem_R;
    Results_PearsonCorrEphys.(animalID).REM.(dataType).meanR = meanREM_R;
    Results_PearsonCorrEphys.(animalID).REM.(dataType).stdR = stdREM_R;
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_PearsonCorrEphys.mat','Results_PearsonCorrEphys')
cd([rootFolder delim 'Data'])

end
