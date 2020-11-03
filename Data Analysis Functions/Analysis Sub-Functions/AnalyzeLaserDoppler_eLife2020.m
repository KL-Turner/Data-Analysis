function [AnalysisResults] = AnalyzeLaserDoppler_eLife2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner%
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the laser Doppler flowmetry during different arousal states (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T108','T109','T110','T111','T119','T120','T121','T122','T123'};
modelType = 'Forest';
params.Offset = 2;
params.minTime.Rest = 10;
params.minTime.Whisk = params.Offset + 5;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLocation)
    % find and load RestData.mat struct
    eventDataFileStruct = dir('*_EventData.mat');
    eventDataFile = {eventDataFileStruct.name}';
    eventDataFileID = char(eventDataFile);
    load(eventDataFileID)
    % find and load manual baseline event information
    manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
    manualBaselineFile = {manualBaselineFileStruct.name}';
    manualBaselineFileID = char(manualBaselineFile);
    load(manualBaselineFileID)
    % find and load RestData.mat struct
    restDataFileStruct = dir('*_RestData.mat');
    restDataFile = {restDataFileStruct.name}';
    restDataFileID = char(restDataFile);
    load(restDataFileID)
    % find and load SleepData.mat strut
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
    % find and load RestingBaselines.mat struct
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    % lowpass filter
    samplingRate = EventData.flow.data.whisk.samplingRate;
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    % criteria for whisking
    WhiskCriteria.Fieldname = {'duration','duration','puffDistance'};
    WhiskCriteria.Comparison = {'gt','lt','gt'};
    WhiskCriteria.Value = {2,5,5};
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
    %% analyze LDF during periods of rest
    % pull data from RestData.mat structure
    [restLogical] = FilterEvents_IOS_eLife2020(RestData.flow.data,RestCriteria);
    [puffLogical] = FilterEvents_IOS_eLife2020(RestData.flow.data,RestPuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.flow.data.fileIDs(combRestLogical,:);
    restFlowData = RestData.flow.data.NormData(combRestLogical,:);
    restEventTimes = RestData.flow.data.eventTimes(combRestLogical,:);
    restDurations = RestData.flow.data.durations(combRestLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [finalRestFlowData,~,~,~] = RemoveInvalidData_IOS_eLife2020(restFlowData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    idx = 1;
    % filter LDF
    for gg = 1:length(finalRestFlowData)
        if isempty(finalRestFlowData{gg,1}) == false
            procRestData{idx,1} = filtfilt(sos,g,finalRestFlowData{gg,1}(1:end)); %#ok<*AGROW>
            idx = idx + 1;
        end
    end
    % tkae mean LDF during resting epochs
    for n = 1:length(procRestData)
        restFlowMean(n,1) = mean(procRestData{n,1})*100;
    end
    % save results
    AnalysisResults.(animalID).LDFlow = [];
    AnalysisResults.(animalID).LDFlow.Rest.mean = restFlowMean;
    AnalysisResults.(animalID).LDFlow.Rest.indData = procRestData;
    %% analyze LDF during periods of moderate whisking (2-5 seconds)
    % pull data from EventData.mat structure
    [whiskLogical] = FilterEvents_IOS_eLife2020(EventData.flow.data.whisk,WhiskCriteria);
    [puffLogical] = FilterEvents_IOS_eLife2020(EventData.flow.data.whisk,WhiskPuffCriteria);
    combWhiskLogical = logical(whiskLogical.*puffLogical);
    whiskFlowData = EventData.flow.data.whisk.NormData(combWhiskLogical,:);
    whiskFileIDs = EventData.flow.data.whisk.fileIDs(combWhiskLogical,:);
    whiskEventTimes = EventData.flow.data.whisk.eventTime(combWhiskLogical,:);
    whiskDurations = EventData.flow.data.whisk.duration(combWhiskLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [finalWhiskData,~,~,~] = RemoveInvalidData_IOS_eLife2020(whiskFlowData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    % filter LDF and mean-subtract 2-seconds prior to whisk
    for gg = 1:size(finalWhiskData,1)
        procWhiskDataA = filtfilt(sos,g,finalWhiskData(gg,:));
        procWhiskDataB = procWhiskDataA - mean(procWhiskDataA(1:params.Offset*samplingRate));
        procWhiskDataC{gg,1} = procWhiskDataB(params.Offset*samplingRate:params.minTime.Whisk*samplingRate)*100;
    end
    % take mean LDF during whisking epochs from onset through 5 seconds
    for n = 1:length(procWhiskDataC)
        whiskFlowMean(n,1) = mean(procWhiskDataC{n,1});
    end
    % save results
    AnalysisResults.(animalID).LDFlow.Whisk.mean = whiskFlowMean;
    AnalysisResults.(animalID).LDFlow.Whisk.indData = procWhiskDataC;
    %% analyze LDF during periods of NREM sleep
    % pull data from SleepData.mat structure
    [nremData,~,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).NREM.data.DopplerFlow,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    % filter LDF and take mean during NREM epochs
    idx = 1;
    for n = 1:length(nremData)
        if sum(isnan(nremData{n,1})) == 0
            nremFlowMean(idx,1) = mean(filtfilt(sos,g,nremData{n,1}))*100;
            nremFlowInd{idx,1} = filtfilt(sos,g,nremData{n,1})*100;
            idx = idx + 1;
        end
    end
    % save results
    AnalysisResults.(animalID).LDFlow.NREM.mean = nremFlowMean;
    AnalysisResults.(animalID).LDFlow.NREM.indData = nremFlowInd;
    %% analyze LDF during periods of REM sleep
    % pull data from SleepData.mat structure
    [remData,~,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).REM.data.DopplerFlow,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    % filter LDF and take mean during NREM epochs
    idx = 1;
    for n = 1:length(remData)
        if sum(isnan(remData{n,1})) == 0
            remFlowMean(idx,1) = mean(filtfilt(sos,g,remData{n,1}))*100;
            remFlowInd{idx,1} = filtfilt(sos,g,remData{n,1})*100;
            idx = idx + 1;
        end
    end
    % save results
    AnalysisResults.(animalID).LDFlow.REM.mean = remFlowMean;
    AnalysisResults.(animalID).LDFlow.REM.indData = remFlowInd; 
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
