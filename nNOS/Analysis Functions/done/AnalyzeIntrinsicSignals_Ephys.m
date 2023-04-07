function [Results_IntSig_Ephys] = AnalyzeIntrinsicSignals_Ephys(animalID,group,set,rootFolder,delim,Results_IntSig_Ephys)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
modelType = 'Forest';
params.minTime.Rest = 10;
params.Offset = 2;
params.minTime.Whisk = params.Offset + 5;
params.minTime.Stim = params.Offset + 2;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% only run analysis for valid animal IDs
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Bilateral Imaging'];
cd(dataLocation)
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
% find and load SleepData.mat strut
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
% criteria for the whisking
WhiskCriteria.Fieldname = {'duration','puffDistance'};
WhiskCriteria.Comparison = {'gt','gt'};
WhiskCriteria.Value = {5,5};
WhiskPuffCriteria.Fieldname = {'puffDistance'};
WhiskPuffCriteria.Comparison = {'gt'};
WhiskPuffCriteria.Value = {5};
% criteria for the resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
% criteria for the stimulation
StimCriteriaA.Value = {'RPadSol'};
StimCriteriaA.Fieldname = {'solenoidName'};
StimCriteriaA.Comparison = {'equal'};
StimCriteriaB.Value = {'LPadSol'};
StimCriteriaB.Fieldname = {'solenoidName'};
StimCriteriaB.Comparison = {'equal'};
% loop variables
hemispheres = {'LH','RH'};
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    % lowpass filter
    samplingRate = RestData.HbT.(hemisphere).CBVCamSamplingRate;
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    %% rest
    [restLogical] = FilterEvents_IOS(RestData.HbT.(hemisphere),RestCriteria);
    [puffLogical] = FilterEvents_IOS(RestData.HbT.(hemisphere),RestPuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.HbT.(hemisphere).fileIDs(combRestLogical,:);
    restEventTimes = RestData.HbT.(hemisphere).eventTimes(combRestLogical,:);
    restDurations = RestData.HbT.(hemisphere).durations(combRestLogical,:);
    restData = RestData.HbT.(hemisphere).data(combRestLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [finalRestData,~,~,~] = RemoveInvalidData_IOS(restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    % filter and average
    for gg = 1:length(finalRestData)
        procRestData{gg,1} = filtfilt(sos,g,finalRestData{gg,1});
        restMean(gg,1) = mean(procRestData{gg,1}(1:end));
    end
    % save results
    Results_IntSig_Ephys.(group).(animalID).(hemisphere).Rest.HbT = restMean;
    %% whisk
    [whiskLogical] = FilterEvents_IOS(EventData.HbT.(hemisphere).whisk,WhiskCriteria);
    [puffLogical] = FilterEvents_IOS(EventData.HbT.(hemisphere).whisk,WhiskPuffCriteria);
    combWhiskLogical = logical(whiskLogical.*puffLogical);
    whiskFileIDs = EventData.HbT.(hemisphere).whisk.fileIDs(combWhiskLogical,:);
    whiskEventTimes = EventData.HbT.(hemisphere).whisk.eventTime(combWhiskLogical,:);
    whiskDurations = EventData.HbT.(hemisphere).whisk.duration(combWhiskLogical,:);
    whiskData = EventData.HbT.(hemisphere).whisk.data(combWhiskLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [finalWhiskData,~,~,~] = RemoveInvalidData_IOS(whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    % filter and average
    for gg = 1:size(finalWhiskData,1)
        procWhiskData_temp = filtfilt(sos,g,finalWhiskData(gg,:));
        procWhiskData(gg,:) = procWhiskData_temp - mean(procWhiskData_temp(1:params.Offset*samplingRate));
        whiskMean{gg,1} = mean(procWhiskData(gg,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
    end
    % save results
    Results_IntSig_Ephys.(group).(animalID).(hemisphere).Whisk.HbT = cell2mat(whiskMean);
    %% stim
    if any(strcmp(hemisphere,{'LH'})) == true
        StimCriteria = StimCriteriaA;
    elseif any(strcmp(hemisphere,{'RH'})) == true
        StimCriteria = StimCriteriaB;
    end
    stimFilter = FilterEvents_IOS(EventData.HbT.(hemisphere).stim,StimCriteria);
    [stimFileIDs] = EventData.HbT.(hemisphere).stim.fileIDs(stimFilter,:);
    [stimEventTimes] = EventData.HbT.(hemisphere).stim.eventTime(stimFilter,:);
    stimDurations = zeros(length(stimEventTimes),1);
    [stimData] = EventData.HbT.(hemisphere).stim.data(stimFilter,:);
    % keep only the data that occurs within the manually-approved awake regions
    [finalStimData,~,~,~] = RemoveInvalidData_IOS(stimData,stimFileIDs,stimDurations,stimEventTimes,ManualDecisions);
    % filter and average
    for gg = 1:size(finalStimData,1)
        procStimData_temp = filtfilt(sos,g,finalStimData(gg,:));
        procStimData(gg,:) = procStimData_temp - mean(procStimData_temp(1:params.Offset*samplingRate));
        stimMean{gg,1} = mean(procStimData(gg,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
    end
    % save results
    Results_IntSig_Ephys.(group).(animalID).(hemisphere).Stim.HbT = cell2mat(stimMean);
    %% NREM
    [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.HbT.(hemisphere),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    % filter and average
    for nn = 1:length(nremData)
        nremMean(nn,1) = mean(filtfilt(sos,g,nremData{nn,1}(1:end)));
    end
    % save results
    Results_IntSig_Ephys.(group).(animalID).(hemisphere).NREM.HbT = nremMean;
    %% REM
    [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.HbT.(hemisphere),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    % filter and average
    for nn = 1:length(remData)
        remMean(nn,1) = mean(filtfilt(sos,g,remData{nn,1}(1:end)));
    end
    % save results
    Results_IntSig_Ephys.(group).(animalID).(hemisphere).REM.HbT = remMean;
    %% isolfurane
    isoDataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Isoflurane Trials'];
    cd(isoDataLocation)
    try
        % pull ProcData.mat file associated with isoflurane administration
        procDataFileStruct = dir('*_ProcData.mat');
        procDataFile = {procDataFileStruct.name}';
        procDataFileID = char(procDataFile);
        load(procDataFileID,'-mat')
        % filter and average
        isoData = ProcData.data.CBV_HbT.(hemisphere)((end - samplingRate*100):end);
        filtIsoData = filtfilt(sos,g,isoData);
        % save results
        Results_IntSig_Ephys.(group).(animalID).(hemisphere).Iso.HbT = mean(filtIsoData);
    catch
        % save results
        Results_IntSig_Ephys.(group).(animalID).(hemisphere).Iso.HbT = [];
    end
    cd(dataLocation)
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_IntSig_Ephys.mat','Results_IntSig_Ephys')
cd([rootFolder delim 'Data'])