function [Results_BehavHbT] = AnalyzeBehavioralHbT(animalID,group,rootFolder,delim,Results_BehavHbT)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the hemodynamic signal [HbT] during different behavioral states (IOS)
%________________________________________________________________________________________________________________________

% function parameters
modelType = 'Forest';
params.minTime.Rest = 10;
params.Offset = 2;
params.minTime.Whisk = params.Offset + 5;
params.minTime.Stim = params.Offset + 2;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% only run analysis for valid animal IDs
dataLocation = [rootFolder delim group delim animalID delim 'Bilateral Imaging'];
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
% determine data types based on experimental group
strBreaks = strfind(group,delim);
groupName = group(strBreaks + 1:end);
% lowpass filter
samplingRate = RestData.CBV_HbT.LH.CBVCamSamplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
%% analyze [HbT] during periods of rest
% pull data from RestData.mat structure
[restLogical] = FilterEvents_IOS(RestData.CBV_HbT.LH,RestCriteria);
[puffLogical] = FilterEvents_IOS(RestData.CBV_HbT.LH,RestPuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFileIDs = RestData.CBV_HbT.LH.fileIDs(combRestLogical,:);
restEventTimes = RestData.CBV_HbT.LH.eventTimes(combRestLogical,:);
restDurations = RestData.CBV_HbT.LH.durations(combRestLogical,:);
LH_RestingData = RestData.CBV_HbT.LH.data(combRestLogical,:);
RH_RestingData = RestData.CBV_HbT.RH.data(combRestLogical,:);
if strcmp(groupName,'IOS GCaMP7s') == true
    fLH_RestingData = RestData.CBV_HbT.frontalLH.data(combRestLogical,:);
    fRH_RestingData = RestData.CBV_HbT.frontalRH.data(combRestLogical,:);
end
% keep only the data that occurs within the manually-approved awake regions
[LH_finalRestData,finalRestFileIDs,~,~] = RemoveInvalidData_IOS(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
[RH_finalRestData,~,~,~] = RemoveInvalidData_IOS(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
if strcmp(groupName,'IOS GCaMP7s') == true
    [fLH_finalRestData,finalRestFileIDs,~,~] = RemoveInvalidData_IOS(fLH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [fRH_finalRestData,~,~,~] = RemoveInvalidData_IOS(fRH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
end
% filter [HbT]
for gg = 1:length(LH_finalRestData)
    LH_ProcRestData{gg,1} = filtfilt(sos,g,LH_finalRestData{gg,1});
    RH_ProcRestData{gg,1} = filtfilt(sos,g,RH_finalRestData{gg,1});
    if strcmp(groupName,'IOS GCaMP7s') == true
        fLH_ProcRestData{gg,1} = filtfilt(sos,g,fLH_finalRestData{gg,1});
        fRH_ProcRestData{gg,1} = filtfilt(sos,g,fRH_finalRestData{gg,1});
    end
end
% take mean [HbT] during resting epochs
for nn = 1:length(LH_ProcRestData)
    LH_restCBVMean(nn,1) = mean(LH_ProcRestData{nn,1}(1:end));
    RH_restCBVMean(nn,1) = mean(RH_ProcRestData{nn,1}(1:end));
    if strcmp(groupName,'IOS GCaMP7s') == true
        fLH_restCBVMean(nn,1) = mean(fLH_ProcRestData{nn,1}(1:end));
        fRH_restCBVMean(nn,1) = mean(fRH_ProcRestData{nn,1}(1:end));
    end
end
% save results
Results_BehavHbT.(animalID).Rest.FileIDs = finalRestFileIDs;
Results_BehavHbT.(animalID).Rest.MeanLH = LH_restCBVMean;
Results_BehavHbT.(animalID).Rest.MeanRH = RH_restCBVMean;
Results_BehavHbT.(animalID).Rest.IndLH = LH_ProcRestData;
Results_BehavHbT.(animalID).Rest.IndRH = RH_ProcRestData;
if strcmp(groupName,'IOS GCaMP7s') == true
    % save results
    Results_BehavHbT.(animalID).Rest.MeanfLH = fLH_restCBVMean;
    Results_BehavHbT.(animalID).Rest.MeanfRH = fRH_restCBVMean;
    Results_BehavHbT.(animalID).Rest.IndfLH = fLH_ProcRestData;
    Results_BehavHbT.(animalID).Rest.IndfRH = fRH_ProcRestData;
end
%% analyze [HbT] during periods of moderate whisking (2-5 seconds)
% pull data from EventData.mat structure
[whiskLogical] = FilterEvents_IOS(EventData.CBV_HbT.LH.whisk,WhiskCriteria);
[puffLogical] = FilterEvents_IOS(EventData.CBV_HbT.LH.whisk,WhiskPuffCriteria);
combWhiskLogical = logical(whiskLogical.*puffLogical);
whiskFileIDs = EventData.CBV_HbT.LH.whisk.fileIDs(combWhiskLogical,:);
whiskEventTimes = EventData.CBV_HbT.LH.whisk.eventTime(combWhiskLogical,:);
whiskDurations = EventData.CBV_HbT.LH.whisk.duration(combWhiskLogical,:);
LH_whiskData = EventData.CBV_HbT.LH.whisk.data(combWhiskLogical,:);
RH_whiskData = EventData.CBV_HbT.RH.whisk.data(combWhiskLogical,:);
if strcmp(groupName,'IOS GCaMP7s') == true
    fLH_whiskData = EventData.CBV_HbT.frontalLH.whisk.data(combWhiskLogical,:);
    fRH_whiskData = EventData.CBV_HbT.frontalRH.whisk.data(combWhiskLogical,:);
end
% keep only the data that occurs within the manually-approved awake regions
[LH_finalWhiskData,finalWhiskFileIDs,~,~] = RemoveInvalidData_IOS(LH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
[RH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
if strcmp(groupName,'IOS GCaMP7s') == true
    [fLH_finalWhiskData,finalWhiskFileIDs,~,~] = RemoveInvalidData_IOS(fLH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    [fRH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(fRH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
end
% filter [HbT] and mean-subtract 2 seconds prior to whisk
for gg = 1:size(LH_finalWhiskData,1)
    LH_ProcWhiskData_temp = filtfilt(sos,g,LH_finalWhiskData(gg,:));
    LH_ProcWhiskData(gg,:) = LH_ProcWhiskData_temp - mean(LH_ProcWhiskData_temp(1:params.Offset*samplingRate));
    RH_ProcWhiskData_temp = filtfilt(sos,g,RH_finalWhiskData(gg,:));
    RH_ProcWhiskData(gg,:) = RH_ProcWhiskData_temp - mean(RH_ProcWhiskData_temp(1:params.Offset*samplingRate));
    if strcmp(groupName,'IOS GCaMP7s') == true
        fLH_ProcWhiskData_temp = filtfilt(sos,g,fLH_finalWhiskData(gg,:));
        fLH_ProcWhiskData(gg,:) = fLH_ProcWhiskData_temp - mean(fLH_ProcWhiskData_temp(1:params.Offset*samplingRate));
        fRH_ProcWhiskData_temp = filtfilt(sos,g,fRH_finalWhiskData(gg,:));
        fRH_ProcWhiskData(gg,:) = fRH_ProcWhiskData_temp - mean(fRH_ProcWhiskData_temp(1:params.Offset*samplingRate));
    end
end
% take mean [HbT] during whisking epochs from onset through 5 seconds
for nn = 1:size(LH_ProcWhiskData,1)
    LH_whiskCBVMean{nn,1} = mean(LH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
    RH_whiskCBVMean{nn,1} = mean(RH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
    LH_whiskCBV{nn,1} = LH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
    RH_whiskCBV{nn,1} = RH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
    if strcmp(groupName,'IOS GCaMP7s') == true
        fLH_whiskCBVMean{nn,1} = mean(fLH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
        fRH_whiskCBVMean{nn,1} = mean(fRH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
        fLH_whiskCBV{nn,1} = fLH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
        fRH_whiskCBV{nn,1} = fRH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
    end
end
% save results
Results_BehavHbT.(animalID).Whisk.FileIDs = finalWhiskFileIDs;
Results_BehavHbT.(animalID).Whisk.MeanLH = cell2mat(LH_whiskCBVMean);
Results_BehavHbT.(animalID).Whisk.MeanRH = cell2mat(RH_whiskCBVMean);
Results_BehavHbT.(animalID).Whisk.IndLH = LH_whiskCBV;
Results_BehavHbT.(animalID).Whisk.IndRH = RH_whiskCBV;
if strcmp(groupName,'IOS GCaMP7s') == true
    Results_BehavHbT.(animalID).Whisk.MeanfLH = cell2mat(fLH_whiskCBVMean);
    Results_BehavHbT.(animalID).Whisk.MeanfRH = cell2mat(fRH_whiskCBVMean);
    Results_BehavHbT.(animalID).Whisk.IndfLH = fLH_whiskCBV;
    Results_BehavHbT.(animalID).Whisk.IndfRH = fRH_whiskCBV;
end
%% analyze [HbT] during periods of stimulation
% pull data from EventData.mat structure
LH_stimFilter = FilterEvents_IOS(EventData.CBV_HbT.LH.stim,StimCriteriaA);
RH_stimFilter = FilterEvents_IOS(EventData.CBV_HbT.RH.stim,StimCriteriaB);
[LH_stimFileIDs] = EventData.CBV_HbT.LH.stim.fileIDs(LH_stimFilter,:);
[RH_stimFileIDs] = EventData.CBV_HbT.RH.stim.fileIDs(RH_stimFilter,:);
[LH_stimEventTimes] = EventData.CBV_HbT.LH.stim.eventTime(LH_stimFilter,:);
[RH_stimEventTimes] = EventData.CBV_HbT.RH.stim.eventTime(RH_stimFilter,:);
LH_stimDurations = zeros(length(LH_stimEventTimes),1);
RH_stimDurations = zeros(length(RH_stimEventTimes),1);
[LH_stimHbTData] = EventData.CBV_HbT.LH.stim.data(LH_stimFilter,:);
[RH_stimHbTData] = EventData.CBV_HbT.RH.stim.data(RH_stimFilter,:);
if strcmp(groupName,'IOS GCaMP7s') == true
    [fLH_stimHbTData] = EventData.CBV_HbT.frontalLH.stim.data(LH_stimFilter,:);
    [fRH_stimHbTData] = EventData.CBV_HbT.frontalRH.stim.data(RH_stimFilter,:);
end
% keep only the data that occurs within the manually-approved awake regions
[LH_finalStimData,LH_finalStimFileIDs,~,~] = RemoveInvalidData_IOS(LH_stimHbTData,LH_stimFileIDs,LH_stimDurations,LH_stimEventTimes,ManualDecisions);
[RH_finalStimData,RH_finalStimFileIDs,~,~] = RemoveInvalidData_IOS(RH_stimHbTData,RH_stimFileIDs,RH_stimDurations,RH_stimEventTimes,ManualDecisions);
if strcmp(groupName,'IOS GCaMP7s') == true
    [fLH_finalStimData,~,~,~] = RemoveInvalidData_IOS(fLH_stimHbTData,LH_stimFileIDs,LH_stimDurations,LH_stimEventTimes,ManualDecisions);
    [fRH_finalStimData,~,~,~] = RemoveInvalidData_IOS(fRH_stimHbTData,RH_stimFileIDs,RH_stimDurations,RH_stimEventTimes,ManualDecisions);
end
% filter [HbT] and mean-subtract 2 seconds prior to stimulus (left hem)
for gg = 1:size(LH_finalStimData,1)
    LH_ProcStimData_temp = filtfilt(sos,g,LH_finalStimData(gg,:));
    LH_ProcStimData(gg,:) = LH_ProcStimData_temp - mean(LH_ProcStimData_temp(1:params.Offset*samplingRate));
    if strcmp(groupName,'IOS GCaMP7s') == true
        fLH_ProcStimData_temp = filtfilt(sos,g,fLH_finalStimData(gg,:));
        fLH_ProcStimData(gg,:) = fLH_ProcStimData_temp - mean(fLH_ProcStimData_temp(1:params.Offset*samplingRate));
    end
end
% filter [HbT] and mean-subtract 2 seconds prior to stimulus (right hem)
for gg = 1:size(RH_finalStimData,1)
    RH_ProcStimData_temp = filtfilt(sos,g,RH_finalStimData(gg,:));
    RH_ProcStimData(gg,:) = RH_ProcStimData_temp - mean(RH_ProcStimData_temp(1:params.Offset*samplingRate));
    if strcmp(groupName,'IOS GCaMP7s') == true
        fRH_ProcStimData_temp = filtfilt(sos,g,fRH_finalStimData(gg,:));
        fRH_ProcStimData(gg,:) = fRH_ProcStimData_temp - mean(fRH_ProcStimData_temp(1:params.Offset*samplingRate));
    end
end
% take mean [HbT] 1-2 seconds after stimulation (left hem)
for nn = 1:size(LH_ProcStimData,1)
    LH_stimCBVMean{nn,1} = mean(LH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
    LH_stimCBV{nn,1} = LH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
    if strcmp(groupName,'IOS GCaMP7s') == true
        fLH_stimCBVMean{nn,1} = mean(fLH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
        fLH_stimCBV{nn,1} = fLH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
    end
end
% take mean [HbT] 1-2 seconds after stimulation (right hem)
for nn = 1:size(RH_ProcStimData,1)
    RH_stimCBVMean{nn,1} = mean(RH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
    RH_stimCBV{nn,1} = RH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
    if strcmp(groupName,'IOS GCaMP7s') == true
        fRH_stimCBVMean{nn,1} = mean(fRH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
        fRH_stimCBV{nn,1} = fRH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
    end
end
% save results
Results_BehavHbT.(animalID).Stim.LH_FileIDs = LH_finalStimFileIDs;
Results_BehavHbT.(animalID).Stim.RH_FileIDs = RH_finalStimFileIDs;
Results_BehavHbT.(animalID).Stim.MeanLH = cell2mat(LH_stimCBVMean);
Results_BehavHbT.(animalID).Stim.MeanRH = cell2mat(RH_stimCBVMean);
Results_BehavHbT.(animalID).Stim.IndLH = LH_stimCBV;
Results_BehavHbT.(animalID).Stim.IndRH = RH_stimCBV;
if strcmp(groupName,'IOS GCaMP7s') == true
    Results_BehavHbT.(animalID).Stim.MeanfLH = cell2mat(fLH_stimCBVMean);
    Results_BehavHbT.(animalID).Stim.MeanfRH = cell2mat(fRH_stimCBVMean);
    Results_BehavHbT.(animalID).Stim.IndfLH = fLH_stimCBV;
    Results_BehavHbT.(animalID).Stim.IndfRH = fRH_stimCBV;
end
%% analyze [HbT] during periods of NREM sleep
% pull data from SleepData.mat structure
[LH_nremData,nremFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV_HbT.LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
[RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV_HbT.RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
if strcmp(groupName,'IOS GCaMP7s') == true
    [fLH_nremData,nremFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV_HbT.frontalLH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [fRH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV_HbT.frontalRH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
end
% filter and take mean [HbT] during NREM epochs
for nn = 1:length(LH_nremData)
    LH_nremCBVMean(nn,1) = mean(filtfilt(sos,g,LH_nremData{nn,1}(1:end)));
    RH_nremCBVMean(nn,1) = mean(filtfilt(sos,g,RH_nremData{nn,1}(1:end)));
    if strcmp(groupName,'IOS GCaMP7s') == true
        fLH_nremCBVMean(nn,1) = mean(filtfilt(sos,g,fLH_nremData{nn,1}(1:end)));
        fRH_nremCBVMean(nn,1) = mean(filtfilt(sos,g,fRH_nremData{nn,1}(1:end)));
    end
end
% save results
Results_BehavHbT.(animalID).NREM.FileIDs = nremFileIDs;
Results_BehavHbT.(animalID).NREM.MeanLH = LH_nremCBVMean;
Results_BehavHbT.(animalID).NREM.MeanRH = RH_nremCBVMean;
Results_BehavHbT.(animalID).NREM.IndLH = LH_nremData;
Results_BehavHbT.(animalID).NREM.IndRH = RH_nremData;
if strcmp(groupName,'IOS GCaMP7s') == true
    Results_BehavHbT.(animalID).NREM.MeanfLH = fLH_nremCBVMean;
    Results_BehavHbT.(animalID).NREM.MeanfRH = fRH_nremCBVMean;
    Results_BehavHbT.(animalID).NREM.IndfLH = fLH_nremData;
    Results_BehavHbT.(animalID).NREM.IndfRH = fRH_nremData;
end
%% analyze [HbT] during periods of REM sleep
% pull data from SleepData.mat structure
[LH_remData,remFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV_HbT.LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
[RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV_HbT.RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
if strcmp(groupName,'IOS GCaMP7s') == true
    [fLH_remData,remFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV_HbT.frontalLH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [fRH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV_HbT.frontalRH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
end
% filter and take mean [HbT] during REM epochs
for nn = 1:length(LH_remData)
    LH_remCBVMean(nn,1) = mean(filtfilt(sos,g,LH_remData{nn,1}(1:end)));
    RH_remCBVMean(nn,1) = mean(filtfilt(sos,g,RH_remData{nn,1}(1:end)));
    if strcmp(groupName,'IOS GCaMP7s') == true
        fLH_remCBVMean(nn,1) = mean(filtfilt(sos,g,fLH_remData{nn,1}(1:end)));
        fRH_remCBVMean(nn,1) = mean(filtfilt(sos,g,fRH_remData{nn,1}(1:end)));
    end
end
% save results
Results_BehavHbT.(animalID).REM.FileIDs = remFileIDs;
Results_BehavHbT.(animalID).REM.MeanLH = LH_remCBVMean;
Results_BehavHbT.(animalID).REM.MeanRH = RH_remCBVMean;
Results_BehavHbT.(animalID).REM.IndLH = LH_remData;
Results_BehavHbT.(animalID).REM.IndRH = RH_remData;
if strcmp(groupName,'IOS GCaMP7s') == true
    Results_BehavHbT.(animalID).REM.MeanfLH = fLH_remCBVMean;
    Results_BehavHbT.(animalID).REM.MeanfRH = fRH_remCBVMean;
    Results_BehavHbT.(animalID).REM.IndfLH = fLH_remData;
    Results_BehavHbT.(animalID).REM.IndfRH = fRH_remData;
end
%% save data
cd(rootFolder)
save('Results_BehavHbT.mat','Results_BehavHbT')

end
