function [Results_BehavDiameter] = AnalyzeBehavioralDiameter(animalID,rootFolder,delim,Results_BehavDiameter)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the hemodynamic signal [] during different behavioral states (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
modelType = 'Forest';
params.minTime.Rest = 10;
params.Offset = 2;
params.minTime.Whisk = params.Offset + 5;
params.minTime.Stim = params.Offset + 2;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs
dataLocation = [rootFolder delim 'Data' delim animalID delim 'Bilateral Imaging'];
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
% lowpass filter
samplingRate = RestData.Pupil.pupilArea.CBVCamSamplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
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
%% analyze [] during periods of rest
% pull data from RestData.mat structure
[restLogical] = FilterEvents_IOS(RestData.Pupil.pupilArea,RestCriteria);
[puffLogical] = FilterEvents_IOS(RestData.Pupil.pupilArea,RestPuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFileIDs = RestData.Pupil.pupilArea.fileIDs(combRestLogical,:);
restEventTimes = RestData.Pupil.pupilArea.eventTimes(combRestLogical,:);
restDurations = RestData.Pupil.pupilArea.durations(combRestLogical,:);
restData = RestData.Pupil.pupilArea.data(combRestLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[finalRestData,finalRestFileIDs,~,~] = RemoveInvalidData_IOS(restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
% filter []
for gg = 1:length(finalRestData)
    procRestData{gg,1} = filtfilt(sos,g,finalRestData{gg,1}); %#ok<*AGROW>
end
% take mean [] during resting epochs
for nn = 1:length(procRestData)
    restMeanData(nn,1) = mean(procRestData{nn,1}(1:end));
end
% save results
Results_BehavDiameter.(animalID).Rest.meanPupilArea = restMeanData;
Results_BehavDiameter.(animalID).Rest.indPupilArea = procRestData;
Results_BehavDiameter.(animalID).Rest.fileIDs = finalRestFileIDs;
%% analyze [] during periods of moderate whisking (2-5 seconds)
% pull data from EventData.mat structure
[whiskLogical] = FilterEvents_IOS(EventData.Pupil.pupilArea.whisk,WhiskCriteria);
[puffLogical] = FilterEvents_IOS(EventData.Pupil.pupilArea.whisk,WhiskPuffCriteria);
combWhiskLogical = logical(whiskLogical.*puffLogical);
whiskFileIDs = EventData.Pupil.pupilArea.whisk.fileIDs(combWhiskLogical,:);
whiskEventTimes = EventData.Pupil.pupilArea.whisk.eventTime(combWhiskLogical,:);
whiskDurations = EventData.Pupil.pupilArea.whisk.duration(combWhiskLogical,:);
whiskData = EventData.Pupil.pupilArea.whisk.data(combWhiskLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[finalWhiskData,finalWhiskFileIDs,~,~] = RemoveInvalidData_IOS(whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
% filter [] and mean-subtract 2 seconds prior to whisk
for gg = 1:size(finalWhiskData,1)
    procWhiskData = filtfilt(sos,g,finalWhiskData(gg,:));
end
% take mean [] during whisking epochs from onset through 5 seconds
for nn = 1:size(procWhiskData,1)
    whiskMean{nn,1} = mean(procWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
    whisknd{nn,1} = procWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
end
% save results
Results_BehavDiameter.(animalID).Whisk.meanPupilArea = cell2mat(whiskMean);
Results_BehavDiameter.(animalID).Whisk.indPupilArea = whisknd;
Results_BehavDiameter.(animalID).Whisk.fileIDs = finalWhiskFileIDs;
%% analyze [] during periods of stimulation
% pull data from EventData.mat structure
LH_stimFilter = FilterEvents_IOS(EventData.Pupil.pupilArea.stim,StimCriteriaA);
RH_stimFilter = FilterEvents_IOS(EventData.Pupil.pupilArea.stim,StimCriteriaB);
[LH_stimData] = EventData.Pupil.pupilArea.stim.data(LH_stimFilter,:);
[RH_stimData] = EventData.Pupil.pupilArea.stim.data(RH_stimFilter,:);
[LH_stimFileIDs] = EventData.Pupil.pupilArea.stim.fileIDs(LH_stimFilter,:);
[RH_stimFileIDs] = EventData.Pupil.pupilArea.stim.fileIDs(RH_stimFilter,:);
[LH_stimEventTimes] = EventData.Pupil.pupilArea.stim.eventTime(LH_stimFilter,:);
[RH_stimEventTimes] = EventData.Pupil.pupilArea.stim.eventTime(RH_stimFilter,:);
LH_stimDurations = zeros(length(LH_stimEventTimes),1);
RH_stimDurations = zeros(length(RH_stimEventTimes),1);
% keep only the data that occurs within the manually-approved awake regions
[LH_finalStimData,LH_finalStimFileIDs,~,~] = RemoveInvalidData_IOS(LH_stimData,LH_stimFileIDs,LH_stimDurations,LH_stimEventTimes,ManualDecisions);
[RH_finalStimData,RH_finalStimFileIDs,~,~] = RemoveInvalidData_IOS(RH_stimData,RH_stimFileIDs,RH_stimDurations,RH_stimEventTimes,ManualDecisions);
% filter [] and mean-subtract 2 seconds prior to stimulus (left hem)
for gg = 1:size(LH_finalStimData,1)
    LH_ProcStimData = filtfilt(sos,g,LH_finalStimData(gg,:));
end
% filter [] and mean-subtract 2 seconds prior to stimulus (right hem)
for gg = 1:size(RH_finalStimData,1)
    RH_ProcStimData = filtfilt(sos,g,RH_finalStimData(gg,:));
end
% take mean [] 1-2 seconds after stimulation (left hem)
for nn = 1:size(LH_ProcStimData,1)
    LH_stimMean{nn,1} = mean(LH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
    LH_stim{nn,1} = LH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
end
% take mean [] 1-2 seconds after stimulation (right hem)
for nn = 1:size(RH_ProcStimData,1)
    RH_stimMean{nn,1} = mean(RH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
    RH_stim{nn,1} = RH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
end
% save results
Results_BehavDiameter.(animalID).Stim.meanPupilArea = cat(1,cell2mat(LH_stimMean),cell2mat(RH_stimMean));
Results_BehavDiameter.(animalID).Stim.indPupilArea = cat(1,LH_stim,RH_stim);
Results_BehavDiameter.(animalID).Stim.fileIDs = cat(1,LH_finalStimFileIDs,RH_finalStimFileIDs);
%% analyze [] during periods of NREM sleep
% pull data from SleepData.mat structure
if isempty(SleepData.(modelType).NREM.data.Pupil) == false
    [nremData,nremFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Pupil.pupilArea,SleepData.(modelType).NREM.data.Pupil.fileIDs,SleepData.(modelType).NREM.data.Pupil.binTimes);
    % filter and take mean [] during NREM epochs
    for nn = 1:length(nremData)
        try
            nremMean(nn,1) = mean(filtfilt(sos,g,nremData{nn,1}(1:end)));
        catch
            nremMean(nn,1) = nanmean(nremData{nn,1}(1:end));
        end
    end
    % save results
    Results_BehavDiameter.(animalID).NREM.meanPupilArea = nremMean;
    Results_BehavDiameter.(animalID).NREM.indPupilArea = nremData;
    Results_BehavDiameter.(animalID).NREM.FileIDs = nremFileIDs;
else
    % save results
    Results_BehavDiameter.(animalID).NREM.meanPupilArea = [];
    Results_BehavDiameter.(animalID).NREM.indPupilArea = [];
    Results_BehavDiameter.(animalID).NREM.FileIDs = [];
end
%% analyze [] during periods of REM sleep
% pull data from SleepData.mat structure
if isempty(SleepData.(modelType).REM.data.Pupil) == false
    [remData,remFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Pupil.pupilArea,SleepData.(modelType).REM.data.Pupil.fileIDs,SleepData.(modelType).REM.data.Pupil.binTimes);
    % filter and take mean [] during REM epochs
    for nn = 1:length(remData)
        try
            remMean(nn,1) = mean(filtfilt(sos,g,remData{nn,1}(1:end)));
        catch
            remMean(nn,1) = nanmean(remData{nn,1}(1:end));
        end
    end
    % save results
    Results_BehavDiameter.(animalID).REM.meanPupilArea = remMean;
    Results_BehavDiameter.(animalID).REM.indPupilArea = remData;
    Results_BehavDiameter.(animalID).REM.FileIDs = remFileIDs;
else
    % save results
    Results_BehavDiameter.(animalID).REM.meanPupilArea = [];
    Results_BehavDiameter.(animalID).REM.indPupilArea = [];
    Results_BehavDiameter.(animalID).REM.FileIDs = [];
end
% save data
cd([rootFolder delim])
save('Results_BehavDiameter.mat','Results_BehavDiameter')

end
