function [Results_Evoked] = AnalyzeEvokedResponses_Pupil(animalID,rootFolder,delim,Results_Evoked)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

%% function parameters
dataLocation = [rootFolder delim 'Data' delim animalID delim 'Bilateral Imaging'];
%% only run analysis for valid animal IDs
cd(dataLocation)
% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID,'-mat')
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
% find and load AllSpecStruct.mat struct
allSpecStructFileStruct = dir('*_AllSpecStructB.mat');
allSpecStructFile = {allSpecStructFileStruct.name}';
allSpecStructFileID = char(allSpecStructFile);
load(allSpecStructFileID,'-mat')
% criteria for whisking
WhiskCriteriaA.Fieldname = {'duration','duration','puffDistance'};
WhiskCriteriaA.Comparison = {'gt','lt','gt'};
WhiskCriteriaA.Value = {0.5,2,5};
WhiskCriteriaB.Fieldname = {'duration','duration','puffDistance'};
WhiskCriteriaB.Comparison = {'gt','lt','gt'};
WhiskCriteriaB.Value = {2,5,5};
WhiskCriteriaC.Fieldname = {'duration','puffDistance'};
WhiskCriteriaC.Comparison = {'gt','gt'};
WhiskCriteriaC.Value = {5,5};
WhiskCriteriaNames = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
% criteria for stimulation
StimCriteriaA.Value = {'LPadSol'};
StimCriteriaA.Fieldname = {'solenoidName'};
StimCriteriaA.Comparison = {'equal'};
StimCriteriaB.Value = {'RPadSol'};
StimCriteriaB.Fieldname = {'solenoidName'};
StimCriteriaB.Comparison = {'equal'};
StimCriteriaC.Value = {'AudSol'};
StimCriteriaC.Fieldname = {'solenoidName'};
StimCriteriaC.Comparison = {'equal'};
stimCriteriaNames = {'stimCriteriaA','stimCriteriaB','stimCriteriaC'};
%% analyze whisking-evoked responses
% pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
samplingRate = EventData.Pupil.pupilArea.whisk.samplingRate;
offset = EventData.Pupil.pupilArea.whisk.epoch.offset;
for bb = 1:length(WhiskCriteriaNames)
    whiskCriteriaName = WhiskCriteriaNames{1,bb};
    if strcmp(whiskCriteriaName,'ShortWhisks') == true
        WhiskCriteria = WhiskCriteriaA;
    elseif strcmp(whiskCriteriaName,'IntermediateWhisks') == true
        WhiskCriteria = WhiskCriteriaB;
    elseif strcmp(whiskCriteriaName,'LongWhisks') == true
        WhiskCriteria = WhiskCriteriaC;
    end
    % pull data from EventData.mat structure
    [whiskLogical] = FilterEvents_IOS(EventData.Pupil.pupilArea.whisk,WhiskCriteria);
    combWhiskLogical = logical(whiskLogical);
    [allWhiskData] = EventData.Pupil.pupilArea.whisk.data(combWhiskLogical,:);
    [allWhiskFileIDs] = EventData.Pupil.pupilArea.whisk.fileIDs(combWhiskLogical,:);
    [allWhiskEventTimes] = EventData.Pupil.pupilArea.whisk.eventTime(combWhiskLogical,:);
    allWhiskDurations = EventData.Pupil.pupilArea.whisk.duration(combWhiskLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [finalWhiskData,~,~,~] = RemoveInvalidData_IOS(allWhiskData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
    % lowpass filter each whisking event and nanmean-subtract by the first 2 seconds
    clear procWhiskData
    dd = 1;
    for cc = 1:size(finalWhiskData,1)
        whiskHbTarray = finalWhiskData(cc,:);
        filtWhiskarray = sgolayfilt(whiskHbTarray,3,17);
        procWhiskData(dd,:) = filtWhiskarray - nanmean(filtWhiskarray(1:(offset*samplingRate))); %#ok<*AGROW>
    end
    meanWhiskData = nanmean(procWhiskData,1);
    stdWhiskData = nanstd(procWhiskData,0,1);
    % save results
    Results_Evoked.(animalID).Whisk.pupilArea.(whiskCriteriaName).pupilArea.mean = meanWhiskData;
    Results_Evoked.(animalID).Whisk.pupilArea.(whiskCriteriaName).pupilArea.stdev = stdWhiskData;
end
%% analyze stimulus-evoked responses
for gg = 1:length(stimCriteriaNames)
    stimCriteriaName = stimCriteriaNames{1,gg};
    if strcmp(stimCriteriaName,'stimCriteriaA') == true
        StimCriteria = StimCriteriaA;
        solenoid = 'LPadSol';
    elseif strcmp(stimCriteriaName,'stimCriteriaB') == true
        StimCriteria = StimCriteriaB;
        solenoid = 'RPadSol';
    elseif strcmp(stimCriteriaName,'stimCriteriaC') == true
        StimCriteria = StimCriteriaC;
        solenoid = 'AudSol';
    end
    % pull data from EventData.mat structure
    allStimFilter = FilterEvents_IOS(EventData.Pupil.pupilArea.stim,StimCriteria);
    [allStimData] = EventData.Pupil.pupilArea.stim.data(allStimFilter,:);
    [allStimFileIDs] = EventData.Pupil.pupilArea.stim.fileIDs(allStimFilter,:);
    [allStimEventTimes] = EventData.Pupil.pupilArea.stim.eventTime(allStimFilter,:);
    allStimDurations = zeros(length(allStimEventTimes),1);
    % keep only the data that occurs within the manually-approved awake regions
    [finalStimData,~,~,~] = RemoveInvalidData_IOS(allStimData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
    % lowpass filter each stim event and nanmean-subtract by the first 2 seconds
    clear procStimData
    for hh = 1:size(finalStimData,1)
        stimHbTarray = finalStimData(hh,:);
        filtStimarray = sgolayfilt(stimHbTarray,3,17);
        procStimData(hh,:) = filtStimarray - nanmean(filtStimarray(1:(offset*samplingRate)));
    end
    meanStimData = nanmean(procStimData,1);
    stdStimData = nanstd(procStimData,0,1);
    % save results
    Results_Evoked.(animalID).Stim.(solenoid).pupilArea.mean = meanStimData;
    Results_Evoked.(animalID).Stim.(solenoid).pupilArea.std = stdStimData;
end
% save data
cd([rootFolder delim])
save('Results_Evoked.mat','Results_Evoked')

end
