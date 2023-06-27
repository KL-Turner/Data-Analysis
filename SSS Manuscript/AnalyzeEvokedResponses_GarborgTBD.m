function [Results_Evoked] = AnalyzeEvokedResponses_GarborgTBD(animalID,rootFolder,delim,Results_Evoked)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Determine whisking and stimulus evoked changes in fluorescence around the edges of the SSS
%________________________________________________________________________________________________________________________

% go to animal's data location
dataLocation = [rootFolder delim 'Data' delim animalID delim 'SSS_B'];
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
% criteria for whisking
WhiskCriteriaA.Fieldname = {'duration','duration','stimDistance'};
WhiskCriteriaA.Comparison = {'gt','lt','gt'};
WhiskCriteriaA.Value = {0.5,2,5};
WhiskCriteriaB.Fieldname = {'duration','duration','stimDistance'};
WhiskCriteriaB.Comparison = {'gt','lt','gt'};
WhiskCriteriaB.Value = {2,5,5};
WhiskCriteriaC.Fieldname = {'duration','stimDistance'};
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
dataTypes = {'SSS','lSSS','rSSS'};
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    %% analyze whisking-evoked responses
    % pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
    samplingRate = EventData.CBV.(dataType).whisk.samplingRate;
    offset = EventData.CBV.(dataType).whisk.epoch.offset;
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
        [whiskLogical] = FilterEvents_IOS(EventData.CBV.(dataType).whisk,WhiskCriteria);
        combWhiskLogical = logical(whiskLogical);
        [allWhiskData] = EventData.CBV.(dataType).whisk.NormData(combWhiskLogical,:);
        [allWhiskFileIDs] = EventData.CBV.(dataType).whisk.fileIDs(combWhiskLogical,:);
        [allWhiskEventTimes] = EventData.CBV.(dataType).whisk.eventTime(combWhiskLogical,:);
        allWhiskDurations = EventData.CBV.(dataType).whisk.duration(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [finalWhiskData,~,~,~] = RemoveInvalidData_IOS(allWhiskData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        % lowpass filter each whisking event and nanmean-subtract by the first 2 seconds
        procWhiskData = [];
        for cc = 1:size(finalWhiskData,1)
            whiskArray = finalWhiskData(cc,:);
            filtWhiskarray = whiskArray;
            % filtWhiskarray = sgolayfilt(whiskArray,3,17);
            procWhiskData(cc,:) = (filtWhiskarray - mean(filtWhiskarray(1:(offset*samplingRate)),'omitnan')).*100;
        end
        meanWhiskData = mean(procWhiskData,1,'omitnan');
        stdWhiskData = std(procWhiskData,0,1,'omitnan');
        % save results
        Results_Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).mean = meanWhiskData;
        Results_Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).stdev = stdWhiskData;
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
        allStimFilter = FilterEvents_IOS(EventData.CBV.(dataType).stim,StimCriteria);
        [allStimData] = EventData.CBV.(dataType).stim.NormData(allStimFilter,:);
        [allStimFileIDs] = EventData.CBV.(dataType).stim.fileIDs(allStimFilter,:);
        [allStimEventTimes] = EventData.CBV.(dataType).stim.eventTime(allStimFilter,:);
        allStimDurations = zeros(length(allStimEventTimes),1);
        % keep only the data that occurs within the manually-approved awake regions
        [finalStimData,~,~,~] = RemoveInvalidData_IOS(allStimData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
        % lowpass filter each stim event and nanmean-subtract by the first 2 seconds
        procStimData = [];
        for kk = 1:size(finalStimData,1)
            stimArray = finalStimData(kk,:);
            filtStimarray = stimArray;
            % filtStimarray = sgolayfilt(stimArray,3,17);
            procStimData(kk,:) = (filtStimarray - mean(filtStimarray(1:(offset*samplingRate)),'omitnan')).*100;
        end
        meanStimData = mean(procStimData,1,'omitnan');
        stdStimData = std(procStimData,0,1,'omitnan');
        % save results
        Results_Evoked.(animalID).Stim.(dataType).(solenoid).mean = meanStimData;
        Results_Evoked.(animalID).Stim.(dataType).(solenoid).std = stdStimData;
    end
end
% save data
cd([rootFolder delim 'Analysis Structures'])
save('Results_Evoked.mat','Results_Evoked')

end
