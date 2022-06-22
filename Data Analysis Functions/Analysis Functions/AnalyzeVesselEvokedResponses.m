function [Results_VesselEvoked] = AnalyzeVesselEvokedResponses(animalID,group,rootFolder,delim,Results_VesselEvoked)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the whisking-evoked arteriole D/D responses (2PLSM)
%________________________________________________________________________________________________________________________

%% only run analysis for valid animal IDs
dataLocation = [rootFolder delim group delim animalID delim 'Combined Imaging'];
cd(dataLocation)
% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID)
% find and load manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID)
% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
% criteria for whisking
whiskCriteriaA.Fieldname = {'duration','duration','puffDistance'};
whiskCriteriaA.Comparison = {'gt','lt','gt'};
whiskCriteriaA.Value = {0.5,2,5};
whiskCriteriaB.Fieldname = {'duration','duration','puffDistance'};
whiskCriteriaB.Comparison = {'gt','lt','gt'};
whiskCriteriaB.Value = {2,5,5};
whiskCriteriaC.Fieldname = {'duration','puffDistance'};
whiskCriteriaC.Comparison = {'gt','gt'};
whiskCriteriaC.Value = {5,5};
whiskCriteriaNames = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
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
for qq = 1:length(whiskCriteriaNames)
    % pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
    samplingRate = 5;
    offset = EventData.vesselDiameter.data.whisk.epoch.offset;
    timeVector = (0:(EventData.vesselDiameter.data.whisk.epoch.duration*samplingRate))/samplingRate - EventData.vesselDiameter.data.whisk.epoch.offset;
    whiskCriteriaName = whiskCriteriaNames{1,qq};
    if strcmp(whiskCriteriaName,'ShortWhisks') == true
        WhiskCriteria = whiskCriteriaA;
    elseif strcmp(whiskCriteriaName,'IntermediateWhisks') == true
        WhiskCriteria = whiskCriteriaB;
    elseif strcmp(whiskCriteriaName,'LongWhisks') == true
        WhiskCriteria = whiskCriteriaC;
    end
    % pull data from EventData.mat structure
    [whiskLogical] = FilterEvents_2P(EventData.vesselDiameter.data.whisk,WhiskCriteria);
    whiskLogical = logical(whiskLogical);
    whiskingData = EventData.vesselDiameter.data.whisk.data(whiskLogical,:);
    whiskFileIDs = EventData.vesselDiameter.data.whisk.fileIDs(whiskLogical,:);
    whiskVesselIDs = EventData.vesselDiameter.data.whisk.vesselIDs(whiskLogical,:);
    whiskEventTimes = EventData.vesselDiameter.data.whisk.eventTime(whiskLogical,:);
    whiskDurations = EventData.vesselDiameter.data.whisk.duration(whiskLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [finalWhiskData,finalWhiskFileIDs,finalWhiskVesselIDs,~,~] = RemoveInvalidData_2P(whiskingData,whiskFileIDs,whiskVesselIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    clear procWhiskData
    % filter and detrend data
    for aa = 1:size(finalWhiskData,1)
        whiskStrDay = ConvertDate_2P(finalWhiskFileIDs{aa,1}(1:6));
        normWhiskData = (finalWhiskData(aa,:) - RestingBaselines.manualSelection.vesselDiameter.data.(finalWhiskVesselIDs{aa,1}).(whiskStrDay))./RestingBaselines.manualSelection.vesselDiameter.data.(finalWhiskVesselIDs{aa,1}).(whiskStrDay);
        filtWhiskData = sgolayfilt(normWhiskData,3,17);
        procWhiskData{aa,1} = filtWhiskData - mean(filtWhiskData(1:(offset*samplingRate))); %#ok<*AGROW>
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    reshapedWhiskData = zeros(length(procWhiskData{1,1}),length(procWhiskData));
    for bb = 1:length(procWhiskData)
        reshapedWhiskData(:,bb) = procWhiskData{bb,1};
    end
    % remove veins from the artery list
    uniqueWhiskVesselIDs = unique(finalWhiskVesselIDs);
    dd = 1;
    clear whiskArterioleIDs
    for cc = 1:length(uniqueWhiskVesselIDs)
        if strcmp(uniqueWhiskVesselIDs{cc,1}(1),'V') == false
            whiskArterioleIDs{dd,1} = uniqueWhiskVesselIDs{cc,1};
            dd = dd + 1;
        end
    end
    % split the data based on different arteries
    clear whiskArterioleEvoked
    for ee = 1:length(whiskArterioleIDs)
        whiskArterioleID = whiskArterioleIDs{ee,1};
        gg = 1;
        for ff = 1:length(finalWhiskVesselIDs)
            if strcmp(whiskArterioleID,finalWhiskVesselIDs{ff,1}) == true
                whiskArterioleEvoked.(whiskArterioleID)(gg,:) = reshapedWhiskData(:,ff);
                gg = gg + 1;
            end
        end
    end
    % take mean/std of each arteriole's whisk-evoked response
    for aa = 1:length(whiskArterioleIDs)
        vID = whiskArterioleIDs{aa,1};
        meanWhiskEvokedDiam.(vID) = mean(whiskArterioleEvoked.(vID),1)*100;
        stdWhiskEvokedDiam.(vID) = std(whiskArterioleEvoked.(vID),0,1)*100;
        baselineDates = fieldnames(RestingBaselines.manualSelection.vesselDiameter.data.(vID));
        baselineDiameters = [];
        for bb = 1:length(baselineDates)
            baselineDiameters = cat(1,baselineDiameters,RestingBaselines.manualSelection.vesselDiameter.data.(vID).(baselineDates{bb,1}));
        end
        % save results
        Results_VesselEvoked.(animalID).Whisk.(whiskCriteriaName).(vID).mean = meanWhiskEvokedDiam.(vID);
        Results_VesselEvoked.(animalID).Whisk.(whiskCriteriaName).(vID).StD = stdWhiskEvokedDiam.(vID);
        Results_VesselEvoked.(animalID).Whisk.(whiskCriteriaName).(vID).baseline = mean(baselineDiameters);
        Results_VesselEvoked.(animalID).Whisk.(whiskCriteriaName).(vID).timeVector = timeVector;
    end
end
%% analyze stimulus-evoked responses
for zz = 1:length(stimCriteriaNames)
    StimCriteriaName = stimCriteriaNames{1,zz};
    if strcmp(StimCriteriaName,'stimCriteriaA') == true
        StimCriteria = StimCriteriaA;
        solenoid = 'LPadSol';
    elseif strcmp(StimCriteriaName,'stimCriteriaB') == true
        StimCriteria = StimCriteriaB;
        solenoid = 'RPadSol';
    elseif strcmp(StimCriteriaName,'stimCriteriaC') == true
        StimCriteria = StimCriteriaC;
        solenoid = 'AudSol';
    end
    % pull data from EventData.mat structure
    [stimLogical] = FilterEvents_IOS(EventData.vesselDiameter.data.stim,StimCriteria);
    stimLogical = logical(stimLogical);
    stimData = EventData.vesselDiameter.data.stim.data(stimLogical,:);
    stimFileIDs = EventData.vesselDiameter.data.stim.fileIDs(stimLogical,:);
    stimVesselIDs = EventData.vesselDiameter.data.stim.vesselIDs(stimLogical,:);
    stimEventTimes = EventData.vesselDiameter.data.stim.eventTime(stimLogical,:);
    stimDurations = zeros(length(stimEventTimes),1);
    % keep only the data that occurs within the manually-approved awake regions
    [finalStimData,finalStimFileIDs,finalStimVesselIDs,~,~] = RemoveInvalidData_2P(stimData,stimFileIDs,stimVesselIDs,stimDurations,stimEventTimes,ManualDecisions);
    clear procStimData
    %% filter and detrend data
    for aa = 1:size(finalStimData,1)
        StimStrDay = ConvertDate_2P(finalStimFileIDs{aa,1}(1:6));
        normStimData = (finalStimData(aa,:) - RestingBaselines.manualSelection.vesselDiameter.data.(finalStimVesselIDs{aa,1}).(StimStrDay))./RestingBaselines.manualSelection.vesselDiameter.data.(finalStimVesselIDs{aa,1}).(StimStrDay);
        filtStimData = sgolayfilt(normStimData,3,17);
        procStimData{aa,1} = filtStimData - mean(filtStimData(1:(offset*samplingRate))); %#ok<*AGROW>
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    reshapedStimData = zeros(length(procStimData{1,1}),length(procStimData));
    for bb = 1:length(procStimData)
        reshapedStimData(:,bb) = procStimData{bb,1};
    end
    % remove veins from the artery list
    uniqueStimVesselIDs = unique(finalStimVesselIDs);
    dd = 1;
    clear StimArterioleIDs
    for cc = 1:length(uniqueStimVesselIDs)
        if strcmp(uniqueStimVesselIDs{cc,1}(1),'V') == false
            StimArterioleIDs{dd,1} = uniqueStimVesselIDs{cc,1};
            dd = dd + 1;
        end
    end
    % split the data based on different arteries
    clear StimArterioleEvoked
    for ee = 1:length(StimArterioleIDs)
        StimArterioleID = StimArterioleIDs{ee,1};
        gg = 1;
        for ff = 1:length(finalStimVesselIDs)
            if strcmp(StimArterioleID,finalStimVesselIDs{ff,1}) == true
                StimArterioleEvoked.(StimArterioleID)(gg,:) = reshapedStimData(:,ff);
                gg = gg + 1;
            end
        end
    end
    % take mean/std of each arteriole's Stim-evoked response
    for aa = 1:length(StimArterioleIDs)
        vID = StimArterioleIDs{aa,1};
        meanStimEvokedDiam.(vID) = mean(StimArterioleEvoked.(vID),1)*100;
        stdStimEvokedDiam.(vID) = std(StimArterioleEvoked.(vID),0,1)*100;
        baselineDates = fieldnames(RestingBaselines.manualSelection.vesselDiameter.data.(vID));
        baselineDiameters = [];
        for bb = 1:length(baselineDates)
            baselineDiameters = cat(1,baselineDiameters,RestingBaselines.manualSelection.vesselDiameter.data.(vID).(baselineDates{bb,1}));
        end
        % save results
        Results_VesselEvoked.(animalID).Stim.(solenoid).(vID).count = size(StimArterioleEvoked,1);
        Results_VesselEvoked.(animalID).Stim.(solenoid).(vID).mean = meanStimEvokedDiam.(vID);
        Results_VesselEvoked.(animalID).Stim.(solenoid).(vID).StD = stdStimEvokedDiam.(vID);
        Results_VesselEvoked.(animalID).Stim.(solenoid).(vID).baseline = mean(baselineDiameters);
        Results_VesselEvoked.(animalID).Stim.(solenoid).(vID).timeVector = timeVector;
    end
end
% save data
cd(rootFolder)
save('Results_VesselEvoked.mat','Results_VesselEvoked')

end
