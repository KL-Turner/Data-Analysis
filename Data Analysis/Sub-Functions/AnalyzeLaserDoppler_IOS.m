function [AnalysisResults] = AnalyzeLaserDoppler_IOS(AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner%
%
%   Purpose:
%________________________________________________________________________________________________________________________

%% find and load RestData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID)

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

% determine the animal's ID use the EventData.mat file's name for the current folder
fileBreaks = strfind(eventDataFileID,'_');
animalID = eventDataFileID(1:fileBreaks(1)-1);

%% Whisking-evoked responses
% criteria for the FilterEvents data struct
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

% filter the EventData.mat structure for whisking events that meet the desired criteria
% pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
samplingRate = EventData.flow.data.whisk.samplingRate;
timeVector = (0:(EventData.flow.data.whisk.epoch.duration*samplingRate))/samplingRate - EventData.flow.data.whisk.epoch.offset;
offset = EventData.flow.data.whisk.epoch.offset;
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);

for b = 1:length(whiskCriteriaNames)
    whiskCriteriaName = whiskCriteriaNames{1,b};
    if strcmp(whiskCriteriaName,'ShortWhisks') == true
        whiskCriteria = whiskCriteriaA;
    elseif strcmp(whiskCriteriaName,'IntermediateWhisks') == true
        whiskCriteria = whiskCriteriaB;
    elseif strcmp(whiskCriteriaName,'LongWhisks') == true
        whiskCriteria = whiskCriteriaC;
    end
    disp(['AnalyzeEvokedResponses: ' whiskCriteriaName ' whisker events']); disp(' ')
    allWhiskFilter = FilterEvents_IOS(EventData.flow.data.whisk,whiskCriteria);
    [allWhiskFlowData] = EventData.flow.data.whisk.NormData(allWhiskFilter,:);
    [allWhiskFileIDs] = EventData.flow.data.whisk.fileIDs(allWhiskFilter,:);
    
    % identify the unique days and the unique number of files from the list of all whisking events
    whiskUniqueDays = GetUniqueDays_IOS(allWhiskFileIDs);
    whiskUniqueFiles = unique(allWhiskFileIDs);
    whiskNumberOfFiles = length(unique(allWhiskFileIDs));
    
    % decimate the file list to only include those files that occur within the desired number of target minutes
    clear whiskFiltLogical
    for c = 1:length(whiskUniqueDays)
        whiskDay = whiskUniqueDays(c);
        d = 1;
        for n = 1:whiskNumberOfFiles
            whiskFile = whiskUniqueFiles(n);
            whiskFileID = whiskFile{1}(1:6);
            if strcmp(whiskDay,whiskFileID) && sum(strcmp(whiskFile,manualFileIDs)) == 1
                whiskFiltLogical{c,1}(n,1) = 1; %#ok<*AGROW>
                d = d + 1;
            else
                whiskFiltLogical{c,1}(n,1) = 0;
            end
        end
    end
    whiskFinalLogical = any(sum(cell2mat(whiskFiltLogical'),2),2);
    
    % extract all the whisking events that correspond to the acceptable file list and the acceptable whisking criteria
    clear whiskFileFilter
    filtWhiskFiles = whiskUniqueFiles(whiskFinalLogical,:);
    for e = 1:length(allWhiskFileIDs)
        whiskLogic = strcmp(allWhiskFileIDs{e},filtWhiskFiles);
        whiskLogicSum = sum(whiskLogic);
        if whiskLogicSum == 1
            whiskFileFilter(e,1) = 1;
        else
            whiskFileFilter(e,1) = 0;
        end
    end
    finalWhiskFileFilter = logical(whiskFileFilter);
    finalWhiskFlowData = allWhiskFlowData(finalWhiskFileFilter,:);
    
    % lowpass filter each whisking event and mean-subtract by the first 2 seconds
    clear procWhiskFlowData
    for f = 1:size(finalWhiskFlowData,1)
        whiskFlowArray = finalWhiskFlowData(f,:);
        filtWhiskFlowArray = sgolayfilt(whiskFlowArray,3,17);
        procWhiskFlowData(f,:) = filtWhiskFlowArray - mean(filtWhiskFlowArray(1:(offset*samplingRate)));
    end
    meanWhiskFlowData = mean(procWhiskFlowData,1)*100;
    stdWhiskFlowData = std(procWhiskFlowData,0,1)*100;
    
    % summary figure
    whiskEvoked = figure;
    plot(timeVector,meanWhiskFlowData,'k')
    hold on
    plot(timeVector,meanWhiskFlowData + stdWhiskFlowData,'color',colors_IOS('battleship grey'))
    plot(timeVector,meanWhiskFlowData - stdWhiskFlowData,'color',colors_IOS('battleship grey'))
    title([animalID ' whisking-evoked flow averages - ' whiskCriteriaName])
    xlabel('Time (sec)')
    ylabel('Perfusion Units (%)')
    axis tight
    axis square
    
    % save results
    AnalysisResults.EvokedAvgs.Whisk.DopplerFlow.(whiskCriteriaName).flowMean = meanWhiskFlowData;
    AnalysisResults.EvokedAvgs.Whisk.DopplerFlow.(whiskCriteriaName).flowStD = stdWhiskFlowData;
    AnalysisResults.EvokedAvgs.Whisk.DopplerFlow.(whiskCriteriaName).timeVector = timeVector;
    
    % save figure
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Combined Imaging/Figures/Stim and Whisk Evoked Averages/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(whiskEvoked, [dirpath animalID '_DopplerFlow_' whiskCriteriaName '_WhiskEvokedAverages']);
    close(whiskEvoked)
end

%% Stimulus-evoked responses
% Criteria for the FilterEvents data struct
stimCriteriaA.Value = {'LPadSol'};
stimCriteriaA.Fieldname = {'solenoidName'};
stimCriteriaA.Comparison = {'equal'};

stimCriteriaB.Value = {'RPadSol'};
stimCriteriaB.Fieldname = {'solenoidName'};
stimCriteriaB.Comparison = {'equal'};

stimCriteriaC.Value = {'AudSol'};
stimCriteriaC.Fieldname = {'solenoidName'};
stimCriteriaC.Comparison = {'equal'};

stimCriteriaNames = {'stimCriteriaA','stimCriteriaB','stimCriteriaC'};

% filter the EventData.mat structure for stimulus events that meet the desired criteria
for j = 1:length(stimCriteriaNames)
    stimCriteriaName = stimCriteriaNames{1,j};
    if strcmp(stimCriteriaName,'stimCriteriaA') == true
        stimCriteria = stimCriteriaA;
        solenoid = 'LPadSol';
    elseif strcmp(stimCriteriaName,'stimCriteriaB') == true
        stimCriteria = stimCriteriaB;
        solenoid = 'RPadSol';
    elseif strcmp(stimCriteriaName,'stimCriteriaC') == true
        stimCriteria = stimCriteriaC;
        solenoid = 'AudSol';
    end
    disp(['AnalyzeEvokedResponses: ' solenoid ' stimulus events']); disp(' ')
    allStimFilter = FilterEvents_IOS(EventData.flow.data.stim,stimCriteria);
    [allStimFlowData] = EventData.flow.data.stim.data(allStimFilter,:);
    [allStimFileIDs] = EventData.flow.data.stim.fileIDs(allStimFilter,:);
    
    % identify the unique days and the unique number of files from the list of all whisking events
    stimUniqueDays = GetUniqueDays_IOS(allStimFileIDs);
    stimUniqueFiles = unique(allStimFileIDs);
    stimNumberOfFiles = length(unique(allStimFileIDs));
    
    % decimate the file list to only include those files that occur within the desired number of target minutes
    clear stimFiltLogical
    for k = 1:length(stimUniqueDays)
        stimDay = stimUniqueDays(k);
        m = 1;
        for n = 1:stimNumberOfFiles
            stimFile = stimUniqueFiles(n);
            stimFileID = stimFile{1}(1:6);
            if strcmp(stimDay,stimFileID) && sum(strcmp(stimFile,manualFileIDs)) == 1
                stimFiltLogical{c,1}(n,1) = 1; %#ok<*AGROW>
                m = m + 1;
            else
                stimFiltLogical{c,1}(n,1) = 0;
            end
        end
    end
    stimFinalLogical = any(sum(cell2mat(stimFiltLogical'),2),2);
    
    % extract all the whisking events that correspond to the acceptable file list and the acceptable whisking criteria
    clear stimFileFilter
    filtStimFiles = stimUniqueFiles(stimFinalLogical,:);
    for o = 1:length(allStimFileIDs)
        stimLogic = strcmp(allStimFileIDs{o},filtStimFiles);
        stimLogicSum = sum(stimLogic);
        if stimLogicSum == 1
            stimFileFilter(o,1) = 1;
        else
            stimFileFilter(o,1) = 0;
        end
    end
    finalStimFileFilter = logical(stimFileFilter);
    finalStimFlowData = allStimFlowData(finalStimFileFilter,:);
    
    % lowpass filter each whisking event and mean-subtract by the first 2 seconds
    clear procStimFlowData
    for p = 1:size(finalStimFlowData,1)
        stimFlowArray = finalStimFlowData(p,:);
        filtStimFlowArray = sgolayfilt(stimFlowArray,3,17);
        procStimFlowData(p,:) = filtStimFlowArray - mean(filtStimFlowArray(1:(offset*samplingRate)));
    end
    meanStimFlowData = mean(procStimFlowData,1)*100;
    stdStimFlowData = std(procStimFlowData,0,1)*100;
   
    % summary figure
    stimEvoked = figure;
    plot(timeVector,meanStimFlowData,'k')
    hold on
    plot(timeVector,meanStimFlowData + stdStimFlowData,'color',colors_IOS('battleship grey'))
    plot(timeVector,meanStimFlowData - stdStimFlowData,'color',colors_IOS('battleship grey'))
    title([animalID ' stimulus-evoked flow averages - ' solenoid])
    xlabel('Time (sec)')
    ylabel('Perfusion Units (%)')
    axis tight
    axis square
    
    % save results
    AnalysisResults.EvokedAvgs.Stim.DopperFlow.(solenoid).flowMean = meanStimFlowData;
    AnalysisResults.EvokedAvgs.Stim.DopperFlow.(solenoid).flowStD = stdStimFlowData;
    AnalysisResults.EvokedAvgs.Stim.DopperFlow.(solenoid).timeVector = timeVector;
    
    % save figure
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Combined Imaging/Figures/Stim and Whisk Evoked Averages/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(stimEvoked,[dirpath animalID '_DopperFlow_' solenoid '_StimEvokedAverages']);
%     close(stimEvoked)
end

%% Average behavior-dependent flow
% Analyze mean CBV during periods of rest
disp('AnalyzeMeanFlow: average laser doppler flow during rest')
[restLogical] = FilterEvents_IOS(RestData.flow.data,RestCriteria);
[puffLogical] = FilterEvents_IOS(RestData.flow.data,RestPuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFiles = RestData.flow.data.fileIDs(combRestLogical,:);
restingData = RestData.flow.data.NormData(combRestLogical,:);

% identify the unique days and the unique number of files from the list of unstim resting events
restUniqueDays = GetUniqueDays_IOS(restFiles);
restUniqueFiles = unique(restFiles);
restNumberOfFiles = length(unique(restFiles));

% decimate the file list to only include those files that occur within the desired number of target minutes
clear restFiltLogical
for c = 1:length(restUniqueDays)
    restDay = restUniqueDays(c);
    d = 1;
    for e = 1:restNumberOfFiles
        restFile = restUniqueFiles(e);
        restFileID = restFile{1}(1:6);
        if strcmp(restDay,restFileID) && sum(strcmp(restFile,manualFileIDs)) == 1
            restFiltLogical{c,1}(e,1) = 1; %#ok<*AGROW>
            d = d + 1;
        else
            restFiltLogical{c,1}(e,1) = 0;
        end
    end
end
restFinalLogical = any(sum(cell2mat(restFiltLogical'),2),2);

% extract unstim the resting events that correspond to the acceptable file list and the acceptable resting criteria
clear restFileFilter
filtRestFiles = restUniqueFiles(restFinalLogical,:);
for f = 1:length(restFiles)
    restLogic = strcmp(restFiles{f},filtRestFiles);
    restLogicSum = sum(restLogic);
    if restLogicSum == 1
        restFileFilter(f,1) = 1;
    else
        restFileFilter(f,1) = 0;
    end
end
restFinalFileFilter = logical(restFileFilter);
finalRestData = restingData(restFinalFileFilter,:);

% only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
% original epoch create so we can add a sample of two back to the end for those just under 10 seconds
% lowpass filter and detrend each segment
[B, A] = butter(3,1/(samplingRate/2),'low');
clear procRestData
for g = 1:length(finalRestData)
    procRestData{g,1} = filtfilt(B,A,finalRestData{g,1});
end

% analyze correlation coefficient between resting epochs
for n = 1:length(procRestData)
    restFlowMean(n,1) = mean(procRestData{n,1});
end

% save results
AnalysisResults.MeanFlow.Rest = restFlowMean;

%% Analyze mean CBV during periods of extended whisking
disp('AnalzeMeanCBV: average hemodynamics during whisking')
[whiskLogical] = FilterEvents_IOS(EventData.flow.data.whisk,WhiskCriteria);
[puffLogical] = FilterEvents_IOS(EventData.flow.data.whisk,WhiskPuffCriteria);
combWhiskLogical = logical(whiskLogical.*puffLogical);
whiskFiles = EventData.flow.data.whisk.fileIDs(combWhiskLogical,:);
whiskData = EventData.flow.data.whisk.NormData(combWhiskLogical,:);

% identify the unique days and the unique number of files from the list of unstim resting events
whiskUniqueDays = GetUniqueDays_IOS(whiskFiles);
whiskUniqueFiles = unique(whiskFiles);
whiskNumberOfFiles = length(unique(whiskFiles));

% decimate the file list to only include those files that occur within the desired number of target minutes
clear whiskFiltLogical
for c = 1:length(whiskUniqueDays)
    whiskDay = whiskUniqueDays(c);
    d = 1;
    for e = 1:whiskNumberOfFiles
        whiskFile = whiskUniqueFiles(e);
        whiskFileID = whiskFile{1}(1:6);
        if strcmp(whiskDay,whiskFileID) && sum(strcmp(whiskFile,manualFileIDs)) == 1
            whiskFiltLogical{c,1}(e,1) = 1; %#ok<*AGROW>
            d = d + 1;
        else
            whiskFiltLogical{c,1}(e,1) = 0;
        end
    end
end
whiskFinalLogical = any(sum(cell2mat(whiskFiltLogical'),2),2);

% extract unstim the resting events that correspond to the acceptable file list and the acceptable resting criteria
clear whiskFileFilter
filtWhiskFiles = whiskUniqueFiles(whiskFinalLogical,:);
for f = 1:length(whiskFiles)
    whiskLogic = strcmp(whiskFiles{f},filtWhiskFiles);
    whiskLogicSum = sum(whiskLogic);
    if whiskLogicSum == 1
        whiskFileFilter(f,1) = 1;
    else
        whiskFileFilter(f,1) = 0;
    end
end
whiskFinalFileFilter = logical(whiskFileFilter);
finalWhiskData = whiskData(whiskFinalFileFilter,:);

% only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
% original epoch create so we can add a sample of two back to the end for those just under 10 seconds
% lowpass filter and detrend each segment
[B, A] = butter(3,1/(samplingRate/2),'low');
clear procWhiskData
for g = 1:size(finalWhiskData,1)
    procWhiskData(g,:) = filtfilt(B,A,finalWhiskData(g,:));
end

% analyze correlation coefficient between resting epochs
for n = 1:size(procWhiskData,1)
    whiskFlowMean{n,1} = mean(procWhiskData(n,samplingRate*2:end),2);
end

% save results
AnalysisResults.MeanCBV.Whisk = cell2mat(whiskFlowMean);

%% Analyze mean CBV during periods of NREM sleep
% pull data from SleepData.mat structure
disp('AnalzeMeanFlow: average laser doppler flow during NREM')
nremData = SleepData.NREM.data.DopplerFlow;

% analyze correlation coefficient between NREM epochs
for n = 1:length(nremData)
    nremFlowMean(n,1) = mean(nremData{n,1});
end

% save results
AnalysisResults.MeanCBV.NREM = nremFlowMean;

%% Analyze mean CBV during periods of REM sleep
% pull data from SleepData.mat structure
disp('AnalzeMeanCBV: average hemodynamics during REM')
remData = SleepData.REM.data.DopplerFlow;

% analyze correlation coefficient between NREM epochs
for n = 1:length(remData)
    remFlowMean(n,1) = mean(remData{n,1});
end

% save results
AnalysisResults.MeanCBV.REM = remFlowMean;

%% save results structure
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

end
