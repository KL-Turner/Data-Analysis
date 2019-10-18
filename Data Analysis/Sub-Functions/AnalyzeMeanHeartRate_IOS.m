function [AnalysisResults] = AnalyzeMeanHeartRate_IOS(params,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%________________________________________________________________________________________________________________________

% list of unstim Procdata.mat files
procDataFileStruct = dir('*_Procdata.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID)

% find and load RestData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID)

% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)

% find and load SleepData.mat strut
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID)

% identify animal's ID and pull important infortmat
fileBreaks = strfind(restDataFileID, '_');
animalID = restDataFileID(1:fileBreaks(1)-1);
trialDuration_min = RestData.CBV.LH.trialDuration_sec/60;   % min
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
fileTarget = params.targetMinutes/trialDuration_min;
filterSets = {'manualSelection','setDuration','entireDuration'};

whiskCriteria.Fieldname = {'duration','puffDistance'};
whiskCriteria.Comparison = {'gt','gt'};
whiskCriteria.Value = {5,5};

RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};

PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};

%% Analyze heart rate during long whisking events
for a = 1:length(filterSets)
    filterSet = filterSets{1,a};
    disp('Extracting average extended-whisking data heart rate'); disp(' ')
    allWhiskFilter = FilterEvents_IOS(EventData.CBV.LH.whisk,whiskCriteria);
    [allWhiskFileIDs] = EventData.CBV.LH.whisk.fileIDs(allWhiskFilter,:);
    [allWhiskEventTimes] = EventData.CBV.LH.whisk.eventTime(allWhiskFilter,:);
    [allWhiskDurations] = EventData.CBV.LH.whisk.duration(allWhiskFilter,:);
    
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
            if strcmp(filterSet,'manualSelection') == true
                if strcmp(whiskDay,whiskFileID) && sum(strcmp(whiskFile,manualFileIDs)) == 1
                    whiskFiltLogical{c,1}(n,1) = 1; %#ok<*AGROW>
                    d = d + 1;
                else
                    whiskFiltLogical{c,1}(n,1) = 0;
                end
            elseif strcmp(filterSet,'setDuration') == true
                if strcmp(whiskDay,whiskFileID) && d <= fileTarget
                    whiskFiltLogical{c,1}(n,1) = 1;
                    d = d + 1;
                else
                    whiskFiltLogical{c,1}(n,1) = 0;
                end
            elseif strcmp(filterSet,'entireDuration') == true
                if strcmp(whiskDay,whiskFileID)
                    whiskFiltLogical{c,1}(n,1) = 1;
                    d = d + 1;
                else
                    whiskFiltLogical{c,1}(n,1) = 0;
                end
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
    finalWhiskFileIDs = allWhiskFileIDs(finalWhiskFileFilter,:);
    finalWhiskEventTimes = allWhiskEventTimes(finalWhiskFileFilter,:);
    finalWhiskDurations = allWhiskDurations(finalWhiskFileFilter,:);
    
    clear whiskingHeartRate
    for a = 1:length(finalWhiskFileIDs)
        whiskFileID = [animalID '_' finalWhiskFileIDs{a,1} '_ProcData.mat'];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            if strcmp(whiskFileID,procDataFileID) == true
                load(whiskFileID)
                heartRate = ProcData.data.heartRate;
                eventTime = floor(finalWhiskEventTimes(a,1));
                duration = floor(finalWhiskDurations(a,1));
                try
                    whiskingHeartRate(a,1) = mean(heartRate(eventTime:eventTime + duration));
                catch
                    whiskingHeartRate(a,1) = mean(heartRate(1:eventTime + duration));
                end
                break
            end
        end
    end
    disp(['Average whisking heart rate (' filterSet '): ' num2str(mean(whiskingHeartRate)) ' Hz.']); disp(' ')
    
    % save results
    AnalysisResults.MeanHR.Whisk.(filterSet) = whiskingHeartRate;
end

%% Analyze heart rate during rest data
for b = 1:length(filterSets)
    filterSet = filterSets{1,b};
    disp(['Extracting average awake heart rate - ' filterSet]); disp(' ')
    %% Analyze coherence during periods of rest
    % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
    [restLogical] = FilterEvents_IOS(RestData.CBV.LH,RestCriteria);
    [puffLogical] = FilterEvents_IOS(RestData.CBV.LH,PuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFiles = RestData.CBV.LH.fileIDs(combRestLogical,:);
    restEventTimes = RestData.CBV.LH.eventTimes(combRestLogical,:);
    restDurations = RestData.CBV.LH.durations(combRestLogical,:);
    
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
            if strcmp(filterSet,'manualSelection') == true
                if strcmp(restDay,restFileID) && sum(strcmp(restFile,manualFileIDs)) == 1
                    restFiltLogical{c,1}(e,1) = 1; %#ok<*AGROW>
                    d = d + 1;
                else
                    restFiltLogical{c,1}(e,1) = 0;
                end
            elseif strcmp(filterSet,'setDuration') == true
                if strcmp(restDay,restFileID) && d <= fileTarget
                    restFiltLogical{c,1}(e,1) = 1;
                    d = d + 1;
                else
                    restFiltLogical{c,1}(e,1) = 0;
                end
            elseif strcmp(filterSet,'entireDuration') == true
                if strcmp(restDay,restFileID)
                    restFiltLogical{c,1}(e,1) = 1;
                    d = d + 1;
                else
                    restFiltLogical{c,1}(e,1) = 0;
                end
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
    finalRestFileList = restFiles(restFinalFileFilter,:);
    finalRestEventTimes = restEventTimes(restFinalFileFilter,:);
    finalRestDurations = restDurations(restFinalFileFilter,:);
    
    clear restingHeartRate
    for a = 1:length(finalRestFileList)
        restFileID = [animalID '_' finalRestFileList{a,1} '_ProcData.mat'];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            if strcmp(restFileID,procDataFileID) == true
                load(restFileID)
                heartRate = ProcData.data.heartRate;
                eventTime = floor(finalRestEventTimes(a,1));
                duration = floor(finalRestDurations(a,1));
                try
                    restingHeartRate(a,1) = mean(heartRate(eventTime:eventTime + duration));
                catch
                    restingHeartRate(a,1) = mean(heartRate(1:eventTime + duration));
                end
                break
            end
        end
    end
    disp(['Average awake heart rate (' filterSet '): ' num2str(mean(restingHeartRate)) ' Hz.']); disp(' ')
    
    % save results
    AnalysisResults.MeanHR.Rest.(filterSet) = restingHeartRate;
end

%% Analyze heart rate during periods of NREM sleep
% pull data from SleepData.mat structure
disp('Extracting average NREM heart rate'); disp(' ')
nremData = SleepData.NREM.data.HeartRate;

% analyze correlation coefficient between NREM epochs
for n = 1:length(nremData)
    nremHRMean(n,1) = mean(nremData{n,1});
end
disp(['Average NREM heart rate: ' num2str(mean(nremHRMean)) ' Hz.']); disp(' ')

% save results
AnalysisResults.MeanHR.NREM = nremHRMean;

%% Analyze heart rate during periods of REM sleep
% pull data from SleepData.mat structure
disp('Extracting average REM heart rate'); disp(' ')
remData = SleepData.REM.data.HeartRate;

% analyze correlation coefficient between REM epochs
for n = 1:length(remData)
    remHRMean(n,1) = mean(remData{n,1});
end
disp(['Average REM heart rate: ' num2str(mean(remHRMean)) ' Hz.']); disp(' ')

% save results
AnalysisResults.MeanHR.REM = remHRMean;

%% Analyze heart rate during unstim unstimulated data
disp('Extracting average unstimulated data heart rate'); disp(' ')
for o = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(o,:);
    load(procDataFileID);
    if isempty(ProcData.data.solenoids.LPadSol) == true
        stimLogical(o,1) = 1;
    else
        stimLogical(o,1) = 0;
    end
end
stimLogical = logical(stimLogical);
unstim_procDataFileIDs = procDataFileIDs(stimLogical,:);
for p = 1:size(unstim_procDataFileIDs)
    unstim_procDataFileID = unstim_procDataFileIDs(p,:);
    load(unstim_procDataFileID)
    % pull data from each file
    unstimData{p,1} = ProcData.data.heartRate;
end

% analyze correlation coefficient between NREM epochs
for n = 1:length(unstimData)
    unstimDataMean(n,1) = mean(unstimData{n,1});
end
disp(['Average unstimulated data heart rate: ' num2str(mean(unstimDataMean)) ' Hz.']); disp(' ')

% save results
AnalysisResults.MeanHR.Unstim = unstimDataMean;

%% Analyze heart rate during all data
disp('Extracting average all data heart rate'); disp(' ')
for p = 1:size(procDataFileIDs)
    procDataFileID = procDataFileIDs(p,:);
    load(procDataFileID)
    % pull data from each file
    allData{p,1} = ProcData.data.heartRate;
end

% analyze correlation coefficient between NREM epochs
for n = 1:length(allData)
    allDataMean(n,1) = mean(allData{n,1});
end
disp(['Average all data heart rate: ' num2str(mean(allDataMean)) ' Hz.']); disp(' ')

% save results
AnalysisResults.MeanHR.All = allDataMean;

%% save results strucure
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

end
