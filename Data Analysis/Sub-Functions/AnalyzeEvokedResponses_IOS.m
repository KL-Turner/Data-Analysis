function [AnalysisResults] = AnalyzeEvokedResponses_IOS(dataType, params, AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Use epochs from the EventData.mat struct to determine the average hemodynamic and neural responses to
%            both volitional whisking and whisker stimuli
%________________________________________________________________________________________________________________________
%
%   Inputs: dataType (character) of the hemisphere to analyze. Typically this is labeled as LH or RH
%           params (struct) - how many minutes of each day to use
%           AnalysisResults (struct) - where to save the data for later between-animal averages and comparisons
%
%   Outputs:  AnalysisResults (struct) - where to save the data for later between-animal averages and comparisons
%________________________________________________________________________________________________________________________

% find and load RestData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID)

% determine the animal's ID use the EventData.mat file's name for the current folder
fileBreaks = strfind(eventDataFileID, '_');
animalID = eventDataFileID(1:fileBreaks(1)-1);

% pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
samplingRate = EventData.CBV.(dataType).whisk.samplingRate;
trialDuration_sec = EventData.CBV.(dataType).whisk.trialDuration_sec;
trialDurationMin = trialDuration_sec/samplingRate;
timeVector = (0:(EventData.CBV.(dataType).whisk.epoch.duration*samplingRate))/samplingRate - EventData.CBV.(dataType).whisk.epoch.offset; % 12 seconds
offset = EventData.CBV.(dataType).whisk.epoch.offset;
eventWindow = EventData.CBV.(dataType).whisk.epoch.duration;

%% Whisking-evoked responses
% criteria for the FilterEvents data struct
whiskCriteriaA.Fieldname = {'duration', 'duration', 'puffDistance'};
whiskCriteriaA.Comparison = {'gt','lt','gt'};
whiskCriteriaA.Value = {0.5, 2, 5};

whiskCriteriaB.Fieldname = {'duration', 'duration', 'puffDistance'};
whiskCriteriaB.Comparison = {'gt','lt','gt'};
whiskCriteriaB.Value = {2, 5, 5};

whiskCriteriaC.Fieldname = {'duration', 'puffDistance'};
whiskCriteriaC.Comparison = {'gt','gt'};
whiskCriteriaC.Value = {5, 5};

whiskCriteriaNames = {'whiskCriteriaA', 'whiskCriteriaB', 'whiskCriteriaC'};

% filter the EventData.mat structure for whisking events that meet the desired criteria
for a = 1:length(whiskCriteriaNames)
    whiskCriteriaName = whiskCriteriaNames{1,a};
    if strcmp(whiskCriteriaName, 'whiskCriteriaA') == true
        whiskCriteria = whiskCriteriaA;
    elseif strcmp(whiskCriteriaName, 'whiskCriteriaB') == true
        whiskCriteria = whiskCriteriaB;
    elseif strcmp(whiskCriteriaName, 'whiskCriteriaC') == true
        whiskCriteria = whiskCriteriaC;
    end
    allWhiskFilter = FilterEvents_IOS(EventData.CBV.(dataType).whisk, whiskCriteria);
    [allWhiskCBVData] = EventData.CBV.(dataType).whisk.NormData(allWhiskFilter, :);
    [allWhiskHbTData] = EventData.CBV_HbT.(dataType).whisk.data(allWhiskFilter, :);
    [allWhiskFileIDs] = EventData.CBV.(dataType).whisk.fileIDs(allWhiskFilter, :);
    [allWhiskEventTimes] = EventData.CBV.(dataType).whisk.eventTime(allWhiskFilter, :);
    
    % identify the unique days and the unique number of files from the list of all whisking events
    whiskUniqueDays = GetUniqueDays_IOS(allWhiskFileIDs);
    whiskUniqueFiles = unique(allWhiskFileIDs);
    whiskNumberOfFiles = length(unique(allWhiskFileIDs));
    fileTarget = params.targetMinutes/trialDurationMin;
    
    % decimate the file list to only include those files that occur within the desired number of target minutes
    % this is done to exclude the data that occurs later in the trial when the animal is sleeping
    for b = 1:length(whiskUniqueDays)
        whiskDay = whiskUniqueDays(b);
        c = 1;
        for nOF = 1:whiskNumberOfFiles
            file = whiskUniqueFiles(nOF);
            whiskFileID = file{1}(1:6);
            if strcmp(whiskDay, whiskFileID) && c <= fileTarget
                whiskFiltLogical{b,1}(nOF,1) = 1; %#ok<*AGROW>
                c = c + 1;
            else
                whiskFiltLogical{b,1}(nOF,1) = 0;
            end
        end
    end
    whiskFinalLogical = any(sum(cell2mat(whiskFiltLogical'), 2), 2);
    
    clear whiskFileFilter
    % extract all the whisking events that correspond to the acceptable file list and the acceptable whisking criteria
    filtWhiskFiles = whiskUniqueFiles(whiskFinalLogical,:);
    for d = 1:length(allWhiskFileIDs)
        logic = strcmp(allWhiskFileIDs{d}, filtWhiskFiles);
        logicSum = sum(logic);
        if logicSum == 1
            whiskFileFilter(d,1) = 1;
        else
            whiskFileFilter(d,1) = 0;
        end
    end
    finalWhiskFileFilter = logical(whiskFileFilter);
    finalWhiskCBVData = allWhiskCBVData(finalWhiskFileFilter,:);
    finalWhiskHbTData = allWhiskHbTData(finalWhiskFileFilter,:);
    finalWhiskFileIDs = allWhiskFileIDs(finalWhiskFileFilter,:);
    finalWhiskFileEventTimes = allWhiskEventTimes(finalWhiskFileFilter,:);
    
    % lowpass filter each whisking event and mean-subtract by the first 2 seconds
    for e = 1:size(finalWhiskCBVData,1)
        whiskCBVarray = finalWhiskCBVData(e,:);
        whiskHbTarray = finalWhiskHbTData(e,:);
        filtWhiskCBVarray = sgolayfilt(whiskCBVarray,3,17)*100;
        filtWhiskHbTarray = sgolayfilt(whiskHbTarray,3,17);
        procWhiskCBVData(e,:) = filtWhiskCBVarray - mean(filtWhiskCBVarray(1:(offset*samplingRate)));
        procWhiskHbTData(e,:) = filtWhiskHbTarray - mean(filtWhiskHbTarray(1:(offset*samplingRate)));
    end
    meanWhiskCBVData = mean(procWhiskCBVData,1);
    stdWhiskCBVData = std(procWhiskCBVData,0,1);
    meanWhiskHbTData = mean(procWhiskHbTData,1);
    stdWhiskHbTData = std(procWhiskHbTData,0,1);

    % extract LFP from spectrograms associated with the whisking indecies
    whiskZhold = [];
    for f = 1:length(finalWhiskFileIDs)
        % load normalized one-second bin data from each file
        whiskFileID = finalWhiskFileIDs{f, 1};
        specDataFileID = [animalID '_' whiskFileID '_SpecData.mat'];
        specField = ['cortical_' dataType];
        load(specDataFileID)
        whiskS_Data = SpecData.(specField).oneSec.normS;
        whiskSLength = size(whiskS_Data,2);
        whiskBinSize = ceil(whiskSLength/trialDuration_sec);
        whiskSamplingDiff = samplingRate/whiskBinSize;
        
        % Find the start time and duration
        whiskDuration = floor(floor(size(meanWhiskCBVData,2))/samplingRate);
        startTime = floor(floor(finalWhiskFileEventTimes(f,1)*samplingRate)/whiskSamplingDiff);
        if startTime == 0
            startTime = 1;
        end
        
        % Take the S_data from the start time throughout the duration
        try
            whiskS_Vals = whiskS_Data(:,(startTime - (offset*whiskBinSize)):(startTime + ((whiskDuration - offset)*whiskBinSize)));
        catch
            whiskS_Vals = whiskS_Data(:,end - (whiskDuration*whiskBinSize):end);
        end
        
        % Mean subtract each row with detrend
        transpWhiskS_Vals = whiskS_Vals';   % Transpose since detrend goes down columns
        dTWhiskS_Vals = detrend(transpWhiskS_Vals);
        whiskZhold = cat(3, whiskZhold, dTWhiskS_Vals');   % transpose back to original orientation
    end
    
    % figure time/frequency axis and average each S data matrix through time
    whiskS = mean(whiskZhold, 3);
    T = ((1:size(whiskS,2))/(eventWindow-offset))-offset;
    F = SpecData.(specField).oneSec.F;
    
    % summary figure
    whiskEvoked = figure;
    subplot(3,1,1);
    imagesc(T, F, whiskS)
    title([animalID ' ' dataType ' Whisking-evoked averages - ' whiskCriteriaName])
    xlabel('Time (sec)')
    ylabel('Freq (Hz)')
    ylim([1 100])
    caxis([-0.5 1])
    set(gca, 'Ticklength', [0 0])
    axis xy
    axis square
    
    subplot(3,1,2);
    plot(timeVector, meanWhiskCBVData, 'k')
    hold on
    plot(timeVector, meanWhiskCBVData + stdWhiskCBVData, 'color', colors_IOS('battleship grey'))
    plot(timeVector, meanWhiskCBVData - stdWhiskCBVData, 'color', colors_IOS('battleship grey'))
    xlabel('Time (sec)')
    ylabel('\DeltaR/R (%)')
    axis tight
    axis square
    
    subplot(3,1,3);
    plot(timeVector, meanWhiskHbTData)
    hold on
    plot(timeVector, meanWhiskHbTData + stdWhiskHbTData, 'color', colors_IOS('battleship grey'))
    plot(timeVector, meanWhiskHbTData - stdWhiskHbTData, 'color', colors_IOS('battleship grey'))
    xlabel('Time (sec)')
    ylabel('\DeltaHbT')
    axis tight
    axis square
    
    % save results
    AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV.Refl = meanWhiskCBVData;
    AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV.ReflStD = stdWhiskCBVData;
    AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV.HbT = meanWhiskHbTData;
    AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV.HbTStD = stdWhiskHbTData;
    AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV.timeVector = timeVector;
    AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.S = whiskS;
    AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.T = T;
    AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.F = F;
    
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Analysis Evoked Averages/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(whiskEvoked, [dirpath animalID '_' dataType '_' whiskCriteriaName '_WhiskEvokedAverages']);
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

stimCriteriaNames = {'stimCriteriaA', 'stimCriteriaB', 'stimCriteriaC'};

% filter the EventData.mat structure for stimulus events that meet the desired criteria
for a = 1:length(stimCriteriaNames)
    stimCriteriaName = stimCriteriaNames{1,a};
    if strcmp(stimCriteriaName, 'stimCriteriaA') == true
        stimCriteria = stimCriteriaA;
        solenoid = 'LPadSol';
    elseif strcmp(stimCriteriaName, 'stimCriteriaB') == true
        stimCriteria = stimCriteriaB;
        solenoid = 'RPadSol';
    elseif strcmp(stimCriteriaName, 'stimCriteriaC') == true
        stimCriteria = stimCriteriaC;
        solenoid = 'AudSol';
    end
    allStimFilter = FilterEvents_IOS(EventData.CBV.(dataType).stim, stimCriteria);
    [allStimCBVData] = EventData.CBV.(dataType).stim.NormData(allStimFilter, :);
    [allStimHbTData] = EventData.CBV_HbT.(dataType).stim.data(allStimFilter, :);
    [allStimFileIDs] = EventData.CBV.(dataType).stim.fileIDs(allStimFilter, :);
    [allStimEventTimes] = EventData.CBV.(dataType).stim.eventTime(allStimFilter, :);
    
    % identify the unique days and the unique number of files from the list of all whisking events
    stimUniqueDays = GetUniqueDays_IOS(allStimFileIDs);
    stimUniqueFiles = unique(allStimFileIDs);
    stimNumberOfFiles = length(unique(allStimFileIDs));
    fileTarget = params.targetMinutes/trialDurationMin;
    
    % decimate the file list to only include those files that occur within the desired number of target minutes
    % this is done to exclude the data that occurs later in the trial when the animal is sleeping
    for g = 1:length(stimUniqueDays)
        stimDay = stimUniqueDays(g);
        h = 1;
        for nOF = 1:stimNumberOfFiles
            file = stimUniqueFiles(nOF);
            stimFileID = file{1}(1:6);
            if strcmp(stimDay, stimFileID) && h <= fileTarget
                stimFiltLogical{g, 1}(nOF, 1) = 1;
                h = h + 1;
            else
                stimFiltLogical{g, 1}(nOF, 1) = 0;
            end
        end
    end
    stimFinalLogical = any(sum(cell2mat(stimFiltLogical'), 2), 2);
        
    clear stimFileFilter
    % extract all the whisking events that correspond to the acceptable file list and the acceptable whisking criteria
    filtStimFiles = stimUniqueFiles(stimFinalLogical, :);
    for j = 1:length(allStimFileIDs)
        logic = strcmp(allStimFileIDs{j}, filtStimFiles);
        logicSum = sum(logic);
        if logicSum == 1
            stimFileFilter(j, 1) = 1;
        else
            stimFileFilter(j, 1) = 0;
        end
    end
    finalStimFileFilter = logical(stimFileFilter);
    finalStimCBVData = allStimCBVData(finalStimFileFilter, :);
    finalStimHbTData = allStimHbTData(finalStimFileFilter, :);
    finalStimFileIDs = allStimFileIDs(finalStimFileFilter, :);
    finalStimFileEventTimes = allStimEventTimes(finalStimFileFilter, :);
    
    % lowpass filter each whisking event and mean-subtract by the first 2 seconds
    for k = 1:size(finalStimCBVData,1)
        stimCBVarray = finalStimCBVData(k,:);
        stimHbTarray = finalStimHbTData(k,:);
        filtStimCBVarray = sgolayfilt(stimCBVarray,3,17)*100;
        filtStimHbTarray = sgolayfilt(stimHbTarray,3,17);
        procStimCBVData(k,:) = filtStimCBVarray - mean(filtStimCBVarray(1:(offset*samplingRate)));
        procStimHbTData(k,:) = filtStimHbTarray - mean(filtStimHbTarray(1:(offset*samplingRate)));
    end
    meanStimCBVData = mean(procStimCBVData,1);
    stdStimCBVData = std(procStimCBVData,0,1);
    meanStimHbTData = mean(procStimHbTData,1);
    stdStimHbTData = std(procStimHbTData,0,1);

    % extract LFP from spectrograms associated with the stimuli indecies
    stimZhold = [];
    for m = 1:length(finalStimFileIDs)
        % load normalized one-second bin data from each file
        stimFileID = finalStimFileIDs{m, 1};
        specDataFileID = [animalID '_' stimFileID '_SpecData.mat'];
        specField = ['cortical_' dataType];
        load(specDataFileID)
        stimS_Data = SpecData.(specField).oneSec.normS;
        stimSLength = size(stimS_Data, 2);
        stimBinSize = ceil(stimSLength/trialDuration_sec);
        stimSamplingDiff = samplingRate/stimBinSize;
        
        % find the start time and duration
        stimDuration = floor(floor(size(meanStimCBVData,2))/samplingRate);
        startTime = floor(floor(finalStimFileEventTimes(m,1)*samplingRate)/stimSamplingDiff);
        if startTime == 0
            startTime = 1;
        end
        
        % take the S_data from the start time throughout the duration
        try
            stimS_Vals = stimS_Data(:,(startTime - (offset*stimBinSize)):(startTime + ((stimDuration - offset)*stimBinSize)));
        catch
            stimS_Vals = stimS_Data(:,end - (stimDuration*stimBinSize):end);
        end
        
        % mean subtract each row with detrend
        transpStimS_Vals = stimS_Vals';   % Transpose since detrend goes down columns
        dTStimS_Vals = detrend(transpStimS_Vals);
        stimZhold = cat(3,stimZhold,dTStimS_Vals');   % transpose back to original orientation
    end
    
    % figure time/frequency axis and average each S data matrix through time
    stimS = mean(stimZhold,3);
    
    % summary figure
    stimEvoked = figure;
    subplot(3,1,1);
    imagesc(T, F, stimS)
    title([animalID ' ' dataType ' ' solenoid ' Stimulus-evoked averages'])
    xlabel('Time (sec)')
    ylabel('Freq (Hz)')
    ylim([1 100])
    caxis([-0.5 1])
    set(gca, 'Ticklength', [0 0])
    axis xy
    axis square
    
    subplot(3,1,2);
    plot(timeVector, meanStimCBVData, 'k')
    hold on
    plot(timeVector, meanStimCBVData + stdStimCBVData, 'color', colors_IOS('battleship grey'))
    plot(timeVector, meanStimCBVData - stdStimCBVData, 'color', colors_IOS('battleship grey'))
    xlabel('Time (sec)')
    ylabel('\DeltaR/R (%)')
    axis tight
    axis square
    
    subplot(3,1,3);
    plot(timeVector, meanStimHbTData)
    hold on
    plot(timeVector, meanStimHbTData + stdStimHbTData, 'color', colors_IOS('battleship grey'))
    plot(timeVector, meanStimHbTData - stdStimHbTData, 'color', colors_IOS('battleship grey'))
    xlabel('Time (sec)')
    ylabel('\DeltaHbT')
    axis tight
    axis square

    % save results
    AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).CBV.Refl = meanStimCBVData;
    AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).CBV.ReflStD = stdStimCBVData;
    AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).CBV.HbT = meanStimHbTData;
    AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).CBV.HbTStD = stdStimHbTData;
    AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).CBV.timeVector = timeVector;
    AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).LFP.S = stimS;
    AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).LFP.T = T;
    AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).LFP.F = F;
    
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Analysis Evoked Averages/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(stimEvoked, [dirpath animalID '_' dataType '_' stimCriteriaName '_StimEvokedAverages']);
end

save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

end

