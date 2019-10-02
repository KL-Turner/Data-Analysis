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
%
%   Last Revised: September 25th, 2019
%________________________________________________________________________________________________________________________

% find and load RestData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID)

% find and load RestingBaselines.mat struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)

% find and load AllSpecStruct.mat struct
allSpecStructFileStruct = dir('*_AllSpecStruct.mat');
allSpecStructFile = {allSpecStructFileStruct.name}';
allSpecStructFileID = char(allSpecStructFile);
load(allSpecStructFileID)

% determine the animal's ID use the EventData.mat file's name for the current folder
fileBreaks = strfind(eventDataFileID,'_');
animalID = eventDataFileID(1:fileBreaks(1)-1);

% pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
samplingRate = EventData.CBV.(dataType).whisk.samplingRate;
trialDuration_sec = EventData.CBV.(dataType).whisk.trialDuration_sec;
trialDuration_min = trialDuration_sec/samplingRate;
timeVector = (0:(EventData.CBV.(dataType).whisk.epoch.duration*samplingRate))/samplingRate - EventData.CBV.(dataType).whisk.epoch.offset;
offset = EventData.CBV.(dataType).whisk.epoch.offset;
eventWindow = EventData.CBV.(dataType).whisk.epoch.duration;
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
fileTarget = params.targetMinutes/trialDuration_min;
filterSets = {'manualSelection','setDuration','entireDuration'};

for a = 1:length(filterSets)
    filterSet = filterSets{1,a};
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
    
    whiskCriteriaNames = {'whiskCriteriaA','whiskCriteriaB','whiskCriteriaC'};
    
    % filter the EventData.mat structure for whisking events that meet the desired criteria
    for b = 1:length(whiskCriteriaNames)
        whiskCriteriaName = whiskCriteriaNames{1,b};
        if strcmp(whiskCriteriaName,'whiskCriteriaA') == true
            whiskCriteria = whiskCriteriaA;
        elseif strcmp(whiskCriteriaName,'whiskCriteriaB') == true
            whiskCriteria = whiskCriteriaB;
        elseif strcmp(whiskCriteriaName,'whiskCriteriaC') == true
            whiskCriteria = whiskCriteriaC;
        end
        allWhiskFilter = FilterEvents_IOS(EventData.CBV.(dataType).whisk,whiskCriteria);
        [allWhiskCBVData] = EventData.CBV.(dataType).whisk.NormData(allWhiskFilter,:);
        [allWhiskHbTData] = EventData.CBV_HbT.(dataType).whisk.data(allWhiskFilter,:);
        [allWhiskCorticalMUAData] = EventData.(['cortical_' dataType]).muaPower.whisk.NormData(allWhiskFilter,:);
        [allWhiskHippocampalMUAData] = EventData.hippocampus.muaPower.whisk.NormData(allWhiskFilter,:);
        [allWhiskFileIDs] = EventData.CBV.(dataType).whisk.fileIDs(allWhiskFilter,:);
        [allWhiskEventTimes] = EventData.CBV.(dataType).whisk.eventTime(allWhiskFilter,:);
        
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
        finalWhiskCBVData = allWhiskCBVData(finalWhiskFileFilter,:);
        finalWhiskHbTData = allWhiskHbTData(finalWhiskFileFilter,:);
        finalWhiskCorticalMUAData = allWhiskCorticalMUAData(finalWhiskFileFilter,:);
        finalWhiskHippocampalMUAData = allWhiskHippocampalMUAData(finalWhiskFileFilter,:);
        finalWhiskFileIDs = allWhiskFileIDs(finalWhiskFileFilter,:);
        finalWhiskFileEventTimes = allWhiskEventTimes(finalWhiskFileFilter,:);
        
        % lowpass filter each whisking event and mean-subtract by the first 2 seconds
        for f = 1:size(finalWhiskCBVData,1)
            whiskCBVarray = finalWhiskCBVData(f,:);
            whiskHbTarray = finalWhiskHbTData(f,:);
            whiskCorticalMUAarray = finalWhiskCorticalMUAData(f,:);
            whiskHippocampalMUAarray = finalWhiskHippocampalMUAData(f,:);
            filtWhiskCBVarray = sgolayfilt(whiskCBVarray,3,17)*100;
            filtWhiskHbTarray = sgolayfilt(whiskHbTarray,3,17);
            filtWhiskCorticalMUAarray = sgolayfilt(whiskCorticalMUAarray,3,17);
            filtWhiskHippocampalMUAarray = sgolayfilt(whiskHippocampalMUAarray,3,17);
            procWhiskCBVData(f,:) = filtWhiskCBVarray - mean(filtWhiskCBVarray(1:(offset*samplingRate)));
            procWhiskHbTData(f,:) = filtWhiskHbTarray - mean(filtWhiskHbTarray(1:(offset*samplingRate)));
            procWhiskCorticalMUAData(f,:) = filtWhiskCorticalMUAarray - mean(filtWhiskCorticalMUAarray(1:(offset*samplingRate)));
            procWhiskHippocampalMUAData(f,:) = filtWhiskHippocampalMUAarray - mean(filtWhiskHippocampalMUAarray(1:(offset*samplingRate)));
        end
        meanWhiskCBVData = mean(procWhiskCBVData,1);
        stdWhiskCBVData = std(procWhiskCBVData,0,1);
        meanWhiskHbTData = mean(procWhiskHbTData,1);
        stdWhiskHbTData = std(procWhiskHbTData,0,1);
        meanWhiskCorticalMUAData = mean(procWhiskCorticalMUAData,1);
        stdWhiskCorticalMUAData = std(procWhiskCorticalMUAData,0,1);
        meanWhiskHippocampalMUAData = mean(procWhiskHippocampalMUAData,1);
        stdWhiskHippocampalMUAData = std(procWhiskHippocampalMUAData,0,1);
        
        % extract LFP from spectrograms associated with the whisking indecies
        whiskCorticalZhold = [];
        whiskHippocampalZhold = [];
        for g = 1:length(finalWhiskFileIDs)
            % load normalized one-second bin data from each file
            whiskFileID = finalWhiskFileIDs{g,1};
            whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecData.mat'];
            whiskSpecField = ['cortical_' dataType];
            for h = 1:length(AllSpecData.(whiskSpecField).fileIDs)
                if strcmp(AllSpecData.(whiskSpecField).fileIDs{h,1},whiskSpecDataFileID) == true
                    whiskCorticalS_Data = AllSpecData.(whiskSpecField).oneSec.normS{h,1};
                    whiskHippocampalS_Data = AllSpecData.hippocampus.oneSec.normS{h,1};
                    F = AllSpecData.(whiskSpecField).oneSec.F{h,1};
                end
            end
            whiskSLength = size(whiskCorticalS_Data,2);
            whiskBinSize = ceil(whiskSLength/trialDuration_sec);
            whiskSamplingDiff = samplingRate/whiskBinSize;
            
            % find the start time and duration
            whiskDuration = floor(floor(size(meanWhiskCBVData,2))/samplingRate);
            whiskStartTime = floor(floor(finalWhiskFileEventTimes(g,1)*samplingRate)/whiskSamplingDiff);
            if whiskStartTime == 0
                whiskStartTime = 1;
            end
            
            % take the S_data from the start time throughout the duration
            try
                whiskCorticalS_Vals = whiskCorticalS_Data(:,(whiskStartTime - (offset*whiskBinSize)):(whiskStartTime + ((whiskDuration - offset)*whiskBinSize)));
                whiskHippocampalS_Vals = whiskHippocampalS_Data(:,(whiskStartTime - (offset*whiskBinSize)):(whiskStartTime + ((whiskDuration - offset)*whiskBinSize)));
            catch
                whiskCorticalS_Vals = whiskCorticalS_Data(:,end - (whiskDuration*whiskBinSize):end);
                whiskHippocampalS_Vals = whiskHippocampalS_Data(:,end - (whiskDuration*whiskBinSize):end);
            end
            
            % mean subtract each row with detrend
            transpWhiskCorticalS_Vals = whiskCorticalS_Vals';   % Transpose since detrend goes down columns
            transpWhiskHippocampalS_Vals = whiskHippocampalS_Vals';
            dTWhiskCorticalS_Vals = detrend(transpWhiskCorticalS_Vals,'constant');
            dTWhiskHippocampalS_Vals = detrend(transpWhiskHippocampalS_Vals,'constant');
            whiskCorticalZhold = cat(3,whiskCorticalZhold,dTWhiskCorticalS_Vals');   % transpose back to original orientation
            whiskHippocampalZhold = cat(3,whiskHippocampalZhold,dTWhiskHippocampalS_Vals');
        end
        
        % figure time/frequency axis and average each S data matrix through time
        meanWhiskCorticalS = mean(whiskCorticalZhold,3);
        meanWhiskHippocampalS = mean(whiskHippocampalZhold,3);
        T = ((1:size(meanWhiskCorticalS,2))/(eventWindow - offset)) - offset;
        
        % summary figure
        whiskEvoked = figure;
        subplot(2,3,1);
        plot(timeVector,meanWhiskCorticalMUAData,'k')
        hold on
        plot(timeVector,meanWhiskCorticalMUAData + stdWhiskCorticalMUAData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanWhiskCorticalMUAData - stdWhiskCorticalMUAData,'color',colors_IOS('battleship grey'))
        title('Cortex')
        title([animalID ' ' dataType ' ' filterSet ' whisking-evoked averages - ' whiskCriteriaName])
        xlabel('Time (sec)')
        ylabel('Fold-change (Norm Power)')
        axis tight
        axis square
        
        subplot(2,3,4);
        plot(timeVector,meanWhiskHippocampalMUAData,'k')
        hold on
        plot(timeVector,meanWhiskHippocampalMUAData + stdWhiskHippocampalMUAData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanWhiskHippocampalMUAData - stdWhiskHippocampalMUAData,'color',colors_IOS('battleship grey'))
        title([animalID ' ' dataType ' ' filterSet ' whisking-evoked averages - ' whiskCriteriaName])
        title('Hippocampus')
        xlabel('Time (sec)')
        ylabel('Fold-change (Norm Power)')
        axis tight
        axis square
        
        subplot(2,3,2);
        imagesc(T,F,meanWhiskCorticalS)
        title('Cortex')
        xlabel('Time (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        caxis([-0.5 1])
        set(gca,'Ticklength',[0 0])
        axis xy
        axis square
        
        subplot(2,3,5);
        imagesc(T,F,meanWhiskHippocampalS)
        title('Hippocampus')
        xlabel('Time (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        caxis([-0.5 1])
        set(gca,'Ticklength',[0 0])
        axis xy
        axis square

        subplot(2,3,3);
        plot(timeVector,meanWhiskCBVData,'k')
        hold on
        plot(timeVector,meanWhiskCBVData + stdWhiskCBVData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanWhiskCBVData - stdWhiskCBVData,'color',colors_IOS('battleship grey'))
        xlabel('Time (sec)')
        ylabel('\DeltaR/R (%)')
        axis tight
        axis square

        subplot(2,3,6);
        plot(timeVector,meanWhiskHbTData,'k')
        hold on
        plot(timeVector,meanWhiskHbTData + stdWhiskHbTData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanWhiskHbTData - stdWhiskHbTData,'color',colors_IOS('battleship grey'))
        xlabel('Time (sec)')
        ylabel('\DeltaHbT')
        axis tight
        axis square

        % save results
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).CBV.Refl = meanWhiskCBVData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).CBV.ReflStD = stdWhiskCBVData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).CBV.HbT = meanWhiskHbTData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).CBV.HbTStD = stdWhiskHbTData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).MUA.corticalData = meanWhiskCorticalMUAData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).MUA.corticalStD = stdWhiskCorticalMUAData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).MUA.hippocampalData = meanWhiskHippocampalMUAData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).MUA.hippocampalStD = stdWhiskHippocampalMUAData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).timeVector = timeVector;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).LFP.corticalS = meanWhiskCorticalS;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).LFP.hippocampalS = meanWhiskHippocampalS;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).LFP.T = T;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(filterSet).(whiskCriteriaName).LFP.F = F;
        
        % save figure
        [pathstr, ~, ~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Analysis Evoked Averages/'];
        if ~exist(dirpath, 'dir')
            mkdir(dirpath);
        end
        savefig(whiskEvoked, [dirpath animalID '_' dataType '_' filterSet '_' whiskCriteriaName '_WhiskEvokedAverages']);
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
        allStimFilter = FilterEvents_IOS(EventData.CBV.(dataType).stim, stimCriteria);
        [allStimCBVData] = EventData.CBV.(dataType).stim.NormData(allStimFilter,:);
        [allStimHbTData] = EventData.CBV_HbT.(dataType).stim.data(allStimFilter,:);
        [allStimMUAData] = EventData.(['cortical_' dataType]).muaPower.stim.NormData(allStimFilter,:);
        [allStimFileIDs] = EventData.CBV.(dataType).stim.fileIDs(allStimFilter,:);
        [allStimEventTimes] = EventData.CBV.(dataType).stim.eventTime(allStimFilter,:);
        
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
                if strcmp(filterSet,'manualSelection') == true
                    if strcmp(stimDay,stimFileID) && sum(strcmp(stimFile,manualFileIDs)) == 1
                        stimFiltLogical{c,1}(n,1) = 1; %#ok<*AGROW>
                        m = m + 1;
                    else
                        stimFiltLogical{c,1}(n,1) = 0;
                    end
                elseif strcmp(filterSet,'setDuration') == true
                    if strcmp(stimDay,stimFileID) && m <= fileTarget
                        stimFiltLogical{c,1}(n,1) = 1;
                        m = m + 1;
                    else
                        stimFiltLogical{c,1}(n,1) = 0;
                    end
                elseif strcmp(filterSet,'entireDuration') == true
                    if strcmp(stimDay,stimFileID)
                        stimFiltLogical{c,1}(n,1) = 1;
                        m = m + 1;
                    else
                        stimFiltLogical{c,1}(n,1) = 0;
                    end
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
        finalStimCBVData = allStimCBVData(finalStimFileFilter, :);
        finalStimHbTData = allStimHbTData(finalStimFileFilter, :);
        finalStimMUAData = allStimMUAData(finalStimFileFilter, :);
        finalStimFileIDs = allStimFileIDs(finalStimFileFilter, :);
        finalStimFileEventTimes = allStimEventTimes(finalStimFileFilter, :);
        
        % lowpass filter each whisking event and mean-subtract by the first 2 seconds
        for p = 1:size(finalStimCBVData,1)
            stimCBVarray = finalStimCBVData(p,:);
            stimHbTarray = finalStimHbTData(p,:);
            stimMUAarray = finalStimMUAData(p,:);
            filtStimCBVarray = sgolayfilt(stimCBVarray,3,17)*100;
            filtStimHbTarray = sgolayfilt(stimHbTarray,3,17);
            filtStimMUAarray = sgolayfilt(stimMUAarray,3,17);
            procStimCBVData(p,:) = filtStimCBVarray - mean(filtStimCBVarray(1:(offset*samplingRate)));
            procStimHbTData(p,:) = filtStimHbTarray - mean(filtStimHbTarray(1:(offset*samplingRate)));
            procStimMUAData(p,:) = filtStimMUAarray - mean(filtStimMUAarray(1:(offset*samplingRate)));
        end
        meanStimCBVData = mean(procStimCBVData,1);
        stdStimCBVData = std(procStimCBVData,0,1);
        meanStimHbTData = mean(procStimHbTData,1);
        stdStimHbTData = std(procStimHbTData,0,1);
        meanStimMUAData = mean(procStimMUAData,1);
        stdStimMUAData = std(procStimMUAData,0,1);
        
        % extract LFP from spectrograms associated with the stimuli indecies
        stimZhold = [];
        for q = 1:length(finalStimFileIDs)
            % load normalized one-second bin data from each file
            stimFileID = finalStimFileIDs{q,1};
            stimSpecDataFileID = [animalID '_' stimFileID '_SpecData.mat'];
            stimSpecField = ['cortical_' dataType];
            for r = 1:length(AllSpecData.(stimSpecField).fileIDs)
                if strcmp(AllSpecData.(stimSpecField).fileIDs{r,1},stimSpecDataFileID) == true
                    stimS_Data = AllSpecData.(stimSpecField).oneSec.normS{r,1};
                end
            end
            stimSLength = size(stimS_Data,2);
            stimBinSize = ceil(stimSLength/trialDuration_sec);
            stimSamplingDiff = samplingRate/stimBinSize;
            
            % find the start time and duration
            stimDuration = floor(floor(size(meanStimCBVData,2))/samplingRate);
            stimStartTime = floor(floor(finalStimFileEventTimes(q,1)*samplingRate)/stimSamplingDiff);
            if stimStartTime == 0
                stimStartTime = 1;
            end
            
            % take the S_data from the start time throughout the duration
            try
                stimS_Vals = stimS_Data(:,(stimStartTime - (offset*stimBinSize)):(stimStartTime + ((stimDuration - offset)*stimBinSize)));
            catch
                stimS_Vals = stimS_Data(:,end - (stimDuration*stimBinSize):end);
            end
            
            % mean subtract each row with detrend
            transpStimS_Vals = stimS_Vals';   % Transpose since detrend goes down columns
            dTStimS_Vals = detrend(transpStimS_Vals,'constant');
            stimZhold = cat(3,stimZhold,dTStimS_Vals');   % transpose back to original orientation
        end
        
        % figure time/frequency axis and average each S data matrix through time
        meanStimS = mean(stimZhold,3);
        
        % summary figure
        stimEvoked = figure;
        subplot(2,2,1);
        plot(timeVector,meanStimMUAData,'k')
        hold on
        plot(timeVector,meanStimMUAData + stdStimMUAData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanStimMUAData - stdStimMUAData,'color',colors_IOS('battleship grey'))
        title([animalID ' ' dataType ' ' filterSet ' ' solenoid ' stimulus-evoked averages'])
        xlabel('Time (sec)')
        ylabel('Fold-change (Norm Power)')
        axis tight
        axis square

        subplot(2,2,3);
        imagesc(T,F,meanStimS)
        xlabel('Time (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        caxis([-0.5 1])
        set(gca,'Ticklength',[0 0])
        axis xy
        axis square

        subplot(2,2,2);
        plot(timeVector,meanStimCBVData,'k')
        hold on
        plot(timeVector,meanStimCBVData + stdStimCBVData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanStimCBVData - stdStimCBVData,'color',colors_IOS('battleship grey'))
        xlabel('Time (sec)')
        ylabel('\DeltaR/R (%)')
        axis tight
        axis square

        subplot(2,2,4);
        plot(timeVector,meanStimHbTData,'k')
        hold on
        plot(timeVector,meanStimHbTData + stdStimHbTData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanStimHbTData - stdStimHbTData,'color',colors_IOS('battleship grey'))
        xlabel('Time (sec)')
        ylabel('\DeltaHbT')
        axis tight
        axis square

        % save results
        AnalysisResults.EvokedAvgs.Stim.(dataType).(filterSet).(solenoid).CBV.Refl = meanStimCBVData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(filterSet).(solenoid).CBV.ReflStD = stdStimCBVData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(filterSet).(solenoid).CBV.HbT = meanStimHbTData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(filterSet).(solenoid).CBV.HbTStD = stdStimHbTData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(filterSet).(solenoid).MUA.data = meanStimMUAData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(filterSet).(solenoid).MUA.StD = stdStimMUAData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(filterSet).(solenoid).timeVector = timeVector;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(filterSet).(solenoid).LFP.S = meanStimS;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(filterSet).(solenoid).LFP.T = T;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(filterSet).(solenoid).LFP.F = F;
        
        % save figure
        [pathstr, ~, ~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Analysis Evoked Averages/'];
        if ~exist(dirpath, 'dir')
            mkdir(dirpath);
        end
        savefig(stimEvoked, [dirpath animalID '_' dataType '_' filterSet '_' stimCriteriaName '_StimEvokedAverages']);
    end
end
% save results structure
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

end
