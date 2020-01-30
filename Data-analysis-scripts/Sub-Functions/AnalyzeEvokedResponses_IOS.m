function [AnalysisResults] = AnalyzeEvokedResponses_IOS(dataTypes,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner%
%
%   Purpose: Use epochs from the EventData.mat struct to determine the average hemodynamic and neural responses to
%            both volitional whisking and whisker stimuli
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
for a = 1:length(dataTypes)
    dataType = dataTypes{1,a};
    neuralDataType = ['cortical_' dataType(4:end)];
    % pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
    samplingRate = EventData.CBV_HbT.(dataType).whisk.samplingRate;
    trialDuration_sec = EventData.CBV_HbT.(dataType).whisk.trialDuration_sec;
    timeVector = (0:(EventData.CBV_HbT.(dataType).whisk.epoch.duration*samplingRate))/samplingRate - EventData.CBV_HbT.(dataType).whisk.epoch.offset;
    offset = EventData.CBV_HbT.(dataType).whisk.epoch.offset;
    eventWindow = EventData.CBV_HbT.(dataType).whisk.epoch.duration;
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
        allWhiskFilter = FilterEvents_IOS(EventData.CBV_HbT.(dataType).whisk,whiskCriteria);
        [allWhiskHbTData] = EventData.CBV_HbT.(dataType).whisk.data(allWhiskFilter,:);
        [allWhiskCBVData] = EventData.CBV.(dataType).whisk.NormData(allWhiskFilter,:);
        [allWhiskCorticalMUAData] = EventData.(neuralDataType).muaPower.whisk.NormData(allWhiskFilter,:);
        [allWhiskHippocampalMUAData] = EventData.hippocampus.muaPower.whisk.NormData(allWhiskFilter,:);
        [allWhiskFileIDs] = EventData.CBV_HbT.(dataType).whisk.fileIDs(allWhiskFilter,:);
        [allWhiskEventTimes] = EventData.CBV_HbT.(dataType).whisk.eventTime(allWhiskFilter,:);
        
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
        finalWhiskHbTData = allWhiskHbTData(finalWhiskFileFilter,:);
        finalWhiskCBVData = allWhiskCBVData(finalWhiskFileFilter,:);
        finalWhiskCorticalMUAData = allWhiskCorticalMUAData(finalWhiskFileFilter,:);
        finalWhiskHippocampalMUAData = allWhiskHippocampalMUAData(finalWhiskFileFilter,:);
        finalWhiskFileIDs = allWhiskFileIDs(finalWhiskFileFilter,:);
        finalWhiskFileEventTimes = allWhiskEventTimes(finalWhiskFileFilter,:);
        
        % lowpass filter each whisking event and mean-subtract by the first 2 seconds
        clear procWhiskHbTData procWhiskCBVData procWhiskCorticalMUAData procWhiskHippocampalMUAData
        for f = 1:size(finalWhiskHbTData,1)
            whiskHbTarray = finalWhiskHbTData(f,:);
            whiskCBVarray = finalWhiskCBVData(f,:);
            whiskCorticalMUAarray = finalWhiskCorticalMUAData(f,:);
            whiskHippocampalMUAarray = finalWhiskHippocampalMUAData(f,:);
            filtWhiskHbTarray = sgolayfilt(whiskHbTarray,3,17);
            filtWhiskCBVarray = sgolayfilt(whiskCBVarray,3,17);
            filtWhiskCorticalMUAarray = sgolayfilt(whiskCorticalMUAarray,3,17);
            filtWhiskHippocampalMUAarray = sgolayfilt(whiskHippocampalMUAarray,3,17);
            procWhiskHbTData(f,:) = filtWhiskHbTarray - mean(filtWhiskHbTarray(1:(offset*samplingRate)));
            procWhiskCBVData(f,:) = filtWhiskCBVarray - mean(filtWhiskCBVarray(1:(offset*samplingRate)));
            procWhiskCorticalMUAData(f,:) = filtWhiskCorticalMUAarray - mean(filtWhiskCorticalMUAarray(1:(offset*samplingRate)));
            procWhiskHippocampalMUAData(f,:) = filtWhiskHippocampalMUAarray - mean(filtWhiskHippocampalMUAarray(1:(offset*samplingRate)));
        end
        meanWhiskHbTData = mean(procWhiskHbTData,1);
        stdWhiskHbTData = std(procWhiskHbTData,0,1);
        meanWhiskCBVData = mean(procWhiskCBVData,1)*100;
        stdWhiskCBVData = std(procWhiskCBVData,0,1)*100;
        meanWhiskCorticalMUAData = mean(procWhiskCorticalMUAData,1)*100;
        stdWhiskCorticalMUAData = std(procWhiskCorticalMUAData,0,1)*100;
        meanWhiskHippocampalMUAData = mean(procWhiskHippocampalMUAData,1)*100;
        stdWhiskHippocampalMUAData = std(procWhiskHippocampalMUAData,0,1)*100;
        
        % extract LFP from spectrograms associated with the whisking indecies
        whiskCorticalZhold = [];
        whiskHippocampalZhold = [];
        for g = 1:length(finalWhiskFileIDs)
            % load normalized one-second bin data from each file
            whiskFileID = finalWhiskFileIDs{g,1};
            whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecData.mat'];
            whiskSpecField = neuralDataType;
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
            whiskDuration = floor(floor(size(meanWhiskHbTData,2))/samplingRate);
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
        sgtitle([animalID ' ' dataType ' whisking-evoked averages - ' whiskCriteriaName])
        subplot(2,3,1);
        plot(timeVector,meanWhiskCorticalMUAData,'k')
        hold on
        plot(timeVector,meanWhiskCorticalMUAData + stdWhiskCorticalMUAData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanWhiskCorticalMUAData - stdWhiskCorticalMUAData,'color',colors_IOS('battleship grey'))
        title('Cortical MUA')
        xlabel('Time (sec)')
        ylabel('Fold-change (Norm Power)')
        axis tight
        axis square
        
        subplot(2,3,4);
        plot(timeVector,meanWhiskHippocampalMUAData,'k')
        hold on
        plot(timeVector,meanWhiskHippocampalMUAData + stdWhiskHippocampalMUAData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanWhiskHippocampalMUAData - stdWhiskHippocampalMUAData,'color',colors_IOS('battleship grey'))
        title([animalID ' ' dataType ' ' whiskCriteriaName ' whisking-evoked averages'])
        title('Hippocampal MUA')
        xlabel('Time (sec)')
        ylabel('Fold-change (Norm Power)')
        axis tight
        axis square
        
        subplot(2,3,2);
        imagesc(T,F,meanWhiskCorticalS)
        title('Cortical LFP')
        xlabel('Time (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        caxis([-0.5 1])
        set(gca,'Ticklength',[0 0])
        axis xy
        axis square
        
        subplot(2,3,5);
        imagesc(T,F,meanWhiskHippocampalS)
        title('Hippocampal LFP')
        xlabel('Time (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        caxis([-0.5 1])
        set(gca,'Ticklength',[0 0])
        axis xy
        axis square
        
        subplot(2,3,[3,6]);
        plot(timeVector,meanWhiskHbTData,'k')
        hold on
        plot(timeVector,meanWhiskHbTData + stdWhiskHbTData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanWhiskHbTData - stdWhiskHbTData,'color',colors_IOS('battleship grey'))
        title('Hemodynamic response')
        xlabel('Time (sec)')
        ylabel('\DeltaHbT (\muM)')
        axis tight
        axis square
        
        % save results
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV_HbT.HbT = meanWhiskHbTData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV_HbT.HbTStD = stdWhiskHbTData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV.CBV = meanWhiskCBVData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV.CBVStD = stdWhiskCBVData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).MUA.corticalData = meanWhiskCorticalMUAData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).MUA.corticalStD = stdWhiskCorticalMUAData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).MUA.hippocampalData = meanWhiskHippocampalMUAData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).MUA.hippocampalStD = stdWhiskHippocampalMUAData;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).timeVector = timeVector;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.corticalS = meanWhiskCorticalS;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.hippocampalS = meanWhiskHippocampalS;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.T = T;
        AnalysisResults.EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.F = F;
        
        % save figure
        [pathstr, ~, ~] = fileparts(cd);
        dirpath = [pathstr '/Combined Imaging/Figures/Stim and Whisk Evoked Averages/'];
        if ~exist(dirpath, 'dir')
            mkdir(dirpath);
        end
        savefig(whiskEvoked, [dirpath animalID '_' dataType '_' whiskCriteriaName '_WhiskEvokedAverages']);
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
        allStimFilter = FilterEvents_IOS(EventData.CBV_HbT.(dataType).stim,stimCriteria);
        [allStimHbTData] = EventData.CBV_HbT.(dataType).stim.data(allStimFilter,:);
        [allStimCBVData] = EventData.CBV.(dataType).stim.NormData(allStimFilter,:);
        [allStimCortMUAData] = EventData.(neuralDataType).muaPower.stim.NormData(allStimFilter,:);
        [allStimHipMUAData] = EventData.hippocampus.muaPower.stim.NormData(allStimFilter,:);
        [allStimFileIDs] = EventData.CBV_HbT.(dataType).stim.fileIDs(allStimFilter,:);
        [allStimEventTimes] = EventData.CBV_HbT.(dataType).stim.eventTime(allStimFilter,:);
        
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
        finalStimHbTData = allStimHbTData(finalStimFileFilter,:);
        finalStimCBVData = allStimCBVData(finalStimFileFilter,:);
        finalStimCortMUAData = allStimCortMUAData(finalStimFileFilter,:);
        finalStimHipMUAData = allStimHipMUAData(finalStimFileFilter,:);
        finalStimFileIDs = allStimFileIDs(finalStimFileFilter,:);
        finalStimFileEventTimes = allStimEventTimes(finalStimFileFilter,:);
        
        % lowpass filter each whisking event and mean-subtract by the first 2 seconds
        clear procStimHbTData procStimCBVData procStimCortMUAData procStimHipMUAData
        for p = 1:size(finalStimHbTData,1)
            stimHbTarray = finalStimHbTData(p,:);
            stimCBVarray = finalStimCBVData(p,:);
            stimCortMUAarray = finalStimCortMUAData(p,:);
            stimHipMUAarray = finalStimHipMUAData(p,:);
            filtStimHbTarray = sgolayfilt(stimHbTarray,3,17);
            filtStimCBVarray = sgolayfilt(stimCBVarray,3,17)*100;
            filtStimCortMUAarray = sgolayfilt(stimCortMUAarray,3,17);
            filtStimHipMUAarray = sgolayfilt(stimHipMUAarray,3,17);
            procStimHbTData(p,:) = filtStimHbTarray - mean(filtStimHbTarray(1:(offset*samplingRate)));
            procStimCBVData(p,:) = filtStimCBVarray - mean(filtStimCBVarray(1:(offset*samplingRate)));
            procStimCortMUAData(p,:) = filtStimCortMUAarray - mean(filtStimCortMUAarray(1:(offset*samplingRate)));
            procStimHipMUAData(p,:) = filtStimHipMUAarray - mean(filtStimHipMUAarray(1:(offset*samplingRate)));
        end
        meanStimHbTData = mean(procStimHbTData,1);
        stdStimHbTData = std(procStimHbTData,0,1);
        meanStimCBVData = mean(procStimCBVData,1)*100;
        stdStimCBVData = std(procStimCBVData,0,1)*100;
        meanStimCortMUAData = mean(procStimCortMUAData,1)*100;
        stdStimCortMUAData = std(procStimCortMUAData,0,1)*100;
        meanStimHipMUAData = mean(procStimHipMUAData,1)*100;
        stdStimHipMUAData = std(procStimHipMUAData,0,1)*100;
        
        % extract LFP from spectrograms associated with the stimuli indecies
        stimCortZhold = [];
        stimHipZhold = [];
        for q = 1:length(finalStimFileIDs)
            % load normalized one-second bin data from each file
            stimFileID = finalStimFileIDs{q,1};
            stimSpecDataFileID = [animalID '_' stimFileID '_SpecData.mat'];
            stimSpecField = neuralDataType;
            for r = 1:length(AllSpecData.(stimSpecField).fileIDs)
                if strcmp(AllSpecData.(stimSpecField).fileIDs{r,1},stimSpecDataFileID) == true
                    stimCortS_Data = AllSpecData.(stimSpecField).oneSec.normS{r,1};
                    stimHipS_Data = AllSpecData.hippocampus.oneSec.normS{r,1};
                end
            end
            stimSLength = size(stimCortS_Data,2);
            stimBinSize = ceil(stimSLength/trialDuration_sec);
            stimSamplingDiff = samplingRate/stimBinSize;
            
            % find the start time and duration
            stimDuration = floor(floor(size(meanStimHbTData,2))/samplingRate);
            stimStartTime = floor(floor(finalStimFileEventTimes(q,1)*samplingRate)/stimSamplingDiff);
            if stimStartTime == 0
                stimStartTime = 1;
            end
            
            % take the S_data from the start time throughout the duration
            try
                stimCortS_Vals = stimCortS_Data(:,(stimStartTime - (offset*stimBinSize)):(stimStartTime + ((stimDuration - offset)*stimBinSize)));
                stimHipS_Vals = stimHipS_Data(:,(stimStartTime - (offset*stimBinSize)):(stimStartTime + ((stimDuration - offset)*stimBinSize)));
            catch
                stimCortS_Vals = stimCortS_Data(:,end - (stimDuration*stimBinSize):end);
                stimHipS_Vals = stimHipS_Data(:,end - (stimDuration*stimBinSize):end);
            end
            
            % mean subtract each row with detrend
            transpStimCortS_Vals = stimCortS_Vals';   % Transpose since detrend goes down columns
            transpStimHipS_Vals = stimHipS_Vals';   % Transpose since detrend goes down columns
            dTStimCortS_Vals = detrend(transpStimCortS_Vals,'constant');
            dTStimHipS_Vals = detrend(transpStimHipS_Vals,'constant');
            stimCortZhold = cat(3,stimCortZhold,dTStimCortS_Vals');   % transpose back to original orientation
            stimHipZhold = cat(3,stimHipZhold,dTStimHipS_Vals');   % transpose back to original orientation
        end
        
        % figure time/frequency axis and average each S data matrix through time
        meanStimCortS = mean(stimCortZhold,3);
        meanStimHipS = mean(stimHipZhold,3);
        
        % summary figure
        stimEvoked = figure;
        sgtitle([animalID ' ' dataType ' ' solenoid ' stimulus-evoked averages'])
        subplot(2,3,1);
        plot(timeVector,meanStimCortMUAData,'k')
        hold on
        plot(timeVector,meanStimCortMUAData + stdStimCortMUAData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanStimCortMUAData - stdStimCortMUAData,'color',colors_IOS('battleship grey'))
        title('Cortical MUA')
        xlabel('Time (sec)')
        ylabel('Fold-change (Norm Power)')
        axis tight
        axis square
        
        subplot(2,3,2);
        imagesc(T,F,meanStimCortS)
        title('Cortical MUA')
        xlabel('Time (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        caxis([-0.5 1])
        set(gca,'Ticklength',[0 0])
        axis xy
        axis square
        
        subplot(2,3,4);
        plot(timeVector,meanStimHipMUAData,'k')
        hold on
        plot(timeVector,meanStimHipMUAData + stdStimHipMUAData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanStimHipMUAData - stdStimHipMUAData,'color',colors_IOS('battleship grey'))
        title('Hippocampal MUA')
        xlabel('Time (sec)')
        ylabel('Fold-change (Norm Power)')
        axis tight
        axis square
        
        subplot(2,3,5);
        imagesc(T,F,meanStimHipS)
        title('Hippocampal MUA')
        xlabel('Time (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        caxis([-0.5 1])
        set(gca,'Ticklength',[0 0])
        axis xy
        axis square
        
        subplot(2,3,[3,6]);
        plot(timeVector,meanStimHbTData,'k')
        hold on
        plot(timeVector,meanStimHbTData + stdStimHbTData,'color',colors_IOS('battleship grey'))
        plot(timeVector,meanStimHbTData - stdStimHbTData,'color',colors_IOS('battleship grey'))
        title('Hemodynamics')
        xlabel('Time (sec)')
        ylabel('\DeltaHbT (\muM)')
        axis tight
        axis square
        
        % save results
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).CBV_HbT.HbT = meanStimHbTData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).CBV_HbT.HbTStD = stdStimHbTData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).CBV.CBV = meanStimCBVData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).CBV.CBVStD = stdStimCBVData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).MUA.corticalData = meanStimCortMUAData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).MUA.corticalStD = stdStimCortMUAData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).MUA.hippocampalData = meanStimHipMUAData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).MUA.hippocampalStD = stdStimHipMUAData;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).timeVector = timeVector;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).LFP.corticalS = meanStimCortS;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).LFP.hippocampalS = meanStimHipS;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).LFP.T = T;
        AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoid).LFP.F = F;
        
        % save figure
        [pathstr, ~, ~] = fileparts(cd);
        dirpath = [pathstr '/Combined Imaging/Figures/Stim and Whisk Evoked Averages/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(stimEvoked,[dirpath animalID '_' dataType '_' solenoid '_StimEvokedAverages']);
        close(stimEvoked)
    end
end

% save results structure
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

end
