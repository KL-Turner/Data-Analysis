function [AnalysisResults] = AnalyzeXCorr_IOS(CBVdataTypes, neuralDataTypes, baselineType, params, AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the cross-correlation between a CBV signal and a spectrogram during different behaviors
%________________________________________________________________________________________________________________________
%
%   Inputs: CBV dataType (char) denotes which hemisphere to look at
%           neuralDataType (char) denotes which neural electrode to look at
%           baselineType (char) denotes which baseline parameter to use for normalization
%           params (struct) - how many minutes of each day to use
%           AnalysisResults (struct) - where to save the data for later between-animal averages and comparisons
%
%   Outputs: AnalysisResults (struct) - where to save the data for later between-animal averages and comparisons
%
%   Last Revised: September 25th, 2019
%________________________________________________________________________________________________________________________

% list of all ProcData.mat files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID)

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

% find and load AllSpecStruct.mat struct
allSpecStructFileStruct = dir('*_AllSpecStruct.mat');
allSpecStructFile = {allSpecStructFileStruct.name}';
allSpecStructFileID = char(allSpecStructFile);
load(allSpecStructFileID)

% determine the animal's ID use the RestData.mat file's name for the current folder
fileBreaks = strfind(restDataFileID, '_');
animalID = restDataFileID(1:fileBreaks(1)-1);

% pull a few necessary numbers from the RestData.mat struct such as trial duration and sampling rate
samplingRate = RestData.CBV.LH.CBVCamSamplingRate;
trialDuration_sec = RestData.CBV.LH.trialDuration_sec;   % sec
trialDuration_min = trialDuration_sec/60;   % min
sleepBinWidth = 5;   % sec
fiveSecSpecFs = 5;   % sec
frequencyDiff = 6;   % Hz
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
fileTarget = params.targetMinutes/trialDuration_min;
filterSets = {'manualSelection','setDuration','entireDuration'};

for fgl = 1:length(CBVdataTypes)
    CBVdataType = CBVdataTypes{fgl};
    neuralDataType = neuralDataTypes{fgl};
    %% Cross-correlation analysis for resting data
    for a = 1:length(filterSets)
        filterSet = filterSets{1,a};
        disp(['AnalyzeXCorr: ' CBVdataType ' vs ' neuralDataType ' during Rest - ' filterSet]); disp(' ')
        % set criteria for rest event filter
        RestCriteria.Fieldname = {'durations'};
        RestCriteria.Comparison = {'gt'};
        RestCriteria.Value = {params.minTime.Rest};
        
        PuffCriteria.Fieldname = {'puffDistances'};
        PuffCriteria.Comparison = {'gt'};
        PuffCriteria.Value = {5};
        
        % filter the RestData structure for events that meet the desired criteria
        [restLogical] = FilterEvents_IOS(RestData.CBV.(CBVdataType), RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.CBV.(CBVdataType), PuffCriteria);
        restCombLogical = logical(restLogical.*puffLogical);
        allRestFiles = RestData.CBV.(CBVdataType).fileIDs(restCombLogical, :); 
        allRestDurations = RestData.CBV.(CBVdataType).durations(restCombLogical, :);
        allRestEventTimes = RestData.CBV.(CBVdataType).eventTimes(restCombLogical, :);
        allRestingCBVData = RestData.CBV.(CBVdataType).NormData(restCombLogical, :);
        allRestingHbTData = RestData.CBV_HbT.(CBVdataType).data(restCombLogical, :);
        allRestingMUAData = RestData.(neuralDataType).muaPower.NormData(restCombLogical, :);

        % identify the unique days and the unique number of files from the list of all resting events
        restUniqueDays = GetUniqueDays_IOS(RestData.CBV.(CBVdataType).fileIDs);
        restUniqueFiles = unique(RestData.CBV.(CBVdataType).fileIDs);
        restNumberOfFiles = length(unique(RestData.CBV.(CBVdataType).fileIDs));
        
        % decimate the file list to only include those files that occur within the desired number of target minutes
        clear restFiltLogical
        for b = 1:length(restUniqueDays)
            restDay = restUniqueDays(b);
            c = 1;
            for d = 1:restNumberOfFiles
                restFile = restUniqueFiles(d);
                restFileID = restFile{1}(1:6);
                if strcmp(filterSet,'manualSelection') == true
                    if strcmp(restDay,restFileID) && sum(strcmp(restFile,manualFileIDs)) == 1
                        restFiltLogical{b,1}(d,1) = 1; %#ok<*AGROW>
                        c = c + 1;
                    else
                        restFiltLogical{b,1}(d,1) = 0;
                    end
                elseif strcmp(filterSet,'setDuration') == true
                    if strcmp(restDay,restFileID) && c <= fileTarget
                        restFiltLogical{b,1}(d,1) = 1;
                        c = c + 1;
                    else
                        restFiltLogical{b,1}(d,1) = 0;
                    end
                elseif strcmp(filterSet,'entireDuration') == true
                    if strcmp(restDay,restFileID)
                        restFiltLogical{b,1}(d,1) = 1;
                        c = c + 1;
                    else
                        restFiltLogical{b,1}(d,1) = 0;
                    end
                end
            end
        end
        restFinalLogical = any(sum(cell2mat(restFiltLogical'),2),2);
        
        % extract all the resting events that correspond to the acceptable file list and the acceptable resting criteria
        clear restFileFilter
        filtRestFiles = restUniqueFiles(restFinalLogical, :);
        for d = 1:length(allRestFiles)
            restLogic = strcmp(allRestFiles{d}, filtRestFiles);
            restLogicSum = sum(restLogic);
            if restLogicSum == 1
                restFileFilter(d, 1) = 1;
            else
                restFileFilter(d, 1) = 0;
            end
        end
        restFinalFileFilter = logical(restFileFilter);
        restFinalFileIDs = allRestFiles(restFinalFileFilter, :);
        restFinalFileDurations = allRestDurations(restFinalFileFilter, :);
        restFinalFileEventTimes = allRestEventTimes(restFinalFileFilter, :);
        restFinalRestCBVData = allRestingCBVData(restFinalFileFilter, :);
        restFinalRestHbTData = allRestingHbTData(restFinalFileFilter, :);
        restFinalRestMUAData = allRestingMUAData(restFinalFileFilter, :);

        for e = 1:length(restFinalFileIDs)
            restFileID = restFinalFileIDs{e, 1};
            %% Load in CBV from rest period
            restCBV = restFinalRestCBVData{e,1};
            restHbT = restFinalRestHbTData{e,1};
            restMUA = restFinalRestMUAData{e,1};
            
            % low pass filter the epoch below 1 Hz
            [B, A] = butter(4, 1/(30/2), 'low');
            restFiltCBV = filtfilt(B, A, restCBV);
            restFiltHbT = filtfilt(B, A, restHbT);
            restFiltMUA = filtfilt(B, A, restMUA);
            
            % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
            % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
            if length(restFiltCBV) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(restFiltCBV);
                restPadCBV = (ones(1,restChunkSampleDiff))*restFiltCBV(end);
                restPadHbT = (ones(1,restChunkSampleDiff))*restFiltHbT(end);
                restPadMUA = (ones(1,restChunkSampleDiff))*restFiltMUA(end);
                restShortCBV = horzcat(restFiltCBVarray,restPadCBV);
                restShortHbT = horzcat(restFiltHbTarray,restPadHbT);
                restShortMUA = horzcat(restFiltMUAarray,restPadMUA);
            else
                restShortCBV = restFiltCBV(1:params.minTime.Rest*samplingRate);
                restShortHbT = restFiltHbT(1:params.minTime.Rest*samplingRate);
                restShortMUA = restFiltMUA(1:params.minTime.Rest*samplingRate);
            end
            
            % downsample the 10 second epoch to 5 Hz
            restDsCBV = downsample(restShortCBV, frequencyDiff);
            restDsHbT = downsample(restShortHbT, frequencyDiff);
            restDsMUA = downsample(restShortMUA, frequencyDiff);

            % mean subtract the downsampled epoch
            restProcData.CBV{e, 1} = detrend(restDsCBV,'constant');
            restProcData.HbT{e, 1} = detrend(restDsHbT,'constant');
            restProcData.MUA{e, 1} = detrend(restDsMUA,'constant');
           
            % extract LFP from spectrograms associated with the whisking indecies
            specDataFileID = [animalID '_' restFileID '_SpecData.mat'];
            clear S_data
            for g = 1:length(AllSpecData.(neuralDataType).fileIDs)
                if strcmp(AllSpecData.(neuralDataType).fileIDs{g,1},specDataFileID) == true
                    rest_S_data = AllSpecData.(neuralDataType).fiveSec.normS{g,1};
                    rest_F = AllSpecData.(neuralDataType).fiveSec.F{g,1};
                end
            end
            restSLength = size(rest_S_data, 2);
            restBinSize = ceil(restSLength/trialDuration_sec);
            restSegSamplingDiff = samplingRate/restBinSize;
            
            % Find the start time and duration
            restDuration = floor(floor(restFinalFileDurations(e, 1)*samplingRate)/restSegSamplingDiff);
            startTime = floor(floor(restFinalFileEventTimes(e, 1)*samplingRate)/restSegSamplingDiff);
            if startTime == 0
                startTime = 1;
            end
            
            % take the S data from the start time throughout the duration
            try
                restS_Vals = rest_S_data(:, (startTime:(startTime + restDuration)));
            catch
                restS_Vals = rest_S_data(:, end - restDuration:end);
            end
            
            % only take the first min rest time seconds
            shortRestS_Vals = restS_Vals(:,1:params.minTime.Rest*fiveSecSpecFs);
            
            % mean subtract with detrend and lowpass filter each column
            restProcData.S{e, 1} = detrend(shortRestS_Vals','constant')';
        end
        
        % set parameters for cross-correlation analysis
        restCBVvLFPzhold = [];
        restHbTvLFPzhold = [];
        restLagTime = 5;   % seconds
        restFrequency = fiveSecSpecFs;   % Hz
        restMaxLag = restLagTime*restFrequency;
        restCBVvLFPxcVals = ones(length(rest_F),2*restMaxLag + 1);
        restHbTvLFPxcVals = ones(length(rest_F),2*restMaxLag + 1);

        % run cross-correlation analysis - average through time
        for e = 1:length(restProcData.CBV)
            for b = 1:size(restProcData.S{e, 1}, 1)
                restCBVarray = restProcData.CBV{e,1};
                restHbTarray = restProcData.HbT{e,1};
                restMUAarray = restProcData.MUA{e,1};
                restNeuralArray = restProcData.S{e,1}(b,:);
                [restCBVvLFPxcVals(b,:),restLFP_lags] = xcorr(restCBVarray, restNeuralArray, restMaxLag, 'coeff');
                [restHbTvLFPxcVals(b,:),~] = xcorr(restHbTarray, restNeuralArray, restMaxLag, 'coeff');
            end
            [restCBVvMUAxcVals(e,:),restMUA_lags] = xcorr(restCBVarray, restMUAarray, restMaxLag, 'coeff');
            [restHbTvMUAxcVals(e,:),~] = xcorr(restHbTarray, restMUAarray, restMaxLag, 'coeff');
            restCBVvLFPzhold = cat(3,restCBVvLFPzhold,restCBVvLFPxcVals);
            restHbTvLFPzhold = cat(3,restHbTvLFPzhold,restHbTvLFPxcVals);
        end
        restMeanCBVvLFPxcVals = mean(restCBVvLFPzhold,3);
        restMeanHbTvLFPxcVals = mean(restHbTvLFPzhold,3);
        restMeanCBVvMUAxcVals = mean(restCBVvMUAxcVals,1);
        restStdCBVvMUAxcVals = std(restCBVvMUAxcVals,0,1);
        restMeanHbTvMUAxcVals = mean(restHbTvMUAxcVals,1);
        restStdHbTvMUAxcVals = std(restHbTvMUAxcVals,0,1);

        % summary figure
        titleID = strrep(CBVdataType, '_', ' ');
        RestingXCorr = figure;
        subplot(2,2,1)
        plot(restMUA_lags,restMeanCBVvMUAxcVals,'k')
        hold on
        plot(restMUA_lags,restMeanCBVvMUAxcVals + restStdCBVvMUAxcVals,'color',colors_IOS('battleship grey'))
        plot(restMUA_lags,restMeanCBVvMUAxcVals - restStdCBVvMUAxcVals,'color',colors_IOS('battleship grey'))
        title([animalID ' ' titleID ' ' filterSet ' CBV Refl resting cross-correlation'])
        xticks([-restMaxLag -restMaxLag/2 0 restMaxLag/2 restMaxLag])
        xticklabels({'-5', '-2.5', '0', '2.5' '5'})
        xlim([-restLagTime*restFrequency restLagTime*restFrequency])
        xlabel('Lags (sec)')
        ylabel('Cross-correlation')
        axis xy
        axis square
        
        subplot(2,2,2)
        plot(restMUA_lags,restMeanHbTvMUAxcVals,'k')
        hold on
        plot(restMUA_lags,restMeanHbTvMUAxcVals + restStdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
        plot(restMUA_lags,restMeanHbTvMUAxcVals - restStdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
        title([animalID ' ' titleID ' ' filterSet ' CBV HbT resting cross-correlation'])
        xticks([-restMaxLag -restMaxLag/2 0 restMaxLag/2 restMaxLag])
        xticklabels({'-5', '-2.5', '0', '2.5' '5'})
        xlim([-restLagTime*restFrequency restLagTime*restFrequency])
        xlabel('Lags (sec)')
        ylabel('Cross-correlation')
        axis xy
        axis square
        
        subplot(2,2,3)
        imagesc(restLFP_lags,rest_F,restMeanCBVvLFPxcVals)
        xticks([-restMaxLag -restMaxLag/2 0 restMaxLag/2 restMaxLag])
        xticklabels({'-5', '-2.5', '0', '2.5' '5'})
        xlim([-restLagTime*restFrequency restLagTime*restFrequency])
        xlabel('Lags (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        colorbar
        axis xy
        axis square
        
        subplot(2,2,4)
        imagesc(restLFP_lags,rest_F,restMeanHbTvLFPxcVals)
        xticks([-restMaxLag -restMaxLag/2 0 restMaxLag/2 restMaxLag])
        xticklabels({'-5', '-2.5', '0', '2.5' '5'})
        xlim([-restLagTime*restFrequency restLagTime*restFrequency])
        xlabel('Lags (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        colorbar
        axis xy
        axis square
        
        % save results
        AnalysisResults.XCorr.Rest.(CBVdataType).(filterSet).LFP_lags = restLFP_lags;
        AnalysisResults.XCorr.Rest.(CBVdataType).(filterSet).MUA_lags = restMUA_lags;
        AnalysisResults.XCorr.Rest.(CBVdataType).(filterSet).F = rest_F;
        AnalysisResults.XCorr.Rest.(CBVdataType).(filterSet).CBVvLFPxcVals = restMeanCBVvLFPxcVals;
        AnalysisResults.XCorr.Rest.(CBVdataType).(filterSet).HbTvLFPxcVals = restMeanHbTvLFPxcVals;
        AnalysisResults.XCorr.Rest.(CBVdataType).(filterSet).CBVvMUAxcVals = restMeanCBVvMUAxcVals;
        AnalysisResults.XCorr.Rest.(CBVdataType).(filterSet).CBVvMUAxcVals_std = restStdCBVvMUAxcVals;
        AnalysisResults.XCorr.Rest.(CBVdataType).(filterSet).HbTvMUAxcVals = restMeanHbTvMUAxcVals;
        AnalysisResults.XCorr.Rest.(CBVdataType).(filterSet).HbTvMUAxcVals_std = restStdHbTvMUAxcVals;

        % save figure
        [pathstr, ~, ~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Analysis XCorr/'];
        if ~exist(dirpath, 'dir')
            mkdir(dirpath);
        end
        savefig(RestingXCorr, [dirpath animalID '_' CBVdataType '_' filterSet  '_RestingXCorr']);
    end
    
    %% Cross-correlation analysis for NREM sleep data
    disp(['AnalyzeXCorr: ' CBVdataType ' vs ' neuralDataType ' during NREM.']); disp(' ')
    NREM_sleepTime = params.minTime.NREM;   % seconds
    NREM_allSleepFileIDs = SleepData.NREM.FileIDs;
    NREM_uniqueSleepFileIDs = unique(SleepData.NREM.FileIDs);
    NREM_sleepBins = NREM_sleepTime/sleepBinWidth;
    k = 1;
    for m = 1:length(NREM_uniqueSleepFileIDs)
        % pull out the bin times (there may be multiple events) in each unique NREM sleep file
        NREM_uniqueSleepFileID = char(NREM_uniqueSleepFileIDs(m));
        n = 1;
        clear NREM_binTimes
        for o = 1:length(NREM_allSleepFileIDs)
            NREM_sleepFileID = char(NREM_allSleepFileIDs(o));
            if strcmp(NREM_uniqueSleepFileID,NREM_sleepFileID)
                NREM_binTimes{n,1} = SleepData.NREM.BinTimes{o,1};
                n = n + 1;
            end
        end
        
        % pull out the Spectrogram data that matches the unique NREM sleep file
        NREM_specDataFileID = [animalID '_' NREM_uniqueSleepFileID '_SpecData.mat'];
        load(NREM_specDataFileID)
        NREM_S_Data = SpecData.(neuralDataType).fiveSec.normS;
        for p = 1:length(NREM_binTimes)
            NREM_Bins = NREM_binTimes{p,1};
            NREM_x_Length = size(NREM_S_Data,2);
            NREM_x_binLength = ceil(NREM_x_Length/trialDuration_sec);
            try
                NREM_sleepNeural_Vals{k,1} = NREM_S_Data(:,(NREM_Bins(1) - sleepBinWidth)*NREM_x_binLength + 1:(NREM_Bins(NREM_sleepBins))*NREM_x_binLength);
            catch
                NREM_sleepNeural_Vals{k,1} = NREM_S_Data(:,end - NREM_sleepTime*NREM_x_binLength + 1:end);
            end
            k = k + 1;
        end
    end
    
    % detrend spectrogram neural values
    for q = 1:length(NREM_sleepNeural_Vals)
        NREM_ind_sleepNeural_Vals = (NREM_sleepNeural_Vals{q,1})';
        NREM_dT_sleepNeural_Vals{q,1} = (detrend(NREM_ind_sleepNeural_Vals,'constant'))';
    end
    
    % lowpass filter and detrend reflectance and HbT during corresponding sleep events
    for r = 1:length(SleepData.NREM.data.CBV.(CBVdataType))
        NREM_CBV_Vals = SleepData.NREM.data.CBV.(CBVdataType){r,1}(1:(NREM_sleepTime*samplingRate));
        NREM_HbT_Vals = SleepData.NREM.data.CBV_HbT.(CBVdataType){r,1}(1:(NREM_sleepTime*samplingRate));
        NREM_MUA_Vals = SleepData.NREM.data.(neuralDataType).muaPower{r,1}(1:(NREM_sleepTime*samplingRate));
        NREM_dsCBV_Vals = downsample(NREM_CBV_Vals,frequencyDiff);
        NREM_dsHbT_Vals = downsample(NREM_HbT_Vals,frequencyDiff);
        NREM_dsMUA_Vals = downsample(NREM_MUA_Vals,frequencyDiff);
        NREM_dT_sleepCBV_Vals{r,1} = detrend(NREM_dsCBV_Vals,'constant');
        NREM_dT_sleepHbT_Vals{r,1} = detrend(NREM_dsHbT_Vals,'constant');
        NREM_dT_sleepMUA_Vals{r,1} = detrend(NREM_dsMUA_Vals,'constant');
    end
    
    % run cross-correlation analysis - average through time
    NREM_F = SpecData.(neuralDataType).fiveSec.F;
    NREM_CBVvLFPzHold = [];
    NREM_HbTvLFPzHold = [];
    NREM_lagTime = 15;   % Seconds
    NREM_frequency = fiveSecSpecFs;   % Hz
    NREM_maxLag = NREM_lagTime*NREM_frequency;
    NREM_CBVvLFPxcVals = ones(size(NREM_ind_sleepNeural_Vals,2),2*NREM_maxLag + 1);
    NREM_HbTvLFPxcVals = ones(size(NREM_ind_sleepNeural_Vals,2),2*NREM_maxLag + 1);
    for s = 1:length(NREM_dT_sleepNeural_Vals)
        for t = 1:size(NREM_dT_sleepNeural_Vals{s,1},1)
            NREM_CBV_array = NREM_dT_sleepCBV_Vals{s,1};
            NREM_HbT_array = NREM_dT_sleepHbT_Vals{s,1};
            NREM_MUA_array = NREM_dT_sleepMUA_Vals{s,1};
            NREM_Neural_array = NREM_dT_sleepNeural_Vals{s,1}(t,:);
            [NREM_CBVvLFPxcVals(t,:),NREM_LFP_lags] = xcorr(NREM_CBV_array,NREM_Neural_array,NREM_maxLag,'coeff');
            [NREM_HbTvLFPxcVals(t,:),~] = xcorr(NREM_HbT_array,NREM_Neural_array,NREM_maxLag,'coeff');
        end
        [NREM_CBVvMUAxcVals(s,:),NREM_MUA_lags] = xcorr(NREM_CBV_array,NREM_MUA_array,NREM_maxLag,'coeff');
        [NREM_HbTvMUAxcVals(s,:),~] = xcorr(NREM_HbT_array,NREM_MUA_array,NREM_maxLag,'coeff');
        NREM_CBVvLFPzHold = cat(3,NREM_CBVvLFPzHold,NREM_CBVvLFPxcVals);
        NREM_HbTvLFPzHold = cat(3,NREM_HbTvLFPzHold,NREM_HbTvLFPxcVals);
    end
    NREM_meanCBVvLFPxcVals = mean(NREM_CBVvLFPzHold,3);
    NREM_meanHbTvLFPxcVals = mean(NREM_HbTvLFPzHold,3);
    NREM_meanCBVvMUAxcVals = mean(NREM_CBVvMUAxcVals,1);
    NREM_stdCBVvMUAxcVals = std(NREM_CBVvMUAxcVals,0,1);
    NREM_meanHbTvMUAxcVals = mean(NREM_HbTvMUAxcVals,1);
    NREM_stdHbTvMUAxcVals = std(NREM_HbTvMUAxcVals,0,1);
    
    % summary figure
    NREMXCorr = figure;
    subplot(2,2,1)
    plot(NREM_MUA_lags,NREM_meanCBVvMUAxcVals,'k')
    hold on
    plot(NREM_MUA_lags,NREM_meanCBVvMUAxcVals + NREM_stdCBVvMUAxcVals,'color',colors_IOS('battleship grey'))
    plot(NREM_MUA_lags,NREM_meanCBVvMUAxcVals - NREM_stdCBVvMUAxcVals,'color',colors_IOS('battleship grey'))
    title([animalID ' ' titleID ' ' filterSet ' CBV Refl NREM cross-correlation'])
    xticks([-NREM_maxLag -NREM_maxLag/2 0 NREM_maxLag/2 NREM_maxLag])
    xticklabels({'-15', '-7.5', '0', '7.5' '15'})
    xlim([-NREM_lagTime*NREM_frequency NREM_lagTime*NREM_frequency])
    xlabel('Lags (sec)')
    ylabel('Cross-correlation')
    axis xy
    axis square
    
    subplot(2,2,2)
    plot(NREM_MUA_lags,NREM_meanHbTvMUAxcVals,'k')
    hold on
    plot(NREM_MUA_lags,NREM_meanHbTvMUAxcVals + NREM_stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
    plot(NREM_MUA_lags,NREM_meanHbTvMUAxcVals - NREM_stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
    title([animalID ' ' titleID ' ' filterSet ' CBV HbT NREM cross-correlation'])
    xticks([-NREM_maxLag -NREM_maxLag/2 0 NREM_maxLag/2 NREM_maxLag])
    xticklabels({'-15', '-7.5', '0', '7.5' '15'})
    xlim([-NREM_lagTime*NREM_frequency NREM_lagTime*NREM_frequency])
    xlabel('Lags (sec)')
    ylabel('Cross-correlation')
    axis xy
    axis square
    
    subplot(2,2,3)
    imagesc(NREM_LFP_lags,NREM_F,NREM_meanCBVvLFPxcVals)
    xticks([-NREM_maxLag -NREM_maxLag/2 0 NREM_maxLag/2 NREM_maxLag])
    xticklabels({'-15', '-7.5', '0', '7.5' '15'})
    xlim([-NREM_lagTime*NREM_frequency NREM_lagTime*NREM_frequency])
    xlabel('Lags (sec)')
    ylabel('Freq (Hz)')
    ylim([1 100])
    colorbar
    axis xy
    axis square
    
    subplot(2,2,4)
    imagesc(NREM_LFP_lags,NREM_F,NREM_meanHbTvLFPxcVals)
    xticks([-NREM_maxLag -NREM_maxLag/2 0 NREM_maxLag/2 NREM_maxLag])
    xticklabels({'-15', '-7.5', '0', '7.5' '15'})
    xlim([-NREM_lagTime*NREM_frequency NREM_lagTime*NREM_frequency])
    xlabel('Lags (sec)')
    ylabel('Freq (Hz)')
    ylim([1 100])
    colorbar
    axis xy
    axis square
    
    % save results
    AnalysisResults.XCorr.NREM.(CBVdataType).LFP_lags = NREM_LFP_lags;
    AnalysisResults.XCorr.NREM.(CBVdataType).MUA_lags = NREM_MUA_lags;
    AnalysisResults.XCorr.NREM.(CBVdataType).F = NREM_F;
    AnalysisResults.XCorr.NREM.(CBVdataType).CBVvLFPxcVals = NREM_meanCBVvLFPxcVals;
    AnalysisResults.XCorr.NREM.(CBVdataType).HbTvLFPxcVals = NREM_meanHbTvLFPxcVals;
    AnalysisResults.XCorr.NREM.(CBVdataType).CBVvMUAxcVals = NREM_meanCBVvMUAxcVals;
    AnalysisResults.XCorr.NREM.(CBVdataType).CBVvMUAxcVals_std = NREM_stdCBVvMUAxcVals;
    AnalysisResults.XCorr.NREM.(CBVdataType).HbTvMUAxcVals = NREM_meanHbTvMUAxcVals;
    AnalysisResults.XCorr.NREM.(CBVdataType).HbTvMUAxcVals_std = NREM_stdHbTvMUAxcVals;
    
    % save figure
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/XCorr/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(NREMXCorr, [dirpath animalID '_' CBVdataType '_NREMXCorr']);
    
    %% Cross-correlation analysis for REM sleep data
    disp(['AnalyzeXCorr: ' CBVdataType ' vs ' neuralDataType ' during REM.']); disp(' ')
    REM_sleepTime = params.minTime.REM;   % seconds
    REM_allSleepFileIDs = SleepData.REM.FileIDs;
    REM_uniqueSleepFileIDs = unique(SleepData.REM.FileIDs);
    REM_sleepBins = REM_sleepTime/sleepBinWidth;
    u = 1;
    for v = 1:length(REM_uniqueSleepFileIDs)
        % pull out the bin times (there may be multiple events) in each unique REM sleep file
        REM_uniqueSleepFileID = char(REM_uniqueSleepFileIDs(v));
        w = 1;
        clear REM_binTimes
        for x = 1:length(REM_allSleepFileIDs)
            REM_sleepFileID = char(REM_allSleepFileIDs(x));
            if strcmp(REM_uniqueSleepFileID,REM_sleepFileID)
                REM_binTimes{w,1} = SleepData.REM.BinTimes{x,1};
                w = w + 1;
            end
        end
        
        % pull out the Spectrogram data that matches the unique REM sleep file
        REM_specDataFileID = [animalID '_' REM_uniqueSleepFileID '_SpecData.mat'];
        load(REM_specDataFileID)
        REM_S_Data = SpecData.(neuralDataType).fiveSec.normS;
        for y = 1:length(REM_binTimes)
            REM_Bins = REM_binTimes{y,1};
            REM_x_Length = size(REM_S_Data,2);
            REM_x_binLength = ceil(REM_x_Length/trialDuration_sec);
            try
                REM_sleepNeural_Vals{u,1} = REM_S_Data(:,(REM_Bins(1) - sleepBinWidth)*REM_x_binLength + 1:(REM_Bins(REM_sleepBins))*REM_x_binLength);
            catch
                REM_sleepNeural_Vals{u,1} = REM_S_Data(:,end - REM_sleepTime*REM_x_binLength + 1:end);
            end
            u = u + 1;
        end
    end
    
    % detrend spectrogram neural values
    for z = 1:length(REM_sleepNeural_Vals)
        REM_ind_sleepNeural_Vals = (REM_sleepNeural_Vals{z,1})';
        REM_dT_sleepNeural_Vals{z,1} = (detrend(REM_ind_sleepNeural_Vals,'constant'))';
    end
    
    % lowpass filter and detrend reflectance and HbT during corresponding sleep events
    for r = 1:length(SleepData.REM.data.CBV.(CBVdataType))
        REM_CBV_Vals = SleepData.REM.data.CBV.(CBVdataType){r,1}(1:(REM_sleepTime*samplingRate));
        REM_HbT_Vals = SleepData.REM.data.CBV_HbT.(CBVdataType){r,1}(1:(REM_sleepTime*samplingRate));
        REM_MUA_Vals = SleepData.REM.data.(neuralDataType).muaPower{r,1}(1:(REM_sleepTime*samplingRate));
        REM_dsCBV_Vals = downsample(REM_CBV_Vals,frequencyDiff);
        REM_dsHbT_Vals = downsample(REM_HbT_Vals,frequencyDiff);
        REM_dsMUA_Vals = downsample(REM_MUA_Vals,frequencyDiff);
        REM_dT_sleepCBV_Vals{r,1} = detrend(REM_dsCBV_Vals,'constant');
        REM_dT_sleepHbT_Vals{r,1} = detrend(REM_dsHbT_Vals,'constant');
        REM_dT_sleepMUA_Vals{r,1} = detrend(REM_dsMUA_Vals,'constant');
    end
    
    % run cross-correlation analysis - average through time
    REM_F = SpecData.(neuralDataType).fiveSec.F;
    REM_CBVvLFPzHold = [];
    REM_HbTvLFPzHold = [];
    REM_lagTime = 15;   % Seconds
    REM_frequency = fiveSecSpecFs;   % Hz
    REM_maxLag = REM_lagTime*REM_frequency;
    REM_CBVvLFPxcVals = ones(size(REM_ind_sleepNeural_Vals,2),2*REM_maxLag + 1);
    REM_HbTvLFPxcVals = ones(size(REM_ind_sleepNeural_Vals,2),2*REM_maxLag + 1);
    for s = 1:length(REM_dT_sleepNeural_Vals)
        for t = 1:size(REM_dT_sleepNeural_Vals{s,1},1)
            REM_CBV_array = REM_dT_sleepCBV_Vals{s,1};
            REM_HbT_array = REM_dT_sleepHbT_Vals{s,1};
            REM_MUA_array = REM_dT_sleepMUA_Vals{s,1};
            REM_Neural_array = REM_dT_sleepNeural_Vals{s,1}(t,:);
            [REM_CBVvLFPxcVals(t,:),REM_LFP_lags] = xcorr(REM_CBV_array,REM_Neural_array,REM_maxLag,'coeff');
            [REM_HbTvLFPxcVals(t,:),~] = xcorr(REM_HbT_array,REM_Neural_array,REM_maxLag,'coeff');
        end
        [REM_CBVvMUAxcVals(s,:),REM_MUA_lags] = xcorr(REM_CBV_array,REM_MUA_array,REM_maxLag,'coeff');
        [REM_HbTvMUAxcVals(s,:),~] = xcorr(REM_HbT_array,REM_MUA_array,REM_maxLag,'coeff');
        REM_CBVvLFPzHold = cat(3,REM_CBVvLFPzHold,REM_CBVvLFPxcVals);
        REM_HbTvLFPzHold = cat(3,REM_HbTvLFPzHold,REM_HbTvLFPxcVals);
    end
    REM_meanCBVvLFPxcVals = mean(REM_CBVvLFPzHold,3);
    REM_meanHbTvLFPxcVals = mean(REM_HbTvLFPzHold,3);
    REM_meanCBVvMUAxcVals = mean(REM_CBVvMUAxcVals,1);
    REM_stdCBVvMUAxcVals = std(REM_CBVvMUAxcVals,0,1);
    REM_meanHbTvMUAxcVals = mean(REM_HbTvMUAxcVals,1);
    REM_stdHbTvMUAxcVals = std(REM_HbTvMUAxcVals,0,1);
    
    % summary figure
    REMXCorr = figure;
    subplot(2,2,1)
    plot(REM_MUA_lags,REM_meanCBVvMUAxcVals,'k')
    hold on
    plot(REM_MUA_lags,REM_meanCBVvMUAxcVals + REM_stdCBVvMUAxcVals,'color',colors_IOS('battleship grey'))
    plot(REM_MUA_lags,REM_meanCBVvMUAxcVals - REM_stdCBVvMUAxcVals,'color',colors_IOS('battleship grey'))
    title([animalID ' ' titleID ' ' filterSet ' CBV Refl REM cross-correlation'])
    xticks([-REM_maxLag -REM_maxLag/2 0 REM_maxLag/2 REM_maxLag])
    xticklabels({'-15', '-7.5', '0', '7.5' '15'})
    xlim([-REM_lagTime*REM_frequency REM_lagTime*REM_frequency])
    xlabel('Lags (sec)')
    ylabel('Cross-correlation')
    axis xy
    axis square
    
    subplot(2,2,2)
    plot(REM_MUA_lags,REM_meanHbTvMUAxcVals,'k')
    hold on
    plot(REM_MUA_lags,REM_meanHbTvMUAxcVals + REM_stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
    plot(REM_MUA_lags,REM_meanHbTvMUAxcVals - REM_stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
    title([animalID ' ' titleID ' ' filterSet ' CBV HbT REM cross-correlation'])
    xticks([-REM_maxLag -REM_maxLag/2 0 REM_maxLag/2 REM_maxLag])
    xticklabels({'-15', '-7.5', '0', '7.5' '15'})
    xlim([-REM_lagTime*REM_frequency REM_lagTime*REM_frequency])
    xlabel('Lags (sec)')
    ylabel('Cross-correlation')
    axis xy
    axis square
    
    subplot(2,2,3)
    imagesc(REM_LFP_lags,REM_F,REM_meanCBVvLFPxcVals)
    xticks([-REM_maxLag -REM_maxLag/2 0 REM_maxLag/2 REM_maxLag])
    xticklabels({'-15', '-7.5', '0', '7.5' '15'})
    xlim([-REM_lagTime*REM_frequency REM_lagTime*REM_frequency])
    xlabel('Lags (sec)')
    ylabel('Freq (Hz)')
    ylim([1 100])
    colorbar
    axis xy
    axis square
    
    subplot(2,2,4)
    imagesc(REM_LFP_lags,REM_F,REM_meanHbTvLFPxcVals)
    xticks([-REM_maxLag -REM_maxLag/2 0 REM_maxLag/2 REM_maxLag])
    xticklabels({'-15', '-7.5', '0', '7.5' '15'})
    xlim([-REM_lagTime*REM_frequency REM_lagTime*REM_frequency])
    xlabel('Lags (sec)')
    ylabel('Freq (Hz)')
    ylim([1 100])
    colorbar
    axis xy
    axis square
    
    % save results
    AnalysisResults.XCorr.REM.(CBVdataType).LFP_lags = REM_LFP_lags;
    AnalysisResults.XCorr.REM.(CBVdataType).MUA_lags = REM_MUA_lags;
    AnalysisResults.XCorr.REM.(CBVdataType).F = REM_F;
    AnalysisResults.XCorr.REM.(CBVdataType).CBVvLFPxcVals = REM_meanCBVvLFPxcVals;
    AnalysisResults.XCorr.REM.(CBVdataType).HbTvLFPxcVals = REM_meanHbTvLFPxcVals;
    AnalysisResults.XCorr.REM.(CBVdataType).CBVvMUAxcVals = REM_meanCBVvMUAxcVals;
    AnalysisResults.XCorr.REM.(CBVdataType).CBVvMUAxcVals_std = REM_stdCBVvMUAxcVals;
    AnalysisResults.XCorr.REM.(CBVdataType).HbTvMUAxcVals = REM_meanHbTvMUAxcVals;
    AnalysisResults.XCorr.REM.(CBVdataType).HbTvMUAxcVals_std = REM_stdHbTvMUAxcVals;
    
    % save figure
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/XCorr/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(REMXCorr, [dirpath animalID '_' CBVdataType '_REMXCorr']);
    
    %% Cross-correlation analysis for all data - no behavioral characterization
    disp(['AnalyzeXCorr: ' CBVdataType ' vs ' neuralDataType ' during all behaviors.']); disp(' ')
    for dd = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(dd, :);
        load(procDataFileID);
        [~,AD_fileDate,~] = GetFileInfo_IOS(procDataFileID);
        AD_Date = ConvertDate_IOS(AD_fileDate);
        
        % extract LFP from spectrograms associated with the whisking indecies
        AD_specDataFileID = [procDataFileID(1:end-12) 'SpecData.mat'];
        for ee = 1:length(AllSpecData.(neuralDataType).fileIDs)
            if strcmp(AllSpecData.(neuralDataType).fileIDs{ee,1},AD_specDataFileID) == true
                AD_S_Data = AllSpecData.(neuralDataType).fiveSec.normS{ee,1};
                AD_F = AllSpecData.(neuralDataType).fiveSec.F{ee,1};
            end
        end
        
        % mean subtract each row with detrend
        transp_AD_S_Vals = AD_S_Data';   % Transpose since detrend goes down columns
        dT_AD_S_Data{dd,1} = detrend(transp_AD_S_Vals,'constant')';
        
        % lowpass filter and detrend reflectance and HbT during corresponding sleep events
        AD_CBV = (ProcData.data.CBV.(CBVdataType) - RestingBaselines.(baselineType).CBV.(CBVdataType).(AD_Date))/RestingBaselines.(baselineType).CBV.(CBVdataType).(AD_Date);
        AD_MUA = (ProcData.data.(neuralDataType).muaPower - RestingBaselines.(baselineType).(neuralDataType).muaPower.(AD_Date))/RestingBaselines.(baselineType).(neuralDataType).muaPower.(AD_Date);
        AD_HbT = ProcData.data.CBV_HbT.(CBVdataType);
        filt_AD_CBV = filtfilt(B,A,AD_CBV);
        filt_AD_HbT = filtfilt(B,A,AD_HbT);
        filt_AD_MUA = filtfilt(B,A,AD_MUA);
        dS_AD_CBV = downsample(filt_AD_CBV,frequencyDiff);
        dS_AD_HbT = downsample(filt_AD_HbT,frequencyDiff);
        dS_AD_MUA = downsample(filt_AD_MUA,frequencyDiff);
        AD_SampleDiff = length(dS_AD_CBV) - size(AD_S_Data,2);
        short_AD_CBV = dS_AD_CBV((AD_SampleDiff/2):(end-(AD_SampleDiff/2))-1);
        short_AD_HbT = dS_AD_HbT((AD_SampleDiff/2):(end-(AD_SampleDiff/2))-1);
        short_AD_MUA = dS_AD_MUA((AD_SampleDiff/2):(end-(AD_SampleDiff/2))-1);
        dT_AD_CBV{dd,1} = detrend(short_AD_CBV,'constant');
        dT_AD_HbT{dd,1} = detrend(short_AD_HbT,'constant');
        dT_AD_MUA{dd,1} = detrend(short_AD_MUA,'constant');
    end
    
    % run cross-correlation analysis - average through time
    AD_CBVvLFPzHold = [];
    AD_HbTvLFPzHold = [];
    AD_lagTime = 5;   % seconds
    AD_frequency = fiveSecSpecFs;   % Hz
    AD_maxLag = AD_lagTime*AD_frequency;
    AD_CBVvLFPxcVals = ones(size(AD_S_Data,1),2*AD_maxLag + 1);
    AD_HbTvLFPxcVals = ones(size(AD_S_Data,1),2*AD_maxLag + 1);
    for ff = 1:length(dT_AD_S_Data)
        for gg = 1:size(dT_AD_S_Data{ff,1},1)
            AD_CBVarray = dT_AD_CBV{ff,1};
            AD_HbTarray = dT_AD_HbT{ff,1};
            AD_MUAarray = dT_AD_MUA{ff,1};
            AD_neuralArray = dT_AD_S_Data{ff,1}(gg,:);
            [AD_CBVvLFPxcVals(gg,:),AD_LFP_lags] = xcorr(AD_CBVarray,AD_neuralArray,AD_maxLag,'coeff');
            [AD_HbTvLFPxcVals(gg,:),~] = xcorr(AD_HbTarray,AD_neuralArray,AD_maxLag,'coeff');
        end
        [AD_CBVvMUAxcVals(ff,:),AD_MUA_lags] = xcorr(AD_CBVarray,AD_MUAarray,AD_maxLag,'coeff');
        [AD_HbTvMUAxcVals(ff,:),~] = xcorr(AD_HbTarray,AD_MUAarray,AD_maxLag,'coeff');
        AD_CBVvLFPzHold = cat(3,AD_CBVvLFPzHold,AD_CBVvLFPxcVals);
        AD_HbTvLFPzHold = cat(3,AD_HbTvLFPzHold,AD_HbTvLFPxcVals);
    end
    AD_meanCBVvLFPxcVals = mean(AD_CBVvLFPzHold,3);
    AD_meanHbTvLFPxcVals = mean(AD_HbTvLFPzHold,3);
    AD_meanCBVvMUAxcVals = mean(AD_CBVvMUAxcVals,1);
    AD_stdCBVvMUAxcVals = std(AD_CBVvMUAxcVals,0,1);
    AD_meanHbTvMUAxcVals = mean(AD_HbTvMUAxcVals,1);
    AD_stdHbTvMUAxcVals = std(AD_HbTvMUAxcVals,0,1);
    
    % summary figure
    AllDataXCorr = figure;
    subplot(2,2,1)
    plot(AD_MUA_lags,AD_meanCBVvMUAxcVals,'k')
    hold on
    plot(AD_MUA_lags,AD_meanCBVvMUAxcVals + AD_stdCBVvMUAxcVals,'color',colors_IOS('battleship grey'))
    plot(AD_MUA_lags,AD_meanCBVvMUAxcVals - AD_stdCBVvMUAxcVals,'color',colors_IOS('battleship grey'))
    title([animalID ' ' titleID ' ' filterSet ' CBV Refl all data cross-correlation'])
    xticks([-AD_maxLag -AD_maxLag/2 0 AD_maxLag/2 AD_maxLag])
    xticklabels({'-5', '-2.5', '0', '2.5' '5'})
    xlim([-AD_lagTime*AD_frequency AD_lagTime*AD_frequency])
    xlabel('Lags (sec)')
    ylabel('Cross-correlation')
    axis xy
    axis square
    
    subplot(2,2,2)
    plot(AD_MUA_lags,AD_meanHbTvMUAxcVals,'k')
    hold on
    plot(AD_MUA_lags,AD_meanHbTvMUAxcVals + AD_stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
    plot(AD_MUA_lags,AD_meanHbTvMUAxcVals - AD_stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
    title([animalID ' ' titleID ' ' filterSet ' CBV HbT all data cross-correlation'])
    xticks([-AD_maxLag -AD_maxLag/2 0 AD_maxLag/2 AD_maxLag])
    xticklabels({'-5', '-2.5', '0', '2.5' '5'})
    xlim([-AD_lagTime*AD_frequency AD_lagTime*AD_frequency])
    xlabel('Lags (sec)')
    ylabel('Cross-correlation')
    axis xy
    axis square
    
    subplot(2,2,3)
    imagesc(AD_LFP_lags,AD_F,AD_meanCBVvLFPxcVals)
    xticks([-AD_maxLag -AD_maxLag/2 0 AD_maxLag/2 AD_maxLag])
    xticklabels({'-5', '-2.5', '0', '2.5' '5'})
    xlim([-AD_lagTime*AD_frequency AD_lagTime*AD_frequency])
    xlabel('Lags (sec)')
    ylabel('Freq (Hz)')
    ylim([1 100])
    colorbar
    axis xy
    axis square
    
    subplot(2,2,4)
    imagesc(AD_LFP_lags,AD_F,AD_meanHbTvLFPxcVals)
    xticks([-AD_maxLag -AD_maxLag/2 0 AD_maxLag/2 AD_maxLag])
    xticklabels({'-5', '-2.5', '0', '2.5' '5'})
    xlim([-AD_lagTime*AD_frequency AD_lagTime*AD_frequency])
    xlabel('Lags (sec)')
    ylabel('Freq (Hz)')
    ylim([1 100])
    colorbar
    axis xy
    axis square
    
    % save results
    AnalysisResults.XCorr.AllData.(CBVdataType).LFP_lags = AD_LFP_lags;
    AnalysisResults.XCorr.AllData.(CBVdataType).MUA_lags = AD_MUA_lags;
    AnalysisResults.XCorr.AllData.(CBVdataType).F = AD_F;
    AnalysisResults.XCorr.AllData.(CBVdataType).CBVvLFPxcVals = AD_meanCBVvLFPxcVals;
    AnalysisResults.XCorr.AllData.(CBVdataType).HbTvLFPxcVals = AD_meanHbTvLFPxcVals;
    AnalysisResults.XCorr.AllData.(CBVdataType).CBVvMUAxcVals = AD_meanCBVvMUAxcVals;
    AnalysisResults.XCorr.AllData.(CBVdataType).CBVvMUAxcVals_std = AD_stdCBVvMUAxcVals;
    AnalysisResults.XCorr.AllData.(CBVdataType).HbTvMUAxcVals = AD_meanHbTvMUAxcVals;
    AnalysisResults.XCorr.AllData.(CBVdataType).HbTvMUAxcVals_std = AD_stdHbTvMUAxcVals;
    
    % save figure
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Analysis XCorr/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(AllDataXCorr, [dirpath animalID '_' CBVdataType '_AllDataXCorr']);
end

%% save final results structure
save([animalID '_AnalysisResults'], 'AnalysisResults')

end
