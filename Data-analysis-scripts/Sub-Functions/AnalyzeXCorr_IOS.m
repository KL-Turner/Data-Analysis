function [AnalysisResults] = AnalyzeXCorr_IOS(dataTypes,params,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
%   Purpose: Analyze the cross-correlation between a CBV_HbT signal and a spectrogram during different behaviors.
%________________________________________________________________________________________________________________________

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
fileBreaks = strfind(restDataFileID,'_');
animalID = restDataFileID(1:fileBreaks(1)-1);

for z = 1:length(dataTypes)
    dataType = dataTypes{1,z};
    neuralDataType = ['cortical_' dataType(4:end)];
    % pull a few necessary numbers from the RestData.mat struct such as trial duration and sampling rate
    samplingRate = RestData.CBV_HbT.LH.CBVCamSamplingRate;
    trialDuration_sec = RestData.CBV_HbT.LH.trialDuration_sec;   % sec
    sleepBinWidth = 5;   % sec
    oneSecSpecFs = 10;   % sec   5 for fiveSec, 10 for oneSec
    frequencyDiff = 3;   % Hz    6 for fiveSec, 3 for oneSec
    manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
    
    %% Cross-correlation analysis for resting data
    disp(['AnalyzeXCorr: ' dataType ' vs ' neuralDataType ' during Rest']); disp(' ')
    % set criteria for rest event filter
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    
    PuffCriteria.Fieldname = {'puffDistances'};
    PuffCriteria.Comparison = {'gt'};
    PuffCriteria.Value = {5};
    
    % filter the RestData structure for events that meet the desired criteria
    [restLogical] = FilterEvents_IOS(RestData.CBV_HbT.(dataType),RestCriteria);
    [puffLogical] = FilterEvents_IOS(RestData.CBV_HbT.(dataType),PuffCriteria);
    restCombLogical = logical(restLogical.*puffLogical);
    allRestFiles = RestData.CBV_HbT.(dataType).fileIDs(restCombLogical,:);
    allRestDurations = RestData.CBV_HbT.(dataType).durations(restCombLogical,:);
    allRestEventTimes = RestData.CBV_HbT.(dataType).eventTimes(restCombLogical,:);
    allRestingHbTData = RestData.CBV_HbT.(dataType).data(restCombLogical,:);
    allRestingMUAData = RestData.(neuralDataType).muaPower.NormData(restCombLogical,:);
    
    % identify the unique days and the unique number of files from the list of all resting events
    restUniqueDays = GetUniqueDays_IOS(RestData.CBV_HbT.(dataType).fileIDs);
    restUniqueFiles = unique(RestData.CBV_HbT.(dataType).fileIDs);
    restNumberOfFiles = length(unique(RestData.CBV_HbT.(dataType).fileIDs));
    
    % decimate the file list to only include those files that occur within the desired number of target minutes
    clear restFiltLogical
    for b = 1:length(restUniqueDays)
        restDay = restUniqueDays(b);
        c = 1;
        for d = 1:restNumberOfFiles
            restFile = restUniqueFiles(d);
            restFileID = restFile{1}(1:6);
            if strcmp(restDay,restFileID) && sum(strcmp(restFile,manualFileIDs)) == 1
                restFiltLogical{b,1}(d,1) = 1; %#ok<*AGROW>
                c = c + 1;
            else
                restFiltLogical{b,1}(d,1) = 0;
            end
        end
    end
    restFinalLogical = any(sum(cell2mat(restFiltLogical'),2),2);
    
    % extract all the resting events that correspond to the acceptable file list and the acceptable resting criteria
    clear restFileFilter
    filtRestFiles = restUniqueFiles(restFinalLogical,:);
    for d = 1:length(allRestFiles)
        restLogic = strcmp(allRestFiles{d},filtRestFiles);
        restLogicSum = sum(restLogic);
        if restLogicSum == 1
            restFileFilter(d, 1) = 1;
        else
            restFileFilter(d, 1) = 0;
        end
    end
    restFinalFileFilter = logical(restFileFilter);
    restFinalFileIDs = allRestFiles(restFinalFileFilter,:);
    restFinalFileDurations = allRestDurations(restFinalFileFilter,:);
    restFinalFileEventTimes = allRestEventTimes(restFinalFileFilter,:);
    restFinalRestHbTData = allRestingHbTData(restFinalFileFilter,:);
    restFinalRestMUAData = allRestingMUAData(restFinalFileFilter,:);
    
    for e = 1:length(restFinalFileIDs)
        restFileID = restFinalFileIDs{e, 1};
        %% Load in CBV_HbT from rest period
        restHbT = restFinalRestHbTData{e,1};
        restMUA = restFinalRestMUAData{e,1};
        
        % low pass filter the epoch below 1 Hz
        [B, A] = butter(3,1/(samplingRate/2),'low');
        restFiltHbT = filtfilt(B,A,restHbT);
        restFiltMUA = filtfilt(B,A,restMUA);
        
        % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
        % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
        if length(restFiltHbT) < params.minTime.Rest*samplingRate
            restChunkSampleDiff = params.minTime.Rest*samplingRate - length(restFiltHbT);
            restPadHbT = (ones(1,restChunkSampleDiff))*restFiltHbT(end);
            restPadMUA = (ones(1,restChunkSampleDiff))*restFiltMUA(end);
            restShortHbT = horzcat(restFiltHbT,restPadHbT);
            restShortMUA = horzcat(restFiltMUA,restPadMUA);
        else
            restShortHbT = restFiltHbT(1:params.minTime.Rest*samplingRate);
            restShortMUA = restFiltMUA(1:params.minTime.Rest*samplingRate);
        end
        
        % downsample the 10 second epoch to 5 Hz
        restDsHbT = downsample(restShortHbT,frequencyDiff);
        restDsMUA = downsample(restShortMUA,frequencyDiff);
        
        % mean subtract the downsampled epoch
        restProcData.HbT{e,1} = detrend(restDsHbT,'constant');
        restProcData.MUA{e,1} = detrend(restDsMUA,'constant');
        
        % extract LFP from spectrograms associated with the whisking indecies
        specDataFileID = [animalID '_' restFileID '_SpecData.mat'];
        clear S_data
        for g = 1:length(AllSpecData.(neuralDataType).fileIDs)
            if strcmp(AllSpecData.(neuralDataType).fileIDs{g,1},specDataFileID) == true
                rest_S_data = AllSpecData.(neuralDataType).oneSec.normS{g,1};
                rest_F = AllSpecData.(neuralDataType).oneSec.F{g,1};
            end
        end
        restSLength = size(rest_S_data,2);
        restBinSize = ceil(restSLength/trialDuration_sec);
        restSegSamplingDiff = samplingRate/restBinSize;
        
        % Find the start time and duration
        restDuration = floor(floor(restFinalFileDurations(e,1)*samplingRate)/restSegSamplingDiff);
        startTime = floor(floor(restFinalFileEventTimes(e,1)*samplingRate)/restSegSamplingDiff);
        if startTime == 0
            startTime = 1;
        end
        
        % take the S data from the start time throughout the duration
        try
            restS_Vals = rest_S_data(:,(startTime:(startTime + restDuration)));
        catch
            restS_Vals = rest_S_data(:,end - restDuration:end);
        end
        
        % only take the first min rest time seconds
        shortRestS_Vals = restS_Vals(:,1:params.minTime.Rest*oneSecSpecFs);
        
        % mean subtract with detrend and lowpass filter each column
        restProcData.S{e,1} = detrend(shortRestS_Vals','constant')';
    end
    
    % set parameters for cross-correlation analysis
    restHbTvLFPzhold = [];
    restLagTime = 5;   % seconds
    restFrequency = oneSecSpecFs;   % Hz
    restMaxLag = restLagTime*restFrequency;
    restHbTvLFPxcVals = ones(length(rest_F),2*restMaxLag + 1);
    
    % run cross-correlation analysis - average through time
    for e = 1:length(restProcData.HbT)
        for b = 1:size(restProcData.S{e, 1}, 1)
            restHbTarray = restProcData.HbT{e,1};
            restMUAarray = restProcData.MUA{e,1};
            restNeuralArray = restProcData.S{e,1}(b,:);
            [restHbTvLFPxcVals(b,:),restLFP_lags] = xcorr(restHbTarray, restNeuralArray, restMaxLag, 'coeff');
        end
        [restHbTvMUAxcVals(e,:),restMUA_lags] = xcorr(restHbTarray, restMUAarray, restMaxLag, 'coeff');
        restHbTvLFPzhold = cat(3,restHbTvLFPzhold,restHbTvLFPxcVals);
    end
    restMeanHbTvLFPxcVals = mean(restHbTvLFPzhold,3);
    restMeanHbTvMUAxcVals = mean(restHbTvMUAxcVals,1);
    restStdHbTvMUAxcVals = std(restHbTvMUAxcVals,0,1);
    
    % summary figure
    titleID = strrep(dataType, '_', ' ');
    RestingXCorr = figure;   
    sgtitle([animalID ' ' titleID ' resting cross-correlation'])
    subplot(2,1,1)
    plot(restMUA_lags,restMeanHbTvMUAxcVals,'k')
    hold on
    plot(restMUA_lags,restMeanHbTvMUAxcVals + restStdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
    plot(restMUA_lags,restMeanHbTvMUAxcVals - restStdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
    title('MUA XCorr')
    xticks([-restMaxLag -restMaxLag/2 0 restMaxLag/2 restMaxLag])
    xticklabels({'-5','-2.5','0','2.5','5'})
    xlim([-restLagTime*restFrequency restLagTime*restFrequency])
    xlabel('Lags (sec)')
    ylabel('Cross-correlation')
    axis xy
    axis square
    
    subplot(2,1,2)
    imagesc(restLFP_lags,rest_F,restMeanHbTvLFPxcVals)
    title('LFP XCorr')
    xticks([-restMaxLag -restMaxLag/2 0 restMaxLag/2 restMaxLag])
    xticklabels({'-5','-2.5','0','2.5','5'})
    xlim([-restLagTime*restFrequency restLagTime*restFrequency])
    xlabel('Lags (sec)')
    ylabel('Freq (Hz)')
    ylim([1 100])
    colorbar
    axis xy
    axis square
    
    % save results
    AnalysisResults.XCorr.Rest.(dataType).LFP_lags = restLFP_lags;
    AnalysisResults.XCorr.Rest.(dataType).MUA_lags = restMUA_lags;
    AnalysisResults.XCorr.Rest.(dataType).F = rest_F;
    AnalysisResults.XCorr.Rest.(dataType).HbTvLFPxcVals = restMeanHbTvLFPxcVals;
    AnalysisResults.XCorr.Rest.(dataType).HbTvMUAxcVals = restMeanHbTvMUAxcVals;
    AnalysisResults.XCorr.Rest.(dataType).HbTvMUAxcVals_std = restStdHbTvMUAxcVals;
    
    % save figure
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Combined Imaging/Figures/Cross Correlation/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(RestingXCorr, [dirpath animalID '_' dataType '_RestingXCorr']);
    close(RestingXCorr)
    
    %% Cross-correlation analysis for NREM sleep data
    disp(['AnalyzeXCorr: ' dataType ' vs ' neuralDataType ' during NREM.']); disp(' ')
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
        for ee = 1:length(NREM_allSleepFileIDs)
            NREM_sleepFileID = char(NREM_allSleepFileIDs(ee));
            if strcmp(NREM_uniqueSleepFileID,NREM_sleepFileID)
                NREM_binTimes{n,1} = SleepData.NREM.BinTimes{ee,1};
                n = n + 1;
            end
        end
        
        % pull out the Spectrogram data that matches the unique NREM sleep file
        NREM_specDataFileID = [animalID '_' NREM_uniqueSleepFileID '_SpecData.mat'];
        load(NREM_specDataFileID)
        NREM_S_Data = SpecData.(neuralDataType).oneSec.normS;
        for q = 1:length(NREM_binTimes)
            NREM_Bins = NREM_binTimes{q,1};
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
    for r = 1:length(NREM_sleepNeural_Vals)
        NREM_ind_sleepNeural_Vals = (NREM_sleepNeural_Vals{r,1})';
        NREM_dT_sleepNeural_Vals{r,1} = (detrend(NREM_ind_sleepNeural_Vals,'constant'))';
    end
    
    % lowpass filter and detrend reflectance and HbT during corresponding sleep events
    for s = 1:length(SleepData.NREM.data.CBV_HbT.(dataType(4:end)))
        NREM_HbT_Vals = SleepData.NREM.data.CBV_HbT.(dataType(4:end)){s,1}(1:(NREM_sleepTime*samplingRate));
        NREM_MUA_Vals = SleepData.NREM.data.(neuralDataType).muaPower{s,1}(1:(NREM_sleepTime*samplingRate));
        NREM_dsHbT_Vals = downsample(NREM_HbT_Vals,frequencyDiff);
        NREM_dsMUA_Vals = downsample(NREM_MUA_Vals,frequencyDiff);
        NREM_dT_sleepHbT_Vals{s,1} = detrend(NREM_dsHbT_Vals,'constant');
        NREM_dT_sleepMUA_Vals{s,1} = detrend(NREM_dsMUA_Vals,'constant');
    end
    
    % run cross-correlation analysis - average through time
    NREM_F = SpecData.(neuralDataType).oneSec.F;
    NREM_HbTvLFPzHold = [];
    NREM_lagTime = 15;   % Seconds
    NREM_frequency = oneSecSpecFs;   % Hz
    NREM_maxLag = NREM_lagTime*NREM_frequency;
    NREM_HbTvLFPxcVals = ones(size(NREM_ind_sleepNeural_Vals,2),2*NREM_maxLag + 1);
    for t = 1:length(NREM_dT_sleepNeural_Vals)
        for u = 1:size(NREM_dT_sleepNeural_Vals{t,1},1)
            NREM_HbT_array = NREM_dT_sleepHbT_Vals{t,1};
            NREM_MUA_array = NREM_dT_sleepMUA_Vals{t,1};
            NREM_Neural_array = NREM_dT_sleepNeural_Vals{t,1}(u,:);
            [NREM_HbTvLFPxcVals(u,:),NREM_LFP_lags] = xcorr(NREM_HbT_array,NREM_Neural_array,NREM_maxLag,'coeff');
        end
        [NREM_HbTvMUAxcVals(t,:),NREM_MUA_lags] = xcorr(NREM_HbT_array,NREM_MUA_array,NREM_maxLag,'coeff');
        NREM_HbTvLFPzHold = cat(3,NREM_HbTvLFPzHold,NREM_HbTvLFPxcVals);
    end
    NREM_meanHbTvLFPxcVals = mean(NREM_HbTvLFPzHold,3);
    NREM_meanHbTvMUAxcVals = mean(NREM_HbTvMUAxcVals,1);
    NREM_stdHbTvMUAxcVals = std(NREM_HbTvMUAxcVals,0,1);
    
    % summary figure
    NREMXCorr = figure;
    subplot(2,1,1)
    sgtitle([animalID ' ' titleID ' NREM cross-correlation'])
    plot(NREM_MUA_lags,NREM_meanHbTvMUAxcVals,'k')
    hold on
    plot(NREM_MUA_lags,NREM_meanHbTvMUAxcVals + NREM_stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
    plot(NREM_MUA_lags,NREM_meanHbTvMUAxcVals - NREM_stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
    title('MUA XCorr')
    xticks([-NREM_maxLag -NREM_maxLag/2 0 NREM_maxLag/2 NREM_maxLag])
    xticklabels({'-15','-7.5','0','7.5','15'})
    xlim([-NREM_lagTime*NREM_frequency NREM_lagTime*NREM_frequency])
    xlabel('Lags (sec)')
    ylabel('Cross-correlation')
    axis xy
    axis square
    
    subplot(2,1,2)
    imagesc(NREM_LFP_lags,NREM_F,NREM_meanHbTvLFPxcVals)
    title('LFP XCorr')
    xticks([-NREM_maxLag -NREM_maxLag/2 0 NREM_maxLag/2 NREM_maxLag])
    xticklabels({'-15','-7.5','0','7.5','15'})
    xlim([-NREM_lagTime*NREM_frequency NREM_lagTime*NREM_frequency])
    xlabel('Lags (sec)')
    ylabel('Freq (Hz)')
    ylim([1 100])
    colorbar
    axis xy
    axis square
    
    % save results
    AnalysisResults.XCorr.NREM.(dataType).LFP_lags = NREM_LFP_lags;
    AnalysisResults.XCorr.NREM.(dataType).MUA_lags = NREM_MUA_lags;
    AnalysisResults.XCorr.NREM.(dataType).F = NREM_F;
    AnalysisResults.XCorr.NREM.(dataType).HbTvLFPxcVals = NREM_meanHbTvLFPxcVals;
    AnalysisResults.XCorr.NREM.(dataType).HbTvMUAxcVals = NREM_meanHbTvMUAxcVals;
    AnalysisResults.XCorr.NREM.(dataType).HbTvMUAxcVals_std = NREM_stdHbTvMUAxcVals;
    
    % save figure
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Combined Imaging/Figures/Cross Correlation/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(NREMXCorr, [dirpath animalID '_' dataType '_NREMXCorr']);
    close(NREMXCorr)
    
    %% Cross-correlation analysis for REM sleep data
    disp(['AnalyzeXCorr: ' dataType ' vs ' neuralDataType ' during REM.']); disp(' ')
    REM_sleepTime = params.minTime.REM;   % seconds
    REM_allSleepFileIDs = SleepData.REM.FileIDs;
    REM_uniqueSleepFileIDs = unique(SleepData.REM.FileIDs);
    REM_sleepBins = REM_sleepTime/sleepBinWidth;
    v = 1;
    for w = 1:length(REM_uniqueSleepFileIDs)
        % pull out the bin times (there may be multiple events) in each unique REM sleep file
        REM_uniqueSleepFileID = char(REM_uniqueSleepFileIDs(w));
        x = 1;
        clear REM_binTimes
        for y = 1:length(REM_allSleepFileIDs)
            REM_sleepFileID = char(REM_allSleepFileIDs(y));
            if strcmp(REM_uniqueSleepFileID,REM_sleepFileID)
                REM_binTimes{x,1} = SleepData.REM.BinTimes{y,1};
                x = x + 1;
            end
        end
        
        % pull out the Spectrogram data that matches the unique REM sleep file
        REM_specDataFileID = [animalID '_' REM_uniqueSleepFileID '_SpecData.mat'];
        load(REM_specDataFileID)
        REM_S_Data = SpecData.(neuralDataType).oneSec.normS;
        for z = 1:length(REM_binTimes)
            REM_Bins = REM_binTimes{z,1};
            REM_x_Length = size(REM_S_Data,2);
            REM_x_binLength = ceil(REM_x_Length/trialDuration_sec);
            try
                REM_sleepNeural_Vals{v,1} = REM_S_Data(:,(REM_Bins(1) - sleepBinWidth)*REM_x_binLength + 1:(REM_Bins(REM_sleepBins))*REM_x_binLength);
            catch
                REM_sleepNeural_Vals{v,1} = REM_S_Data(:,end - REM_sleepTime*REM_x_binLength + 1:end);
            end
            v = v + 1;
        end
    end
    
    % detrend spectrogram neural values
    for aa = 1:length(REM_sleepNeural_Vals)
        REM_ind_sleepNeural_Vals = (REM_sleepNeural_Vals{aa,1})';
        REM_dT_sleepNeural_Vals{aa,1} = (detrend(REM_ind_sleepNeural_Vals,'constant'))';
    end
    
    % lowpass filter and detrend reflectance and HbT during corresponding sleep events
    for bb = 1:length(SleepData.REM.data.CBV_HbT.(dataType(4:end)))
        REM_HbT_Vals = SleepData.REM.data.CBV_HbT.(dataType(4:end)){bb,1}(1:(REM_sleepTime*samplingRate));
        REM_MUA_Vals = SleepData.REM.data.(neuralDataType).muaPower{bb,1}(1:(REM_sleepTime*samplingRate));
        REM_dsHbT_Vals = downsample(REM_HbT_Vals,frequencyDiff);
        REM_dsMUA_Vals = downsample(REM_MUA_Vals,frequencyDiff);
        REM_dT_sleepHbT_Vals{bb,1} = detrend(REM_dsHbT_Vals,'constant');
        REM_dT_sleepMUA_Vals{bb,1} = detrend(REM_dsMUA_Vals,'constant');
    end
    
    % run cross-correlation analysis - average through time
    REM_F = SpecData.(neuralDataType).oneSec.F;
    REM_HbTvLFPzHold = [];
    REM_lagTime = 15;   % Seconds
    REM_frequency = oneSecSpecFs;   % Hz
    REM_maxLag = REM_lagTime*REM_frequency;
    REM_HbTvLFPxcVals = ones(size(REM_ind_sleepNeural_Vals,2),2*REM_maxLag + 1);
    for cc = 1:length(REM_dT_sleepNeural_Vals)
        for dd = 1:size(REM_dT_sleepNeural_Vals{cc,1},1)
            REM_HbT_array = REM_dT_sleepHbT_Vals{cc,1};
            REM_MUA_array = REM_dT_sleepMUA_Vals{cc,1};
            REM_Neural_array = REM_dT_sleepNeural_Vals{cc,1}(dd,:);
            [REM_HbTvLFPxcVals(dd,:),REM_LFP_lags] = xcorr(REM_HbT_array,REM_Neural_array,REM_maxLag,'coeff');
        end
        [REM_HbTvMUAxcVals(cc,:),REM_MUA_lags] = xcorr(REM_HbT_array,REM_MUA_array,REM_maxLag,'coeff');
        REM_HbTvLFPzHold = cat(3,REM_HbTvLFPzHold,REM_HbTvLFPxcVals);
    end
    REM_meanHbTvLFPxcVals = mean(REM_HbTvLFPzHold,3);
    REM_meanHbTvMUAxcVals = mean(REM_HbTvMUAxcVals,1);
    REM_stdHbTvMUAxcVals = std(REM_HbTvMUAxcVals,0,1);
    
    % summary figure
    REMXCorr = figure;
    subplot(2,1,1)
    sgtitle([animalID ' ' titleID ' REM cross-correlation'])
    plot(REM_MUA_lags,REM_meanHbTvMUAxcVals,'k')
    hold on
    plot(REM_MUA_lags,REM_meanHbTvMUAxcVals + REM_stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
    plot(REM_MUA_lags,REM_meanHbTvMUAxcVals - REM_stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
    title('MUA XCorr')
    xticks([-REM_maxLag -REM_maxLag/2 0 REM_maxLag/2 REM_maxLag])
    xticklabels({'-15','-7.5','0','7.5','15'})
    xlim([-REM_lagTime*REM_frequency REM_lagTime*REM_frequency])
    xlabel('Lags (sec)')
    ylabel('Cross-correlation')
    axis xy
    axis square
    
    subplot(2,1,2)
    imagesc(REM_LFP_lags,REM_F,REM_meanHbTvLFPxcVals)
    title('LFP XCorr')
    xticks([-REM_maxLag -REM_maxLag/2 0 REM_maxLag/2 REM_maxLag])
    xticklabels({'-15','-7.5','0','7.5','15'})
    xlim([-REM_lagTime*REM_frequency REM_lagTime*REM_frequency])
    xlabel('Lags (sec)')
    ylabel('Freq (Hz)')
    ylim([1 100])
    colorbar
    axis xy
    axis square
    
    % save results
    AnalysisResults.XCorr.REM.(dataType).LFP_lags = REM_LFP_lags;
    AnalysisResults.XCorr.REM.(dataType).MUA_lags = REM_MUA_lags;
    AnalysisResults.XCorr.REM.(dataType).F = REM_F;
    AnalysisResults.XCorr.REM.(dataType).HbTvLFPxcVals = REM_meanHbTvLFPxcVals;
    AnalysisResults.XCorr.REM.(dataType).HbTvMUAxcVals = REM_meanHbTvMUAxcVals;
    AnalysisResults.XCorr.REM.(dataType).HbTvMUAxcVals_std = REM_stdHbTvMUAxcVals;
    
    % save figure
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Combined Imaging/Figures/Cross Correlation/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(REMXCorr, [dirpath animalID '_' dataType '_REMXCorr']);
    close(REMXCorr)
end

%% save final results structure
save([animalID '_AnalysisResults'],'AnalysisResults')

end
