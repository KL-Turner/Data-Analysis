function [AnalysisResults] = AnalyzeXCorr_IOS2(CBVdataTypes, neuralDataTypes, baselineType, params, AnalysisResults)
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
filterSets = {'setDuration', 'manualSelection', 'entireDuration'};

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
                    F = AllSpecData.(neuralDataType).fiveSec.F{g,1};
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
            shortRestS_Vals = restS_Vals(:, 1:params.minTime.Rest*fiveSecSpecFs);
            
            % mean subtract with detrend and lowpass filter each column
            transpRestS_Vals = shortRestS_Vals';        % Transpose since detrend goes down columns
            dtRestS_Vals = detrend(transpRestS_Vals);
            restProcData.S{e, 1} = dtRestS_Vals';              % transpose back to original orientation
        end
        
        % set parameters for cross-correlation analysis
        % F = SpecData.(neuralDataType).fiveSec.F;
        restCBVvLFPzhold = [];
        restHbTvLFPzhold = [];
        restCBVvMUAzhold = [];
        restHbTvMUAzhold = [];
        restLagTime = 5;   % seconds
        restFrequency = fiveSecSpecFs;   % Hz
        restMaxLag = restLagTime*restFrequency;
        restCBVvLFPxcVals = ones(length(F),2*restMaxLag + 1);
        restHbTvLFPxcVals = ones(length(F),2*restMaxLag + 1);
%         restCBVvMUAxcVals = ones(length(F),2*restMaxLag + 1);
%         restHbTvMUAxcVals = ones(length(F),2*restMaxLag + 1);

        % run cross-correlation analysis - average through time
        for e = 1:length(restProcData.CBV)
            for b = 1:size(restProcData.S{e, 1}, 1)
                CBV_array = restProcData.CBV{e, 1};
                HbT_array = restProcData.HbT{e, 1};
                Neural_array = restProcData.S{e, 1}(b, :);
                [restCBVvLFPxcVals(b,:),restLags] = xcorr(CBV_array, Neural_array, restMaxLag, 'coeff');
                [restHbTvLFPxcVals(b,:),~] = xcorr(HbT_array, Neural_array, restMaxLag, 'coeff');
            end
            restCBVvLFPzhold = cat(3,restCBVvLFPzhold,restCBVvLFPxcVals);
            restHbTvLFPzhold = cat(3,restHbTvLFPzhold,restHbTvLFPxcVals);
        end
        restMeanCBVvLFPxcVals = mean(restCBVvLFPzhold,3);
        restMeanHbTvLFPxcVals = mean(restHbTvLFPzhold,3);

        % summary figure
        titleID = strrep(CBVdataType, '_', ' ');
        RestingXCorr = figure;
        subplot(1,2,1)
        imagesc(restLags, F, restMeanCBVvLFPxcVals)
        title([animalID ' ' titleID ' ' filterSet ' Resting Cross Correlation'])
        xticks([-restMaxLag -restMaxLag/2 0 restMaxLag/2 restMaxLag])
        xticklabels({'-5', '-2.5', '0', '2.5' '5'})
        xlim([-restLagTime*restFrequency restLagTime*restFrequency])
        xlabel('Lags (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        colorbar
        axis xy
        axis square
        
           subplot(1,2,2)
        imagesc(restLags, F, restMeanHbTvLFPxcVals)
        title([animalID ' ' titleID ' ' filterSet ' Resting Cross Correlation'])
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
        AnalysisResults.XCorr.Rest.(CBVdataType).lags = lags;
        AnalysisResults.XCorr.Rest.(CBVdataType).F = F;
        AnalysisResults.XCorr.Rest.(CBVdataType).CBVvLFPxcVals = meanXC_Vals;
        AnalysisResults.XCorr.Rest.(CBVdataType).HbTvLFPxcVals = meanXC_Vals;
        AnalysisResults.XCorr.Rest.(CBVdataType).CBVvMUAxcVals = meanXC_Vals;
        AnalysisResults.XCorr.Rest.(CBVdataType).HbTvMUAxcVals = meanXC_Vals;

        % save figure
        [pathstr, ~, ~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Analysis XCorr/'];
        if ~exist(dirpath, 'dir')
            mkdir(dirpath);
        end
        savefig(RestingXCorr, [dirpath animalID '_' CBVdataType '_' filterSet  '_RestingXCorr']);
    end
end

%% Cross-correlation analysis for NREM sleep data
disp(['AnalyzeXCorr: ' CBVdataType ' vs ' neuralDataType ' during NREM.']); disp(' ')
NREM_sleepTime = params.minTime.NREM;   % seconds
NREM_allSleepFileIDs = SleepData.NREM.FileIDs;
NREM_uniqueSleepFileIDs = unique(SleepData.NREM.FileIDs);
NREM_sleepBins = NREM_sleepTime/sleepBinWidth;
x = 1;
for uID = 1:length(NREM_uniqueSleepFileIDs)
    
    % Pull out the bin times (there may be multiple events) in each unique NREM sleep file
    NREM_uniqueSleepFileID = char(NREM_uniqueSleepFileIDs(uID));
    y = 1;
    clear NREM_binTimes
    for iID = 1:length(NREM_allSleepFileIDs)
        NREM_sleepFileID = char(NREM_allSleepFileIDs(iID));
        if strcmp(NREM_uniqueSleepFileID, NREM_sleepFileID)
            NREM_binTimes{y, 1} = SleepData.NREM.BinTimes{iID, 1};
            y = y + 1;
        end
    end
    
    % Pull out the Spectrogram data that matches the unique NREM sleep file
    NREM_specDataFileID = [animalID '_' NREM_uniqueSleepFileID '_SpecData.mat'];
            load(NREM_specDataFileID)
            NREM_S_Data = SpecData.(neuralDataType).fiveSec.normS;

            
    for rBN = 1:length(NREM_binTimes)
        NREM_Bins = NREM_binTimes{rBN, 1};
        NREM_x_Length = size(NREM_S_Data, 2);
        NREM_x_binLength = ceil(NREM_x_Length/900);
        try
            NREM_sleepNeural_Vals{x, 1} = NREM_S_Data(:, (NREM_Bins(1) - sleepBinWidth)*NREM_x_binLength + 1:(NREM_Bins(NREM_sleepBins))*NREM_x_binLength);
        catch
            NREM_sleepNeural_Vals{x, 1} = NREM_S_Data(:, end - NREM_sleepTime*NREM_x_binLength + 1:end);
        end
        x = x + 1;
    end
end

for sNV = 1:length(NREM_sleepNeural_Vals)
    NREM_ind_sleepNeural_Vals = (NREM_sleepNeural_Vals{sNV, 1})';
    NREM_dT_sleepNeural_Vals{sNV, 1} = (detrend(NREM_ind_sleepNeural_Vals, 'constant'))';
end

% CBV
for restCBV = 1:length(SleepData.NREM.data.CBV.(CBVdataType))
    NREM_CBV_Vals = SleepData.NREM.data.CBV.(CBVdataType){restCBV, 1}(1:(NREM_sleepTime*samplingRate));
    NREM_dsCBV_Vals = downsample(NREM_CBV_Vals, 6);
    NREM_dT_sleepCBV_Vals{restCBV, 1} = detrend(NREM_dsCBV_Vals, 'constant');
end

% Cross Correlation
NREM_F = SpecData.(neuralDataType).fiveSec.F;
NREM_z_hold = [];
NREM_lagTime = 15;   % Seconds
NREM_frequency = fiveSecSpecFs;   % Hz
NREM_maxLag = NREM_lagTime*NREM_frequency;
NREM_XC_Vals = ones(size(NREM_ind_sleepNeural_Vals, 2), 2*NREM_maxLag + 1);

for x = 1:length(NREM_dT_sleepNeural_Vals)
    for y = 1:size(NREM_dT_sleepNeural_Vals{x, 1}, 1)
        NREM_CBV_array = NREM_dT_sleepCBV_Vals{x, 1};
        NREM_Neural_array = NREM_dT_sleepNeural_Vals{x, 1}(y, :);
        [NREM_XC_Vals(y, :), NREM_lags] = xcorr(NREM_CBV_array, NREM_Neural_array, NREM_maxLag, 'coeff');
    end
    NREM_z_hold = cat(3, NREM_z_hold, NREM_XC_Vals);
end
NREM_meanXC_Vals = mean(NREM_z_hold, 3);

NREMXCorr = figure;
imagesc(NREM_lags, NREM_F, NREM_meanXC_Vals)
title([animalID ' ' titleID ' NREM Cross Correlation'])
xticks([-NREM_maxLag -NREM_maxLag/2 0 NREM_maxLag/2 NREM_maxLag])
xticklabels({'-15', '-7.5', '0', '7.5' '15'})
xlabel('Lags (sec)')
ylabel('Freq (Hz)')
xlim([-NREM_lagTime*NREM_frequency NREM_lagTime*NREM_frequency])
ylim([1 100])
colorbar
axis xy

AnalysisResults.XCorr.NREM.(CBVdataType).lags = NREM_lags;
AnalysisResults.XCorr.NREM.(CBVdataType).F = NREM_F;
AnalysisResults.XCorr.NREM.(CBVdataType).XC_Vals = NREM_meanXC_Vals;
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

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
x = 1;
for uID = 1:length(REM_uniqueSleepFileIDs)
    
    % Pull out the bin times (there may be multiple events) in each unique NREM sleep file
    REM_uniqueSleepFileID = char(REM_uniqueSleepFileIDs(uID));
    y = 1;
    clear REM_binTimes
    for iID = 1:length(REM_allSleepFileIDs)
        REM_sleepFileID = char(REM_allSleepFileIDs(iID));
        if strcmp(REM_uniqueSleepFileID, REM_sleepFileID)
            REM_binTimes{y, 1} = SleepData.REM.BinTimes{iID, 1};
            y = y + 1;
        end
    end
    
    % Pull out the Spectrogram data that matches the unique NREM sleep file
    REM_specDataFileID = [animalID '_' REM_uniqueSleepFileID '_SpecData.mat'];
    load(REM_specDataFileID)
    REM_S_Data = SpecData.(neuralDataType).fiveSec.normS;
    
    for rBN = 1:length(REM_binTimes)
        REM_Bins = REM_binTimes{rBN, 1};
        REM_x_Length = size(REM_S_Data, 2);
        REM_x_binLength = ceil(REM_x_Length/900);
        try
            REM_sleepNeural_Vals{x, 1} = REM_S_Data(:, (REM_Bins(1) - sleepBinWidth)*REM_x_binLength + 1:(REM_Bins(REM_sleepBins))*REM_x_binLength);
        catch
            REM_sleepNeural_Vals{x, 1} = REM_S_Data(:, end - REM_sleepTime*REM_x_binLength + 1:end);
        end
        x = x + 1;
    end
end

for sNV = 1:length(REM_sleepNeural_Vals)
    REM_ind_sleepNeural_Vals = (REM_sleepNeural_Vals{sNV, 1})';
    REM_dT_sleepNeural_Vals{sNV, 1} = (detrend(REM_ind_sleepNeural_Vals, 'constant'))';
end

% CBV
for restCBV = 1:length(SleepData.REM.data.CBV.(CBVdataType))
    REM_CBV_Vals = SleepData.REM.data.CBV.(CBVdataType){restCBV, 1}(1:(REM_sleepTime*samplingRate));
    REM_dsCBV_Vals = downsample(REM_CBV_Vals, 6);
    REM_dT_sleepCBV_Vals{restCBV, 1} = detrend(REM_dsCBV_Vals, 'constant');
end

% Cross Correlation
REM_F = SpecData.(neuralDataType).fiveSec.F;
REM_z_hold = [];
REM_lagTime = 15;   % Seconds
REM_frequency = fiveSecSpecFs;   % Hz
REM_maxLag = REM_lagTime*REM_frequency;
REM_XC_Vals = ones(size(REM_ind_sleepNeural_Vals, 2), 2*REM_maxLag + 1);

for x = 1:length(REM_dT_sleepNeural_Vals)
    for y = 1:size(REM_dT_sleepNeural_Vals{x,1},1)
        REM_CBV_array = REM_dT_sleepCBV_Vals{x, 1};
        REM_Neural_array = REM_dT_sleepNeural_Vals{x, 1}(y, :);
        [REM_XC_Vals(y, :), REM_lags] = xcorr(REM_CBV_array, REM_Neural_array, REM_maxLag, 'coeff');
    end
    REM_z_hold = cat(3, REM_z_hold, REM_XC_Vals);
end
REM_meanXC_Vals = mean(REM_z_hold, 3);

REMXCorr = figure;
imagesc(REM_lags, REM_F, REM_meanXC_Vals)
title([animalID ' ' titleID ' REM Cross Correlation'])
xticks([-REM_maxLag -REM_maxLag/2 0 REM_maxLag/2 REM_maxLag])
xticklabels({'-15', '-7.5', '0', '7.5' '15'})
xlabel('Lags (sec)')
ylabel('Freq (Hz)')
xlim([-REM_lagTime*REM_frequency REM_lagTime*REM_frequency])
ylim([1 100])
colorbar
axis xy

AnalysisResults.XCorr.REM.(CBVdataType).lags = REM_lags;
AnalysisResults.XCorr.REM.(CBVdataType).F = REM_F;
AnalysisResults.XCorr.REM.(CBVdataType).XC_Vals = REM_meanXC_Vals;
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/XCorr/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(NREMXCorr, [dirpath animalID '_' CBVdataType '_REMXCorr']);

%% BLOCK PURPOSE: All data - no behavioral characterization
% disp(['AnalyzeXCorr: ' CBVdataType ' vs ' neuralDataType ' during all behaviors.']); disp(' ')
% for pDF = 1:size(procDataFileIDs,1)
%     procDataFileID = procDataFileIDs(pDF, :);
%     load(procDataFileID);
%     [~, fileDate, allData_FileID] = GetFileInfo_IOS(procDataFileID);
%     strDay = ConvertDate_IOS(fileDate);
%     
%     %% Neural Data associated with fileID
%     specDataFileID = [procDataFileID(1:end-12) 'SpecData.mat'];
%     load(specDataFileID)
%     allData_S_Data = SpecData.(neuralDataType).fiveSec.normS';
%     dT_allData_S_Data{pDF, 1} = (detrend(allData_S_Data, 'constant'))';
%     
%     %% CBV Data associated with fileID
%     allData_CBV = (ProcData.data.CBV.(CBVdataType) - RestingBaselines.(baselineType).CBV.(CBVdataType).(strDay))/RestingBaselines.(baselineType).CBV.(CBVdataType).(strDay);
%     [B, A] = butter(4, 1/(30/2), 'low');
%     filt_allData_CBV = filtfilt(B, A, allData_CBV);
%     dS_allData_CBV = downsample(filt_allData_CBV, 6);
%     dT_allData_CBV{pDF, 1} = detrend(dS_allData_CBV, 'constant');
% end
% 
% allData_F = SpecData.(neuralDataType).fiveSec.F;
% allData_z_hold = [];
% allData_lagTime = 5;       % Seconds
% allData_frequency = 5;     % Hz
% allData_maxLag = allData_lagTime*allData_frequency;    % Number of points
% allData_XC_Vals = ones(size(allData_S_Data, 2), 2*allData_maxLag + 1);   % Pre-allocate size of cross-corr matrix
% 
% for x = 1:length(dT_allData_S_Data)
%     for y = 1:size(dT_allData_S_Data{x,1},2)
%         allData_CBV_array = dT_allData_CBV{x,1};
%         allData_Neural_array = dT_allData_S_Data{x,1}(y,:);
%         [allData_XC_Vals(y,:), allData_lags] = xcorr(allData_CBV_array, allData_Neural_array, allData_maxLag, 'coeff');
%     end
%     allData_z_hold = cat(3, allData_z_hold, allData_XC_Vals);
% end
% allData_meanXC_Vals = mean(allData_z_hold, 3);
% 
% AllDataXCorr = figure;
% imagesc(allData_lags, allData_F, allData_meanXC_Vals)
% title([animalID ' ' titleID ' No Behavior (All Data) Cross Correlation, One Sec Bins'])
% xticks([-allData_maxLag -allData_maxLag/2 0 allData_maxLag/2 allData_maxLag])
% xticklabels({'-5', '-2.5', '0', '2.5' '5'})
% xlabel('Lags (sec)')
% ylabel('Freq (Hz)')
% xlim([-allData_lagTime*allData_frequency allData_lagTime*allData_frequency])
% ylim([1 100])
% colorbar
% axis xy
% 
% AnalysisResults.XCorr.AllData.(CBVdataType).lags = allData_lags;
% AnalysisResults.XCorr.AllData.(CBVdataType).F = allData_F;
% AnalysisResults.XCorr.AllData.(CBVdataType).XC_Vals = allData_meanXC_Vals;
% save([animalID '_AnalysisResults.mat'], 'AnalysisResults');
% 
% [pathstr, ~, ~] = fileparts(cd);
% dirpath = [pathstr '/Figures/Analysis XCorr/'];
% 
% if ~exist(dirpath, 'dir')
%     mkdir(dirpath);
% end
% 
% savefig(AllDataXCorr, [dirpath animalID '_' CBVdataType '_AllDataXCorr']);

%% save results
end