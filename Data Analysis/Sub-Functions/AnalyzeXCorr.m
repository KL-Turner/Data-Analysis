function [ComparisonData] = AnalyzeXCorr(animal, CBVdataType, neuralDataType, params, RestData, RestingBaselines, SpectrogramData, SleepData, ComparisonData)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: //
%________________________________________________________________________________________________________________________
%
%   Inputs: //
%
%   Outputs: //
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: RestData
% The RestData.mat struct has all resting events, regardless of duration. We want to set the threshold for rest as anything
% that is greater than 10 seconds.
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};

PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};

% Use the RestCriteria we specified earlier to find all resting events that are greater than the criteria
[restLogical] = FilterEvents(RestData.CBV.(CBVdataType), RestCriteria);   % Output is a logical
[puffLogical] = FilterEvents(RestData.CBV.(CBVdataType), PuffCriteria);   % Output is a logical
combRestLogical = logical(restLogical.*puffLogical);
allRestFiles = RestData.CBV.(CBVdataType).fileIDs(combRestLogical, :);   % Overall logical for all resting file names that meet criteria
allRestDurations = RestData.CBV.(CBVdataType).durations(combRestLogical, :);
allRestEventTimes = RestData.CBV.(CBVdataType).eventTimes(combRestLogical, :);
allRestingData = RestData.CBV.(CBVdataType).data(combRestLogical, :);   % Pull out data from all those resting files that meet criteria

uniqueDays = GetUniqueDays(RestData.CBV.(CBVdataType).fileIDs);   % Find the unique days of imaging
uniqueFiles = unique(RestData.CBV.(CBVdataType).fileIDs);   % Find the unique files from the filelist. This removes duplicates
% since most files have more than one resting event
numberOfFiles = length(unique(RestData.CBV.(CBVdataType).fileIDs));   % Find the number of unique files
fileTarget = params.targetMinutes / 5;   % Divide that number of unique files by 5 (minutes) to get the number of files that
% corresponds to the desired targetMinutes

% Loop through each unique day in order to create a logical to filter the file list so that it only includes the first
% x number of files that fall within the targetMinutes requirement
for uD = 1:length(uniqueDays)
    day = uniqueDays(uD);
    x = 1;
    for nOF = 1:numberOfFiles
        file = uniqueFiles(nOF);
        fileID = file{1}(1:6);
        if strcmp(params.Infusion, 'y')
            if strcmp(day, fileID) && x > fileTarget
                filtLogical{uD, 1}(nOF, 1) = 1;
            else
                filtLogical{uD, 1}(nOF, 1) = 0;
                x = x + 1;
            end
        else
            if strcmp(day, fileID) && x <= fileTarget
                filtLogical{uD, 1}(nOF, 1) = 1;
                x = x + 1;
            else
                filtLogical{uD, 1}(nOF, 1) = 0;
            end
        end
    end
end

% Combine the 3 logicals so that it reflects the first "x" number of files from each day
finalLogical = any(sum(cell2mat(filtLogical'), 2), 2);

% Now that the appropriate files from each day are identified, loop through each file name with respect to the original
% list of ALL resting files, only keeping the ones that fall within the first targetMinutes of each day.
filtRestFiles = uniqueFiles(finalLogical, :);
for rF = 1:length(allRestFiles)
    logic = strcmp(allRestFiles{rF}, filtRestFiles);
    logicSum = sum(logic);
    if logicSum == 1
        fileFilter(rF, 1) = 1;
    else
        fileFilter(rF, 1) = 0;
    end
end

finalFileFilter = logical(fileFilter);
finalFileIDs = allRestFiles(finalFileFilter, :);
finalFileDurations = allRestDurations(finalFileFilter, :);
finalFileEventTimes = allRestEventTimes(finalFileFilter, :);
finalRestData = allRestingData(finalFileFilter, :);

disp('Add lowpass to neural data'); disp(' ')
keyboard

for a = 1:length(finalFileIDs)   % Loop through each non-unique file
    fileID = finalFileIDs{a, 1};   % Specify the fileID from the list
    strDay = fileID(1:6);   % First 6 characters denote the date
    date = ConvertDate(strDay);   % Convert the numeric date string to a month/day
    
    %% Load in CBV from rest period
    CBVSamplingRate = 30;   % CBV sampling rate is typically always 30 Hz
    CBV = (finalRestData{a, 1} - RestingBaselines.CBV.(CBVdataType).(date)) /  RestingBaselines.CBV.(CBVdataType).(date);   % Loop through the list of rest CBV data - the index is consistent with the file number
    
    % Low pass filter the epoch below 2 Hz
    [B, A] = butter(4, 2 / (30 / 2), 'low');
    filtCBV = filtfilt(B, A, CBV);
    
    % Only take the first 15 seconds of the epoch
    shortCBV = filtCBV(1:10*CBVSamplingRate);
    
    % Downsample the 10 second epoch to 5 Hz
    dsCBV = downsample(shortCBV, 6);
    
    % Mean subtract the downsampled epoch
    Data.CBV{a, 1} = detrend(dsCBV, 'constant');
    
    
    %% Load in Neural Data from rest period
    for s = 1:length(SpectrogramData.(neuralDataType).FileIDs)
        if strcmp(fileID, SpectrogramData.(neuralDataType).FileIDs{s, 1})
            S_data = SpectrogramData.(neuralDataType).OneSec.S_Norm{s, 1};  % S data for this specific file
        end
    end
    sLength = size(S_data, 2);                                  % Length of the data across time (number of samples)
    binSize = ceil(sLength / 300);                              % Find the number of bins needed to divide this into 300 seconds
    samplingDiff = CBVSamplingRate / binSize;                   % Number to divide by for data to be at 5 Hz
    
    % Find the start time and duration
    restDuration = floor(floor(finalFileDurations(a, 1)*CBVSamplingRate) / samplingDiff);
    startTime = floor(floor(finalFileEventTimes(a, 1)*CBVSamplingRate) / samplingDiff);
    if startTime == 0
        startTime = 1;
    end
    
    % Take the S_data from the start time throughout the duration
    try
        restS_Vals = S_data(:, (startTime:(startTime + restDuration)));
    catch
        restS_Vals = S_data(:, end - restDuration:end);
    end
    
    % Only take the first 15 seconds
    shortRestS_Vals = restS_Vals(:, 1:50);
    
    % Mean subtract each row with detrend
    transpRestS_Vals = shortRestS_Vals';        % Transpose since detrend goes down columns
    dtRestS_Vals = detrend(transpRestS_Vals);   % detrend
    Data.S{a, 1} = dtRestS_Vals';              % transpose back to original orientation and save into Data.S
end

F = SpectrogramData.(neuralDataType).OneSec.F{1,1};
z_hold = [];
lagTime = 5;       % Seconds
frequency = 5;     % Hz
maxLag = lagTime*frequency;    % Number of points
XC_Vals = ones(length(F), 2*maxLag + 1);   % Pre-allocate size of cross-corr matrix

for a = 1:length(Data.CBV)
    for b = 1:size(Data.S{a, 1}, 1)
        CBV_array = Data.CBV{a, 1};
        Neural_array = Data.S{a, 1}(b, :);
        [XC_Vals(b, :), lags] = xcorr(CBV_array, Neural_array, maxLag, 'coeff');
    end
    z_hold = cat(3, z_hold, XC_Vals);
end

meanXC_Vals = mean(z_hold, 3);

RestingXCorr = figure;
imagesc(lags, F, meanXC_Vals)
title([animal ' ' neuralDataType ' Resting Cross Correlation'])
xticks([-maxLag -maxLag/2 0 maxLag/2 maxLag])
xticklabels({'-5', '-2.5', '0', '2.5' '5'})
xlim([-lagTime*frequency lagTime*frequency])
xlabel('Lags (sec)')
ylabel('Freq (Hz)')
ylim([1 100])
colorbar
axis xy

boostrapResults = [];
for a = 1:length(Data.CBV)
    disp(['Boostrapping CBV event ' num2str(a) ' of ' num2str(length(Data.CBV))]); disp(' ')
    for boot = 1:1000
        shuffledCBV_array_Inds = randperm(length(Data.CBV{a, 1}));
        shuffledCBV_array = Data.CBV{a,1}(shuffledCBV_array_Inds);
        shufXC_Vals = ones(length(F), 2*maxLag + 1);   % Pre-allocate size of cross-corr matrix
        for b = 1:size(Data.S{a, 1}, 1)
            Neural_array = Data.S{a, 1}(b, :);
            [shufXC_Vals(b, :), lags] = xcorr(shuffledCBV_array, Neural_array, maxLag, 'coeff');
        end
        boostrapResults = cat(3, boostrapResults, shufXC_Vals);
    end
end

test = prctile(boostrapResults, [2.5 97.5], 3);
test2 = test(:,:,1) - test(:,:,2);

ComparisonData.XCorr.Rest.(neuralDataType).lags = lags;
ComparisonData.XCorr.Rest.(neuralDataType).F = F;
ComparisonData.XCorr.Rest.(neuralDataType).XC_Vals = meanXC_Vals;
save([animal '_ComparisonData.mat'], 'ComparisonData');

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/XCorr/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(RestingXCorr, [dirpath animal '_' neuralDataType '_RestingXCorr']);

% %% BLOCK PURPOSE: All data - no behavioral characterization
% procDataFiles = ls('*_ProcData.mat');
% for pDF = 1:size(procDataFiles, 1)
%     procDataFile = procDataFiles(pDF, :);
%     load(procDataFile);
%     [~, ~, fileDate, allData_FileID] = GetFileInfo_IOS(procDataFile);
%     strDay = ConvertDate(fileDate);
%     
%     %% Neural Data associated with fileID
%     allData_spectrogramFileIDs = SpectrogramData.(neuralDataType).FileIDs;
%     for sID = 1:length(allData_spectrogramFileIDs)
%         allData_specFileID = char(allData_spectrogramFileIDs(sID));
%         if strcmp(allData_FileID, allData_specFileID)
%             allData_S_Data = (SpectrogramData.(neuralDataType).OneSec.S_Norm{sID, 1}(:, 1:1495))';
%         end
%     end
%     
%     dT_allData_S_Data{pDF, 1} = (detrend(allData_S_Data, 'constant'))';
%     
%     %% CBV Data associated with fileID
%     allData_CBV = (ProcData.Data.CBV.(CBVdataType)(1:8970) - RestingBaselines.CBV.(CBVdataType).(strDay)) / RestingBaselines.CBV.(CBVdataType).(strDay);
%     [B, A] = butter(4, 2 / (30 / 2), 'low');
%     filt_allData_CBV = filtfilt(B, A, allData_CBV);
%     dS_allData_CBV = downsample(filt_allData_CBV, 6);
%     dT_allData_CBV{pDF, 1} = detrend(dS_allData_CBV, 'constant');
% end
% 
% allData_F = SpectrogramData.(neuralDataType).OneSec.F{pDF, 1};
% allData_z_hold = [];
% allData_lagTime = 5;       % Seconds
% allData_frequency = 5;     % Hz
% allData_maxLag = allData_lagTime*allData_frequency;    % Number of points
% allData_XC_Vals = ones(size(allData_S_Data, 2), 2*allData_maxLag + 1);   % Pre-allocate size of cross-corr matrix
% 
% for x = 1:length(dT_allData_S_Data)
%     for y = 1:size(dT_allData_S_Data{x, 1}, 1)
%         allData_CBV_array = dT_allData_CBV{x, 1};
%         allData_Neural_array = dT_allData_S_Data{x, 1}(y, :);
%         [allData_XC_Vals(y, :), allData_lags] = xcorr(allData_CBV_array, allData_Neural_array, allData_maxLag, 'coeff');
%     end
%     allData_z_hold = cat(3, allData_z_hold, allData_XC_Vals);
% end
% allData_meanXC_Vals = mean(allData_z_hold, 3);
% 
% AllDataXCorr = figure;
% imagesc(allData_lags, allData_F, allData_meanXC_Vals)
% title([neuralDataType ' No Behavior Cross Correlation, One Sec Bins'])
% xticks([-allData_maxLag -allData_maxLag/2 0 allData_maxLag/2 allData_maxLag])
% xticklabels({'-5', '-2.5', '0', '2.5' '5'})
% xlabel('Lags (sec)')
% ylabel('Freq (Hz)')
% xlim([-allData_lagTime*allData_frequency allData_lagTime*allData_frequency])
% ylim([1 100])
% colorbar
% axis xy
% 
% ComparisonData.XCorr.AllData.(neuralDataType).lags = allData_lags;
% ComparisonData.XCorr.AllData.(neuralDataType).F = allData_F;
% ComparisonData.XCorr.AllData.(neuralDataType).XC_Vals = allData_meanXC_Vals;
% save([animal '_ComparisonData.mat'], 'ComparisonData');
% 
% [pathstr, ~, ~] = fileparts(cd);
% dirpath = [pathstr '/Figures/XCorr/'];
% 
% if ~exist(dirpath, 'dir')
%     mkdir(dirpath);
% end
% 
% savefig(AllDataXCorr, [dirpath animal '_' neuralDataType '_AllDataXCorr']);
% 
%% BLOCK PURPOSE: NREM SleepData
if ~isempty(SleepData)
    NREM_sleepTime = params.minTime.NREM;   % seconds
    NREM_allSleepFileIDs = SleepData.NREM.FileIDs;
    NREM_uniqueSleepFileIDs = unique(SleepData.NREM.FileIDs);
    NREM_spectrogramFileIDs = SpectrogramData.(neuralDataType).FileIDs;
    NREM_sleepBins = NREM_sleepTime / 5;
    x = 1;
    for uID = 1:length(NREM_uniqueSleepFileIDs)
        
        %% Pull out the bin times (there may be multiple events) in each unique NREM sleep file
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
        
        %% Pull out the Spectrogram data that matches the unique NREM sleep file
        for sID = 1:length(NREM_spectrogramFileIDs)
            NREM_spectrogramFileID = char(NREM_spectrogramFileIDs(sID));
            if strcmp(NREM_uniqueSleepFileID, NREM_spectrogramFileID)
                NREM_S_Data = SpectrogramData.(neuralDataType).OneSec.S_Norm{sID, 1};
            end
        end
        
        for rBN = 1:length(NREM_binTimes)
            NREM_Bins = NREM_binTimes{rBN, 1};
            NREM_x_Length = size(NREM_S_Data, 2);
            NREM_x_binLength = ceil(NREM_x_Length / 300);
            try
                NREM_sleepNeural_Vals{x, 1} = NREM_S_Data(:, (NREM_Bins(1) - 5)*NREM_x_binLength + 1:(NREM_Bins(NREM_sleepBins))*NREM_x_binLength);
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
    
    %% CBV
    CBVSamplingRate = 30;
    for CBV = 1:length(SleepData.NREM.Data.CBV.(CBVdataType))
        NREM_CBV_Vals = SleepData.NREM.Data.CBV.(CBVdataType){CBV, 1}(1:(NREM_sleepTime*CBVSamplingRate));
        NREM_dsCBV_Vals = downsample(NREM_CBV_Vals, 6);
        NREM_dT_sleepCBV_Vals{CBV, 1} = detrend(NREM_dsCBV_Vals, 'constant');
    end
    
    %% Cross Correlation
    NREM_F = SpectrogramData.(neuralDataType).OneSec.F{uID, 1};
    NREM_z_hold = [];
    NREM_lagTime = 5;       % Seconds
    NREM_frequency = 5;     % Hz
    NREM_maxLag = NREM_lagTime*NREM_frequency;    % Number of points
    NREM_XC_Vals = ones(size(NREM_ind_sleepNeural_Vals, 2), 2*NREM_maxLag + 1);   % Pre-allocate size of cross-corr matrix
    
    for x = 1:length(NREM_dT_sleepNeural_Vals)
        for y = 1:size(NREM_dT_sleepNeural_Vals{x, 1}, 1)
            NREM_CBV_array = NREM_dT_sleepCBV_Vals{x, 1};
            NREM_Neural_array = NREM_dT_sleepNeural_Vals{x, 1}(y, :);
            [NREM_XC_Vals(y, :), NREM_lags] = xcorr(NREM_CBV_array, NREM_Neural_array, NREM_maxLag, 'coeff');
        end
        NREM_z_hold = cat(3, NREM_z_hold, NREM_XC_Vals);
    end
    NREM_meanXC_Vals = mean(NREM_z_hold, 3);
    
%     nremboostrapResults = [];
%     for a = 1:length(Data.CBV)
%         disp(['Boostrapping CBV event ' num2str(a) ' of ' num2str(length(Data.CBV))]); disp(' ')
%         for boot = 1:1000
%             shuffledCBV_array_Inds = randperm(length(NREM_dT_sleepCBV_Vals{a, 1}));
%             shuffledCBV_array = NREM_dT_sleepCBV_Vals{a,1}(shuffledCBV_array_Inds);
%             nremshufXC_Vals = ones(length(F), 2*maxLag + 1);   % Pre-allocate size of cross-corr matrix
%             for b = 1:size(NREM_S_Data{a, 1}, 1)
%                 NREM_Neural_array = NREM_S_Data{a, 1}(b, :);
%                 [nremshufXC_Vals(b, :), lags] = xcorr(shuffledCBV_array, NREM_Neural_array, maxLag, 'coeff');
%             end
%             nremboostrapResults = cat(3, nremboostrapResults, nremshufXC_Vals);
%         end
%     end
%     
%     test3 = prctile(nremboostrapResults, [2.5 97.5], 3);
%     test4 = test3(:,:,1) - test3(:,:,2);

    NREMXCorr = figure;
    imagesc(NREM_lags, NREM_F, NREM_meanXC_Vals)
    title([neuralDataType ' NREM Cross Correlation, One Sec Bins'])
    xticks([-NREM_maxLag -NREM_maxLag/2 0 NREM_maxLag/2 NREM_maxLag])
    xticklabels({'-5', '-2.5', '0', '2.5' '5'})
    xlabel('Lags (sec)')
    ylabel('Freq (Hz)')
    xlim([-NREM_lagTime*NREM_frequency NREM_lagTime*NREM_frequency])
    ylim([1 100])
    colorbar
    axis xy
    
    ComparisonData.XCorr.NREM.(neuralDataType).lags = NREM_lags;
    ComparisonData.XCorr.NREM.(neuralDataType).F = NREM_F;
    ComparisonData.XCorr.NREM.(neuralDataType).XC_Vals = NREM_meanXC_Vals;
    save([animal '_ComparisonData.mat'], 'ComparisonData');
    
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/XCorr/'];
    
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    
    savefig(NREMXCorr, [dirpath animal '_' neuralDataType '_NREMXCorr']);
    
    %% BLOCK PURPOSE: REM SleepData
    if ~isempty(SleepData.REM)
        REM_sleepTime = params.minTime.REM;   % seconds
        REM_allSleepFileIDs = SleepData.REM.FileIDs;
        REM_uniqueSleepFileIDs = unique(SleepData.REM.FileIDs);
        REM_spectrogramFileIDs = SpectrogramData.(neuralDataType).FileIDs;
        REM_sleepBins = REM_sleepTime / 5;
        x = 1;
        for uID = 1:length(REM_uniqueSleepFileIDs)
            
            %% Pull out the bin times (there may be multiple events) in each unique REM sleep file
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
            
            %% Pull out the Spectrogram data that matches the unique REM sleep file
            for sID = 1:length(REM_spectrogramFileIDs)
                REM_spectrogramFileID = char(REM_spectrogramFileIDs(sID));
                if strcmp(REM_uniqueSleepFileID, REM_spectrogramFileID)
                    REM_S_Data = SpectrogramData.(neuralDataType).OneSec.S_Norm{sID, 1};
                end
            end
            
            for rBN = 1:length(REM_binTimes)
                REM_Bins = REM_binTimes{rBN, 1};
                REM_x_Length = size(REM_S_Data, 2);
                REM_x_binLength = ceil(REM_x_Length / 300);
                try
                    REM_sleepNeural_Vals{x, 1} = REM_S_Data(:, (REM_Bins(1) - 5)*REM_x_binLength + 1:(REM_Bins(REM_sleepBins))*REM_x_binLength);
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
        
        %% CBV
        CBVSamplingRate = 30;
        for CBV = 1:length(SleepData.REM.Data.CBV.(CBVdataType))
            REM_CBV_Vals = SleepData.REM.Data.CBV.(CBVdataType){CBV, 1}(1:(REM_sleepTime*CBVSamplingRate));
            REM_dsCBV_Vals = downsample(REM_CBV_Vals, 6);
            REM_dT_sleepCBV_Vals{CBV, 1} = detrend(REM_dsCBV_Vals, 'constant');
        end
        
        %% Cross Correlation
        REM_F = SpectrogramData.(neuralDataType).OneSec.F{uID, 1};
        REM_z_hold = [];
        REM_lagTime = 5;       % Seconds
        REM_frequency = 5;     % Hz
        REM_maxLag = REM_lagTime*REM_frequency;    % Number of points
        REM_XC_Vals = ones(size(REM_ind_sleepNeural_Vals, 2), 2*REM_maxLag + 1);   % Pre-allocate size of cross-corr matrix
        
        for x = 1:length(REM_dT_sleepNeural_Vals)
            for y = 1:size(REM_dT_sleepNeural_Vals{x, 1}, 1)
                REM_CBV_array = REM_dT_sleepCBV_Vals{x, 1};
                REM_Neural_array = REM_dT_sleepNeural_Vals{x, 1}(y, :);
                [REM_XC_Vals(y, :), REM_lags] = xcorr(REM_CBV_array, REM_Neural_array, REM_maxLag, 'coeff');
            end
            REM_z_hold = cat(3, REM_z_hold, REM_XC_Vals);
        end
        REM_meanXC_Vals = mean(REM_z_hold, 3);
        
        REMXCorr = figure;
        imagesc(REM_lags, REM_F, REM_meanXC_Vals)
        title([neuralDataType ' REM Cross Correlation, One Sec Bins'])
        xticks([-REM_maxLag -REM_maxLag/2 0 REM_maxLag/2 REM_maxLag])
        xticklabels({'-5', '-2.5', '0', '2.5' '5'})
        xlabel('Lags (sec)')
        ylabel('Freq (Hz)')
        xlim([-REM_lagTime*REM_frequency REM_lagTime*REM_frequency])
        ylim([1 100])
        colorbar
        axis xy
        
        ComparisonData.XCorr.REM.(neuralDataType).lags = REM_lags;
        ComparisonData.XCorr.REM.(neuralDataType).F = REM_F;
        ComparisonData.XCorr.REM.(neuralDataType).XC_Vals = REM_meanXC_Vals;
        save([animal '_ComparisonData.mat'], 'ComparisonData');
        
        [pathstr, ~, ~] = fileparts(cd);
        dirpath = [pathstr '/Figures/XCorr/'];
        
        if ~exist(dirpath, 'dir')
            mkdir(dirpath);
        end
        
        savefig(REMXCorr, [dirpath animal '_' neuralDataType '_REMXCorr']);
    end
    
end

end
