function [AnalysisResults] = AnalyzeXCorr_IOS(CBVdataType, neuralDataType, baselineType, params, AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the cross-correlation between a CBV signal and a spectrogram during different behaviors
%________________________________________________________________________________________________________________________
%
%   Inputs: //
%
%   Outputs: //
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

fileBreaks = strfind(restDataFileID, '_');
animalID = restDataFileID(1:fileBreaks(1)-1);

trialDuration = 900;   % sec
trialDurationMin = 15;   % min
sleepBinWidth = 5;   % sec
fiveSecSpecFs = 5;   % sec
CBVsamplingRate = RestData.data.CBV.LH.CBVCamSamplingRate;

%% Cross-correlation analysis for resting data
disp(['AnalyzeXCorr: ' CBVdataType ' vs ' neuralDataType ' during Rest.']); disp(' ')
% set criteria for rest event filter
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};

PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};

% Use the RestCriteria we specified earlier to find all resting events that are greater than the criteria
[restLogical] = FilterEvents_IOS(RestData.CBV.(CBVdataType), RestCriteria);   % Output is e logical
[puffLogical] = FilterEvents_IOS(RestData.CBV.(CBVdataType), PuffCriteria);   % Output is e logical
combRestLogical = logical(restLogical.*puffLogical);
allRestFiles = RestData.CBV.(CBVdataType).fileIDs(combRestLogical, :);   % Overall logical for all resting file names that meet criteria
allRestDurations = RestData.CBV.(CBVdataType).durations(combRestLogical, :);
allRestEventTimes = RestData.CBV.(CBVdataType).eventTimes(combRestLogical, :);
allRestingData = RestData.CBV.(CBVdataType).data(combRestLogical, :);   % Pull out data from all those resting files that meet criteria

uniqueDays = GetUniqueDays_IOS(RestData.CBV.(CBVdataType).fileIDs);   % Find the unique days of imaging
uniqueFiles = unique(RestData.CBV.(CBVdataType).fileIDs);   
numberOfFiles = length(unique(RestData.CBV.(CBVdataType).fileIDs));   % Find the number of unique files
fileTarget = params.targetMinutes/trialDurationMin;  

% Loop through each unique day in order to create e logical to filter the file list so that it only includes the first
% x number of files that fall within the targetMinutes requirement
for e = 1:length(uniqueDays)
    day = uniqueDays(e);
    b = 1;
    for c = 1:numberOfFiles
        file = uniqueFiles(c);
        fileID = file{1}(1:6);
        if strcmp(day, fileID) && b <= fileTarget
            filtLogical{e, 1}(c, 1) = 1;
            b = b + 1;
        else
            filtLogical{e, 1}(c, 1) = 0;
        end
    end
end

% Combine the 3 logicals so that it reflects the first "x" number of files from each day
finalLogical = any(sum(cell2mat(filtLogical'), 2), 2);

% Now that the appropriate files from each day are identified, loop through each file name with respect to the original
% list of ALL resting files, only keeping the ones that fall within the first targetMinutes of each day.
filtRestFiles = uniqueFiles(finalLogical, :);
for d = 1:length(allRestFiles)
    logic = strcmp(allRestFiles{d}, filtRestFiles);
    logicSum = sum(logic);
    if logicSum == 1
        fileFilter(d, 1) = 1;
    else
        fileFilter(d, 1) = 0;
    end
end

finalFileFilter = logical(fileFilter);
finalFileIDs = allRestFiles(finalFileFilter, :);
finalFileDurations = allRestDurations(finalFileFilter, :);
finalFileEventTimes = allRestEventTimes(finalFileFilter, :);
finalRestData = allRestingData(finalFileFilter, :);

for e = 1:length(finalFileIDs)   
    fileID = finalFileIDs{e, 1}; 
    strDay = fileID(1:6); 
    date = ConvertDate_IOS(strDay);
    
    %% Load in CBV from rest period
    CBV = (finalRestData{e, 1} - RestingBaselines.(baselineType).CBV.(CBVdataType).(date))/RestingBaselines.(baselineType).CBV.(CBVdataType).(date);   % Loop through the list of rest CBV data - the index is consistent with the file number
    
    % Low pass filter the epoch below 1 Hz
    [B, A] = butter(4, 1/(30/2), 'low');
    filtCBV = filtfilt(B, A, CBV);
    
    % Only take the first 10 seconds of the epoch
    shortCBV = filtCBV(1:params.minTime.Rest*CBVSamplingRate);
    
    % Downsample the 10 second epoch to 5 Hz
    dsCBV = downsample(shortCBV, 6);
    
    % Mean subtract the downsampled epoch
    Data.CBV{e, 1} = detrend(dsCBV, 'constant');
    
    
    %% Load in Neural Data from rest period
    specDataFileID = [animalID '_' fileID '_SpecData.mat'];
    load(specDataFileID);
    S_data = SpecData.(neuralDataType).fiveSec.normS;
    sLength = size(S_data, 2);                                  
    binSize = ceil(sLength/trialDuration); 
    samplingDiff = CBVSamplingRate/binSize;
    
    % Find the start time and duration
    restDuration = floor(floor(finalFileDurations(e, 1)*CBVSamplingRate)/samplingDiff);
    startTime = floor(floor(finalFileEventTimes(e, 1)*CBVSamplingRate)/samplingDiff);
    if startTime == 0
        startTime = 1;
    end
    
    % Take the S_data from the start time throughout the duration
    try
        restS_Vals = S_data(:, (startTime:(startTime + restDuration)));
    catch
        restS_Vals = S_data(:, end - restDuration:end);
    end
    
    % Only take the first min rest time seconds
    shortRestS_Vals = restS_Vals(:, 1:params.minTime.Rest*fiveSecSpecFs);
    
    % Mean subtract each row with detrend
    transpRestS_Vals = shortRestS_Vals';        % Transpose since detrend goes down columns
    dtRestS_Vals = detrend(transpRestS_Vals);   % detrend
    Data.S{e, 1} = dtRestS_Vals';              % transpose back to original orientation and save into Data.S
end

F = SpecData.(neuralDataType).fiveSec.F;
z_hold = [];
lagTime = 5;       % Seconds
frequency = fiveSecSpecFs;     % Hz
maxLag = lagTime*frequency;    % Number of points
XC_Vals = ones(length(F), 2*maxLag + 1);   % Pre-allocate size of cross-corr matrix

for e = 1:length(Data.CBV)
    for b = 1:size(Data.S{e, 1}, 1)
        CBV_array = Data.CBV{e, 1};
        Neural_array = Data.S{e, 1}(b, :);
        [XC_Vals(b, :), lags] = xcorr(CBV_array, Neural_array, maxLag, 'coeff');
    end
    z_hold = cat(3, z_hold, XC_Vals);
end

meanXC_Vals = mean(z_hold, 3);

titleID = strrep(CBVdataType, '_', ' ');
RestingXCorr = figure;
imagesc(lags, F, meanXC_Vals)
title([animalID ' ' titleID ' Resting Cross Correlation'])
xticks([-maxLag -maxLag/2 0 maxLag/2 maxLag])
xticklabels({'-5', '-2.5', '0', '2.5' '5'})
xlim([-lagTime*frequency lagTime*frequency])
xlabel('Lags (sec)')
ylabel('Freq (Hz)')
ylim([1 100])
colorbar
axis xy

AnalysisResults.XCorr.Rest.(CBVdataType).lags = lags;
AnalysisResults.XCorr.Rest.(CBVdataType).F = F;
AnalysisResults.XCorr.Rest.(CBVdataType).XC_Vals = meanXC_Vals;
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Analysis XCorr/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(RestingXCorr, [dirpath animalID '_' CBVdataType '_RestingXCorr']);

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
for CBV = 1:length(SleepData.NREM.data.CBV.(CBVdataType))
    NREM_CBV_Vals = SleepData.NREM.data.CBV.(CBVdataType){CBV, 1}(1:(NREM_sleepTime*CBVSamplingRate));
    NREM_dsCBV_Vals = downsample(NREM_CBV_Vals, 6);
    NREM_dT_sleepCBV_Vals{CBV, 1} = detrend(NREM_dsCBV_Vals, 'constant');
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
for CBV = 1:length(SleepData.REM.data.CBV.(CBVdataType))
    REM_CBV_Vals = SleepData.REM.data.CBV.(CBVdataType){CBV, 1}(1:(REM_sleepTime*CBVSamplingRate));
    REM_dsCBV_Vals = downsample(REM_CBV_Vals, 6);
    REM_dT_sleepCBV_Vals{CBV, 1} = detrend(REM_dsCBV_Vals, 'constant');
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
save([animalID '_AnalysisResults'], 'AnalysisResults')

end
