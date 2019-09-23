function [AnalysisResults] = AnalyzeEvokedResponses_IOS(dataType, params, AnalysisResults)
%___________________________________________________________________________________________________
% Written by Kevin L. Turner, Jr.
% Adapted from codes credited to Dr. Patrick J. Drew and Aaron T. Winder
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%___________________________________________________________________________________________________
%
%   Purpose:
%___________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%___________________________________________________________________________________________________

% find and load RestData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID)

fileBreaks = strfind(eventDataFileID, '_');
animalID = eventDataFileID(1:fileBreaks(1)-1);
trialDuration = 900;
trialDurationMin = 15;
samplingRate = EventData.CBV.(dataType).whisk.samplingRate;      % Can be LH or RH, does not matter.
timeVector = (0:(EventData.CBV.(dataType).whisk.epoch.duration*samplingRate))/samplingRate - EventData.CBV.(dataType).whisk.epoch.offset; % 12 seconds

%% Whisking-evoked responses
% Criteria for the FilterEvents data struct
whiskCriteria.Fieldname = {'duration', 'duration', 'puffDistance'};
whiskCriteria.Comparison = {'gt','lt','gt'};
whiskCriteria.Value = {0.3, 2, 5};

allWhiskFilter = FilterEvents_IOS(EventData.CBV.(dataType).whisk, whiskCriteria);

% CBV whisk
[allWhiskCBVData] = EventData.CBV.(dataType).whisk.NormData(allWhiskFilter, :);
[allWhiskFileIDs] = EventData.CBV.(dataType).whisk.fileIDs(allWhiskFilter, :);
[allWhiskEventTimes] = EventData.CBV.(dataType).whisk.eventTime(allWhiskFilter, :);

whiskUniqueDays = GetUniqueDays_IOS(allWhiskFileIDs);
whiskUniqueFiles = unique(allWhiskFileIDs);
whiskNumberOfFiles = length(unique(allWhiskFileIDs));
fileTarget = params.targetMinutes/trialDurationMin;

for uD = 1:length(whiskUniqueDays)
    whiskDay = whiskUniqueDays(uD);
    x = 1;
    for nOF = 1:whiskNumberOfFiles
        file = whiskUniqueFiles(nOF);
        whiskFileID = file{1}(1:6);
        if strcmp(whiskDay, whiskFileID) && x <= fileTarget
            whiskFiltLogical{uD, 1}(nOF, 1) = 1;
            x = x + 1;
        else
            whiskFiltLogical{uD, 1}(nOF, 1) = 0;
        end
    end
end

% Combine the 3 logicals so that it reflects the first "x" number of files from each day
whiskFinalLogical = any(sum(cell2mat(whiskFiltLogical'), 2), 2);

filtWhiskFiles = whiskUniqueFiles(whiskFinalLogical, :);
for rF = 1:length(allWhiskFileIDs)
    logic = strcmp(allWhiskFileIDs{rF}, filtWhiskFiles);
    logicSum = sum(logic);
    if logicSum == 1
        whiskFileFilter(rF, 1) = 1;
    else
        whiskFileFilter(rF, 1) = 0;
    end
end

finalWhiskFileFilter = logical(whiskFileFilter);
finalWhiskCBVData = allWhiskCBVData(finalWhiskFileFilter, :);
finalWhiskFileIDs = allWhiskFileIDs(finalWhiskFileFilter, :);
finalWhiskFileEventTimes = allWhiskEventTimes(finalWhiskFileFilter, :);

offset = 2;
for a = 1:size(finalWhiskCBVData,1)
    cbvArray = finalWhiskCBVData(a,:);
%          filtCBVarray = sgolayfilt(cbvArray, 3, 17)*100;
    procCBVarray = cbvArray - mean(cbvArray(1:(offset*samplingRate)));
    cbvData(a,:) = procCBVarray;
end

meanWhiskCBVData = mean(cbvData, 1);

% LFP
whiskZhold = [];
for ii = 1:length(finalWhiskFileIDs)   % Loop through each non-unique file
    whiskFileID = finalWhiskFileIDs{ii, 1};
    whiskStrDay = whiskFileID(1:6);
    whiskDate = ConvertDate_IOS(whiskStrDay);
    
    % Load in Neural Data from rest period
    specDataFileID = [animalID '_' whiskFileID '_SpecData.mat'];
    specField = ['cortical_' dataType];
    load(specDataFileID)
    whiskS_Data = SpecData.(specField).oneSec.normS;
    whiskSLength = size(whiskS_Data, 2);
    whiskBinSize = ceil(whiskSLength/trialDuration);
    whiskSamplingDiff = samplingRate/whiskBinSize;
    
    % Find the start time and duration
    whiskDuration = floor(floor(size(meanWhiskCBVData, 2)) / samplingRate);
    startTime = floor(floor(finalWhiskFileEventTimes(ii, 1)*samplingRate) / whiskSamplingDiff);
    if startTime == 0
        startTime = 1;
    end
    
    % Take the S_data from the start time throughout the duration
    try
        whiskS_Vals = whiskS_Data(:, (startTime - (offset*whiskBinSize)):(startTime + ((whiskDuration - offset)*whiskBinSize)));
    catch
        whiskS_Vals = whiskS_Data(:, end - (whiskDuration*whiskBinSize):end);
    end
    
    % Mean subtract each row with detrend
    transpWhiskS_Vals = whiskS_Vals';   % Transpose since detrend goes down columns
    dTWhiskS_Vals = detrend(transpWhiskS_Vals);   % detrend
    whiskZhold = cat(3, whiskZhold, dTWhiskS_Vals');   % transpose back to original orientation and save into Data.S
end

T = 1:whiskDuration;
F = SpecData.(specField).oneSec.F;
whiskS = mean(whiskZhold, 3);

% Figure
whiskEvoked = figure;
subplot(2,1,1);
plot(timeVector, meanWhiskCBVData*100)
xticklabels('')
ylabel('Reflectance (%)')
axis tight
title([animalID ' ' dataType ' Whisking-evoked averages'])

subplot(2,1,2);
imagesc(T, F, whiskS)
set(gca, 'Ticklength', [0 0])
xticklabels('')
ylim([1 100])
ylabel('Freq (Hz)')
colorbar
caxis([-0.5 1])
axis xy

AnalysisResults.EvokedAvgs.Whisk.(dataType).CBV.data = meanWhiskCBVData;
AnalysisResults.EvokedAvgs.Whisk.(dataType).CBV.timeVector = timeVector;
AnalysisResults.EvokedAvgs.Whisk.(dataType).LFP.data = whiskS;
AnalysisResults.EvokedAvgs.Whisk.(dataType).LFP.T = T;
AnalysisResults.EvokedAvgs.Whisk.(dataType).LFP.F = F;

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Analysis Evoked Averages/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(whiskEvoked, [dirpath animalID '_' dataType '_WhiskEvokedAverages']);


%% Stimulus-evoked responses
if ~isempty(EventData.CBV.(dataType).stim.fileIDs)
    
    % Criteria for the FilterEvents data struct
    stimCriteria.Fieldname = {'solenoidName'};
    stimCriteria.Comparison = {'equal'};
    if strcmp(dataType, 'LH')
        stimCriteria.Value = {'RPadSol'};
    else
        stimCriteria.Value = {'LPadSol'};
    end
    
    allStimFilter = FilterEvents_IOS(EventData.CBV.(dataType).stim, stimCriteria);
    
    % CBV whisk
    [allStimCBVData] = EventData.CBV.(dataType).stim.NormData(allStimFilter, :);
    [allStimFileIDs] = EventData.CBV.(dataType).stim.fileIDs(allStimFilter, :);
    [allStimEventTimes] = EventData.CBV.(dataType).stim.eventTime(allStimFilter, :);
    
    stimUniqueDays = GetUniqueDays_IOS(allStimFileIDs);
    stimUniqueFiles = unique(allStimFileIDs);
    stimNumberOfFiles = length(unique(allStimFileIDs));
    fileTarget = params.targetMinutes/trialDurationMin;
    
    for uD = 1:length(stimUniqueDays)
        stimDay = stimUniqueDays(uD);
        x = 1;
        for nOF = 1:stimNumberOfFiles
            file = stimUniqueFiles(nOF);
            stimFileID = file{1}(1:6);
            if strcmp(stimDay, stimFileID) && x <= fileTarget
                stimFiltLogical{uD, 1}(nOF, 1) = 1;
                x = x + 1;
            else
                stimFiltLogical{uD, 1}(nOF, 1) = 0;
            end
        end
    end
    
    % Combine the 3 logicals so that it reflects the first "x" number of files from each day
    stimFinalLogical = any(sum(cell2mat(stimFiltLogical'), 2), 2);
    
    filtStimFiles = stimUniqueFiles(stimFinalLogical, :);
    for rF = 1:length(allStimFileIDs)
        logic = strcmp(allStimFileIDs{rF}, filtStimFiles);
        logicSum = sum(logic);
        if logicSum == 1
            stimFileFilter(rF, 1) = 1;
        else
            stimFileFilter(rF, 1) = 0;
        end
    end
    
    finalStimFileFilter = logical(stimFileFilter);
    
    finalStimCBVData = allStimCBVData(finalStimFileFilter, :);
    finalStimFileIDs = allStimFileIDs(finalStimFileFilter, :);
    finalStimFileEventTimes = allStimEventTimes(finalStimFileFilter, :);
    
    [B, A] = butter(4, 1/(30/2), 'low');
    for a = 1:size(finalStimCBVData,1)
        cbvArray = finalStimCBVData(a,:);
        %          filtCBVarray = sgolayfilt(cbvArray, 3, 17)*100;
        procCBVarray = cbvArray - mean(cbvArray(1:(offset*samplingRate)));
        stimCBVdata(a,:) = procCBVarray;
    end
    
    meanStimCBVData = mean(stimCBVdata, 1);
    
    % LFP
    stimZhold = [];
    for ii = 1:length(finalStimFileIDs)   % Loop through each non-unique file
        stimFileID = finalStimFileIDs{ii, 1};
        stimStrDay = stimFileID(1:6);
        stimDate = ConvertDate(stimStrDay);
        
        % Load in Neural Data from rest period
        specDataFileID = [animalID '_' stimFileID '_SpecData.mat'];
        specField = ['cortical_' dataType];
        load(specDataFileID)
        stimSLength = size(stimS_Data, 2);
        stimBinSize = ceil(stimSLength/trialDuration);
        stimSamplingDiff = samplingRate/stimBinSize;
        
        % Find the start time and duration
        stimDuration = floor(floor(size(dTStimMUAData, 1))/samplingRate);
        startTime = floor(floor(finalStimFileEventTimes(ii,1)*samplingRate)/stimSamplingDiff);
        if startTime == 0
            startTime = 1;
        end
        
        % Take the S_data from the start time throughout the duration
        try
            stimS_Vals = stimS_Data(:, (startTime - (offset*stimBinSize)):(startTime + ((stimDuration - offset)*stimBinSize)));
        catch
            stimS_Vals = stimS_Data(:, end - (stimDuration*stimBinSize):end);
        end
    
        % Mean subtract each row with detrend
        transpStimS_Vals = stimS_Vals';   % Transpose since detrend goes down columns
        dTStimS_Vals = detrend(transpStimS_Vals);   % detrend
        stimZhold = cat(3, stimZhold, dTStimS_Vals');   % transpose back to original orientation and save into Data.S
    end
    
    stimS = mean(stimZhold, 3);
    
    % Figure
    stimEvoked = figure;
    subplot(2,1,1);
    plot(timeVector, meanStimCBVData*100)
    xticklabels('')
    ylabel('Reflectance (%)')
    axis tight
    title([animalID ' ' dataType ' (Contralateral) Stimulus-evoked averages'])
    
    subplot(2,1,2);
    imagesc(T, F, stimS)
    set(gca, 'Ticklength', [0 0])
    xticklabels('')
    ylim([1 100])
    ylabel('Freq (Hz)')
    colorbar
    caxis([-0.5 1])
    axis xy
    
    AnalysisResults.EvokedAvgs.Stim.(dataType).CBV.data = meanStimCBVData;
    AnalysisResults.EvokedAvgs.Stim.(dataType).CBV.timeVector = timeVector;
    AnalysisResults.EvokedAvgs.Stim.(dataType).LFP.data = stimS;
    AnalysisResults.EvokedAvgs.Stim.(dataType).LFP.T = T;
    AnalysisResults.EvokedAvgs.Stim.(dataType).LFP.F = F;
        
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Analysis Evoked Averages/'];
    
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    
    savefig(stimEvoked, [dirpath animalID '_' dataType '_StimEvokedAverages']);
end

save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

end



