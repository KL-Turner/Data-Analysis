function [ComparisonData] = AnalyzeEvokedResponses_2P(RestingBaselines, EventData, SpectrogramData)
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

samplingRate = 30;      % Can be LH or RH, does not matter.

%% Whisking-evoked responses
% Criteria for the FilterEvents data struct
whiskCriteria1.Fieldname = {'duration', 'duration', 'puffDistance'};
whiskCriteria1.Comparison = {'gt','lt','gt'};
whiskCriteria1.Value = {0.5, 2, 5};

whiskCriteria2.Fieldname = {'duration', 'duration', 'puffDistance'};
whiskCriteria2.Comparison = {'gt','lt','gt'};
whiskCriteria2.Value = {2, 5, 5};

whiskCriteria3.Fieldname = {'duration', 'puffDistance'};
whiskCriteria3.Comparison = {'gt','gt'};
whiskCriteria3.Value = {5, 5};

allWhiskFilter1 = FilterEvents(EventData.Vessel_Diameter.whisk, whiskCriteria1);
allWhiskFilter2 = FilterEvents(EventData.Vessel_Diameter.whisk, whiskCriteria2);
allWhiskFilter3 = FilterEvents(EventData.Vessel_Diameter.whisk, whiskCriteria3);

[allWhiskData1] = EventData.Vessel_Diameter.whisk.data(allWhiskFilter1, :);
[allWhiskVesselIDs1] = EventData.Vessel_Diameter.whisk.vesselIDs(allWhiskFilter1, :);
[allWhiskEventTimes1] = EventData.Vessel_Diameter.whisk.eventTime(allWhiskFilter1, :);

[allWhiskData2] = EventData.Vessel_Diameter.whisk.data(allWhiskFilter2, :);
[allWhiskVesselIDs2] = EventData.Vessel_Diameter.whisk.vesselIDs(allWhiskFilter2, :);
[allWhiskEventTimes2] = EventData.Vessel_Diameter.whisk.eventTime(allWhiskFilter2, :);

[allWhiskData3] = EventData.Vessel_Diameter.whisk.data(allWhiskFilter3, :);
[allWhiskVesselIDs3] = EventData.Vessel_Diameter.whisk.vesselIDs(allWhiskFilter3, :);
[allWhiskEventTimes3] = EventData.Vessel_Diameter.whisk.eventTime(allWhiskFilter3, :);

whiskUniqueVessels = unique(allWhiskVesselIDs1);
data = [];
for x = 1:length(whiskUniqueVessels)
    vID = whiskUniqueVessels(x);
    a = 1;
    for y = 1:length(allWhiskVesselIDs2)
        vesselID = allWhiskVesselIDs2(y);
        if strcmp(vID, vesselID)
            data{x,1}(a,:) = allWhiskData2(y, :);
            a = a+1;
        end
    end
end










filtWhiskFiles = whiskUniqueVessels(whiskFinalLogical, :);
for rF = 1:length(allWhiskVesselIDs1)
    logic = strcmp(allWhiskVesselIDs1{rF}, filtWhiskFiles);
    logicSum = sum(logic);
    if logicSum == 1
        whiskFileFilter(rF, 1) = 1;
    else
        whiskFileFilter(rF, 1) = 0;
    end
end

finalWhiskFileFilter = logical(whiskFileFilter);

finalWhiskCBVData = allWhiskData1(finalWhiskFileFilter, :);
finalWhiskFileIDs = allWhiskVesselIDs1(finalWhiskFileFilter, :);
finalWhiskFileEventTimes = allWhiskEventTimes1(finalWhiskFileFilter, :);

[B, A] = butter(4, 2 / (30 / 2), 'low');
filtWhiskCBVData = filtfilt(B, A, finalWhiskCBVData')';
transposedWhiskCBVData = filtWhiskCBVData';
dTWhiskCBVData = detrend(transposedWhiskCBVData, 'constant');
meanWhiskCBVData = mean(dTWhiskCBVData', 1);

% MUA whisk
[allWhiskMUAData] = EventData.MUA_Power.(dataType).whisk.NormData(allWhiskFilter1, :);
[allWhiskVesselIDs1] = EventData.MUA_Power.(dataType).whisk.fileIDs(allWhiskFilter1, :);
[allWhiskEventTimes1] = EventData.MUA_Power.(dataType).whisk.eventTime(allWhiskFilter1, :);

finalWhiskMUAData = allWhiskMUAData(finalWhiskFileFilter, :);

% filtWhiskMUAData = filt(B, A, finalWhiskMUAData')';
transposedWhiskMUAData = finalWhiskMUAData';
dTWhiskMUAData = detrend(transposedWhiskMUAData, 'constant');
meanWhiskMUAData = mean(dTWhiskMUAData', 1);

% LFP
whiskZhold = [];
for ii = 1:length(finalWhiskFileIDs)   % Loop through each non-unique file
    whiskFileID = finalWhiskFileIDs{ii, 1};
    whiskStrDay = whiskFileID(1:6);
    whiskDate = ConvertDate(whiskStrDay);
    
    % Load in Neural Data from rest period
    for s = 1:length(SpectrogramData.(dataType).FileIDs)
        if strcmp(whiskFileID, SpectrogramData.(dataType).FileIDs{s, 1})
            whiskS_Data = SpectrogramData.(dataType).OneSec.S_Norm{s, 1};  % S data for this specific file
        end
    end
    whiskSLength = size(whiskS_Data, 2);                                  % Length of the data across time (number of samples)
    whiskBinSize = ceil(whiskSLength / 300);                              % Find the number of bins needed to divide this into 300 seconds
    whiskSamplingDiff = samplingRate / whiskBinSize;                   % Number to divide by for data to be at 5 Hz
    
    % Find the start time and duration
    whiskDuration = floor(floor(size(dTWhiskMUAData, 1)) / samplingRate);
    startTime = floor(floor(finalWhiskFileEventTimes(ii, 1)*samplingRate) / whiskSamplingDiff);
    if startTime == 0
        startTime = 1;
    end
    
    % Take the S_data from the start time throughout the duration
    try
        whiskS_Vals = whiskS_Data(:, (startTime - (2*whiskBinSize)):(startTime + ((whiskDuration - 2)*whiskBinSize)));
    catch
        whiskS_Vals = whiskS_Data(:, end - (whiskDuration*whiskBinSize):end);
    end
    
    % Mean subtract each row with detrend
    transpWhiskS_Vals = whiskS_Vals';   % Transpose since detrend goes down columns
    dTWhiskS_Vals = detrend(transpWhiskS_Vals);   % detrend
    whiskZhold = cat(3, whiskZhold, dTWhiskS_Vals');   % transpose back to original orientation and save into Data.S
end

T = 1:whiskDuration;
F = SpectrogramData.(dataType).OneSec.F{1, 1};
whiskS = mean(whiskZhold, 3);

% Figure
whiskEvoked = figure;
ax1 = subplot(3,1,1);
plot(timeVector, meanWhiskCBVData*100)
xticklabels('')
ylabel('Reflectance (%)')
axis tight
title([animal ' ' dataType ' Whisking-evoked averages'])

ax2 = subplot(3,1,2);
imagesc(T, F, whiskS)
set(gca, 'Ticklength', [0 0])
xticklabels('')
ylim([1 100])
ylabel('Freq (Hz)')
colorbar
caxis([-0.5 1])
axis xy

ax3 =subplot(3,1,3);
plot(timeVector, meanWhiskMUAData)
axis tight
xlabel('Peristimulus time (sec)')
ylabel('Normalized Power')

ComparisonData.Evoked.Whisk.(dataType).CBV.data = meanWhiskCBVData;
ComparisonData.Evoked.Whisk.(dataType).CBV.timeVector = timeVector;
ComparisonData.Evoked.Whisk.(dataType).MUA.data = meanWhiskMUAData;
ComparisonData.Evoked.Whisk.(dataType).MUA.timeVector = timeVector;
ComparisonData.Evoked.Whisk.(dataType).LFP.data = whiskS;
ComparisonData.Evoked.Whisk.(dataType).LFP.T = T;
ComparisonData.Evoked.Whisk.(dataType).LFP.F = F;

save([animal '_ComparisonData.mat'], 'ComparisonData');

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Evoked Averages/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(whiskEvoked, [dirpath animal '_' dataType '_WhiskEvokedAverages']);


%% Stimulus-evoked responses
if ~isempty(EventData.CBV.(dataType).stim.fileIDs)
    
    % Criteria for the FilterEvents data struct
    stimCriteria.Fieldname = {'solenoidName'};
    stimCriteria.Comparison = {'equal'};
    if strcmp(dataType, 'LH')
        stimCriteria.Value = {'solenoidRightPad'};
    else
        stimCriteria.Value = {'solenoidLeftPad'};
    end
    
    allStimFilter = FilterEvents(EventData.CBV.(dataType).stim, stimCriteria);
    
    % CBV whisk
    [allStimCBVData] = EventData.CBV.(dataType).stim.NormData(allStimFilter, :);
    [allStimFileIDs] = EventData.CBV.(dataType).stim.fileIDs(allStimFilter, :);
    [allStimEventTimes] = EventData.CBV.(dataType).stim.eventTime(allStimFilter, :);
    
    stimUniqueDays = GetUniqueDays(allStimFileIDs);
    stimUniqueFiles = unique(allStimFileIDs);
    stimNumberOfFiles = length(unique(allStimFileIDs));
    fileTarget = params.targetMinutes / 5;
    
    for uD = 1:length(stimUniqueDays)
        stimDay = stimUniqueDays(uD);
        x = 1;
        for nOF = 1:stimNumberOfFiles
            file = stimUniqueFiles(nOF);
            stimFileID = file{1}(1:6);
            if strcmp(params.Infusion, 'y')
                if strcmp(stimDay, stimFileID) && x > fileTarget
                    stimFiltLogical{uD, 1}(nOF, 1) = 1;
                else
                    stimFiltLogical{uD, 1}(nOF, 1) = 0;
                    x = x + 1;
                end
            else
                if strcmp(stimDay, stimFileID) && x <= fileTarget
                    stimFiltLogical{uD, 1}(nOF, 1) = 1;
                    x = x + 1;
                else
                    stimFiltLogical{uD, 1}(nOF, 1) = 0;
                end
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
    
    [B, A] = butter(4, 2 / (30 / 2), 'low');
    filtStimCBVData = filtfilt(B, A, finalStimCBVData')';
    transposedStimCBVData = filtStimCBVData';
    dTStimCBVData = detrend(transposedStimCBVData, 'constant');
    meanStimCBVData = mean(dTStimCBVData', 1);
    
    % MUA whisk
    [allStimMUAData] = EventData.MUA_Power.(dataType).stim.NormData(allStimFilter, :);
    [allStimFileIDs] = EventData.MUA_Power.(dataType).stim.fileIDs(allStimFilter, :);
    [allStimEventTimes] = EventData.MUA_Power.(dataType).stim.eventTime(allStimFilter, :);
    
    finalStimMUAData = allStimMUAData(finalStimFileFilter, :);
    
%     filtStimMUAData = filtfilt(B, A, finalStimMUAData')';
    transposedStimMUAData = finalStimMUAData';
    dTStimMUAData = detrend(transposedStimMUAData, 'constant');
    meanStimMUAData = mean(dTStimMUAData', 1);
    
    % LFP
    stimZhold = [];
    for ii = 1:length(finalStimFileIDs)   % Loop through each non-unique file
        stimFileID = finalStimFileIDs{ii, 1};
        stimStrDay = stimFileID(1:6);
        stimDate = ConvertDate(stimStrDay);
        
        % Load in Neural Data from rest period
        for s = 1:length(SpectrogramData.(dataType).FileIDs)
            if strcmp(stimFileID, SpectrogramData.(dataType).FileIDs{s, 1})
                stimS_Data = SpectrogramData.(dataType).OneSec.S_Norm{s, 1};  % S data for this specific file
            end
        end
        stimSLength = size(stimS_Data, 2);                                  % Length of the data across time (number of samples)
        stimBinSize = ceil(stimSLength / 300);                              % Find the number of bins needed to divide this into 300 seconds
        stimSamplingDiff = samplingRate / stimBinSize;                   % Number to divide by for data to be at 5 Hz
        
        % Find the start time and duration
        stimDuration = floor(floor(size(dTStimMUAData, 1)) / samplingRate);
        startTime = floor(floor(finalStimFileEventTimes(ii, 1)*samplingRate) / stimSamplingDiff);
        if startTime == 0
            startTime = 1;
        end
        
        % Take the S_data from the start time throughout the duration
        try
            stimS_Vals = stimS_Data(:, (startTime - (2*stimBinSize)):(startTime + ((stimDuration - 2)*stimBinSize)));
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
    ax1 = subplot(3,1,1);
    plot(timeVector, meanStimCBVData*100)
    xticklabels('')
    ylabel('Reflectance (%)')
    axis tight
    title([animal ' ' dataType ' (Contralateral) Stimulus-evoked averages'])
    
    ax2 = subplot(3,1,2);
    imagesc(T, F, stimS)
    set(gca, 'Ticklength', [0 0])
    xticklabels('')
    ylim([1 100])
    ylabel('Freq (Hz)')
    colorbar
    caxis([-0.5 1])
    axis xy
    
    ax3 =subplot(3,1,3);
    plot(timeVector, meanStimMUAData)
    axis tight
    xlabel('Peristimulus time (sec)')
    ylabel('Normalized Power')
    
    ComparisonData.Evoked.Stim.(dataType).CBV.data = meanStimCBVData;
    ComparisonData.Evoked.Stim.(dataType).CBV.timeVector = timeVector;
    ComparisonData.Evoked.Stim.(dataType).MUA.data = meanStimMUAData;
    ComparisonData.Evoked.Stim.(dataType).MUA.timeVector = timeVector;
    ComparisonData.Evoked.Stim.(dataType).LFP.data = stimS;
    ComparisonData.Evoked.Stim.(dataType).LFP.T = T;
    ComparisonData.Evoked.Stim.(dataType).LFP.F = F;
    
    save([animal '_ComparisonData.mat'], 'ComparisonData');
    
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Evoked Averages/'];
    
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    
    savefig(stimEvoked, [dirpath animal '_' dataType '_StimEvokedAverages']);
end

end



