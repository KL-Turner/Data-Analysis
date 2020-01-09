function [] = CreateAnimalHypnograms_IOS()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Loads sleep scoring results structure and creates a hypnogram on the scored data across the days
%________________________________________________________________________________________________________________________

%% Identify the sleep scoring results structure in the current folder
sleepScoringResultsFileStruct = dir('*_SVM_SleepScoringResults.mat'); 
sleepScoringResultsFile = {sleepScoringResultsFileStruct.name}';
sleepScoringResultsFileID = char(sleepScoringResultsFile);
load(sleepScoringResultsFileID)
% Identify the unique file IDs, unique imaging days, and scoring labels from the file list
allScoringLabels = SVMResults.labels;
allFileIDs = SVMResults.fileIDs;
uniqueFileIDs = unique(allFileIDs);
for a = 1:size(uniqueFileIDs,1)
    [animalID,allFileDates{a,1},~] = GetFileInfo_IOS(uniqueFileIDs{a,1});
    allFileDays{a,1} = ConvertDate_IOS(allFileDates{a,1});
end
data.uniqueDates = unique(allFileDates);
data.uniqueDays = unique(allFileDays);
% Determine how many 5 second bins there are for each individual day
for b = 1:size(data.uniqueDays,1)
    data.dayBinLengths{b,1} = find(~cellfun('isempty',strfind(allFileIDs,data.uniqueDates{b,1})));
end
% Extract the file IDs and scoring labels that correspond to each day's indeces
for c = 1:size(data.dayBinLengths,1)
    dayInds = data.dayBinLengths{c,1};
    data.dayScoreLabels{c,1} = allScoringLabels(dayInds);
    data.dayScoreFileIDs{c,1} = allFileIDs(dayInds);
end
trialDuration = 15;   % minutes
binTime = 5;   % seconds
fileBinLength = (trialDuration*60)/binTime;
% Further break down each day's scores into the scores for each individual file
for d = 1:size(data.uniqueDays,1)
    uniqueDay = data.uniqueDays{d,1};
    uniqueDayFileIDs = unique(data.dayScoreFileIDs{d,1});
    for e = 1:size(uniqueDayFileIDs,1)
        if e == 1
            data.(uniqueDay).indFileData{e,1} = data.dayScoreLabels{d,1}(1:fileBinLength);
        else
            data.(uniqueDay).indFileData{e,1} = data.dayScoreLabels{d,1}((e - 1)*fileBinLength + 1:e*fileBinLength);
        end
    end
end
% Calculate the time difference between every file to append padding 'Time Pad' to the end of
% the leading file's score labels
for f = 1:size(data.uniqueDays,1)
    uniqueDay = data.uniqueDays{f,1};
    uniqueDayFileIDs = unique(data.dayScoreFileIDs{f,1});
    % start with file 2 to focus on the differences between each file
    for g = 2:size(uniqueDayFileIDs,1)
        leadFileID = uniqueDayFileIDs{g - 1,1};
        lagFileID = uniqueDayFileIDs{g,1};
        [~,~,leadFileInfo] = GetFileInfo_IOS(leadFileID);
        [~,~,lagFileInfo] = GetFileInfo_IOS(lagFileID);
        leadFileStr = ConvertDateTime_IOS(leadFileInfo);
        lagFileStr = ConvertDateTime_IOS(lagFileInfo);
        leadFileTime = datevec(leadFileStr);
        lagFileTime = datevec(lagFileStr);
        timeDifference = etime(lagFileTime,leadFileTime) - (trialDuration*60);   % seconds
        timePadBins.(uniqueDay){g - 1,1} = cell(floor(timeDifference/binTime),1);
        timePadBins.(uniqueDay){g - 1,1}(:) = {'Time Pad'};
    end
end
% Apply the time padding to the end of each file
for h = 1:size(data.uniqueDays,1)
    uniqueDay = data.uniqueDays{h,1};
    for j = 1:size(data.(uniqueDay).indFileData,1) - 1
        data.(uniqueDay).indFileData{j,1} = vertcat(data.(uniqueDay).indFileData{j,1},timePadBins.(uniqueDay){j,1});
    end
end
% Concatendate the data for each day now that the padding is added at the end of each file
for k = 1:size(data.uniqueDays,1)
    uniqueDay = data.uniqueDays{k,1};
    data.(uniqueDay).catData = [];
    for m = 1:size(data.(uniqueDay).indFileData,1)
        data.(uniqueDay).catData = vertcat(data.(uniqueDay).catData,data.(uniqueDay).indFileData{m,1});
    end
end

%% Hypnogram figure
% Prepare indeces for each behavior
for n = 1:size(data.uniqueDays,1)
    uniqueDay = data.uniqueDays{n,1};
    data.(uniqueDay).NotSleep_inds = NaN(1,size(data.(uniqueDay).catData,1));
    data.(uniqueDay).NREM_inds = NaN(1,size(data.(uniqueDay).catData,1));
    data.(uniqueDay).REM_inds = NaN(1,size(data.(uniqueDay).catData,1));
    data.(uniqueDay).TimePad_inds = NaN(1,size(data.(uniqueDay).catData,1));
    for o = 1:size(data.(uniqueDay).catData,1)
        if strcmp(data.(uniqueDay).catData{o,1},'Not Sleep') == true
            data.(uniqueDay).NotSleep_inds(1,o) = 1;
            data.(uniqueDay).NREM_inds(1,o) = NaN;
            data.(uniqueDay).REM_inds(1,o) = NaN;
            data.(uniqueDay).TimePad_inds(1,o) = NaN;
        elseif strcmp(data.(uniqueDay).catData{o,1},'NREM Sleep') == true
            data.(uniqueDay).NotSleep_inds(1,o) = NaN;
            data.(uniqueDay).NREM_inds(1,o) = 1;
            data.(uniqueDay).REM_inds(1,o) = NaN;
            data.(uniqueDay).TimePad_inds(1,o) = NaN;
        elseif strcmp(data.(uniqueDay).catData{o,1},'REM Sleep') == true
            data.(uniqueDay).NotSleep_inds(1,o) = NaN;
            data.(uniqueDay).NREM_inds(1,o) = NaN;
            data.(uniqueDay).REM_inds(1,o) = 1;
            data.(uniqueDay).TimePad_inds(1,o) = NaN;
        elseif strcmp(data.(uniqueDay).catData{o,1},'Day Pad') == true
            data.(uniqueDay).NotSleep_inds(1,o) = NaN;
            data.(uniqueDay).NREM_inds(1,o) = NaN;
            data.(uniqueDay).REM_inds(1,o) = NaN;
            data.(uniqueDay).TimePad_inds(1,o) = 1;
        end
    end
end
% Generate figure
hypFig = figure;
sgtitle([animalID ' Hyponogram'])
timeConv = 60*(60/binTime);
for p = 1:size(data.uniqueDays,1)
    uniqueDay = data.uniqueDays{p,1};
    % Create new subplot for each day
    ax(p) = subplot(size(data.uniqueDays,1),1,p);
    b1 = bar((1:length(data.(uniqueDay).NotSleep_inds))/timeConv,data.(uniqueDay).NotSleep_inds,'k','BarWidth',1);
    hold on
    b2 = bar((1:length(data.(uniqueDay).NREM_inds))/timeConv,data.(uniqueDay).NREM_inds,'b','BarWidth',1);
    b3 = bar((1:length(data.(uniqueDay).REM_inds))/timeConv,data.(uniqueDay).REM_inds,'r','BarWidth',1);
    if p == 1
        legend([b1,b2,b3],'Not Sleep','NREM Sleep','REM Sleep')
    elseif p == size(data.uniqueDays,1)
        xlabel('Time (hrs)')
    end 
    ylabel(data.uniqueDays{p,1})
    set(gca,'YTickLabel',[]);
    set(gca,'Ticklength',[0,0])
    set(gca,'box','off')
end
linkaxes(ax(1:p),'x')

%% Save the figure to directory.
[pathstr,~,~] = fileparts(cd);
dirpath = [pathstr '/Combined Imaging/Figures/Hypnogram/'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(hypFig,[dirpath animalID '_Hypnogram']);
close(hypFig)

end
