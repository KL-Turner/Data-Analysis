function [] = CreateAnimalHypnograms_IOS()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Loads sleep scoring results structure and creates a hypnogram on the scored data across the days
%________________________________________________________________________________________________________________________

%%  Identify the sleep scoring results structure in the current folder
sleepScoringResultsFileStruct = dir('*_SVM_SleepScoringResults.mat'); 
sleepScoringResultsFile = {sleepScoringResultsFileStruct.name}';
sleepScoringResultsFileID = char(sleepScoringResultsFile);
load(sleepScoringResultsFileID)
% Identify the unique file IDs, unique imaging days, and scoring labels from the file list
allScoringLabels = SVMResults.labels;
allFileIDs = SVMResults.fileIDs;
uniqueFileIDs = unique(allFileIDs);
for a = 1:size(uniqueFileIDs,1)
    [animalID,allFileDates{a,1},~] = GetFileInfo_IOS(uniqueFileIDs{a,1}); %#ok<AGROW>
end
uniqueDays = unique(allFileDates);
%  
for b = 1:size(uniqueDays,1)
    dayLengths{b,1} = find(~cellfun('isempty',strfind(allFileIDs,uniqueDays{b,1}))); %#ok<AGROW>
end
%
for c = 1:size(dayLengths,1)
    dayInds = dayLengths{c,1};
    dayScoreLabels{c,1} = allScoringLabels(dayInds);
    dayScoreFileIDs{c,1} = allFileIDs(dayInds);
end
% 
for d = 1:size(uniqueDays,1)
    day = ConvertDate_IOS(uniqueDays{d,1});
    uniqueDayFileIDs = unique(dayScoreFileIDs{d,1});
    for e = 1:size(uniqueDayFileIDs,1)
        if e == 1
            data.(day){e,1} = dayScoreLabels{d,1}(1:fileBinLength);
        else
            data.(day){e,1} = dayScoreLabels{d,1}((e - 1)*fileBinLength + 1:e*fileBinLength);
        end
    end
end

trialDuration = 15;   % minutes
binTime = 5;   % seconds
fileBinLength = (trialDuration*60)/binTime;
for d = 1:size(dayScoreLabels,1)
    indDayLabels = dayScoreLabels{d,1};
    day = ConvertDate_IOS(uniqueDays{d,1});
    uniqueDayFileIDs = unique(dayScoreFileIDs{d,1});
    for e = 2:size(uniqueDayFileIDs,1)
        leadFileID = uniqueDayFileIDs{e-1,1};
        lagFileID = uniqueDayFileIDs{e,1};
        [~,~,leadFileInfo] = GetFileInfo_IOS(leadFileID);
        [~,~,lagFileInfo] = GetFileInfo_IOS(lagFileID);
        leadFileStr = ConvertDateTime_IOS(leadFileInfo);
        lagFileStr = ConvertDateTime_IOS(lagFileInfo);
        leadFileTime = datevec(leadFileStr);
        lagFileTime = datevec(lagFileStr);
        timeDifference = etime(lagFileTime,leadFileTime) - (trialDuration*60);   % seconds
        timePadBins.(day){e-1,1} = cell(floor(timeDifference/binTime),1);
        timePadBins.(day){e-1,1}(:) = {'Time Pad'};
    end
end
%
for f = 1:size(uniqueDays,1)
    day = ConvertDate_IOS(uniqueDays{f,1});
    for g = 1:size(data.(day),1) - 1
        data.(day){g,1} = vertcat(data.(day){g,1},timePadBins.(day){g,1});
    end
end
%
for x = 1:size(uniqueDays,1)
    day = ConvertDate_IOS(uniqueDays{x,1});
    catData.(day) = [];
    for y = 1:size(data.(day),1)
        catData.(day) = vertcat(catData.(day),data.(day){y,1});
    end
end
% 
for x = 1:size(uniqueDays,1)
    day = ConvertDate_IOS(uniqueDays{x,1});
    allCatData = [];
    for y = 1:size(catData,1)
        allCatData = vertcat(allCatData,catData.(day));
    end
end

end
