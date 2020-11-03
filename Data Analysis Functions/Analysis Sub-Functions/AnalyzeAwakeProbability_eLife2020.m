function [AnalysisResults] = AnalyzeAwakeProbability_eLife2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the arousal-state probability of trial duration and resting events (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
modelType = 'Forest';
%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '\' animalID '\Bilateral Imaging\'];
    cd(dataLocation)
    % find and load RestData.mat struct
    restDataFileStruct = dir('*_RestData.mat');
    restDataFile = {restDataFileStruct.name}';
    restDataFileID = char(restDataFile);
    load(restDataFileID)
    % scoring results
    modelScoringResults = [modelType '_ScoringResults.mat'];
    load(modelScoringResults)
    %% determine probabilty of a single "resting event" being awake or asleep based on time
    % criteria for resting
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {0};
    RestPuffCriteria.Fieldname = {'puffDistances'};
    RestPuffCriteria.Comparison = {'gt'};
    RestPuffCriteria.Value = {5};
    % pull data from RestData.mat structure
    [restLogical] = FilterEvents_IOS_eLife2020(RestData.CBV_HbT.adjLH,RestCriteria);
    [puffLogical] = FilterEvents_IOS_eLife2020(RestData.CBV_HbT.adjLH,RestPuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.CBV_HbT.adjLH.fileIDs(combRestLogical,:);
    restEventTimes = RestData.CBV_HbT.adjLH.eventTimes(combRestLogical,:);
    restDurations = RestData.CBV_HbT.adjLH.durations(combRestLogical,:);
    % 5 second events
    bins = {'five','ten','fifteen','twenty','twentyfive','thirty','thirtyfive','forty','fortyfive','fifty','fiftyfive','sixty','sixtyplus'};
    for a = 1:length(bins)
        data.(bins{1,a}) = [];
    end
    % group data based on resting duration
    a5 = 1; a10 = 1; a15 = 1; a20 = 1; a25 = 1; a30 = 1; a35 = 1; a40 = 1; a45 = 1; a50 = 1; a55 = 1; a60 = 1; a61 = 1;
    for a = 1:length(restDurations)
        duration = restDurations(a,1);
        if duration <= 5
            data.five.fileIDs{a5,1} = restFileIDs{a,1};
            data.five.durations(a5,1) = restDurations(a,1);
            data.five.eventTimes(a5,1) = restEventTimes(a,1);
            a5 = a5 + 1;
        elseif duration > 5 && duration <= 10
            data.ten.fileIDs{a10,1} = restFileIDs{a,1};
            data.ten.durations(a10,1) = restDurations(a,1);
            data.ten.eventTimes(a10,1) = restEventTimes(a,1);
            a10 = a10 + 1;
        elseif duration > 10 && duration <= 15
            data.fifteen.fileIDs{a15,1} = restFileIDs{a,1};
            data.fifteen.durations(a15,1) = restDurations(a,1);
            data.fifteen.eventTimes(a15,1) = restEventTimes(a,1);
            a15 = a15 + 1;
        elseif duration > 15 && duration <= 20
            data.twenty.fileIDs{a20,1} = restFileIDs{a,1};
            data.twenty.durations(a20,1) = restDurations(a,1);
            data.twenty.eventTimes(a20,1) = restEventTimes(a,1);
            a20 = a20 + 1;
        elseif duration > 20 && duration <= 25
            data.twentyfive.fileIDs{a25,1} = restFileIDs{a,1};
            data.twentyfive.durations(a25,1) = restDurations(a,1);
            data.twentyfive.eventTimes(a25,1) = restEventTimes(a,1);
            a25 = a25 + 1;
        elseif duration > 25 && duration <= 30
            data.thirty.fileIDs{a30,1} = restFileIDs{a,1};
            data.thirty.durations(a30,1) = restDurations(a,1);
            data.thirty.eventTimes(a30,1) = restEventTimes(a,1);
            a30 = a30 + 1;
        elseif duration > 30 && duration <= 35
            data.thirtyfive.fileIDs{a35,1} = restFileIDs{a,1};
            data.thirtyfive.durations(a35,1) = restDurations(a,1);
            data.thirtyfive.eventTimes(a35,1) = restEventTimes(a,1);
            a35 = a35 + 1;
        elseif duration > 35 && duration <= 40
            data.forty.fileIDs{a40,1} = restFileIDs{a,1};
            data.forty.durations(a40,1) = restDurations(a,1);
            data.forty.eventTimes(a40,1) = restEventTimes(a,1);
            a40 = a40 + 1;
        elseif duration > 40 && duration <= 45
            data.fortyfive.fileIDs{a45,1} = restFileIDs{a,1};
            data.fortyfive.durations(a45,1) = restDurations(a,1);
            data.fortyfive.eventTimes(a45,1) = restEventTimes(a,1);
            a45 = a45 + 1;
        elseif duration > 45 && duration <= 50
            data.fifty.fileIDs{a50,1} = restFileIDs{a,1};
            data.fifty.durations(a50,1) = restDurations(a,1);
            data.fifty.eventTimes(a50,1) = restEventTimes(a,1);
            a50 = a50 + 1;
        elseif duration > 50 && duration <= 55
            data.fiftyfive.fileIDs{a55,1} = restFileIDs{a,1};
            data.fiftyfive.durations(a55,1) = restDurations(a,1);
            data.fiftyfive.eventTimes(a55,1) = restEventTimes(a,1);
            a55 = a55 + 1;
        elseif duration > 55 && duration <= 60
            data.sixty.fileIDs{a60,1} = restFileIDs{a,1};
            data.sixty.durations(a60,1) = restDurations(a,1);
            data.sixty.eventTimes(a60,1) = restEventTimes(a,1);
            a60 = a60 + 1;
        elseif duration > 60
            data.sixtyplus.fileIDs{a61,1} = restFileIDs{a,1};
            data.sixtyplus.durations(a61,1) = restDurations(a,1);
            data.sixtyplus.eventTimes(a61,1) = restEventTimes(a,1);
            a61 = a61 + 1;
        end
    end
    % go through each event and determine the probabilty of awake
    for b = 1:length(bins)
        bin = bins{1,b};
        for c = 1:length(data.(bin).fileIDs)
            fileID = data.(bin).fileIDs{c,1};
            eventTime = floor(data.(bin).eventTimes(c,1));
            duration = floor(data.(bin).durations(c,1));
            indFileScores = [];
            for e = 1:length(ScoringResults.fileIDs)
                if strcmp(fileID,ScoringResults.fileIDs{e,1})
                    indFileScores = ScoringResults.labels{e,1};
                end
            end
            % determine starting place of eventTime
            startRemTime = rem(eventTime,5);
            binStartTime = eventTime - startRemTime;
            binStartIndex = (binStartTime/5) + 1;
            % determine ending place of eventTime
            binEndIndex = floor(duration/5);
            % extract data from index
            indexScores = indFileScores(binStartIndex:binStartIndex + binEndIndex);
            indexLogical = strcmp(indexScores,'Not Sleep');
            if sum(indexLogical) == length(indexLogical)
                data.(bin).awakeLogical(c,1) = 1;
            else
                data.(bin).awakeLogical(c,1) = 0;
            end
        end
    end
    % save results
    for f = 1:length(bins)
        AnalysisResults.(animalID).SleepProbability.(bins{1,f}).awakeLogical = data.(bins{1,f}).awakeLogical;
    end
    %% analyze trial hypogram and awake probability based on trial time
    % identify the unique file IDs, unique imaging days, and scoring labels from the file list
    allScoringLabels = ScoringResults.alllabels;
    allFileIDs = ScoringResults.allfileIDs;
    uniqueFileIDs = unique(allFileIDs);
    for a = 1:size(uniqueFileIDs,1)
        [animalID,allFileDates{a,1},~] = GetFileInfo_IOS_eLife2020(uniqueFileIDs{a,1}); %#ok<*AGROW>
        allFileDays{a,1} = ConvertDate_IOS_eLife2020(allFileDates{a,1});
    end
    data.uniqueDates = unique(allFileDates);
    data.uniqueDays = unique(allFileDays);
    % determine how many 5 second bins there are for each individual day
    for b = 1:size(data.uniqueDays,1)
        data.dayBinLengths{b,1} = find(~cellfun('isempty',strfind(allFileIDs,data.uniqueDates{b,1})));
    end
    % extract the file IDs and scoring labels that correspond to each day's indeces
    for c = 1:size(data.dayBinLengths,1)
        dayInds = data.dayBinLengths{c,1};
        data.dayScoreLabels{c,1} = allScoringLabels(dayInds);
        data.dayScoreFileIDs{c,1} = allFileIDs(dayInds);
    end
    trialDuration = 15;   % minutes
    binTime = 5;   % seconds
    fileBinLength = (trialDuration*60)/binTime;
    % further break down each day's scores into the scores for each individual file
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
    % calculate the time difference between every file to append padding 'Time Pad' to the end of the leading file's score labels
    for f = 1:size(data.uniqueDays,1)
        uniqueDay = data.uniqueDays{f,1};
        uniqueDayFileIDs = unique(data.dayScoreFileIDs{f,1});
        % start with file 2 to focus on the differences between each file
        for g = 2:size(uniqueDayFileIDs,1)
            leadFileID = uniqueDayFileIDs{g - 1,1};
            lagFileID = uniqueDayFileIDs{g,1};
            [~,~,leadFileInfo] = GetFileInfo_IOS_eLife2020(leadFileID);
            [~,~,lagFileInfo] = GetFileInfo_IOS_eLife2020(lagFileID);
            leadFileStr = ConvertDateTime_IOS_eLife2020(leadFileInfo);
            lagFileStr = ConvertDateTime_IOS_eLife2020(lagFileInfo);
            leadFileTime = datevec(leadFileStr);
            lagFileTime = datevec(lagFileStr);
            timeDifference = etime(lagFileTime,leadFileTime) - (trialDuration*60);   % seconds
            timePadBins.(uniqueDay){g - 1,1} = cell(floor(timeDifference/binTime),1);
            timePadBins.(uniqueDay){g - 1,1}(:) = {'Time Pad'};
        end
    end
    % apply the time padding to the end of each file
    for h = 1:size(data.uniqueDays,1)
        uniqueDay = data.uniqueDays{h,1};
        for j = 1:size(data.(uniqueDay).indFileData,1) - 1
            data.(uniqueDay).indFileData{j,1} = vertcat(data.(uniqueDay).indFileData{j,1},timePadBins.(uniqueDay){j,1});
        end
    end
    % concatendate the data for each day now that the padding is added at the end of each file
    for k = 1:size(data.uniqueDays,1)
        uniqueDay = data.uniqueDays{k,1};
        data.(uniqueDay).catData = [];
        for m = 1:size(data.(uniqueDay).indFileData,1)
            data.(uniqueDay).catData = vertcat(data.(uniqueDay).catData,data.(uniqueDay).indFileData{m,1});
        end
    end
    % prepare indeces for each behavior
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
                data.(uniqueDay).AwakeProb_inds(1,o) = 1;
                data.(uniqueDay).NREMProb_inds(1,o) = 0;
                data.(uniqueDay).REMProb_inds(1,o) = 0;
            elseif strcmp(data.(uniqueDay).catData{o,1},'NREM Sleep') == true
                data.(uniqueDay).NotSleep_inds(1,o) = NaN;
                data.(uniqueDay).NREM_inds(1,o) = 1;
                data.(uniqueDay).REM_inds(1,o) = NaN;
                data.(uniqueDay).TimePad_inds(1,o) = NaN;
                data.(uniqueDay).AwakeProb_inds(1,o) = 0;
                data.(uniqueDay).NREMProb_inds(1,o) = 1;
                data.(uniqueDay).REMProb_inds(1,o) = 0;
            elseif strcmp(data.(uniqueDay).catData{o,1},'REM Sleep') == true
                data.(uniqueDay).NotSleep_inds(1,o) = NaN;
                data.(uniqueDay).NREM_inds(1,o) = NaN;
                data.(uniqueDay).REM_inds(1,o) = 1;
                data.(uniqueDay).TimePad_inds(1,o) = NaN;
                data.(uniqueDay).AwakeProb_inds(1,o) = 0;
                data.(uniqueDay).NREMProb_inds(1,o) = 0;
                data.(uniqueDay).REMProb_inds(1,o) = 1;
            elseif strcmp(data.(uniqueDay).catData{o,1},'Time Pad') == true
                data.(uniqueDay).NotSleep_inds(1,o) = NaN;
                data.(uniqueDay).NREM_inds(1,o) = NaN;
                data.(uniqueDay).REM_inds(1,o) = NaN;
                data.(uniqueDay).TimePad_inds(1,o) = 1;
                data.(uniqueDay).AwakeProb_inds(1,o) = NaN;
                data.(uniqueDay).NREMProb_inds(1,o) = NaN;
                data.(uniqueDay).REMProb_inds(1,o) = NaN;
            end
        end
    end
    % save figures if desired
    if strcmp(saveFigs,'y') == true
        hypFig = figure;
        sgtitle([animalID ' Hyponogram'])
        timeConv = 60*(60/binTime);
        for p = 1:size(data.uniqueDays,1)
            uniqueDay = data.uniqueDays{p,1};
            % create new subplot for each day
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
        % save figure
        [pathstr,~,~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Hyponogram/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(hypFig,[dirpath animalID '_Hypnogram']);
        close(hypFig)
    end
    % save results
    for a = 1:size(data.uniqueDays,1)
        AnalysisResults.(animalID).SleepProbability.Hypnogram.(data.uniqueDays{a,1}) = data.(data.uniqueDays{a,1});
    end
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
