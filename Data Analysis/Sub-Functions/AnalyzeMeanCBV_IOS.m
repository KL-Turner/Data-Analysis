function [AnalysisResults] = AnalyzeMeanCBV_IOS(params,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID)

% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID)

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

% identify animal's ID and pull important infortmat
fileBreaks = strfind(restDataFileID,'_');
animalID = restDataFileID(1:fileBreaks(1)-1);
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
samplingRate = RestData.CBV.adjLH.CBVCamSamplingRate;

WhiskCriteria.Fieldname = {'duration','puffDistance'};
WhiskCriteria.Comparison = {'gt','gt'};
WhiskCriteria.Value = {5,5};

RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};

RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};

WhiskPuffCriteria.Fieldname = {'puffDistance'};
WhiskPuffCriteria.Comparison = {'gt'};
WhiskPuffCriteria.Value = {5};

%% Analyze mean CBV during periods of rest
disp('AnalzeMeanCBV: average hemodynamics during rest')
[restLogical] = FilterEvents_IOS(RestData.CBV_HbT.adjLH,RestCriteria);
[puffLogical] = FilterEvents_IOS(RestData.CBV_HbT.adjLH,RestPuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFiles = RestData.CBV_HbT.adjLH.fileIDs(combRestLogical,:);
LH_RestingData = RestData.CBV_HbT.adjLH.data(combRestLogical,:);
RH_RestingData = RestData.CBV_HbT.adjRH.data(combRestLogical,:);

% identify the unique days and the unique number of files from the list of unstim resting events
restUniqueDays = GetUniqueDays_IOS(restFiles);
restUniqueFiles = unique(restFiles);
restNumberOfFiles = length(unique(restFiles));

% decimate the file list to only include those files that occur within the desired number of target minutes
clear restFiltLogical
for c = 1:length(restUniqueDays)
    restDay = restUniqueDays(c);
    d = 1;
    for e = 1:restNumberOfFiles
        restFile = restUniqueFiles(e);
        restFileID = restFile{1}(1:6);
        if strcmp(restDay,restFileID) && sum(strcmp(restFile,manualFileIDs)) == 1
            restFiltLogical{c,1}(e,1) = 1; %#ok<*AGROW>
            d = d + 1;
        else
            restFiltLogical{c,1}(e,1) = 0;
        end
    end
end
restFinalLogical = any(sum(cell2mat(restFiltLogical'),2),2);

% extract unstim the resting events that correspond to the acceptable file list and the acceptable resting criteria
clear restFileFilter
filtRestFiles = restUniqueFiles(restFinalLogical,:);
for f = 1:length(restFiles)
    restLogic = strcmp(restFiles{f},filtRestFiles);
    restLogicSum = sum(restLogic);
    if restLogicSum == 1
        restFileFilter(f,1) = 1;
    else
        restFileFilter(f,1) = 0;
    end
end
restFinalFileFilter = logical(restFileFilter);
LH_finalRestData = LH_RestingData(restFinalFileFilter,:);
RH_finalRestData = RH_RestingData(restFinalFileFilter,:);

% only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
% original epoch create so we can add a sample of two back to the end for those just under 10 seconds
% lowpass filter and detrend each segment
[B, A] = butter(3,1/(samplingRate/2),'low');
clear LH_ProcRestData
clear RH_ProcRestData
for g = 1:length(LH_finalRestData)
    LH_ProcRestData{g,1} = filtfilt(B,A,LH_finalRestData{g,1});
    RH_ProcRestData{g,1} = filtfilt(B,A,RH_finalRestData{g,1});
end

% analyze correlation coefficient between resting epochs
for n = 1:length(LH_ProcRestData)
    LH_restCBVMean(n,1) = mean(LH_ProcRestData{n,1});
    RH_restCBVMean(n,1) = mean(RH_ProcRestData{n,1});
end

% save results
AnalysisResults.MeanCBV.Rest.CBV_HbT.adjLH = LH_restCBVMean;
AnalysisResults.MeanCBV.Rest.CBV_HbT.adjRH = RH_restCBVMean;

%% Analyze mean CBV during periods of extended whisking
disp('AnalzeMeanCBV: average hemodynamics during whisking')
[whiskLogical] = FilterEvents_IOS(EventData.CBV_HbT.adjLH.whisk,WhiskCriteria);
[puffLogical] = FilterEvents_IOS(EventData.CBV_HbT.adjLH.whisk,WhiskPuffCriteria);
combWhiskLogical = logical(whiskLogical.*puffLogical);
whiskFiles = EventData.CBV_HbT.adjLH.whisk.fileIDs(combWhiskLogical,:);
LH_whiskData = EventData.CBV_HbT.adjLH.whisk.data(combWhiskLogical,:);
RH_whiskData = EventData.CBV_HbT.adjRH.whisk.data(combWhiskLogical,:);

% identify the unique days and the unique number of files from the list of unstim resting events
whiskUniqueDays = GetUniqueDays_IOS(whiskFiles);
whiskUniqueFiles = unique(whiskFiles);
whiskNumberOfFiles = length(unique(whiskFiles));

% decimate the file list to only include those files that occur within the desired number of target minutes
clear whiskFiltLogical
for c = 1:length(whiskUniqueDays)
    whiskDay = whiskUniqueDays(c);
    d = 1;
    for e = 1:whiskNumberOfFiles
        whiskFile = whiskUniqueFiles(e);
        whiskFileID = whiskFile{1}(1:6);
        if strcmp(whiskDay,whiskFileID) && sum(strcmp(whiskFile,manualFileIDs)) == 1
            whiskFiltLogical{c,1}(e,1) = 1; %#ok<*AGROW>
            d = d + 1;
        else
            whiskFiltLogical{c,1}(e,1) = 0;
        end
    end
end
whiskFinalLogical = any(sum(cell2mat(whiskFiltLogical'),2),2);

% extract unstim the resting events that correspond to the acceptable file list and the acceptable resting criteria
clear whiskFileFilter
filtWhiskFiles = whiskUniqueFiles(whiskFinalLogical,:);
for f = 1:length(whiskFiles)
    whiskLogic = strcmp(whiskFiles{f},filtWhiskFiles);
    whiskLogicSum = sum(whiskLogic);
    if whiskLogicSum == 1
        whiskFileFilter(f,1) = 1;
    else
        whiskFileFilter(f,1) = 0;
    end
end
whiskFinalFileFilter = logical(whiskFileFilter);
LH_finalWhiskData = LH_whiskData(whiskFinalFileFilter,:);
RH_finalWhiskData = RH_whiskData(whiskFinalFileFilter,:);

% only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
% original epoch create so we can add a sample of two back to the end for those just under 10 seconds
% lowpass filter and detrend each segment
[B, A] = butter(3,1/(samplingRate/2),'low');
clear LH_ProcWhiskData
clear RH_ProcWhiskData
for g = 1:size(LH_finalWhiskData,1)
    LH_ProcWhiskData(g,:) = filtfilt(B,A,LH_finalWhiskData(g,:));
    RH_ProcWhiskData(g,:) = filtfilt(B,A,RH_finalWhiskData(g,:));
end

% analyze correlation coefficient between resting epochs
for n = 1:size(LH_ProcWhiskData,1)
    LH_whiskCBVMean{n,1} = mean(LH_ProcWhiskData(n,samplingRate*2:end),2);
    RH_whiskCBVMean{n,1} = mean(RH_ProcWhiskData(n,samplingRate*2:end),2);
end

% save results
AnalysisResults.MeanCBV.Whisk.CBV_HbT.adjLH = cell2mat(LH_whiskCBVMean);
AnalysisResults.MeanCBV.Whisk.CBV_HbT.adjRH = cell2mat(RH_whiskCBVMean);

%% Analyze mean CBV during periods of NREM sleep
% pull data from SleepData.mat structure
disp('AnalzeMeanCBV: average hemodynamics during NREM')
LH_nremData = SleepData.NREM.data.CBV_HbT.LH;
RH_nremData = SleepData.NREM.data.CBV_HbT.RH;

% analyze correlation coefficient between NREM epochs
for n = 1:length(LH_nremData)
    LH_nremCBVMean(n,1) = mean(LH_nremData{n,1});
    RH_nremCBVMean(n,1) = mean(RH_nremData{n,1});
end

% save results
AnalysisResults.MeanCBV.NREM.CBV_HbT.adjLH = LH_nremCBVMean;
AnalysisResults.MeanCBV.NREM.CBV_HbT.adjRH = RH_nremCBVMean;

%% Analyze mean CBV during periods of REM sleep
% pull data from SleepData.mat structure
disp('AnalzeMeanCBV: average hemodynamics during REM')
LH_remData = SleepData.REM.data.CBV_HbT.LH;
RH_remData = SleepData.REM.data.CBV_HbT.RH;

% analyze correlation coefficient between NREM epochs
for n = 1:length(LH_remData)
    LH_remCBVMean(n,1) = mean(LH_remData{n,1});
    RH_remCBVMean(n,1) = mean(RH_remData{n,1});
end

% save results
AnalysisResults.MeanCBV.REM.CBV_HbT.adjLH = LH_remCBVMean;
AnalysisResults.MeanCBV.REM.CBV_HbT.adjRH = RH_remCBVMean;

%% save results strucure
save([animalID '_AnalysisResults.mat'],'AnalysisResults');

end
