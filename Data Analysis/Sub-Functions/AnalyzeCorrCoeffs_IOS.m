function [AnalysisResults] = AnalyzeCorrCoeffs_IOS(dataTypes,params,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Pearson's correlation coefficient between hemodynamic or neural envelope signals
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
fileBreaks = strfind(restDataFileID, '_');
animalID = restDataFileID(1:fileBreaks(1)-1);
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
samplingRate = RestData.CBV_HbT.adjLH.CBVCamSamplingRate;

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

for a = 1:length(dataTypes)
    dataType = dataTypes{1,a};
    %% Analyze Pearson's correlation coefficient during periods of rest
    % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
    if strcmp(dataType,'CBV_HbT') == true
        [restLogical] = FilterEvents_IOS(RestData.(dataType).adjLH,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.(dataType).adjLH,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        unstimRestFiles = RestData.(dataType).adjLH.fileIDs(combRestLogical,:);
        LH_unstimRestingData = RestData.(dataType).adjLH.data(combRestLogical,:);
        RH_unstimRestingData = RestData.(dataType).adjRH.data(combRestLogical,:);
    else
        [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        unstimRestFiles = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
        LH_unstimRestingData = RestData.cortical_LH.(dataType).NormData(combRestLogical,:);
        RH_unstimRestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical,:);
    end
    
    % identify the unique days and the unique number of files from the list of unstim resting events
    restUniqueDays = GetUniqueDays_IOS(unstimRestFiles);
    restUniqueFiles = unique(unstimRestFiles);
    restNumberOfFiles = length(unique(unstimRestFiles));
    
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
    for f = 1:length(unstimRestFiles)
        restLogic = strcmp(unstimRestFiles{f},filtRestFiles);
        restLogicSum = sum(restLogic);
        if restLogicSum == 1
            restFileFilter(f,1) = 1;
        else
            restFileFilter(f,1) = 0;
        end
    end
    restFinalFileFilter = logical(restFileFilter);
    LH_finalRestData = LH_unstimRestingData(restFinalFileFilter,:);
    RH_finalRestData = RH_unstimRestingData(restFinalFileFilter,:);
    
    % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    % lowpass filter and detrend each segment
    [B, A] = butter(3,1/(samplingRate/2),'low');
    clear LH_ProcRestData
    clear RH_ProcRestData
    for g = 1:length(LH_finalRestData)
        LH_ProcRestData{g,1} = detrend(filtfilt(B,A,LH_finalRestData{g,1}),'constant');
        RH_ProcRestData{g,1} = detrend(filtfilt(B,A,RH_finalRestData{g,1}),'constant');
    end
    
    % analyze correlation coefficient between resting epochs
    for n = 1:length(LH_ProcRestData)
        rest_CC = corrcoef(LH_ProcRestData{n,1},RH_ProcRestData{n,1});
        rest_R(n,1) = rest_CC(2,1);
    end
    meanRest_R = mean(rest_R);
    stdRest_R = std(rest_R,0,1);
    
    disp([dataType ' Rest Corr Coeff: ' num2str(meanRest_R) ' +/- ' num2str(stdRest_R)]); disp(' ')
    % save results
    AnalysisResults.CorrCoeff.Rest.(dataType).R = rest_R;
    AnalysisResults.CorrCoeff.Rest.(dataType).meanR = meanRest_R;
    AnalysisResults.CorrCoeff.Rest.(dataType).stdR = stdRest_R;
    
    
    %% Analyze Pearson's correlation coefficient during periods of extended whisking
    % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
    if strcmp(dataType,'CBV_HbT') == true
        [whiskLogical] = FilterEvents_IOS(EventData.(dataType).adjLH.whisk,WhiskCriteria);
        [puffLogical] = FilterEvents_IOS(EventData.(dataType).adjLH.whisk,WhiskPuffCriteria);
        combWhiskLogical = logical(whiskLogical.*puffLogical);
        whiskFiles = EventData.(dataType).adjLH.whisk.fileIDs(combWhiskLogical,:);
        LH_whiskData = EventData.(dataType).adjLH.whisk.data(combWhiskLogical,:);
        RH_whiskData = EventData.(dataType).adjRH.whisk.data(combWhiskLogical,:);
    else
        [whiskLogical] = FilterEvents_IOS(EventData.cortical_LH.(dataType).whisk,WhiskCriteria);
        [puffLogical] = FilterEvents_IOS(EventData.cortical_LH.(dataType).whisk,WhiskPuffCriteria);
        combWhiskLogical = logical(whiskLogical.*puffLogical);
        whiskFiles = EventData.cortical_LH.(dataType).whisk.fileIDs(combWhiskLogical,:);
        LH_whiskData = EventData.cortical_LH.(dataType).whisk.NormData(combWhiskLogical,:);
        RH_whiskData = EventData.cortical_RH.(dataType).whisk.NormData(combWhiskLogical,:);
    end
    
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
        LH_ProcWhiskData(g,:) = detrend(filtfilt(B,A,LH_finalWhiskData(g,:)),'constant');
        RH_ProcWhiskData(g,:) = detrend(filtfilt(B,A,RH_finalWhiskData(g,:)),'constant');
    end
    
    % analyze correlation coefficient between resting epochs
    for n = 1:size(LH_ProcWhiskData,1)
        whisk_CC = corrcoef(LH_ProcWhiskData(n,samplingRate*2:end),RH_ProcWhiskData(n,samplingRate*2:end));
        whisk_R(n,1) = whisk_CC(2,1);
    end
    meanWhisk_R = mean(whisk_R);
    stdWhisk_R = std(whisk_R,0,1);
    
    disp([dataType ' Whisk Corr Coeff: ' num2str(meanWhisk_R) ' +/- ' num2str(stdWhisk_R)]); disp(' ')
    % save results
    AnalysisResults.CorrCoeff.Whisk.(dataType).R = whisk_R;
    AnalysisResults.CorrCoeff.Whisk.(dataType).meanR = meanWhisk_R;
    AnalysisResults.CorrCoeff.Whisk.(dataType).stdR = stdWhisk_R;
    
    %% Analyze Pearson's correlation coefficient during periods of NREM sleep
    % pull data from SleepData.mat structure
    if strcmp(dataType,'CBV_HbT') == true
        LH_nremData = SleepData.NREM.data.(dataType).LH;
        RH_nremData = SleepData.NREM.data.(dataType).RH;
    else
        LH_nremData = SleepData.NREM.data.cortical_LH.(dataType);
        RH_nremData = SleepData.NREM.data.cortical_RH.(dataType);
    end
    
    % detrend - data is already lowpass filtered
    for j = 1:length(LH_nremData)
        LH_nremData{j,1} = detrend(LH_nremData{j,1},'constant');
        RH_nremData{j,1} = detrend(RH_nremData{j,1},'constant');
    end
    
    % analyze correlation coefficient between NREM epochs
    for n = 1:length(LH_nremData)
        nrem_CC = corrcoef(LH_nremData{n,1},RH_nremData{n,1});
        nrem_R(n,1) = nrem_CC(2,1);
    end
    meanNREM_R = mean(nrem_R);
    stdNREM_R = std(nrem_R,0,1);
    
    disp([dataType ' NREM Corr Coeff: ' num2str(meanNREM_R) ' +/- ' num2str(stdNREM_R)]); disp(' ')
    % save results
    AnalysisResults.CorrCoeff.NREM.(dataType).R = nrem_R;
    AnalysisResults.CorrCoeff.NREM.(dataType).meanR = meanNREM_R;
    AnalysisResults.CorrCoeff.NREM.(dataType).stdR = stdNREM_R;
    
    %% Analyze Pearson's correlation coefficient during periods of REM sleep
    % pull data from SleepData.mat structure
    if strcmp(dataType,'CBV_HbT') == true
        LH_remData = SleepData.REM.data.(dataType).LH;
        RH_remData = SleepData.REM.data.(dataType).RH;
    else
        LH_remData = SleepData.REM.data.cortical_LH.(dataType);
        RH_remData = SleepData.REM.data.cortical_RH.(dataType);
    end
    
    % detrend - data is already lowpass filtered
    for m = 1:length(LH_remData)
        LH_remData{m,1} = detrend(LH_remData{m,1},'constant');
        RH_remData{m,1} = detrend(RH_remData{m,1},'constant');
    end
    
    % analyze correlation coefficient between NREM epochs
    for n = 1:length(LH_remData)
        rem_CC = corrcoef(LH_remData{n,1},RH_remData{n,1});
        rem_R(n,1) = rem_CC(2,1);
    end
    meanREM_R = mean(rem_R);
    stdREM_R = std(rem_R,0,1);
    
    disp([dataType ' REM Corr Coeff: ' num2str(meanREM_R) ' +/- ' num2str(stdREM_R)]); disp(' ')
    % save results
    AnalysisResults.CorrCoeff.REM.(dataType).R = rem_R;
    AnalysisResults.CorrCoeff.REM.(dataType).meanR = meanREM_R;
    AnalysisResults.CorrCoeff.REM.(dataType).stdR = stdREM_R;
end

%% save results strucure
save([animalID '_AnalysisResults.mat'],'AnalysisResults');

end
