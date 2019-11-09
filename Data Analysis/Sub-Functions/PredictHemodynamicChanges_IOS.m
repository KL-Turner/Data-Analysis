function [AnalysisResults] = PredictHemodynamicChanges_IOS(params,fileSet,CBVdataType,hemDataType,neuralBand,baselineType,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

% list of unstim Procdata.mat files
procDataFileStruct = dir('*_Procdata.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

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
trialDuration_min = RestData.CBV.adjLH.trialDuration_sec/60;   % min
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
fileTarget = params.targetMinutes/trialDuration_min;
filterSets = {'manualSelection','setDuration','entireDuration'};
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

bestKernel = AnalysisResults.HRFs.(fileSet).(CBVdataType).(hemDataType).(neuralBand).bestKernel;

%% Analyze mean CBV during periods of rest
for b = 1:length(filterSets)
    filterSet = filterSets{1,b};
    % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
        [restLogical] = FilterEvents_IOS(RestData.(CBVdataType).(['adj' hemDataType]),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.(CBVdataType).(['adj' hemDataType]),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFiles = RestData.(CBVdataType).(['adj' hemDataType]).fileIDs(combRestLogical,:);
        if strcmp(CBVdataType,'CBV') == true
            restingCBVData = RestData.(CBVdataType).(['adj' hemDataType]).NormData(combRestLogical,:);
        elseif strcmp(CBVdataType,'CBV_HbT') == true
            restingCBVData = RestData.(CBVdataType).(['adj' hemDataType]).data(combRestLogical,:);
        end
    restingNeuralData = RestData.(['cortical_' hemDataType]).(neuralBand).NormData(combRestLogical,:);
    
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
            if strcmp(filterSet,'manualSelection') == true
                if strcmp(restDay,restFileID) && sum(strcmp(restFile,manualFileIDs)) == 1
                    restFiltLogical{c,1}(e,1) = 1; %#ok<*AGROW>
                    d = d + 1;
                else
                    restFiltLogical{c,1}(e,1) = 0;
                end
            elseif strcmp(filterSet,'setDuration') == true
                if strcmp(restDay,restFileID) && d <= fileTarget
                    restFiltLogical{c,1}(e,1) = 1;
                    d = d + 1;
                else
                    restFiltLogical{c,1}(e,1) = 0;
                end
            elseif strcmp(filterSet,'entireDuration') == true
                if strcmp(restDay,restFileID)
                    restFiltLogical{c,1}(e,1) = 1;
                    d = d + 1;
                else
                    restFiltLogical{c,1}(e,1) = 0;
                end
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
    finalRestCBVData = restingCBVData(restFinalFileFilter,:);
    finalRestNeuralData = restingNeuralData(restFinalFileFilter,:);
    
    % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    % lowpass filter and detrend each segment
    [B, A] = butter(3,1/(samplingRate/2),'low');
    clear LH_ProcRestData
    clear RH_ProcRestData
    for g = 1:length(finalRestCBVData)
        procRestCBVData{g,1} = detrend(filtfilt(B,A,finalRestCBVData{g,1}),'constant');
        procRestNeuralData{g,1} = detrend(filtfilt(B,A,finalRestNeuralData{g,1}),'constant');
    end
    
    % analyze correlation coefficient between resting epochs
    for h = 1:length(procRestCBVData)
        restPredictionConvolution = conv(procRestNeuralData{h,1},bestKernel,'full');
        restPredictionVals{h,:} = restPredictionConvolution(1:length(procRestCBVData{h,:}));
        restCorrCoeff = corrcoef(procRestCBVData{h,:},restPredictionVals{h,:});
        restCorrCoeffR(h,1) = restCorrCoeff(2,1);
        % Error Variance
        restSSE(h,1) = sum((procRestCBVData{h,:} - restPredictionVals{h,:}).^2);
        % Total Variance
        restSST(h,1) = sum((procRestCBVData{h,:} - (ones(size(procRestCBVData{h,:},1),1)*mean(procRestCBVData{h,:}))).^2);
        % Check that the sum of the residuals is small compared to SSE + SSR4
        restCorrCoeffR2(h,1) = ones(size(restSSE(h,1))) - restSSE(h,1)./restSST(h,1);
    end
    AnalysisResults.HRF_Predictions.Rest.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).R = restCorrCoeffR;
    AnalysisResults.HRF_Predictions.Rest.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).R2 = restCorrCoeffR2;
    AnalysisResults.HRF_Predictions.Rest.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).SSE = restSSE;
    AnalysisResults.HRF_Predictions.Rest.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).SST = restSST;
end

%% Analyze mean CBV during periods of extended whisking
for b = 1:length(filterSets)
    filterSet = filterSets{1,b};
    % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
    [whiskLogical] = FilterEvents_IOS(EventData.(CBVdataType).(['adj' hemDataType]).whisk,WhiskCriteria);
    [puffLogical] = FilterEvents_IOS(EventData.(CBVdataType).(['adj' hemDataType]).whisk,WhiskPuffCriteria);
    combWhiskLogical = logical(whiskLogical.*puffLogical);
    whiskFiles = EventData.(CBVdataType).(['adj' hemDataType]).whisk.fileIDs(combWhiskLogical,:);
    if strcmp(CBVdataType,'CBV') == true
        whiskCBVData = EventData.(CBVdataType).(['adj' hemDataType]).whisk.NormData(combWhiskLogical,:);
    elseif strcmp(CBVdataType,'CBV_HbT') == true
        whiskCBVData = EventData.(CBVdataType).(['adj' hemDataType]).whisk.data(combWhiskLogical,:);
    end
    whiskNeuralData = EventData.(['cortical_' hemDataType]).(neuralBand).whisk.NormData(combWhiskLogical,:);
    
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
            if strcmp(filterSet,'manualSelection') == true
                if strcmp(whiskDay,whiskFileID) && sum(strcmp(whiskFile,manualFileIDs)) == 1
                    whiskFiltLogical{c,1}(e,1) = 1; %#ok<*AGROW>
                    d = d + 1;
                else
                    whiskFiltLogical{c,1}(e,1) = 0;
                end
            elseif strcmp(filterSet,'setDuration') == true
                if strcmp(whiskDay,whiskFileID) && d <= fileTarget
                    whiskFiltLogical{c,1}(e,1) = 1;
                    d = d + 1;
                else
                    whiskFiltLogical{c,1}(e,1) = 0;
                end
            elseif strcmp(filterSet,'entireDuration') == true
                if strcmp(whiskDay,whiskFileID)
                    whiskFiltLogical{c,1}(e,1) = 1;
                    d = d + 1;
                else
                    whiskFiltLogical{c,1}(e,1) = 0;
                end
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
    finalWhiskCBVData = whiskCBVData(whiskFinalFileFilter,:);
    finalWhiskNeuralData = whiskNeuralData(whiskFinalFileFilter,:);
    
    % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    % lowpass filter and detrend each segment
    [B, A] = butter(3,1/(samplingRate/2),'low');
    clear LH_ProcRestData
    clear RH_ProcRestData
    for g = 1:size(finalWhiskCBVData,1)
        procWhiskCBVData{g,1} = detrend(filtfilt(B,A,finalWhiskCBVData(g,samplingRate*2:end)),'constant');
        procWhiskNeuralData{g,1} = detrend(filtfilt(B,A,finalWhiskNeuralData(g,samplingRate*2:end)),'constant');
    end
    
    % analyze correlation coefficient between resting epochs
    for h = 1:length(procWhiskCBVData)
        whiskPredictionConvolution = conv(procWhiskNeuralData{h,1},bestKernel,'full');
        whiskPredictionVals{h,:} = whiskPredictionConvolution(1:length(procWhiskCBVData{h,:}));
        whiskCorrCoeff = corrcoef(procWhiskCBVData{h,:},whiskPredictionVals{h,:});
        whiskCorrCoeffR(h,1) = whiskCorrCoeff(2,1);
        % Error Variance
        whiskSSE(h,1) = sum((procWhiskCBVData{h,:} - whiskPredictionVals{h,:}).^2);
        % Total Variance
        whiskSST(h,1) = sum((procWhiskCBVData{h,:} - (ones(size(procWhiskCBVData{h,:},1),1)*mean(procWhiskCBVData{h,:}))).^2);
        % Check that the sum of the residuals is small compared to SSE + SSR4
        whiskCorrCoeffR2(h,1) = ones(size(whiskSSE(h,1))) - whiskSSE(h,1)./whiskSST(h,1);
    end
    AnalysisResults.HRF_Predictions.Whisk.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).R = whiskCorrCoeffR;
    AnalysisResults.HRF_Predictions.Whisk.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).R2 = whiskCorrCoeffR2;
    AnalysisResults.HRF_Predictions.Whisk.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).SSE = whiskSSE;
    AnalysisResults.HRF_Predictions.Whisk.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).SST = whiskSST;
end

%% Analyze mean CBV during periods of NREM sleep
% pull data from SleepData.mat structure
nremCBVData = SleepData.NREM.data.(CBVdataType).(hemDataType);
nremNeuralData = SleepData.NREM.data.(['cortical_' hemDataType]).(neuralBand);

% analyze correlation coefficient between resting epochs
for h = 1:length(nremCBVData)
    nremPredictionConvolution = conv(nremNeuralData{h,1},bestKernel,'full');
    nremPredictionVals{h,:} = nremPredictionConvolution(1:length(nremCBVData{h,:}));
    nremCorrCoeff = corrcoef(nremCBVData{h,:},nremPredictionVals{h,:});
    nremCorrCoeffR(h,1) = nremCorrCoeff(2,1);
    % Error Variance
    nremSSE(h,1) = sum((nremCBVData{h,:} - nremPredictionVals{h,:}).^2);
    % Total Variance
    nremSST(h,1) = sum((nremCBVData{h,:} - (ones(size(nremCBVData{h,:},1),1)*mean(nremCBVData{h,:}))).^2);
    % Check that the sum of the residuals is small compared to SSE + SSR4
    nremCorrCoeffR2(h,1) = ones(size(nremSSE(h,1))) - nremSSE(h,1)./nremSST(h,1);
end
AnalysisResults.HRF_Predictions.NREM.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).R = nremCorrCoeffR;
AnalysisResults.HRF_Predictions.NREM.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).R2 = nremCorrCoeffR2;
AnalysisResults.HRF_Predictions.NREM.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).SSE = nremSSE;
AnalysisResults.HRF_Predictions.NREM.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).SST = nremSST;

%% Analyze mean CBV during periods of REM sleep
remCBVData = SleepData.REM.data.(CBVdataType).(hemDataType);
remNeuralData = SleepData.REM.data.(['cortical_' hemDataType]).(neuralBand);

% analyze correlation coefficient between resting epochs
for h = 1:length(remCBVData)
    remPredictionConvolution = conv(remNeuralData{h,1},bestKernel,'full');
    remPredictionVals{h,:} = remPredictionConvolution(1:length(remCBVData{h,:}));
    remCorrCoeff = corrcoef(remCBVData{h,:},remPredictionVals{h,:});
    remCorrCoeffR(h,1) = remCorrCoeff(2,1);
    % Error Variance
    remSSE(h,1) = sum((remCBVData{h,:} - remPredictionVals{h,:}).^2);
    % Total Variance
    remSST(h,1) = sum((remCBVData{h,:} - (ones(size(remCBVData{h,:},1),1)*mean(remCBVData{h,:}))).^2);
    % Check that the sum of the residuals is small compared to SSE + SSR4
    remCorrCoeffR2(h,1) = ones(size(remSSE(h,1))) - remSSE(h,1)./remSST(h,1);
end
AnalysisResults.HRF_Predictions.REM.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).R = remCorrCoeffR;
AnalysisResults.HRF_Predictions.REM.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).R2 = remCorrCoeffR2;
AnalysisResults.HRF_Predictions.REM.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).SSE = remSSE;
AnalysisResults.HRF_Predictions.REM.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).SST = remSST;

%% Analyze mean CBV during periods of unstimulated data
for o = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(o,:);
    load(procDataFileID);
    if isempty(ProcData.data.solenoids.LPadSol) == true
        stimLogical(o,1) = 1;
    else
        stimLogical(o,1) = 0;
    end
end
stimLogical = logical(stimLogical);
unstim_procDataFileIDs = procDataFileIDs(stimLogical,:);
for p = 1:size(unstim_procDataFileIDs)
    unstim_procDataFileID = unstim_procDataFileIDs(p,:);
    load(unstim_procDataFileID)
    [~,fileDate,~] = GetFileInfo_IOS(unstim_procDataFileID);
    US_strDay = ConvertDate_IOS(fileDate);
    % pull data from each file
    if strcmp(CBVdataType,'CBV') == true
        unstimCBVData{p,1} = (ProcData.data.(CBVdataType).(['adj' hemDataType]) - RestingBaselines.(baselineType).(CBVdataType).(['adj' hemDataType]).(US_strDay))/RestingBaselines.(baselineType).(CBVdataType).(['adj' hemDataType]).(US_strDay);
    elseif strcmp(CBVdataType,'CBV_HbT') == true
        unstimCBVData{p,1} = ProcData.data.(CBVdataType).(['adj' hemDataType]);
    end
    unstimNeuralData{p,1} = (ProcData.data.(['cortical_' hemDataType]).(neuralBand) - RestingBaselines.(baselineType).(['cortical_' hemDataType]).(neuralBand).(US_strDay))/RestingBaselines.(baselineType).(['cortical_' hemDataType]).(neuralBand).(US_strDay);
end

% detend and lowpass filter each signal
for q = 1:length(unstimCBVData)
    procUnstimCBVData{q,1} = detrend(filtfilt(B,A,unstimCBVData{q,1}),'constant');
    procUnstimNeuralData{q,1} = detrend(filtfilt(B,A,unstimNeuralData{q,1}),'constant');
end

% analyze correlation coefficient between resting epochs
for h = 1:length(procUnstimCBVData)
    unstimPredictionConvolution = conv(procUnstimNeuralData{h,1},bestKernel,'full');
    unstimPredictionVals{h,:} = unstimPredictionConvolution(1:length(procUnstimCBVData{h,:}));
    unstimCorrCoeff = corrcoef(procUnstimCBVData{h,:},unstimPredictionVals{h,:});
    unstimCorrCoeffR(h,1) = unstimCorrCoeff(2,1);
    % Error Variance
    unstimSSE(h,1) = sum((procUnstimCBVData{h,:} - unstimPredictionVals{h,:}).^2);
    % Total Variance
    unstimSST(h,1) = sum((procUnstimCBVData{h,:} - (ones(size(procUnstimCBVData{h,:},1),1)*mean(procUnstimCBVData{h,:}))).^2);
    % Check that the sum of the residuals is small compared to SSE + SSR4
    unstimCorrCoeffR2(h,1) = ones(size(unstimSSE(h,1))) - unstimSSE(h,1)./unstimSST(h,1);
end
AnalysisResults.HRF_Predictions.Unstim.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).R = unstimCorrCoeffR;
AnalysisResults.HRF_Predictions.Unstim.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).R2 = unstimCorrCoeffR2;
AnalysisResults.HRF_Predictions.Unstim.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).SSE = unstimSSE;
AnalysisResults.HRF_Predictions.Unstim.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).SST = unstimSST;

%% Analyze mean CBV during periods of all data
for p = 1:size(procDataFileIDs)
    procDataFileID = procDataFileIDs(p,:);
    load(procDataFileID)
    [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
    AD_strDay = ConvertDate_IOS(fileDate);
    % pull data from each file
    if strcmp(CBVdataType,'CBV') == true
        allDataCBVData{p,1} = (ProcData.data.(CBVdataType).(['adj' hemDataType]) - RestingBaselines.(baselineType).(CBVdataType).(['adj' hemDataType]).(AD_strDay))/RestingBaselines.(baselineType).(CBVdataType).(['adj' hemDataType]).(AD_strDay);
    elseif strcmp(CBVdataType,'CBV_HbT') == true
        allDataCBVData{p,1} = ProcData.data.(CBVdataType).(['adj' hemDataType]);
    end
    allDataNeuralData{p,1} = (ProcData.data.(['cortical_' hemDataType]).(neuralBand) - RestingBaselines.(baselineType).(['cortical_' hemDataType]).(neuralBand).(AD_strDay))/RestingBaselines.(baselineType).(['cortical_' hemDataType]).(neuralBand).(AD_strDay);
end

% detend and lowpass filter each signal
for q = 1:length(unstimCBVData)
    procAllDataCBVData{q,1} = detrend(filtfilt(B,A,allDataCBVData{q,1}),'constant');
    procAllDataNeuralData{q,1} = detrend(filtfilt(B,A,allDataNeuralData{q,1}),'constant');
end

% analyze correlation coefficient between resting epochs
for h = 1:length(procAllDataCBVData)
    allDataPredictionConvolution = conv(procAllDataNeuralData{h,1},bestKernel,'full');
    allDataPredictionVals{h,:} = allDataPredictionConvolution(1:length(procAllDataCBVData{h,:}));
    allDataCorrCoeff = corrcoef(procAllDataCBVData{h,:},allDataPredictionVals{h,:});
    allDataCorrCoeffR(h,1) = allDataCorrCoeff(2,1);
    % Error Variance
    allDataSSE(h,1) = sum((procAllDataCBVData{h,:} - allDataPredictionVals{h,:}).^2);
    % Total Variance
    allDataSST(h,1) = sum((procAllDataCBVData{h,:} - (ones(size(procAllDataCBVData{h,:},1),1)*mean(procAllDataCBVData{h,:}))).^2);
    % Check that the sum of the residuals is small compared to SSE + SSR4
    allDataCorrCoeffR2(h,1) = ones(size(allDataSSE(h,1))) - allDataSSE(h,1)./allDataSST(h,1);
end
AnalysisResults.HRF_Predictions.AllData.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).R = allDataCorrCoeffR;
AnalysisResults.HRF_Predictions.AllData.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).R2 = allDataCorrCoeffR2;
AnalysisResults.HRF_Predictions.AllData.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).SSE = allDataSSE;
AnalysisResults.HRF_Predictions.AllData.(fileSet).(CBVdataType).(hemDataType).(neuralBand).(filterSet).SST = allDataSST;

%% save results strucure
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

end