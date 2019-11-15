function [AnalysisResults] = PredictHemodynamicChanges_IOS(params,hemDataType,neuralBand,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
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

bestKernel = AnalysisResults.HRFs.(neuralBand).(hemDataType).bestKernel;
figure;
plot(bestKernel)

%% Analyze mean CBV during periods of rest
% use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
[restLogical] = FilterEvents_IOS(RestData.CBV_HbT.(hemDataType),RestCriteria);
[puffLogical] = FilterEvents_IOS(RestData.CBV_HbT.(hemDataType),RestPuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFiles = RestData.CBV_HbT.(hemDataType).fileIDs(combRestLogical,:);
restingCBVData = RestData.CBV_HbT.(hemDataType).data(combRestLogical,:);
restingNeuralData = RestData.(['cortical_' hemDataType(4:end)]).(neuralBand).NormData(combRestLogical,:);

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
finalRestCBVData = restingCBVData(restFinalFileFilter,:);
finalRestNeuralData = restingNeuralData(restFinalFileFilter,:);

% only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
% original epoch create so we can add a sample of two back to the end for those just under 10 seconds
% lowpass filter and detrend each segment
[B, A] = butter(3,1/(samplingRate/2),'low');
clear LH_ProcRestData RH_ProcRestData
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
AnalysisResults.HRF_Predictions.Rest.(neuralBand).(hemDataType).R = restCorrCoeffR;
AnalysisResults.HRF_Predictions.Rest.(neuralBand).(hemDataType).R2 = restCorrCoeffR2;
AnalysisResults.HRF_Predictions.Rest.(neuralBand).(hemDataType).SSE = restSSE;
AnalysisResults.HRF_Predictions.Rest.(neuralBand).(hemDataType).SST = restSST;

%% Analyze mean CBV during periods of extended whisking
% use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
[whiskLogical] = FilterEvents_IOS(EventData.CBV_HbT.(['adj' hemDataType]).whisk,WhiskCriteria);
[puffLogical] = FilterEvents_IOS(EventData.CBV_HbT.(['adj' hemDataType]).whisk,WhiskPuffCriteria);
combWhiskLogical = logical(whiskLogical.*puffLogical);
whiskFiles = EventData.CBV_HbT.(['adj' hemDataType]).whisk.fileIDs(combWhiskLogical,:);
if strcmp(CBVdataType,'CBV') == true
    whiskCBVData = EventData.CBV_HbT.(['adj' hemDataType]).whisk.NormData(combWhiskLogical,:);
elseif strcmp(CBVdataType,'CBV_HbT') == true
    whiskCBVData = EventData.CBV_HbT.(['adj' hemDataType]).whisk.data(combWhiskLogical,:);
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
AnalysisResults.HRF_Predictions.Whisk.(neuralBand).(hemDataType).R = whiskCorrCoeffR;
AnalysisResults.HRF_Predictions.Whisk.(neuralBand).(hemDataType).R2 = whiskCorrCoeffR2;
AnalysisResults.HRF_Predictions.Whisk.(neuralBand).(hemDataType).SSE = whiskSSE;
AnalysisResults.HRF_Predictions.Whisk.(neuralBand).(hemDataType).SST = whiskSST;

%% Analyze mean CBV during periods of NREM sleep
% pull data from SleepData.mat structure
nremCBVData = SleepData.NREM.data.CBV_HbT.(hemDataType);
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
AnalysisResults.HRF_Predictions.NREM.(neuralBand).(hemDataType).R = nremCorrCoeffR;
AnalysisResults.HRF_Predictions.NREM.(neuralBand).(hemDataType).R2 = nremCorrCoeffR2;
AnalysisResults.HRF_Predictions.NREM.(neuralBand).(hemDataType).SSE = nremSSE;
AnalysisResults.HRF_Predictions.NREM.(neuralBand).(hemDataType).SST = nremSST;

%% Analyze mean CBV during periods of REM sleep
remCBVData = SleepData.REM.data.CBV_HbT.(hemDataType);
remNeuralData = SleepData.REM.data.(['cortical_' hemDataType(4:end)]).(neuralBand);

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
AnalysisResults.HRF_Predictions.REM.(neuralBand).(hemDataType).R = remCorrCoeffR;
AnalysisResults.HRF_Predictions.REM.(neuralBand).(hemDataType).R2 = remCorrCoeffR2;
AnalysisResults.HRF_Predictions.REM.(neuralBand).(hemDataType).SSE = remSSE;
AnalysisResults.HRF_Predictions.REM.(neuralBand).(hemDataType).SST = remSST;

%% save results strucure
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

end