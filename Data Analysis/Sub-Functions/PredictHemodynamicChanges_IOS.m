function [AnalysisResults] = PredictHemodynamicChanges(params,fileSets,CBVDataTypes,hemDataTypes,neuralBands,AnalysisResults)
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

% identify animal's ID and pull important infortmat
fileBreaks = strfind(restDataFileID, '_');
animalID = restDataFileID(1:fileBreaks(1)-1);
trialDuration_min = RestData.CBV.LH.trialDuration_sec/60;   % min
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
fileTarget = params.targetMinutes/trialDuration_min;
filterSets = {'manualSelection','setDuration','entireDuration'};
samplingRate = RestData.CBV.LH.CBVCamSamplingRate;

RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};

PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};



for loop through 

    
end

AnalysisResults.HRFpredictions.Rest.(fileSet).(CBVDataType).(hemDataType).(neuralBand).specBand
AnalysisResults.HRFpredictions.Rest.(fileSet).(CBVDataType).(hemDataType).(neuralBand).filtBand




[B, A] = butter(4,1/(30/2),'low');
CBVData = ProcData.data.(CBVdataType).(hemDataType);
if strcmp(CBVdataType,'CBV') == true
    normCBV = (CBVData - RestingBaselines.manualSelection.(CBVdataType).(hemDataType).(strDay))/RestingBaselines.manualSelection.(CBVdataType).(hemDataType).(strDay);
    responseFuncData.procCBV(g,:) = detrend(filtfilt(B,A,normCBV),'constant');
elseif strcmp(CBVdataType,'CBV_HbT') == true
    responseFuncData.procCBV(g,:) = detrend(filtfilt(B,A,CBVData),'constant');
end
%% Neural band power data - normalize by rest, lowpass filer, detrend
neuralData = ProcData.data.(['cortical_' hemDataType]).(neuralBand);
normNeuro = (neuralData - RestingBaselines.manualSelection.(['cortical_' hemDataType]).(neuralBand).(strDay))/RestingBaselines.manualSelection.(['cortical_' hemDataType]).(neuralBand).(strDay);
responseFuncData.procNeuro(g,:) = detrend(filtfilt(B,A,normNeuro),'constant');

%% Evaluate kernel prediction success
for h = 1:size(responseFuncData.procNeuro,1)
    predictionConvolution = conv(responseFuncData.procNeuro(h,:),AnalysisResults.HRFs.(fileSet).(CBVdataType).(hemDataType).(neuralBand).bestKernel,'full');
    predictionVals(h,:) = predictionConvolution(1:length(responseFuncData.procCBV(h,:)));
    corrCoeff = corrcoef(responseFuncData.procCBV(h,:),predictionVals(h,:));
    corrCoeffR(h,1) = corrCoeff(2,1);
    % Error Variance
    SSE(h,1) = sum((responseFuncData.procCBV(h,:) - predictionVals(h,:)).^2);
    % Total Variance
    SST(h,1) = sum((responseFuncData.procCBV(h,:) - (ones(size(responseFuncData.procCBV(h,:),1),1)*mean(responseFuncData.procCBV(h,:)))).^2);
    % Check that the sum of the residuals is small compared to SSE + SSR4
    corrCoeffR2(h,1) = ones(size(SSE(h,1))) - SSE(h,1)./SST(h,1);
end
AnalysisResults.HRFs.(fileSet).(CBVdataType).(hemDataType).(neuralBand).R = corrCoeffR;
AnalysisResults.HRFs.(fileSet).(CBVdataType).(hemDataType).(neuralBand).R2 = corrCoeffR2;
AnalysisResults.HRFs.(fileSet).(CBVdataType).(hemDataType).(neuralBand).SSE = SSE;
AnalysisResults.HRFs.(fileSet).(CBVdataType).(hemDataType).(neuralBand).SST = SST;


end