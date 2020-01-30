function [AnalysisResults] = AnalyzeAwakeHRF_IOS(hemDataType,neuralBand,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Davis Haselden & Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs: //
%________________________________________________________________________________________________________________________

% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)

clear global responseFuncData
global responseFuncData
fileBreaks = strfind(baselineDataFileID,'_');
animalID = baselineDataFileID(1:fileBreaks(1)-1);
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
samplingRate = 30;
for g = 1:length(manualFileIDs)
    % Handle file number, load each successive file
    fileID = [animalID '_' manualFileIDs{g,1} '_ProcData.mat'];
    load(fileID);
    [animalID,fileDate,~] = GetFileInfo_IOS(fileID);
    strDay = ConvertDate_IOS(fileDate); 
    %% CBV data - normalize by rest, lowpass filer, detrend
    [B, A] = butter(3,1/(samplingRate/2),'low');
    CBVData = ProcData.data.CBV_HbT.(hemDataType);
    responseFuncData.procCBV(g,:) = detrend(filtfilt(B,A,CBVData),'constant');
    %% Neural band power data - normalize by rest, lowpass filer, detrend
    neuralData = ProcData.data.(['cortical_' hemDataType(4:end)]).(neuralBand);
    normNeuro = (neuralData - RestingBaselines.manualSelection.(['cortical_' hemDataType(4:end)]).(neuralBand).(strDay))/RestingBaselines.manualSelection.(['cortical_' hemDataType(4:end)]).(neuralBand).(strDay);
    responseFuncData.procNeuro(g,:) = detrend(filtfilt(B,A,normNeuro),'constant');
end

%% Gamma function and minimizing error using fminsearch
disp(['Minimizing gamma-function between the ' hemDataType ' HbT and ' neuralBand]); disp(' ')
responseFuncData.dt = 1/samplingRate;
% Initialize option parameters for fminsearch
options = optimset('Display','iter','FunValCheck','on','MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-8,'TolX',1e-8);
fitFunction = @fitGammaFunctionError;
% Initial conditions from Winder et al. 2017 [A0,B0,amplitude]
x0 = [4,.75,1e-2];
% lb = [0.25,0.1,-inf];
% ub = [5,5,inf];
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% nonlcon = [];
fitGammaFunctionError(x0);
% [bestFitParams,fVal] = fmincon(fitFunction,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
[bestFitParams,fVal] = fminsearch(fitFunction,x0,options);

% report optimized values
responseFuncData.bestX = bestFitParams;
W = responseFuncData.bestX(1);
T = responseFuncData.bestX(2);
A = responseFuncData.bestX(3);
alpha = ((T/W).^2)*8.0*log(2.0);
beta = (W.^2)./(T*8.0*log(2.0));
bestKernel = real(A*((responseFuncData.kernelT/T).^alpha).*exp((responseFuncData.kernelT-T)/(-beta)));
responseFuncData.bestKernel = bestKernel;
responseFuncData.fval = fVal;

responseFuncData.halfMax = min(responseFuncData.bestKernel)/2;
responseFuncData.FWHM = responseFuncData.dt*length(find(responseFuncData.bestKernel > responseFuncData.halfMax));
[~,peakPoint] = max(abs(responseFuncData.bestKernel));    % find the peak of the HRF
responseFuncData.kernelPeak = peakPoint*responseFuncData.dt;  % time of the peak
AnalysisResults.HRFs.(neuralBand).(hemDataType) = responseFuncData;

%% Evaluate kernel prediction success
for h = 1:size(responseFuncData.procNeuro,1)
    predictionConvolution = conv(responseFuncData.procNeuro(h,:),AnalysisResults.HRFs.(neuralBand).(hemDataType).bestKernel,'full');
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
AnalysisResults.HRFs.(neuralBand).(hemDataType).R = corrCoeffR;
AnalysisResults.HRFs.(neuralBand).(hemDataType).R2 = corrCoeffR2;
AnalysisResults.HRFs.(neuralBand).(hemDataType).SSE = SSE;
AnalysisResults.HRFs.(neuralBand).(hemDataType).SST = SST;
save([animalID '_AnalysisResults.mat'],'AnalysisResults');
end

function [unexplainedVariance] = fitGammaFunctionError(X)
% this function returns the unexplained variance
global responseFuncData
kernelMaxT = 5;
kernelT = 0:responseFuncData.dt:kernelMaxT;
responseFuncData.kernelT = kernelT;
for j = 1:size(responseFuncData.procNeuro,1)
    middle = length(kernelT):(size(responseFuncData.procNeuro,2) - length(kernelT));
    W = X(1);
    T = X(2);
    A = X(3);
    alpha = ((T/W).^2)*8.0*log(2.0);
    beta = (W.^2)./(T*8.0*log(2.0));
    kernel = real(A*((kernelT/T).^alpha).*exp((kernelT-T)/(-beta)));
%     checkNAN = isnan(kernel);
%     kernel(checkNAN) = 0;
    predictCBV = conv(responseFuncData.procNeuro(j,:),(kernel));
    responseFuncData.CBVprediction(j,:) = predictCBV(1:size(responseFuncData.procNeuro,2));
    % r^2 on on test data
    unexplainedVarianceAll(j) = var(responseFuncData.CBVprediction(j,middle) - responseFuncData.procCBV(j,middle))/var(responseFuncData.procCBV(j,middle));
end
unexplainedVariance = sum(unexplainedVarianceAll);
if isnan(unexplainedVariance) == true
    keyboard
end
responseFuncData.kernelT = kernelT;
end
