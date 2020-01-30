function [AnalysisResults] = AnalyzeAwakeHRF(params,fileSet,CBVdataType,hemDataType,neuralBand,AnalysisResults)
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

% list of unstim Procdata.mat files
procDataFileStruct = dir('*_Procdata.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)

fileTarget = params.targetMinutes/15;
for a = 1:size(procDataFileIDs,1)
    fileBreaks = strfind(procDataFileIDs(a,:),'_');
    procDataFileName = procDataFileIDs(a,:);
    short_procDataFileList{a,1} = procDataFileName(fileBreaks(1)+1:fileBreaks(end)-1); %#ok<*AGROW>
    procDataFileList{a,1} = procDataFileName; %#ok<*AGROW>
end
uniqueDays = GetUniqueDays_IOS(short_procDataFileList);

for b = 1:length(uniqueDays)
    uniqueDay = uniqueDays(b);
    c = 1;
    for d = 1:size(procDataFileIDs,1)
        fullFileID = procDataFileIDs(d,:);
        [~,fileDate,~] = GetFileInfo_IOS(fullFileID);
        if strcmp(uniqueDay,fileDate) && c <= fileTarget
            fileFilterLogical{b,1}(d,1) = 1;
            c = c + 1;
        else
            fileFilterLogical{b,1}(d,1) = 0;
        end
    end
end
% Combine the 3 logicals so that it reflects the first "x" number of files from each day
finalFileFilterLogical = any(sum(cell2mat(fileFilterLogical'),2),2);
validProcDataFileList = procDataFileList(finalFileFilterLogical,:);

for e = 1:length(uniqueDays)
    uniqueDay = uniqueDays{e,1};
    f = 1;
    for g = 1:length(validProcDataFileList)
        validProcDataFile = validProcDataFileList{g,1};
        [~,validProcDataFileDate,~] = GetFileInfo_IOS(validProcDataFile);
        if strcmp(uniqueDay,validProcDataFileDate) == true
            if rem(e,2) == 1
                if f == 1
                    fileSetALogical(g,1) = 1;
                    f = f + 1;
                elseif f > 1
                    fileSetALogical(g,1) = 0;
                end
            elseif rem(e,2) == 0
                if f == 1
                    fileSetALogical(g,1) = 0;
                    f = f + 1;
                elseif f > 1
                    fileSetALogical(g,1) = 1;
                end
            end
        end
    end
end
fileSetALogical = logical(fileSetALogical);
fileSetBLogical = ~fileSetALogical;
fileSetAList = validProcDataFileList(fileSetALogical,:);
fileSetBList = validProcDataFileList(fileSetBLogical,:);         
if strcmp(fileSet,'fileSetA')
    fileList = fileSetAList;
else
    fileList = fileSetBList;
end
clear global responseFuncData
global responseFuncData

for e = 1:length(fileList)
    % Handle file number, load each successive file
    fileID = fileList{e,1};
    load(fileID);
    [animalID,fileDate,~] = GetFileInfo_IOS(fileID);
    strDay = ConvertDate_IOS(fileDate);
    
    %% CBV data - normalize by rest, lowpass filer, detrend
    [B, A] = butter(4,1/(30/2),'low');
    CBVData = ProcData.data.(CBVdataType).(hemDataType);
    if strcmp(CBVdataType,'CBV') == true
        normCBV = (CBVData - RestingBaselines.manualSelection.(CBVdataType).(hemDataType).(strDay))/RestingBaselines.manualSelection.(CBVdataType).(hemDataType).(strDay);
        responseFuncData.procCBV(e,:) = detrend(filtfilt(B,A,normCBV),'constant');
    elseif strcmp(CBVdataType,'CBV_HbT') == true
        responseFuncData.procCBV(e,:) = detrend(filtfilt(B,A,CBVData),'constant');
    end
    %% Neural band power data - normalize by rest, lowpass filer, detrend
    neuralData = ProcData.data.(['cortical_' hemDataType]).(neuralBand);
    normNeuro = (neuralData - RestingBaselines.manualSelection.(['cortical_' hemDataType]).(neuralBand).(strDay))/RestingBaselines.manualSelection.(['cortical_' hemDataType]).(neuralBand).(strDay);
    responseFuncData.procNeuro(e,:) = detrend(filtfilt(B,A,normNeuro),'constant');
end

%% Gamma function and minimizing error using fminsearch
responseFuncData.dt = 1/30;   % 30 Hz sampling rate
% Initialize option parameters for fminsearch
options = optimset('Display','iter','FunValCheck','on','MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-8,'TolX',1e-8);
fitFunction = @fitGammaFunctionError;
% Initial conditions from Winder et al. 2017 [A0,B0,amplitude]
if strcmp(CBVdataType,'CBV') == true
    x0 = [4,.75,-1e-3];
elseif strcmp(CBVdataType,'CBV_HbT') == true
    x0 = [4,.75,1e-3];
end
lb = [2,0.25,-inf];
ub = [4,0.75,inf];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
fitGammaFunctionError(x0);
[bestFitParams,fVal] = fmincon(fitFunction,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
% [bestFitParams,fVal] = fminsearch(fitFunction,x0,options);

% report optimized values
responseFuncData.bestX = bestFitParams;
responseFuncData.bestKernel = bestFitParams(3)*gampdf(responseFuncData.kernelT,responseFuncData.bestX(1),responseFuncData.bestX(2));
responseFuncData.fval = fVal;

responseFuncData.halfMax = min(responseFuncData.bestKernel)/2;
responseFuncData.FWHM = responseFuncData.dt*length(find(responseFuncData.bestKernel > responseFuncData.halfMax));
[~,peakPoint] = max(abs(responseFuncData.bestKernel));    % find the peak of the HRF
responseFuncData.kernelPeak = peakPoint*responseFuncData.dt;  % time of the peak
AnalysisResults.HRFs.(fileSet).(CBVdataType).(hemDataType).(neuralBand) = responseFuncData;

%% Evaluate kernel prediction success
for f = 1:size(responseFuncData.procNeuro,1)
    predictionConvolution = conv(responseFuncData.procNeuro(f,:),AnalysisResults.HRFs.(fileSet).(CBVdataType).(hemDataType).(neuralBand).bestKernel,'full');
    predictionVals(f,:) = predictionConvolution(1:length(responseFuncData.procCBV(f,:)));
    corrCoeff = corrcoef(responseFuncData.procCBV(f,:),predictionVals(f,:));
    corrCoeffR(f,1) = corrCoeff(2,1);
    % Error Variance
    SSE(f,1) = sum((responseFuncData.procCBV(f,:) - predictionVals(f,:)).^2);
    % Total Variance
    SST(f,1) = sum((responseFuncData.procCBV(f,:) - (ones(size(responseFuncData.procCBV(f,:),1),1)*mean(responseFuncData.procCBV(f,:)))).^2);
    % Check that the sum of the residuals is small compared to SSE + SSR
    corrCoeffR2(f,1) = ones(size(SSE(f,1))) - SSE(f,1)./SST(f,1);
end
AnalysisResults.HRFs.(fileSet).(CBVdataType).(hemDataType).(neuralBand).R = corrCoeffR;
AnalysisResults.HRFs.(fileSet).(CBVdataType).(hemDataType).(neuralBand).R2 = corrCoeffR2;
AnalysisResults.HRFs.(fileSet).(CBVdataType).(hemDataType).(neuralBand).SSE = SSE;
AnalysisResults.HRFs.(fileSet).(CBVdataType).(hemDataType).(neuralBand).SST = SST;
save([animalID '_AnalysisResults.mat'],'AnalysisResults');
end

function [unexplainedVariance] = fitGammaFunctionError(X)
% this function returns the unexplained variance
global responseFuncData

kernelMaxT = 10;
kernelT = 0:responseFuncData.dt:kernelMaxT;
responseFuncData.kernelT = kernelT;

for g = 1:size(responseFuncData.procNeuro,1)
    middle = length(kernelT):(size(responseFuncData.procNeuro,2) - length(kernelT));
    kernel = X(3)*gampdf(kernelT,X(1),X(2));
    predictCBV = conv(responseFuncData.procNeuro(g,:),(kernel));
    responseFuncData.CBVprediction(g,:) = predictCBV(1:size(responseFuncData.procNeuro,2));
    % r^2 on on test data
    unexplainedVarianceAll(g) = var(responseFuncData.CBVprediction(g,middle) - responseFuncData.procCBV(g,middle))/var(responseFuncData.procCBV(g,middle));
end

unexplainedVariance = sum(unexplainedVarianceAll);
responseFuncData.kernelT = kernelT;
end
