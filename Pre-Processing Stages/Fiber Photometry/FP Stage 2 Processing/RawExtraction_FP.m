function RawExtraction_FP(csvFiles)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Kyle W. Gheres and Md Shakhawat Hossain
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

for qq = 1:size(csvFiles,1)
    fileName = csvFiles(qq,:);
    correctionFlag = 'y'; % If GCaMP signal need to be corrected
    opticalChannelNames = {'RHControl','RHGCaMP7s','RHBloodVolume','LHControl','LHGCaMP7s','LHBloodVolume','LH560Raw'};
    % identify the extension
    extInd = strfind(fileName(1,:),'.');
    % identify the underscores
    fileBreaks = strfind(fileName(1,:),'_');
    % file parameters
    FiberData.notes.animalID = fileName(1:fileBreaks(1) - 1);
    FiberData.notes.date = fileName(fileBreaks(1) + 1:fileBreaks(2) - 1);
    FiberData.notes.sessionNumber = fileName(fileBreaks(2) + 1:extInd - 1);
    %% read.CSV
    csvData = csvread(fileName,2,0);
    FiberData.notes.sessionDuration_sec = round(csvData(end,1));
    % 2 = RH 405nm, 3 = RH 465nm, 5 = RH 560nm, 7 = LH 405nm, 8 = LH 465nm, 10 = RL 560nm
    dataChannels = [2,3,5,7,8,10];
    channelData = csvData(:,dataChannels);
    syncData = csvData(:,11);
    %% parameters
    FiberData.notes.channels = opticalChannelNames;
    FiberData.notes.decimation = 10;
    FiberData.notes.samplingRate = round(length(channelData)/(FiberData.notes.sessionDuration_sec));
    FiberData.notes.fitFreq = 0.01;
    FiberData.notes.lowFreq = 1;
    FiberData.notes.finalFreq = [0.01,1];
    FiberData.notes.trialLength = 15.50; % minutes
    FiberData.notes.idleTime = 4.50; % minutes
    %% determine sampling rate
    timeArray = csvData(:,1);
    dTimeArray = diff(timeArray);
    fileBreaks = dTimeArray > 1;
    fileBreaksIndex = find(fileBreaks == 1);
    FiberData.notes.samplingRate = round(fileBreaksIndex(1)/timeArray(fileBreaksIndex(1)));
    %% filter characteristics
    [z1,p1,k1] = butter(3,FiberData.notes.fitFreq/(0.5*FiberData.notes.samplingRate),'low'); % design lowpass filter for hemodynamic correction
    [sosFit,gFit] = zp2sos(z1,p1,k1);
    [z2,p2,k2] = butter(3,FiberData.notes.lowFreq/(0.5*FiberData.notes.samplingRate),'low'); % low pass for optical data to physiologically relevant range
    [sosLow,gLow] = zp2sos(z2,p2,k2);
    [z3,p3,k3] = butter(3,FiberData.notes.finalFreq/(0.5*FiberData.notes.samplingRate),'bandpass'); % low pass filter for locomotion data
    [sosFinal,gFinal] = zp2sos(z3,p3,k3);
    %% Remove photobleaching/metaloism baseline drift
    spacing = 1:1:length(channelData(:,3));
    % low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
    filtData = filtfilt(sosFit,gFit,channelData);
    %% RH correct TRITC blood volume
    [fitVals] = fit(spacing',filtData(:,3),'exp2');
    coeffVals = coeffvalues(fitVals);
    RH_predictedCBV = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
    RH_correctedCBV = channelData(:,3) - RH_predictedCBV';
    FiberData.FitStruct.RH_CBV = fitVals;
    %% RH correct Ca2+ dependent GCaMP
    [fitVals] = fit(spacing',filtData(:,2),'exp2');
    coeffVals = coeffvalues(fitVals);
    RH_predictedGCaMP = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
    RH_correctedGCaMP = channelData(:,2) - RH_predictedGCaMP';
    FiberData.FitStruct.RH_GCaMP = fitVals;
    %% RH correct Ca2+ independent GCaMP
    [fitVals] = fit(spacing',filtData(:,1),'exp2');
    coeffVals = coeffvalues(fitVals);
    RH_predictedGFP = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
    RH_correctedGFP = channelData(:,1) - RH_predictedGFP';
    FiberData.FitStruct.RH_GFP = fitVals;
    %% LH correct TRITC blood volume
    [fitVals] = fit(spacing',filtData(:,6),'exp2');
    coeffVals = coeffvalues(fitVals);
    LH_predictedCBV = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
    LH_correctedCBV = channelData(:,6) - LH_predictedCBV';
    FiberData.FitStruct.LH_CBV = fitVals;
    %% LH correct Ca2+ dependent GCaMP
    [fitVals] = fit(spacing',filtData(:,5),'exp2');
    coeffVals = coeffvalues(fitVals);
    LH_predictedGCaMP = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
    LH_correctedGCaMP = channelData(:,5) - LH_predictedGCaMP';
    FiberData.FitStruct.LH_GCaMP = fitVals;
    %% LH correct Ca2+ independent GCaMP
    [fitVals] = fit(spacing',filtData(:,4),'exp2');
    coeffVals = coeffvalues(fitVals);
    LH_predictedGFP = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
    LH_correctedGFP = channelData(:,4) - LH_predictedGFP';
    FiberData.FitStruct.LH_GFP = fitVals;
    %% detrend data
    detrendData = [RH_correctedGFP,RH_correctedGCaMP,RH_correctedCBV,LH_correctedGFP,LH_correctedGCaMP,LH_correctedCBV];
    lowPassData = filtfilt(sosLow,gLow,detrendData);
    for q = 1:size(lowPassData,2)
        rescaledData(:,q) = rescale(lowPassData(:,q),0,1); % rescale all data between 0 to 1
    end
    %% Z-score data
    avgData = mean(rescaledData,1);
    stdData = std(rescaledData,0,1);
    avgMatrix = repmat(avgData,length(rescaledData),1);
    stdMatrix = repmat(stdData,length(rescaledData),1);
    zScoredFiberData = (rescaledData - avgMatrix)./stdMatrix;
    if strcmpi(correctionFlag,'y')
        %% Hemodynamic response correction code using Mle
        binEdges = (0:0.001:1);
        figure;
        dataHist = histogram2(rescaledData(:,3),rescaledData(:,2),binEdges,binEdges);
        binCounts = dataHist.BinCounts;
        allCounts = sum(binCounts,'all');
        allHist = (binCounts/allCounts)*100;
        xCounts = sum(binCounts,2);
        mleFit.RH.allHist = allHist;
        for binNum = 2:size(binCounts,1)
            if binNum == 2
                percData(1) = (xCounts(1)/sum(xCounts))*100;
            end
            percData(binNum) = (sum(xCounts(1:binNum))/sum(xCounts))*100;
        end
        startInd = find(percData <= 2,1,'last') + 1; % exclude lower 2% of data points
        endInd = find(percData >= 98,1,'first') - 1; % exclude upper 2% of data points
        mleFit.RH.startInd = startInd;
        mleFit.RH.endInd = endInd;
        for binNum = startInd:endInd
            index = (binNum - startInd) + 1;
            leftEdge = binEdges(binNum);
            rightEdge = binEdges(binNum + 1);
            dataInds = rescaledData(:,3) >= leftEdge & rescaledData(:,3) <= rightEdge;
            theData = rescaledData(dataInds,1);
            [phat,pci] = mle(theData,'distrbution','normal');
            figure(101);
            normHist = histogram(theData,binEdges);
            normCounts = normHist.BinCounts./sum(normHist.BinCounts);
            theFit = pdf('normal',normHist.BinEdges(1:(end - 1)),phat(1),phat(2));
            theFit = theFit./sum(theFit);
            rsqr = 1 - (sum((normCounts - theFit).^2)/sum((normCounts - mean(normCounts)).^2));
            mleFit.RH.fitnotes.phat(index,:) = phat;
            mleFit.RH.fitnotes.pci(:,:,index) = pci;
            mleFit.RH.goodness.rsqr(index) = rsqr;
            mleFit.RH.fitData.eGFPData{index} = theData;
            mleFit.RH.fitData.BinEdges(index,:) = [leftEdge,rightEdge];
            mleFit.RH.fitData.normHist(index,:) = normCounts;
            mleFit.RH.fitData.HistFit(index,:) = theFit;
            mleFit.RH.fitData.fitHistEdges(index,:) = normHist.BinEdges;
        end
        [RH_linFit] = fitlm(binEdges(startInd:endInd)',mleFit.RH.fitnotes.phat(:,1),'RobustOpts','on');
        mleFit.RH.linFit = RH_linFit;
        %% Plot the fit
        RH_yInt = table2array(mleFit.RH.linFit.Coefficients(1,1));
        RH_slope = table2array(mleFit.RH.linFit.Coefficients(2,1));
        RH_finalGCaMP = rescaledData(:,2) - (RH_slope*rescaledData(:,3) + RH_yInt);
        Corrected465.RH = RH_finalGCaMP;
        %% LH
        binEdges = 0:0.001:1;
        figure;
        dataHist = histogram2(rescaledData(:,6),rescaledData(:,5),binEdges,binEdges);
        binCounts = dataHist.BinCounts;
        allCounts = sum(binCounts,'all');
        allHist = (binCounts/allCounts)*100;
        xCounts = sum(binCounts,2);
        mleFit.LH.allHist = allHist;
        for binNum = 2:size(binCounts,1)
            if binNum == 2
                percData(1) = (xCounts(1)/sum(xCounts))*100;
            end
            percData(binNum) = (sum(xCounts(1:binNum))/sum(xCounts))*100;
        end
        startInd = find(percData <= 2,1,'last') + 1; % exclude lower 2% of data points
        endInd = find(percData >= 98,1,'first') - 1; % exclude upper 2% of data points
        mleFit.LH.startInd = startInd;
        mleFit.LH.endInd = endInd;
        for binNum = startInd:endInd
            index = (binNum - startInd) + 1;
            leftEdge = binEdges(binNum);
            rightEdge = binEdges(binNum + 1);
            dataInds = rescaledData(:,6) >= leftEdge & rescaledData(:,6) <= rightEdge;
            theData = rescaledData(dataInds,1);
            [phat,pci] = mle(theData,'distrbution','normal');
            figure(102);
            normHist = histogram(theData,binEdges);
            normCounts = normHist.BinCounts./sum(normHist.BinCounts);
            theFit = pdf('normal',normHist.BinEdges(1:(end - 1)),phat(1),phat(2));
            theFit = theFit./sum(theFit);
            rsqr = 1 - (sum((normCounts - theFit).^2)/sum((normCounts - mean(normCounts)).^2));
            mleFit.LH.fitnotes.phat(index,:) = phat;
            mleFit.LH.fitnotes.pci(:,:,index) = pci;
            mleFit.LH.goodness.rsqr(index) = rsqr;
            mleFit.LH.fitData.eGFPData{index} = theData;
            mleFit.LH.fitData.BinEdges(index,:) = [leftEdge,rightEdge];
            mleFit.LH.fitData.normHist(index,:) = normCounts;
            mleFit.LH.fitData.HistFit(index,:) = theFit;
            mleFit.LH.fitData.fitHistEdges(index,:) = normHist.BinEdges;
        end
        [LH_linFit] = fitlm(binEdges(startInd:endInd)',mleFit.LH.fitnotes.phat(:,1),'RobustOpts','on');
        mleFit.LH.linFit = LH_linFit;
        %% Plot the fit
        LH_yInt = table2array(mleFit.LH.linFit.Coefficients(1,1));
        LH_slope = table2array(mleFit.LH.linFit.Coefficients(2,1));
        LH_finalGCaMP = rescaledData(:,5) - (LH_slope*rescaledData(:,6) + LH_yInt);
        Corrected465.LH = LH_finalGCaMP;
    else
        Corrected465.RH = rescaledData(:,2);
        Corrected465.LH = rescaledData(:,5);
    end
    close all
    Smooth465.RH = filtfilt(sosLow,gLow,Corrected465.RH); % bandpass filter data between [0.01 and 1] Hz
    Z465.RH = (Smooth465.RH - mean(Smooth465.RH))/std(Smooth465.RH); % Z score GCaMP data
    Smooth465.LH = filtfilt(sosLow,gLow,Corrected465.LH); % bandpass filter data between [0.01 and 1] Hz
    Z465.LH = (Smooth465.LH - mean(Smooth465.LH))/std(Smooth465.LH); % Z score GCaMP data
    %% correct GCaMP channel for hemodynamic attenuation
    uncorrectedzScoredGCaMPData.LH = filtfilt(sosFinal,gFinal,zScoredFiberData(:,2));
    uncorrectedzScoredGCaMPData.RH = filtfilt(sosFinal,gFinal,zScoredFiberData(:,5));
    zScoredFiberData(:,2) = Z465.RH;
    zScoredFiberData(:,5) = Z465.LH;
    %% split the data into trials
    fields = {'RH_405','RH_465','RH_560','LH_405','LH_465','LH_560'};
    for aa = 1:length(fileBreaksIndex) + 1
        if aa == 1
            for bb = 1:length(fields)
                FiberData.(fields{1,bb}).raw = channelData(1:fileBreaksIndex(1):end,bb);
                FiberData.(fields{1,bb}).zScored = zScoredFiberData(fileBreaksIndex(1):end,bb);
            end
            FiberData.syncData = syncData(1:fileBreaksIndex(1));
            % Photobleaching removed and HbT corrected
            FiberData.RH_405.rescaledData = rescaledData(1:fileBreaksIndex(1),1);
            FiberData.RH_465.rescaledData = Corrected465.RH(1:fileBreaksIndex(1));
            FiberData.RH_560.rescaledData = rescaledData(1:fileBreaksIndex(1),3);
            FiberData.LH_405.rescaledData = rescaledData(1:fileBreaksIndex(1),4);
            FiberData.LH_465.rescaledData = Corrected465.LH(1:fileBreaksIndex(1));
            FiberData.LH_560.rescaledData = rescaledData(1:fileBreaksIndex(1),6);
            % HbT uncorrected corrected
            FiberData.RH_465.uncorrected = rescaledData(1:fileBreaksIndex(1),2);
            FiberData.LH_465.uncorrected = rescaledData(1:fileBreaksIndex(1),5);
            FiberData.RH_465.zScoredUncorrected = uncorrectedzScoredGCaMPData.RH(1:fileBreaksIndex(1));
            FiberData.LH_465.zScoredUncorrected = uncorrectedzScoredGCaMPData.LH(1:fileBreaksIndex(1));
            save([FiberData.notes.animalID '_' FiberData.notes.date '_' FiberData.notes.sessionNumber '_Trial'  num2str(aa) '_FiberData.mat'],'FiberData','SessionData.FitStruct','-v7.3');
        elseif aa == length(fileBreaksIndex) + 1
            for bb = 1:length(fields)
                FiberData.(fields{1,bb}).raw = channelData(fileBreaksIndex(aa - 1) + 1:end,bb);
                FiberData.(fields{1,bb}).zScored = zScoredFiberData(fileBreaksIndex(aa - 1) + 1:end,bb);
            end
            FiberData.syncData = syncData(fileBreaksIndex(aa - 1) + 1:end);
            % Photobleaching removed and HbT corrected
            FiberData.RH_405.rescaledData = rescaledData(fileBreaksIndex(aa - 1) + 1:end,1);
            FiberData.RH_465.rescaledData = Corrected465.RH(fileBreaksIndex(aa - 1) + 1:end);
            FiberData.RH_560.rescaledData = rescaledData(fileBreaksIndex(aa - 1) + 1:end,3);
            FiberData.LH_405.rescaledData = rescaledData(fileBreaksIndex(aa - 1) + 1:end,4);
            FiberData.LH_465.rescaledData = Corrected465.LH(fileBreaksIndex(aa - 1) + 1:end);
            FiberData.LH_560.rescaledData = rescaledData(fileBreaksIndex(aa - 1) + 1:end,6);
            % HbT uncorrected corrected
            FiberData.RH_465.uncorrected = rescaledData(fileBreaksIndex(aa - 1) + 1:end,2);
            FiberData.LH_465.uncorrected = rescaledData(fileBreaksIndex(aa - 1) + 1:end,5);
            FiberData.RH_465.zScoredUncorrected = uncorrectedzScoredGCaMPData.RH(fileBreaksIndex(aa - 1) + 1:end);
            FiberData.LH_465.zScoredUncorrected = uncorrectedzScoredGCaMPData.LH(fileBreaksIndex(aa - 1) + 1:end);
            save([FiberData.notes.animalID '_' FiberData.notes.date '_' FiberData.notes.sessionNumber '_Trial'  num2str(aa) '_FiberData.mat'],'FiberData','SessionData.FitStruct','-v7.3');
        else
            for bb = 1:length(fields)
                FiberData.(fields{1,bb}).raw = channelData(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa),bb);
                FiberData.(fields{1,bb}).zScored = zScoredFiberData(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa),bb);
            end
            FiberData.syncData = syncData(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa));
            % Photobleaching removed and HbT corrected
            FiberData.RH_405.rescaledData = rescaledData(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa),1);
            FiberData.RH_465.rescaledData = Corrected465.RH(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa));
            FiberData.RH_560.rescaledData = rescaledData(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa),3);
            FiberData.LH_405.rescaledData = rescaledData(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa),4);
            FiberData.LH_465.rescaledData = Corrected465.LH(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa));
            FiberData.LH_560.rescaledData = rescaledData(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa),6);
            % HbT uncorrected corrected
            FiberData.RH_465.uncorrected = rescaledData(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa),2);
            FiberData.LH_465.uncorrected = rescaledData(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa),5);
            FiberData.RH_465.zScoredUncorrected = uncorrectedzScoredGCaMPData.RH(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa));
            FiberData.LH_465.zScoredUncorrected = uncorrectedzScoredGCaMPData.LH(fileBreaksIndex(aa - 1) + 1:fileBreaksIndex(aa));
            save([FiberData.notes.animalID '_' FiberData.notes.date '_' FiberData.notes.sessionNumber '_Trial'  num2str(aa) '_FiberData.mat'],'FiberData','SessionData.FitStruct','-v7.3');
        end
    end
end

end
