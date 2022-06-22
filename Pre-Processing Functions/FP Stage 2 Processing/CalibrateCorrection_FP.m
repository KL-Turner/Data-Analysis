function [RH_Coefficients,LH_Coefficients] = CalibrateCorrection_FP(fileName)
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

% identify the extension
extInd = strfind(fileName(1,:),'.');
% identify the underscores
fileBreaks = strfind(fileName(1,:),'_');
% file parameters
TempData.notes.animalname = fileName(1:fileBreaks(1) - 1);
TempData.notes.date = fileName(fileBreaks(1) + 1:fileBreaks(2) - 1);
TempData.notes.sessionnumber = fileName(fileBreaks(2) + 1:extInd - 1);
FiberData = csvread(fileName,2,0);
TempData.notes.DataSeconds_Data = FiberData(end,1); % total time length
% 2 = RH 405nm, 3 = RH 465nm, 5 = RH 560nm, 7 = LH 405nm, 8 = LH 465nm, 10 = LH 560nm
DataChannels = [2,3,5,7,8,10];
TempData.SyncData.LH560Raw = FiberData(:,11); % extract column 11 = LH 560 Raw for syncing the data with LabView
FiberData = FiberData(:,DataChannels);
TempData.notes.Acquisition_Fs = 1.2e4;
TempData.notes.Decimation = 10;
TempData.notes.DataAcquired = 55.50;
TempData.notes.DataSeconds = TempData.notes.DataAcquired*(60);
TempData.notes.DataFs = round(length(FiberData(:,1))/TempData.notes.DataSeconds_Data);
TempData.notes.VelocityChannel = 10;
TempData.notes.StartPad = 5*TempData.notes.DataFs;
TempData.notes.FollowPad = 15*TempData.notes.DataFs;
TempData.notes.Fit_Freq = 0.01;
TempData.notes.low_Freq = 1;
TempData.notes.Final_Freq = [0.01,1];
[z1,p1,k1] = butter(3,TempData.notes.Fit_Freq/(0.5*TempData.notes.DataFs),'low');
[sos_Fit,g_Fit] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(3,TempData.notes.low_Freq/(0.5*TempData.notes.DataFs),'low');
[sos_Low,g_Low] = zp2sos(z2,p2,k2);
TempData.notes.TrialLength = 15.50; % minutes
TempData.notes.IdleTime = 4.50; % minutes
TempData.notes.TrialDataSeconds = TempData.notes.DataFs*TempData.notes.TrialLength*60;
TempData.notes.IdleDataSeconds = TempData.notes.DataFs*TempData.notes.IdleTime*60;
TempData.notes.SingleTrialData = 20*60;
TrialsNumber = TempData.notes.DataSeconds_Data/TempData.notes.SingleTrialData; % how many 15.5 minutes recording
idleRemovedData = [];
for aa = 1:ceil(TrialsNumber)
   if aa == 1
       idleRemovedData = FiberData(1:TempData.notes.TrialDataSeconds,:);
   else
       idleRemovedData = cat(1,idleRemovedData,FiberData((aa - 1)*TempData.notes.TrialDataSeconds + (aa - 1)*TempData.notes.IdleDataSeconds:...
           ((aa)*TempData.notes.TrialDataSeconds + (aa - 1)*TempData.notes.IdleDataSeconds),:));
   end
end
FiberData =  idleRemovedData; % raw data after the idle time is removed
% remove photobleaching/metaloism baseline drift
Spacing = 1:1:length(FiberData(:,3));
FiltData = filtfilt(sos_Fit,g_Fit,FiberData);
%% RH correct TRITC blood volume
[fitVals] = fit(Spacing',FiltData(:,3),'exp2');
coeffVals = coeffvalues(fitVals);
predictedRHCBV = (coeffVals(1)*exp((coeffVals(2).*Spacing))) + (coeffVals(3)*exp((coeffVals(4).*Spacing)));
figTime = (1:length(FiberData(:,3)))/TempData.notes.DataFs;
figA = figure;
plot(figTime,FiberData(:,3),'r');
hold on;
plot(figTime,predictedRHCBV,'m');
plot(figTime,FiltData(:,3),'b');
title('Exponential fit of TRITC metabolism (RH)');
xlabel('Time (sec)');
legend({'Raw TRITC brightness','Exponential Fit','Low pass filtered data fit'});
xlim([0 figTime(end)]);
xticks(1:900:figTime(length(figTime)));
CorrectedRHCBV = FiberData(:,3) - predictedRHCBV';
%% RH correct Ca2+ dependent GCaMP
[fitVals]=fit(Spacing',FiltData(:,2),'exp2');
coeffVals=coeffvalues(fitVals);
predictedRHGCaMP=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing)));
figTime=(1:length(FiberData(:,2)))/TempData.notes.DataFs;
figB = figure;
plot(figTime,FiberData(:,2),'g');
hold on;
plot(figTime,predictedRHGCaMP,'k');
plot(figTime,FiltData(:,2),'r');
title('Exponential fit of Ca dependent GCaMP (RH)'); xlabel('Time (sec)'); legend({'Raw Ca dependent GCaMP brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
xticks(1:900:figTime(length(figTime)));
CorrectedRHGCaMP = FiberData(:,2) - predictedRHGCaMP';
%% RH correct Ca2+ independent GCaMP
[fitVals] = fit(Spacing',FiltData(:,1),'exp2');
coeffVals = coeffvalues(fitVals);
predictedRHGFP = (coeffVals(1)*exp((coeffVals(2).*Spacing))) + (coeffVals(3)*exp((coeffVals(4).*Spacing)));
figTime = (1:length(FiberData(:,1)))/TempData.notes.DataFs;
figC = figure;
plot(figTime,FiberData(:,1),'g');
hold on;
plot(figTime,predictedRHGFP,'k');
plot(figTime,FiltData(:,1),'r');
title('Exponential fit of GFP (RH)');
xlabel('Time (sec)');
legend({'Raw GFP brightness','Exponential Fit','Low pass filtered data fit'});
xlim([0,figTime(end)]);
xticks(1:900:figTime(length(figTime)));
CorrectedRHGFP = FiberData(:,1) - predictedRHGFP';
%% LH correct TRITC blood volume
[fitVals]=fit(Spacing',FiltData(:,6),'exp2');
coeffVals=coeffvalues(fitVals);
predictedLHCBV = (coeffVals(1)*exp((coeffVals(2).*Spacing))) + (coeffVals(3)*exp((coeffVals(4).*Spacing)));
figTime = (1:length(FiberData(:,6)))/TempData.notes.DataFs;
figD = figure;
plot(figTime,FiberData(:,6),'r');
hold on;
plot(figTime,predictedLHCBV,'m');
plot(figTime,FiltData(:,6),'b');
title('Exponential fit of TRITC metabolism (LH)');
xlabel('Time (sec)');
legend({'Raw TRITC brightness','Exponential Fit','Low pass filtered data fit'});
xlim([0 figTime(end)]);
xticks(1:900:figTime(length(figTime)));
CorrectedLHCBV = FiberData(:,6) -  predictedLHCBV';
%% LH correct Ca2+ dependent GCaMP
[fitVals]  = fit(Spacing',FiltData(:,5),'exp2');
coeffVals = coeffvalues(fitVals);
predictedLHGCaMP = (coeffVals(1)*exp((coeffVals(2).*Spacing))) + (coeffVals(3)*exp((coeffVals(4).*Spacing)));
figTime = (1:length(FiberData(:,5)))/TempData.notes.DataFs;
figE = figure;
plot(figTime,FiberData(:,5),'g');
hold on;
plot(figTime,predictedLHGCaMP,'k');
plot(figTime,FiltData(:,5),'r');
title('Exponential fit of Ca dependent GCaMP (LH)');
xlabel('Time (sec)');
legend({'Raw Ca dependent GCaMP brightness','Exponential Fit','Low pass filtered data fit'});
xlim([0,figTime(end)]);
xticks(1:900:figTime(length(figTime)));
CorrectedLHGCaMP = FiberData(:,5) - predictedLHGCaMP';
%% correct Ca2+ independent GCaMP
[fitVals] = fit(Spacing',FiltData(:,4),'exp2');
coeffVals = coeffvalues(fitVals);
predictedLHGFP = (coeffVals(1)*exp((coeffVals(2).*Spacing))) + (coeffVals(3)*exp((coeffVals(4).*Spacing)));
figTime = (1:length(FiberData(:,4)))/TempData.notes.DataFs;
figF = figure;
plot(figTime,FiberData(:,4),'g');
hold on;
plot(figTime,predictedLHGFP,'k');
plot(figTime,FiltData(:,4),'r');
title('Exponential fit of GFP (LH)');
xlabel('Time (sec)');
legend({'Raw GFP brightness','Exponential Fit','Low pass filtered data fit'});
xlim([0 figTime(end)]);
xticks(1:900:figTime(length(figTime)));
CorrectedLHGFP = FiberData(:,4) - predictedLHGFP';
DetrendData = [CorrectedRHGFP,CorrectedRHGCaMP,CorrectedRHCBV,CorrectedLHGFP,CorrectedLHGCaMP,CorrectedLHCBV];
LowPassData = filtfilt(sos_Low,g_Low,DetrendData);
for q = 1:size(LowPassData,2)
    RescaleData(:,q) = rescale(LowPassData(:,q),0,1); % rescale all data between 0 to 1
end
%% z-score data
AvgData = mean(RescaleData,1);
StdData = std(RescaleData,0,1);
AvgMatrix = repmat(AvgData,length(RescaleData),1);
StdMatrix = repmat(StdData,length(RescaleData),1);
ZscoredFiberData = (RescaleData - AvgMatrix)./StdMatrix;
%% hemodynamic response correction code using Mle
BinEdges = (0:0.001:1);
figG = figure;
DataHist = histogram2(RescaleData(:,3),RescaleData(:,2),BinEdges,BinEdges);
BinCounts = DataHist.BinCounts;
allCounts = sum(BinCounts,'all');
allHist = (BinCounts/allCounts)*100;
XCounts = sum(BinCounts,2);
mleFit.RH.allHist = allHist;
for binNum = 2:size(BinCounts,1)
    if binNum == 2
        percData(1) = (XCounts(1)/sum(XCounts))*100;
    end
    percData(binNum) = (sum(XCounts(1:binNum))/sum(XCounts))*100;
end
startInd = find(percData <= 2,1,'last') + 1; % exclude lower 2% of data points
endInd = find(percData >= 98,1,'first') - 1; % exclude upper 2% of data points
mleFit.RH.startInd = startInd;
mleFit.RH.endInd = endInd;
for binNum = startInd:endInd
    index = (binNum - startInd) + 1;
    leftEdge = BinEdges(binNum);
    rightEdge = BinEdges(binNum + 1);
    dataInds = RescaleData(:,3) >= leftEdge & RescaleData(:,3) <= rightEdge;
    theData = RescaleData(dataInds,1);
    [phat,pci] = mle(theData,'distrbution','normal');
    figure(101);
    normHist = histogram(theData,BinEdges);
    normCounts = normHist.BinCounts./sum(normHist.BinCounts);
    theFit = pdf('normal',normHist.BinEdges(1:(end - 1)),phat(1),phat(2));
    theFit =theFit./sum(theFit);
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
fitPeaks.RH(1:length(BinEdges)) = NaN;
fitposStan.RH(1:length(BinEdges)) = NaN;
fitnegStan.RH(1:length(BinEdges)) = NaN;
fitPeaks.RH(startInd:endInd) = mleFit.RH.fitnotes.phat(:,1);
fitposStan.RH(startInd:endInd) = mleFit.RH.fitnotes.phat(:,1) + mleFit.RH.fitnotes.phat(:,2);
fitnegStan.RH(startInd:endInd) = mleFit.RH.fitnotes.phat(:,1) - mleFit.RH.fitnotes.phat(:,2);
[RH_linFit] = fitlm(BinEdges(startInd:endInd)',mleFit.RH.fitnotes.phat(:,1),'RobustOpts','on');
mleFit.RH.linFit = RH_linFit;
%% Plot the fit
linSlope.RH = table2array(RH_linFit.Coefficients(2,1));
linInt.RH = table2array(RH_linFit.Coefficients(1,1));
linPlot.RH = linSlope.RH*BinEdges+linInt.RH;
mleFit.RH.linFit = RH_linFit;
figH = figure;
imagesc(BinEdges,BinEdges,allHist');
hold on;
caxis([0,0.002]);
colorbar('eastoutside');
axis xy;
ylabel('\Delta GCaMP (%)');
xlabel('\Delta TRITC (%)');
title(' \Delta TRITC vs \Delta GCaMP');
plot(BinEdges,fitPeaks.RH,'r','LineWidth',2);
plot(BinEdges,fitposStan.RH,'--r','LineWidth',1);
plot(BinEdges,fitnegStan.RH,'--r','LineWidth',1,'HandleVisibility','off');
plot(BinEdges,linPlot.RH,'w','LineWidth',2);
xlim([0,1]);
ylim([0,1]);
legend({'MLE Peaks','+/-1Std','linear fit of all data'},'TextColor','white','Location','northeast');
legend('boxoff');
y_int_RH = table2array(mleFit.RH.linFit.Coefficients(1,1));
the_slope_RH = table2array(mleFit.RH.linFit.Coefficients(2,1));
FinalGCaMP_RH = RescaleData(:,2) - (the_slope_RH*RescaleData(:,3) + y_int_RH);
FinalGCaMP_Z_RH = (FinalGCaMP_RH - mean(FinalGCaMP_RH))/std(FinalGCaMP_RH);
figI = figure;
plot(ZscoredFiberData(:,5));
hold on;
plot(FinalGCaMP_Z_RH);
plot(ZscoredFiberData(:,6));
legend({'Raw GCaMP','HbT Corrected GCaMP','TRITC blood plasma'});
%% LH
BinEdges = (0:0.001:1);
figJ = figure;
DataHist = histogram2(RescaleData(:,6),RescaleData(:,5),BinEdges,BinEdges);
BinCounts = DataHist.BinCounts;
allCounts = sum(BinCounts,'all');
allHist = (BinCounts/allCounts)*100;
XCounts = sum(BinCounts,2);
mleFit.LH.allHist = allHist;
for binNum = 2:size(BinCounts,1)
    if binNum == 2
        percData(1) = (XCounts(1)/sum(XCounts))*100;
    end
    percData(binNum) = (sum(XCounts(1:binNum))/sum(XCounts))*100;
end
startInd = find(percData <= 2,1,'last') + 1; % exclude lower 2% of data points
endInd = find(percData >= 98,1,'first') - 1; % exclude upper 2% of data points
mleFit.LH.startInd=startInd;
mleFit.LH.endInd=endInd;
for binNum = startInd:endInd
    index = (binNum-startInd) + 1;
    leftEdge = BinEdges(binNum);
    rightEdge = BinEdges(binNum + 1);
    dataInds = RescaleData(:,6) >= leftEdge & RescaleData(:,6)<=rightEdge;
    theData = RescaleData(dataInds,1);
    [phat,pci] = mle(theData,'distrbution','normal');
    figure(102);
    normHist = histogram(theData,BinEdges);
    normCounts = normHist.BinCounts./sum(normHist.BinCounts);
    theFit = pdf('normal',normHist.BinEdges(1:(end - 1)),phat(1),phat(2));
    theFit = theFit./sum(theFit);
    rsqr = 1 - (sum((normCounts - theFit).^2)/sum((normCounts - mean(normCounts)).^2));
    mleFit.LH.fitnotes.phat(index,:) = phat;
    mleFit.LH.fitnotes.pci(:,:,index) = pci;
    mleFit.LH.goodness.rsqr(index)  = rsqr;
    mleFit.LH.fitData.eGFPData{index} = theData;
    mleFit.LH.fitData.BinEdges(index,:) = [leftEdge,rightEdge];
    mleFit.LH.fitData.normHist(index,:) = normCounts;
    mleFit.LH.fitData.HistFit(index,:) = theFit;
    mleFit.LH.fitData.fitHistEdges(index,:) = normHist.BinEdges;
end
fitPeaks.LH(1:length(BinEdges)) = NaN;
fitposStan.LH(1:length(BinEdges)) = NaN;
fitnegStan.LH(1:length(BinEdges)) = NaN;
fitPeaks.LH(startInd:endInd) = mleFit.LH.fitnotes.phat(:,1);
fitposStan.LH(startInd:endInd) = mleFit.LH.fitnotes.phat(:,1) + mleFit.LH.fitnotes.phat(:,2);
fitnegStan.LH(startInd:endInd) = mleFit.LH.fitnotes.phat(:,1) - mleFit.LH.fitnotes.phat(:,2);
[LH_linFit] = fitlm(BinEdges(startInd:endInd)',mleFit.LH.fitnotes.phat(:,1),'RobustOpts','on');
mleFit.LH.linFit = LH_linFit;
%% Plot the fit
linSlope.LH = table2array(LH_linFit.Coefficients(2,1));
linInt.LH = table2array(LH_linFit.Coefficients(1,1));
linPlot.LH = linSlope.LH*BinEdges+linInt.LH;
mleFit.LH.linFit = LH_linFit;
figK = figure;
imagesc(BinEdges,BinEdges,allHist');
hold on;
caxis([0,0.002]);
colorbar('eastoutside');
axis xy;
ylabel('\Delta GCaMP (%)');
xlabel('\Delta TRITC (%)');
title(' \Delta TRITC vs \Delta GCaMP');
plot(BinEdges,fitPeaks.LH,'r','LineWidth',2);
plot(BinEdges,fitposStan.LH,'--r','LineWidth',1);
plot(BinEdges,fitnegStan.LH,'--r','LineWidth',1,'HandleVisibility','off');
plot(BinEdges,linPlot.LH,'w','LineWidth',2);
xlim([0,1]);
ylim([0,1]);
legend({'MLE Peaks','+/-1Std','linear fit of all data'},'TextColor','white','Location','northeast');
legend('boxoff');
y_int_LH = table2array(mleFit.LH.linFit.Coefficients(1,1));
the_slope_LH  = table2array(mleFit.LH.linFit.Coefficients(2,1));
FinalGCaMP_LH = RescaleData(:,5)  - (the_slope_LH*RescaleData(:,6)+y_int_LH);
FinalGCaMP_Z_LH = (FinalGCaMP_LH - mean( FinalGCaMP_LH))/std( FinalGCaMP_LH);
figL = figure;
plot(ZscoredFiberData(:,5));
hold on;
plot(FinalGCaMP_Z_LH);
plot(ZscoredFiberData(:,6));
legend({'Raw GCaMP','HbT Corrected GCaMP','TRITC blood plasma'});
figM = figure;
subplot(2,1,1)
plot(FinalGCaMP_Z_RH);
hold on;
plot(FinalGCaMP_Z_LH);
subplot(2,1,2)
plot(RescaleData(:,2));
hold on;
plot(RescaleData(:,5));
RH_Coefficients = table2array(mleFit.RH.linFit.Coefficients(:,1));
LH_Coefficients = table2array(mleFit.LH.linFit.Coefficients(:,1));
%% save figures
[pathstr,~,~] = fileparts(cd);
dirpath = [pathstr '/Figures/Calibration Correction/'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(figA,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_RH_Corrected_TRITC.fig']);
savefig(figB,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_RH_Dependent_GCaMP.fig']);
savefig(figC,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_RH_Independent_GCaMP.fig']);
savefig(figD,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_LH_Corrected_TRITC.fig']);
savefig(figE,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_LH_Dependent_GCaMP.fig']);
savefig(figF,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_LH_Independent_GCaMP.fig']);
savefig(figG,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_RH_3D_Histogram.fig']);
savefig(figH,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_RH_TRITC_GCaMP_Fit.fig']);
savefig(figI,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_RH_Normalized_Histogram_Fit.fig']);
savefig(figJ,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_LH_3D_Histogram.fig']);
savefig(figK,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_LH_TRITC_GCaMP_Fit.fig']);
savefig(figL,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_LH_Normalized_Histogram_Fit.fig']);
savefig(figM,[dirpath TempData.notes.animalname '_' TempData.notes.date '_' TempData.notes.sessionnumber '_Bilateral_GCaMP_Comparison.fig']);
close all

end
