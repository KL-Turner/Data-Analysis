%% RH
BinEdges=(0:0.001:1);
figure(100);
DataHist=histogram2(greenChannel,blueChannel,BinEdges,BinEdges);

BinCounts=DataHist.BinCounts;

allCounts=sum(BinCounts,'all');

allHist=(BinCounts/allCounts)*100;
XCounts=sum(BinCounts,2);
mleFit.RH.allHist=allHist;
for binNum=2:size(BinCounts,1)
    if binNum==2
        percData(1)=(XCounts(1)/sum(XCounts))*100;
    end
    percData(binNum)=(sum(XCounts(1:binNum))/sum(XCounts))*100;
end
startInd=find(percData<=2,1,'last')+1; %exclude lower 2% of data points
endInd=find(percData>=98,1,'first')-1; %exclude upper 2% of data points
mleFit.RH.startInd=startInd;
mleFit.RH.endInd=endInd;

for binNum=startInd:endInd
    index=(binNum-startInd)+1;
    leftEdge=BinEdges(binNum);
    rightEdge=BinEdges(binNum+1);
    
    dataInds=greenChannel>=leftEdge...
        & greenChannel<=rightEdge;
    theData = blueChannel(dataInds);
    [phat,pci]=mle(theData,'distrbution','normal');
    figure(605); normHist=histogram(theData,BinEdges);
    normCounts=normHist.BinCounts./sum(normHist.BinCounts);
    theFit=pdf('normal',normHist.BinEdges(1:(end-1)),phat(1),phat(2));
    theFit=theFit./sum(theFit);
    
    rsqr=1-(sum((normCounts-theFit).^2)/sum((normCounts-mean(normCounts)).^2));
    
    mleFit.RH.fitParams.phat(index,:)=phat;
    mleFit.RH.fitParams.pci(:,:,index)=pci;
    mleFit.RH.goodness.rsqr(index)=rsqr;
    mleFit.RH.fitData.eGFPData{index}=theData;
    mleFit.RH.fitData.BinEdges(index,:)=[leftEdge,rightEdge];
    mleFit.RH.fitData.normHist(index,:)=normCounts;
    mleFit.RH.fitData.HistFit(index,:)=theFit;
    mleFit.RH.fitData.fitHistEdges(index,:)=normHist.BinEdges;
end
fitPeaks.RH(1:length(BinEdges))=NaN;
fitposStan.RH(1:length(BinEdges))=NaN;
fitnegStan.RH(1:length(BinEdges))=NaN;
fitPeaks.RH(startInd:endInd)=mleFit.RH.fitParams.phat(:,1);
fitposStan.RH(startInd:endInd)=mleFit.RH.fitParams.phat(:,1)+mleFit.RH.fitParams.phat(:,2);
fitnegStan.RH(startInd:endInd)=mleFit.RH.fitParams.phat(:,1)-mleFit.RH.fitParams.phat(:,2);


[RH_linFit]=fitlm(BinEdges(startInd:endInd)',mleFit.RH.fitParams.phat(:,1),'RobustOpts','on');

mleFit.RH.linFit=RH_linFit;

%% Plot the fit
linSlope.RH=table2array(RH_linFit.Coefficients(2,1));
linInt.RH=table2array(RH_linFit.Coefficients(1,1));
slope_pVal.RH=table2array(RH_linFit.Coefficients(2,4));
linPlot.RH=linSlope.RH*BinEdges+linInt.RH;
mleFit.RH.linFit=RH_linFit;
figure(203); hold on;
imagesc(BinEdges,BinEdges,allHist');
caxis([0 0.002]);
h=colorbar('eastoutside');
axis xy;
ylabel('GFP ');
xlabel('HbT');
title(' GFP vs HbT');
plot(BinEdges,fitPeaks.RH,'r','LineWidth',2);
plot(BinEdges,fitposStan.RH,'--r','LineWidth',1);
plot(BinEdges,fitnegStan.RH,'--r','LineWidth',1,'HandleVisibility','off');
plot(BinEdges,linPlot.RH,'w','LineWidth',2);
xlim([0 1]);
ylim([0 1]);
legend({'MLE Peaks','+/-1Std','linear fit of all data'},'TextColor','white','Location','northeast');
legend('boxoff');
slopeTxt=num2str(linSlope.RH);
intTxt=num2str(linInt.RH);
pvalTxt=num2str(slope_pVal.RH);
rsqrTxt=num2str(RH_linFit.Rsquared.Adjusted);
all_an=annotation('textbox',[0.13, 0.64, 0.1, 0.1],'String',['Linear fit y=' slopeTxt 'x+' intTxt ' rsqr=' rsqrTxt],'FitBoxToText','on','LineStyle','none','Color',[1 1 1]);

coeffVals=table2array(mleFit.RH.linFit.Coefficients);

y_int_RH=table2array(mleFit.RH.linFit.Coefficients(1,1));
the_slope_RH=table2array(mleFit.RH.linFit.Coefficients(2,1));
FinalGCaMP_RH=blueChannel-(the_slope_RH*greenChannel+y_int_RH);
% FinalGCaMP_Z_RH=( FinalGCaMP_RH-mean( FinalGCaMP_RH))/std( FinalGCaMP_RH);

figure(49);plot(blueChannel);hold on; plot(FinalGCaMP_RH);plot(greenChannel);
legend({'Raw GFP','HbT Corrected GFP','530nm Reflectance'});