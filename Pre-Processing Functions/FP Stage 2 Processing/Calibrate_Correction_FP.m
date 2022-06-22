function [RH_Coefficients,LH_Coefficients]=Calibrate_Correction_FP(filename)
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

% close all
% directory=cd;
% folderBreaks=strfind(directory,'\');
ChunkData.Params.animalname = 'T325';
ChunkData.Params.sessionnumber = 'Session1';
% if numel(folderBreaks)==5
%     ChunkData.Params.date=directory((folderBreaks(5)+1):end);
%     ChunkData.Params.GFP_Type=directory((folderBreaks(2)+1):(folderBreaks(3)-1));
% else
    ChunkData.Params.date = '210917';
%     ChunkData.Params.GFP_Type=directory((folderBreaks(2)+1):(folderBreaks(3)-1));
% end

%% Read .CSV
fileInf=detectImportOptions(filename);
RawData=csvread(filename,2,0);

ChunkData.Params.DataSeconds_Data = RawData(end,1); % total time length
DataChannels=[2,3,5,7,8,10];%  these are the three demodulated optical channels [2,3,5,7,8,10,11]
% 2 = RH 405nm, 3 = RH 465nm,5 = RH 560nm, 7 = LH 405nm, 8 = LH 465nm,10 = RH 560nm 

% WheelData=filtfilt(sos_ball,g_ball,detrend(RawData(:,ChunkData.Params.VelocityChannel))); %(5*ChunkData.Params.DataFs):(end-(5*ChunkData.Params.DataFs))
ChunkData.SyncData.LH560Raw = RawData(:,11); % Extract column 11 = LH 560 Raw for syncing the data with LabView
RawData=RawData(:,DataChannels); %(5*ChunkData.Params.DataFs):(end-(5*ChunkData.Params.DataFs))
%%
ChunkData.Params.Acquisition_Fs=1.2e4;
ChunkData.Params.Decimation=10;
ChunkData.Params.DataAcquired=55.50;%time in minutes of data acquired

ChunkData.Params.DataSeconds=ChunkData.Params.DataAcquired*(60);
ChunkData.Params.DataFs=round(length(RawData(:,1))/ChunkData.Params.DataSeconds_Data);%ChunkData.Params.Acquisition_Fs/ChunkData.Params.Decimation;
ChunkData.Params.VelocityChannel=10;
ChunkData.Params.StartPad=5*ChunkData.Params.DataFs;%Time to collect in front of running start
ChunkData.Params.FollowPad=15*ChunkData.Params.DataFs;%Time to follow running start
ChunkData.Params.Fit_Freq=0.01;
ChunkData.Params.low_Freq=1;
ChunkData.Params.Final_Freq=[0.01 1];
params.Fs=ChunkData.Params.DataFs; %multitaper estimation parameters
params.tapers=[3 5];%multitaper estimation parameters

[z,p,k]=butter(3,ChunkData.Params.Fit_Freq/(0.5*ChunkData.Params.DataFs),'low'); %design lowpass filter for hemodynamic correction
[sos_Fit,g_Fit]=zp2sos(z,p,k);

[z,p,k]=butter(3,ChunkData.Params.low_Freq/(0.5*ChunkData.Params.DataFs),'low'); %Low pass for optical data to physiologically relevant range
[sos_Low,g_Low]=zp2sos(z,p,k);

[z,p,k]=butter(3,ChunkData.Params.Final_Freq/(0.5*ChunkData.Params.DataFs),'bandpass'); %Low pass filter for locomotion data
[sos_final,g_final]=zp2sos(z,p,k);

[z,p,k]=butter(3,10/(0.5*ChunkData.Params.DataFs),'low'); %Low pass filter for locomotion data
[sos_ball,g_ball]=zp2sos(z,p,k);

%% Find Locomotion points to exclude from baseline calculations
%{
[imp_bin]=velocity_binarize_fiberphotometry(WheelData,ChunkData.Params.DataFs,ChunkData.Params.DataFs,1e-3);
FuseGaps=15*ChunkData.Params.DataFs;
bin_Run=double(imp_bin);
RunInds=find(bin_Run==1);
IndGap=diff(RunInds);
for IndNum=1:length(IndGap)
    if IndGap(IndNum)>1
        if IndGap(IndNum)<=FuseGaps
            bin_Run(RunInds(IndNum):RunInds(IndNum+1))=1;
        else
            bin_Run(RunInds(IndNum):(RunInds(IndNum)+FuseGaps))=1;
        end
    end
end
ExcludeVals=find(bin_Run==1);
%}

%% Remove idle period time of 4:30 minutes between two 15:30 minutes trials
%% Remove idle period time of 4:30 minutes between two 15:30 minutes trials
ChunkData.Params.TrialLength = 15.50; %minutes
ChunkData.Params.IdleTime = 4.50; %minutes
ChunkData.Params.TrialDataSeconds = ChunkData.Params.DataFs*ChunkData.Params.TrialLength*60;
ChunkData.Params.IdleDataSeconds = ChunkData.Params.DataFs*ChunkData.Params.IdleTime*60;
ChunkData.Params.SingleTrialData = 20*60;



TrialsNumber = ChunkData.Params.DataSeconds_Data / ChunkData.Params.SingleTrialData; % how many 15.5 minutes recording
    if TrialsNumber > 2 && TrialsNumber < 3 
    IdleremovedData = [RawData(1:ChunkData.Params.TrialDataSeconds,:);...
                       RawData((ChunkData.Params.TrialDataSeconds+ChunkData.Params.IdleDataSeconds+1):(ChunkData.Params.TrialDataSeconds*2+ChunkData.Params.IdleDataSeconds),:);...
                       RawData((ChunkData.Params.TrialDataSeconds*2+ChunkData.Params.IdleDataSeconds*2)+1:end,:)];
    elseif TrialsNumber > 3  && TrialsNumber < 4               

    IdleremovedData = [RawData(1:ChunkData.Params.TrialDataSeconds,:);...
                       RawData((ChunkData.Params.TrialDataSeconds+ChunkData.Params.IdleDataSeconds+1):(ChunkData.Params.TrialDataSeconds*2+ChunkData.Params.IdleDataSeconds),:);...
                       RawData((ChunkData.Params.TrialDataSeconds*2+ChunkData.Params.IdleDataSeconds*2)+1:(ChunkData.Params.TrialDataSeconds*3+ChunkData.Params.IdleDataSeconds*2),:)
                       RawData((ChunkData.Params.TrialDataSeconds*3+ChunkData.Params.IdleDataSeconds*3)+1:end,:)];    

    elseif TrialsNumber > 1  && TrialsNumber < 2 

    IdleremovedData = [RawData(1:ChunkData.Params.TrialDataSeconds,:);...
                       RawData((ChunkData.Params.TrialDataSeconds+ChunkData.Params.IdleDataSeconds+1):end,:)];
    end
    
CollectedRawData = RawData;% Raw Data without any manipulation
    
RawData =  IdleremovedData; % raw data after the idle time is removed
%% Remove photobleaching/metaloism baseline drift
Spacing=1:1:length(RawData(:,3));
FiltData=filtfilt(sos_Fit,g_Fit,RawData); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends

%% RH

    %%Correct TRITC blood volume
    [fitVals]=fit(Spacing',FiltData(:,3),'exp2');
    coeffVals=coeffvalues(fitVals);
    predictedRHCBV=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing))); 

        figTime=(1:length(RawData(:,3)))/ChunkData.Params.DataFs;
        figure(1);plot(figTime,RawData(:,3),'r'); hold on; plot(figTime,predictedRHCBV,'m'); plot(figTime,FiltData(:,3),'b');
        title('Exponential fit of TRITC metabolism (RH)'); xlabel('Time (sec)'); legend({'Raw TRITC brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
        xticks(1:900:figTime(length(figTime)));

    CorrectedRHCBV=RawData(:,3)-predictedRHCBV';
    FitStruct.RHCBV=fitVals;
    %% Correct Ca2+ dependent GCaMP
    [fitVals]=fit(Spacing',FiltData(:,2),'exp2');
    coeffVals=coeffvalues(fitVals);
    predictedRHGCaMP=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing))); 
        figTime=(1:length(RawData(:,2)))/ChunkData.Params.DataFs;
        figure(2);plot(figTime,RawData(:,2),'g'); hold on; plot(figTime,predictedRHGCaMP,'k'); plot(figTime,FiltData(:,2),'r');
        title('Exponential fit of Ca dependent GCaMP (RH)'); xlabel('Time (sec)'); legend({'Raw Ca dependent GCaMP brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
        xticks(1:900:figTime(length(figTime)));
    CorrectedRHGCaMP=RawData(:,2)-predictedRHGCaMP';
    FitStruct.RHGCaMP=fitVals;

    %% Correct Ca2+ independent GCaMP
    [fitVals]=fit(Spacing',FiltData(:,1),'exp2');
    coeffVals=coeffvalues(fitVals);
    predictedRHGFP=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing))); 
        figTime=(1:length(RawData(:,1)))/ChunkData.Params.DataFs;
        figure(3);plot(figTime,RawData(:,1),'g'); hold on; plot(figTime,predictedRHGFP,'k'); plot(figTime,FiltData(:,1),'r');
        title('Exponential fit of GFP (RH)'); xlabel('Time (sec)'); legend({'Raw GFP brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
        xticks(1:900:figTime(length(figTime)));
    CorrectedRHGFP=RawData(:,1)-predictedRHGFP';
    FitStruct.RHGFP=fitVals;

%% LH

    %Correct TRITC blood volume
    [fitVals]=fit(Spacing',FiltData(:,6),'exp2');
    coeffVals=coeffvalues(fitVals);
    predictedLHCBV=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing))); 
        figTime=(1:length(RawData(:,6)))/ChunkData.Params.DataFs;
        figure(4);plot(figTime,RawData(:,6),'r'); hold on; plot(figTime,predictedLHCBV,'m'); plot(figTime,FiltData(:,6),'b');
        title('Exponential fit of TRITC metabolism (LH)'); xlabel('Time (sec)'); legend({'Raw TRITC brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
        xticks(1:900:figTime(length(figTime)));
    CorrectedLHCBV=RawData(:,6)-predictedLHCBV';
    FitStruct.LHCBV=fitVals;
    %% Correct Ca2+ dependent GCaMP
    [fitVals]=fit(Spacing',FiltData(:,5),'exp2');
    coeffVals=coeffvalues(fitVals);
    predictedLHGCaMP=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing))); 
        figTime=(1:length(RawData(:,5)))/ChunkData.Params.DataFs;
        figure(5);plot(figTime,RawData(:,5),'g'); hold on; plot(figTime,predictedLHGCaMP,'k'); plot(figTime,FiltData(:,5),'r');
        title('Exponential fit of Ca dependent GCaMP (LH)'); xlabel('Time (sec)'); legend({'Raw Ca dependent GCaMP brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
        xticks(1:900:figTime(length(figTime)));
    CorrectedLHGCaMP=RawData(:,5)-predictedLHGCaMP';
    FitStruct.LHGCaMP=fitVals;
    %% Correct Ca2+ independent GCaMP
    [fitVals]=fit(Spacing',FiltData(:,4),'exp2');
    coeffVals=coeffvalues(fitVals);
    predictedLHGFP=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing))); 
        figTime=(1:length(RawData(:,4)))/ChunkData.Params.DataFs;
        figure(6);plot(figTime,RawData(:,4),'g'); hold on; plot(figTime,predictedLHGFP,'k'); plot(figTime,FiltData(:,4),'r');
        title('Exponential fit of GFP (LH)'); xlabel('Time (sec)'); legend({'Raw GFP brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
        xticks(1:900:figTime(length(figTime)));
    CorrectedLHGFP=RawData(:,4)-predictedLHGFP';
    FitStruct.LHGFP=fitVals;
%%
DetrendData=[CorrectedRHGFP,CorrectedRHGCaMP,CorrectedRHCBV,CorrectedLHGFP,CorrectedLHGCaMP,CorrectedLHCBV];
LowPassData=filtfilt(sos_Low,g_Low,DetrendData);

for q=1:size(LowPassData,2)
RescaleData(:,q)=rescale(LowPassData(:,q),0,1); %rescale all data between 0 to 1
end
%% Z-score data
AvgData=mean(RescaleData,1);
StdData=std(RescaleData,0,1);
AvgMatrix=repmat(AvgData,length(RescaleData),1);
StdMatrix=repmat(StdData,length(RescaleData),1);
ZscoredFiberData=(RescaleData-AvgMatrix)./StdMatrix;
%% Hemodynamic response correction code using Mle
    %% RH
            BinEdges=(0:0.001:1);
            figure(100);
            DataHist=histogram2(RescaleData(:,3),RescaleData(:,2),BinEdges,BinEdges);

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

                    dataInds=RescaleData(:,3)>=leftEdge...
                        & RescaleData(:,3)<=rightEdge;
                    theData=RescaleData(dataInds,1);
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
                figure(203);subplot(3,3,1); hold on;
                imagesc(BinEdges,BinEdges,allHist');
                caxis([0 0.002]);
                h=colorbar('eastoutside');
                axis xy;
                ylabel('\Delta GCaMP (%)');
                xlabel('\Delta TRITC (%)');
                title(' \Delta TRITC vs \Delta GCaMP');
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
                FinalGCaMP_RH=RescaleData(:,2)-(the_slope_RH*RescaleData(:,3)+y_int_RH);
                FinalGCaMP_Z_RH=( FinalGCaMP_RH-mean( FinalGCaMP_RH))/std( FinalGCaMP_RH);

                figure(49);plot(ZscoredFiberData(:,5));hold on; plot( FinalGCaMP_Z_RH);plot(ZscoredFiberData(:,6));
                legend({'Raw GCaMP','HbT Corrected GCaMP','TRITC blood plasma'});
                savefig([ChunkData.Params.animalname '_' ChunkData.Params.date '_' ChunkData.Params.sessionnumber '_RH' '_NormalizedHistogramFit.fig']);
    %% LH

            BinEdges=(0:0.001:1);
            figure(101);
            DataHist=histogram2(RescaleData(:,6),RescaleData(:,5),BinEdges,BinEdges);

            BinCounts=DataHist.BinCounts;

            allCounts=sum(BinCounts,'all');

            allHist=(BinCounts/allCounts)*100;
            XCounts=sum(BinCounts,2);
            mleFit.LH.allHist=allHist;
            for binNum=2:size(BinCounts,1)
                if binNum==2
                    percData(1)=(XCounts(1)/sum(XCounts))*100;
                end
                percData(binNum)=(sum(XCounts(1:binNum))/sum(XCounts))*100;
            end
            startInd=find(percData<=2,1,'last')+1; %exclude lower 2% of data points
            endInd=find(percData>=98,1,'first')-1; %exclude upper 2% of data points
            mleFit.LH.startInd=startInd;
            mleFit.LH.endInd=endInd;

            for binNum=startInd:endInd
                index=(binNum-startInd)+1;
                leftEdge=BinEdges(binNum);
                rightEdge=BinEdges(binNum+1);

                    dataInds=RescaleData(:,6)>=leftEdge...
                        & RescaleData(:,6)<=rightEdge;
                    theData=RescaleData(dataInds,1);
                [phat,pci]=mle(theData,'distrbution','normal');
                figure(605); normHist=histogram(theData,BinEdges);
                normCounts=normHist.BinCounts./sum(normHist.BinCounts);
                theFit=pdf('normal',normHist.BinEdges(1:(end-1)),phat(1),phat(2));
                theFit=theFit./sum(theFit);

                rsqr=1-(sum((normCounts-theFit).^2)/sum((normCounts-mean(normCounts)).^2));

                mleFit.LH.fitParams.phat(index,:)=phat;
                mleFit.LH.fitParams.pci(:,:,index)=pci;
                mleFit.LH.goodness.rsqr(index)=rsqr;
                mleFit.LH.fitData.eGFPData{index}=theData;
                mleFit.LH.fitData.BinEdges(index,:)=[leftEdge,rightEdge];
                mleFit.LH.fitData.normHist(index,:)=normCounts;
                mleFit.LH.fitData.HistFit(index,:)=theFit;
                mleFit.LH.fitData.fitHistEdges(index,:)=normHist.BinEdges;
            end
                fitPeaks.LH(1:length(BinEdges))=NaN;
                fitposStan.LH(1:length(BinEdges))=NaN;
                fitnegStan.LH(1:length(BinEdges))=NaN;
                fitPeaks.LH(startInd:endInd)=mleFit.LH.fitParams.phat(:,1);
                fitposStan.LH(startInd:endInd)=mleFit.LH.fitParams.phat(:,1)+mleFit.LH.fitParams.phat(:,2);
                fitnegStan.LH(startInd:endInd)=mleFit.LH.fitParams.phat(:,1)-mleFit.LH.fitParams.phat(:,2);


            [LH_linFit]=fitlm(BinEdges(startInd:endInd)',mleFit.LH.fitParams.phat(:,1),'RobustOpts','on');

            mleFit.LH.linFit=LH_linFit;
            
            %% Plot the fit
                linSlope.LH=table2array(LH_linFit.Coefficients(2,1));
                linInt.LH=table2array(LH_linFit.Coefficients(1,1));
                slope_pVal.LH=table2array(LH_linFit.Coefficients(2,4));
                linPlot.LH=linSlope.LH*BinEdges+linInt.LH;
                mleFit.LH.linFit=LH_linFit;
                figure(204);subplot(3,3,1); hold on;
                imagesc(BinEdges,BinEdges,allHist');
                caxis([0 0.002]);
                h=colorbar('eastoutside');
                axis xy;
                ylabel('\Delta GCaMP (%)');
                xlabel('\Delta TRITC (%)');
                title(' \Delta TRITC vs \Delta GCaMP');
                plot(BinEdges,fitPeaks.LH,'r','LineWidth',2);
                plot(BinEdges,fitposStan.LH,'--r','LineWidth',1);
                plot(BinEdges,fitnegStan.LH,'--r','LineWidth',1,'HandleVisibility','off');
                plot(BinEdges,linPlot.LH,'w','LineWidth',2);
                xlim([0 1]);
                ylim([0 1]);
                legend({'MLE Peaks','+/-1Std','linear fit of all data'},'TextColor','white','Location','northeast');
                legend('boxoff');
                slopeTxt=num2str(linSlope.LH);
                intTxt=num2str(linInt.LH);
                pvalTxt=num2str(slope_pVal.LH);
                rsqrTxt=num2str(LH_linFit.Rsquared.Adjusted);
                all_an=annotation('textbox',[0.13, 0.64, 0.1, 0.1],'String',['Linear fit y=' slopeTxt 'x+' intTxt ' rsqr=' rsqrTxt],'FitBoxToText','on','LineStyle','none','Color',[1 1 1]);

%                 coeffVals=table2array(mleFit.LH.linFit.Coefficients);

                y_int_LH=table2array(mleFit.LH.linFit.Coefficients(1,1));
                the_slope_LH=table2array(mleFit.LH.linFit.Coefficients(2,1));
                FinalGCaMP_LH=RescaleData(:,5)-(the_slope_LH*RescaleData(:,6)+y_int_LH);
                FinalGCaMP_Z_LH=( FinalGCaMP_LH-mean( FinalGCaMP_LH))/std( FinalGCaMP_LH);
                
                
                figure(50);plot(ZscoredFiberData(:,5));hold on; plot( FinalGCaMP_Z_LH);plot(ZscoredFiberData(:,6));
                legend({'Raw GCaMP','HbT Corrected GCaMP','TRITC blood plasma'});
                savefig([ChunkData.Params.animalname '_' ChunkData.Params.date '_' ChunkData.Params.sessionnumber '_LH' '_NormalizedHistogramFit.fig']);
                
                figure(55);plot(FinalGCaMP_Z_RH);hold on; plot( FinalGCaMP_Z_LH);
                
                 figure(56);plot(RescaleData(:,2));hold on; plot( RescaleData(:,5));
                
                
                RH_Coefficients = table2array(mleFit.RH.linFit.Coefficients(:,1));
                LH_Coefficients = table2array(mleFit.LH.linFit.Coefficients(:,1));
                close all
end