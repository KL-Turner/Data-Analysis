function [AnalysisResults] = Fig3_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel 3 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
colorRfcAwake = [(0/256),(64/256),(64/256)];
colorRfcNREM = [(0/256),(174/256),(239/256)];
colorRfcREM = [(190/256),(30/256),(45/256)];
% colorRest = [(0/256),(166/256),(81/256)];
% colorWhisk = [(31/256),(120/256),(179/256)];
% colorStim = [(255/256),(28/256),(206/256)];
% colorNREM = [(191/256),(0/256),(255/256)];
% colorREM = [(254/256),(139/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% information and data for example
if isfield(AnalysisResults,'ExampleTrials') == false
    AnalysisResults.ExampleTrials = [];
end
if isfield(AnalysisResults.ExampleTrials,'T118A') == true
    dsFs = AnalysisResults.ExampleTrials.T118A.dsFs;
    p2Fs = AnalysisResults.ExampleTrials.T118A.p2Fs;
    filtEMG = AnalysisResults.ExampleTrials.T118A.filtEMG;
    filtForceSensor = AnalysisResults.ExampleTrials.T118A.filtForceSensor;
    filtWhiskerAngle = AnalysisResults.ExampleTrials.T118A.filtWhiskerAngle;
    filtVesselDiameter = AnalysisResults.ExampleTrials.T118A.filtVesselDiameter;
    T = AnalysisResults.ExampleTrials.T118A.T;
    F = AnalysisResults.ExampleTrials.T118A.F;
    cortNormS = AnalysisResults.ExampleTrials.T118A.cortNormS;
    hipNormS = AnalysisResults.ExampleTrials.T118A.hipNormS;
else
    animalID = 'T118';
    dataLocation = [rootFolder '\' animalID '\2P Data\'];
    cd(dataLocation)
    exampleMergedFileID = 'T118_RH_191211_14_44_20_015_A1_MergedData.mat';
    load(exampleMergedFileID,'-mat')
    exampleSpecFileID = 'T118_RH_191211_14_44_20_015_A1_SpecData.mat';
    load(exampleSpecFileID,'-mat')
    exampleBaselinesFileID = 'T118_RestingBaselines.mat';
    load(exampleBaselinesFileID,'-mat')
    [~,~,fileDate,~,~,vesselID] = GetFileInfo2_2P_eLife2020(exampleMergedFileID);
    strDay = ConvertDate_2P_eLife2020(fileDate);
    dsFs = MergedData.notes.dsFs;
    p2Fs = MergedData.notes.p2Fs;
    % setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
    [z1,p1,k1] = butter(4,10/(MergedData.notes.dsFs/2),'low');
    [sos1,g1] = zp2sos(z1,p1,k1);
    [z2,p2,k2] = butter(4,0.5/(MergedData.notes.p2Fs/2),'low');
    [sos2,g2] = zp2sos(z2,p2,k2);
    % whisker angle
    filtWhiskerAngle = filtfilt(sos1,g1,MergedData.data.whiskerAngle);
    % pressure sensor
    filtForceSensor = filtfilt(sos1,g1,abs(MergedData.data.forceSensorL));
    % EMG
    EMG = MergedData.data.EMG.data;
    normEMG = EMG - RestingBaselines.manualSelection.EMG.data.(strDay);
    filtEMG = filtfilt(sos1,g1,normEMG);
    % vessel diameter
    vesselDiameter = MergedData.data.vesselDiameter.data;
    normVesselDiameter = (vesselDiameter - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay))./(RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay));
    filtVesselDiameter = filtfilt(sos2,g2,normVesselDiameter)*100;
    % cortical and hippocampal spectrograms
    cortNormS = SpecData.corticalNeural.fiveSec.normS.*100;
    hipNormS = SpecData.hippocampalNeural.fiveSec.normS.*100;
    T = SpecData.corticalNeural.fiveSec.T;
    F = SpecData.corticalNeural.fiveSec.F;
    % update analysis structure
    AnalysisResults.ExampleTrials.T118A.dsFs = dsFs;
    AnalysisResults.ExampleTrials.T118A.p2Fs = p2Fs;
    AnalysisResults.ExampleTrials.T118A.filtEMG = filtEMG;
    AnalysisResults.ExampleTrials.T118A.filtForceSensor = filtForceSensor;
    AnalysisResults.ExampleTrials.T118A.filtWhiskerAngle = filtWhiskerAngle;
    AnalysisResults.ExampleTrials.T118A.filtVesselDiameter = filtVesselDiameter;
    AnalysisResults.ExampleTrials.T118A.T = T;
    AnalysisResults.ExampleTrials.T118A.F = F;
    AnalysisResults.ExampleTrials.T118A.cortNormS = cortNormS;
    AnalysisResults.ExampleTrials.T118A.hipNormS = hipNormS;
    % save results
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end
%% Fig. 3
summaryFigure = figure('Name','Fig3 (e-i)');
sgtitle('Figure 3 - Turner et al. 2020')
%% [3e-i] 2PLSM sleep example
% EMG and force sensor
ax1 = subplot(6,1,1);
p1 = plot((1:length(filtEMG))/dsFs,filtEMG,'color',colors_eLife2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','power (a.u.)'})
ylim([-2,2.5])
yyaxis right
p2 = plot((1:length(filtForceSensor))/dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'EMG','Pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
xlim([15,615])
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_eLife2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% whisker angle
ax2 = subplot(6,1,2);
plot((1:length(filtWhiskerAngle))/dsFs,-filtWhiskerAngle,'color',colors_eLife2020('rich black'),'LineWidth',0.5)
ylabel({'Whisker','angle (deg)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([15,75,135,195,255,315,375,435,495,555,615])
ax2.TickLength = [0.01,0.01];
xlim([15,615])
ylim([-10,40])
% vessel diameter
ax34 = subplot(6,1,[3,4]);
p3 = plot((1:length(filtVesselDiameter))/p2Fs,filtVesselDiameter,'color',colors_eLife2020('dark candy apple red'),'LineWidth',1);
hold on
xline(15,'color',colorRfcAwake,'LineWidth',2);
x2 = xline(55,'color',colorRfcNREM,'LineWidth',2);
x1 = xline(97,'color',colorRfcAwake,'LineWidth',2);
xline(105,'color',colorRfcNREM,'LineWidth',2);
x3 = xline(156,'color',colorRfcREM,'LineWidth',2);
xline(224,'color',colorRfcAwake,'LineWidth',2);
xline(248,'color',colorRfcNREM,'LineWidth',2);
xline(342,'color',colorRfcAwake,'LineWidth',2);
xline(360,'color',colorRfcNREM,'LineWidth',2);
xline(450,'color',colorRfcREM,'LineWidth',2);
xline(537,'color',colorRfcAwake,'LineWidth',2);
ylabel('\DeltaD/D (%)')
legend([p3,x1,x2,x3],'Arteriole diameter','Awake','NREM','REM')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
ax34.TickLength = [0.01,0.01];
xlim([15,615])
ylim([-30,60])
% cortical LFP
ax5 = subplot(6,1,5);
semilog_imagesc_eLife2020(T,F,cortNormS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'Cort LFP','Freq (Hz)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
ax5.TickLength = [0.01,0.01];
xlim([15,615])
% hippocampal LFP
ax6 = subplot(6,1,6);
semilog_imagesc_eLife2020(T,F,hipNormS,'y')
axis xy
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hipp LFP','Freq (Hz)'})
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
ax6.TickLength = [0.01,0.01];
xlim([15,615])
% axes properties
ax1Pos = get(ax1,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig3']);
    % remove surface subplots because they take forever to render
    cla(ax5);
    set(ax5,'YLim',[1,99]);
    cla(ax6);
    set(ax6,'YLim',[1,99]);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig3'])
    close(summaryFigure)
    %% subplot figures
    subplotImgs = figure;
    % example 1 cortical LFP
    subplot(2,1,1);
    semilog_imagesc_eLife2020(T,F,cortNormS,'y')
    caxis([-100,200])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([15,615])
    % example 1 hippocampal LFP
    subplot(2,1,2);
    semilog_imagesc_eLife2020(T,F,hipNormS,'y')
    caxis([-100,100])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([15,615])
    print('-painters','-dtiffn',[dirpath 'Fig3_SpecImages'])
    close(subplotImgs)
    %% Fig. S3
    figure('Name','Fig3 (e-i)');
    sgtitle('Figure Panel 3 - Turner et al. 2020')
    %% [3e-i] 2PLSM sleep example
    % EMG and force sensor
    ax1 = subplot(6,1,1);
    p1 = plot((1:length(filtEMG))/dsFs,filtEMG,'color',colors_eLife2020('rich black'),'LineWidth',0.5);
    ylabel({'EMG','power (a.u.)'})
    ylim([-2,2.5])
    yyaxis right
    p2 = plot((1:length(filtForceSensor))/dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
    ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
    legend([p1,p2],'EMG','Pressure')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    xticks([15,75,135,195,255,315,375,435,495,555,615])
    xlim([15,615])
    ylim([-0.1,2.5])
    ax1.TickLength = [0.01,0.01];
    ax1.YAxis(1).Color = colors_eLife2020('rich black');
    ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
    % whisker angle
    ax2 = subplot(6,1,2);
    plot((1:length(filtWhiskerAngle))/dsFs,-filtWhiskerAngle,'color',colors_eLife2020('rich black'),'LineWidth',0.5)
    ylabel({'Whisker','angle (deg)'})
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([15,75,135,195,255,315,375,435,495,555,615])
    ax2.TickLength = [0.01,0.01];
    xlim([15,615])
    ylim([-10,40])
    % vessel diameter
    ax34 = subplot(6,1,[3,4]);
    p3 = plot((1:length(filtVesselDiameter))/p2Fs,filtVesselDiameter,'color',colors_eLife2020('dark candy apple red'),'LineWidth',1);
    hold on
    xline(15,'color',colorRfcAwake,'LineWidth',2);
    x2 = xline(55,'color',colorRfcNREM,'LineWidth',2);
    x1 = xline(97,'color',colorRfcAwake,'LineWidth',2);
    xline(105,'color',colorRfcNREM,'LineWidth',2);
    x3 = xline(156,'color',colorRfcREM,'LineWidth',2);
    xline(224,'color',colorRfcAwake,'LineWidth',2);
    xline(248,'color',colorRfcNREM,'LineWidth',2);
    xline(342,'color',colorRfcAwake,'LineWidth',2);
    xline(360,'color',colorRfcNREM,'LineWidth',2);
    xline(450,'color',colorRfcREM,'LineWidth',2);
    xline(537,'color',colorRfcAwake,'LineWidth',2);
    ylabel('\DeltaD/D (%)')
    legend([p3,x1,x2,x3],'Arteriole diameter','Awake','NREM','REM')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    xticks([15,75,135,195,255,315,375,435,495,555,615])
    ax34.TickLength = [0.01,0.01];
    xlim([15,615])
    ylim([-30,60])
    % cortical LFP
    ax5 = subplot(6,1,5);
    semilog_imagesc_eLife2020(T,F,cortNormS,'y')
    axis xy
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,200])
    ylabel({'Cort LFP','Freq (Hz)'})
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    xticks([15,75,135,195,255,315,375,435,495,555,615])
    ax5.TickLength = [0.01,0.01];
    xlim([15,615])
    % hippocampal LFP
    ax6 = subplot(6,1,6);
    semilog_imagesc_eLife2020(T,F,hipNormS,'y')
    axis xy
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,100])
    xlabel('Time (min)')
    ylabel({'Hipp LFP','Freq (Hz)'})
    set(gca,'box','off')
    axis tight
    xticks([15,75,135,195,255,315,375,435,495,555,615])
    xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    ax6.TickLength = [0.01,0.01];
    xlim([15,615])
    % axes properties
    ax1Pos = get(ax1,'position');
    ax5Pos = get(ax5,'position');
    ax6Pos = get(ax6,'position');
    ax5Pos(3:4) = ax1Pos(3:4);
    ax6Pos(3:4) = ax1Pos(3:4);
    set(ax5,'position',ax5Pos);
    set(ax6,'position',ax6Pos);
end

end
