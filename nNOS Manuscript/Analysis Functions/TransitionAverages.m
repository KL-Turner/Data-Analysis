function [AnalysisResults] = Fig4_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
IOSanimalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
TwoP_animalIDs = {'T115','T117','T118','T125','T126'};   % T116 has no REM events
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
%% IOS mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,aa};
    for bb = 1:length(transitions)
        transition = transitions{1,bb};
        data.(transition).EMG(aa,:) = AnalysisResults.(animalID).Transitions.(transition).EMG;
        data.(transition).Hip(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).Hip;
        data.(transition).T = AnalysisResults.(animalID).Transitions.(transition).T;
        data.(transition).F = AnalysisResults.(animalID).Transitions.(transition).F;
        data.(transition).Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).Cort;
        data.(transition).HbT(aa,:) = AnalysisResults.(animalID).Transitions.(transition).HbT;
    end
end
% take average for each behavioral transition
for cc = 1:length(transitions)
    transition = transitions{1,cc};
    data.(transition).meanEMG = mean(data.(transition).EMG,1);
    data.(transition).stdEMG = std(data.(transition).EMG,0,1);
    data.(transition).meanHbT = mean(data.(transition).HbT,1);
    data.(transition).stdHbT = std(data.(transition).HbT,0,1);
    data.(transition).meanHip = mean(data.(transition).Hip,3)*100;
    data.(transition).meanCort = mean(data.(transition).Cort,3)*100;
end
T1 = -30 + (1/30):(1/30):30;
T2 = -30 + (1/10):(1/10):30;
%% two photon REM dilation and transition REM to Awake
% cd through each animal's directory and extract the appropriate analysis results
data.VesselTransitions.NREMtoREM.data = []; data.VesselTransitions.REMtoAwake.data = [];
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    evokedBehavFields = fieldnames(AnalysisResults.(animalID).Transitions);
    for bb = 1:length(evokedBehavFields)
        behavField = evokedBehavFields{bb,1};
        vesselIDs = fieldnames(AnalysisResults.(animalID).Transitions.(behavField));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            data.VesselTransitions.(behavField).data = vertcat(data.VesselTransitions.(behavField).data,medfilt1(AnalysisResults.(animalID).Transitions.(behavField).(vesselID).mean,10,'truncate'));
            data.VesselTransitions.(behavField).timeVector = AnalysisResults.(animalID).Transitions.(behavField).(vesselID).timeVector;
        end
    end
end
% take the average of the vessels for each behavior
evokedBehavFields = {'NREMtoREM','REMtoAwake'};
for dd = 1:length(evokedBehavFields)
    behavField = evokedBehavFields{1,dd};
    data.VesselTransitions.(behavField).mean = mean(data.VesselTransitions.(behavField).data,1);
    data.VesselTransitions.(behavField).StD = std(data.VesselTransitions.(behavField).data,0,1);
end
%% Fig. 4 (part two)
summaryFigure_B = figure('Name','Fig4 (f-g)');
sgtitle('Figure panel 4 - Turner et al. 2020')
%% [4f] NREM to REM transition
ax1 = subplot(1,2,1);
p1 = plot(T1,data.NREMtoREM.meanHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',2);
ylabel('\Delta[HbT] (\muM)')
ylim([40,80])
yyaxis right
p2 = plot(data.VesselTransitions.NREMtoREM.timeVector,data.VesselTransitions.NREMtoREM.mean,'color',colors_eLife2020('rich black'),'LineWidth',2);
title('[4f] NREM to REM transition')
xlabel('Time (s)')
ylabel('\DeltaD/D (%)','rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'IOS','2PLSM','Location','NorthWest')
xlim([-30,30])
ylim([5,42])
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ax1.YAxis(1).Color = colors_eLife2020('dark candy apple red');
ax1.YAxis(2).Color = colors_eLife2020('rich black');
%% [4g] REM to Awake transition
ax2 = subplot(1,2,2);
plot(T1,data.REMtoAWAKE.meanHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',2);
ylabel('\Delta[HbT] (\muM)')
ylim([0,80])
yyaxis right
plot(data.VesselTransitions.REMtoAwake.timeVector,data.VesselTransitions.REMtoAwake.mean,'color',colors_eLife2020('rich black'),'LineWidth',2)
title('[4g] REM to Awake transition')
xlabel('Time (s)')
ylabel('\DeltaD/D (%)','rotation',-90,'VerticalAlignment','bottom')
xlim([-30,30])
ylim([0,42])
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
ax2.YAxis(1).Color = colors_eLife2020('dark candy apple red');
ax2.YAxis(2).Color = colors_eLife2020('rich black');
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder '\Summary Figures and Structures\MATLAB Analysis Figures\'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_B,[dirpath 'Fig4_B']);
    set(summaryFigure_B,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig4_B'])
end
%% Fig. 4
summaryFigure_A = figure('Name','Fig4 (b-e)');
sgtitle('Figure 4 - Turner et al. 2020')
%% [4b] Awake to NREM
ax1 = subplot(6,2,1);
% HbT and EMG
p1 = plot(T1,data.AWAKEtoNREM.meanHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanHbT + data.AWAKEtoNREM.stdHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.AWAKEtoNREM.meanHbT - data.AWAKEtoNREM.stdHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',0.5);
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
ylim([-5,50])
yyaxis right
p2 = plot(T1,data.AWAKEtoNREM.meanEMG,'-','color',colors_eLife2020('rich black'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color',colors_eLife2020('rich black'),'LineWidth',0.5);
plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color',colors_eLife2020('rich black'),'LineWidth',0.5);
title('[4b] Awake to NREM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'HbT','EMG')
ax1.YAxis(1).Color = colors_eLife2020('dark candy apple red');
ax1.YAxis(2).Color = colors_eLife2020('rich black');
ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];
% cort neural
ax2 = subplot(6,2,3);
semilog_imagesc_eLife2020(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanCort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
% hippocampal neural
ax3 = subplot(6,2,5);
semilog_imagesc_eLife2020(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [4c] NREM to Awake
ax4 = subplot(6,2,2);
% HbT and EMG
plot(T1,data.NREMtoAWAKE.meanHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanHbT + data.NREMtoAWAKE.stdHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.NREMtoAWAKE.meanHbT - data.NREMtoAWAKE.stdHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',0.5);
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
ylim([-5,50])
yyaxis right
plot(T1,data.NREMtoAWAKE.meanEMG ,'-','color',colors_eLife2020('rich black'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color',colors_eLife2020('rich black'),'LineWidth',0.5);
plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color',colors_eLife2020('rich black'),'LineWidth',0.5);
title('[4c] NREM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax4.YAxis(1).Color = colors_eLife2020('dark candy apple red');
ax4.YAxis(2).Color = colors_eLife2020('rich black');
ylim([-1,0.5])
ax4.TickLength = [0.03,0.03];
% cort neural
ax5 = subplot(6,2,4);
semilog_imagesc_eLife2020(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanCort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
% hippocampal neural
ax6 = subplot(6,2,6);
semilog_imagesc_eLife2020(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [4d] NREM to REM
ax7 = subplot(6,2,7);
% HbT and EMG
plot(T1,data.NREMtoREM.meanHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.meanHbT + data.NREMtoREM.stdHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.NREMtoREM.meanHbT - data.NREMtoREM.stdHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',0.5);
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
ylim([35,90])
yyaxis right
plot(T1,data.NREMtoREM.meanEMG ,'-','color',colors_eLife2020('rich black'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.meanEMG + data.NREMtoREM.stdEMG,'-','color',colors_eLife2020('rich black'),'LineWidth',0.5);
plot(T1,data.NREMtoREM.meanEMG - data.NREMtoREM.stdEMG,'-','color',colors_eLife2020('rich black'),'LineWidth',0.5);
title('[4d] NREM to REM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax7.YAxis(1).Color = colors_eLife2020('dark candy apple red');
ax7.YAxis(2).Color = colors_eLife2020('rich black');
ylim([-2,-0.5])
ax7.TickLength = [0.03,0.03];
% cort neural
ax8 = subplot(6,2,9);
semilog_imagesc_eLife2020(T2,data.NREMtoREM.F,data.NREMtoREM.meanCort,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,300])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
% hippocampal neural
ax9 = subplot(6,2,11);
semilog_imagesc_eLife2020(T2,data.NREMtoREM.F,data.NREMtoREM.meanHip,'y')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,300])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [4e] REM to Awake
ax10 = subplot(6,2,8);
plot(T1,data.REMtoAWAKE.meanHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.meanHbT + data.REMtoAWAKE.stdHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.REMtoAWAKE.meanHbT - data.REMtoAWAKE.stdHbT,'-','color',colors_eLife2020('dark candy apple red'),'LineWidth',0.5);
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
ylim([0,90])
yyaxis right
plot(T1,data.REMtoAWAKE.meanEMG ,'-','color',colors_eLife2020('rich black'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.meanEMG + data.REMtoAWAKE.stdEMG,'-','color',colors_eLife2020('rich black'),'LineWidth',0.5);
plot(T1,data.REMtoAWAKE.meanEMG - data.REMtoAWAKE.stdEMG,'-','color',colors_eLife2020('rich black'),'LineWidth',0.5);
title('[4e] REM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax10.YAxis(1).Color = colors_eLife2020('dark candy apple red');
ax10.YAxis(2).Color = colors_eLife2020('rich black');
ylim([-2,1])
ax10.TickLength = [0.03,0.03];
% cort neural
ax11 = subplot(6,2,10);
semilog_imagesc_eLife2020(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanCort,'y')
axis xy
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,300])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
% hippocampal neural
ax12 = subplot(6,2,12);
semilog_imagesc_eLife2020(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanHip,'y')
c8 = colorbar;
ylabel(c8,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,300])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];
%% axes positionns
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax8Pos = get(ax8,'position');
ax9Pos = get(ax9,'position');
ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
ax12Pos = get(ax12,'position');
ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax4Pos(3:4);
ax6Pos(3:4) = ax4Pos(3:4);
ax8Pos(3:4) = ax7Pos(3:4);
ax9Pos(3:4) = ax7Pos(3:4);
ax11Pos(3:4) = ax10Pos(3:4);
ax12Pos(3:4) = ax10Pos(3:4);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax8,'position',ax8Pos);
set(ax9,'position',ax9Pos);
set(ax11,'position',ax11Pos);
set(ax12,'position',ax12Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_A,[dirpath 'Fig4_A']);
    set(summaryFigure_A,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4_A'])
end

end
