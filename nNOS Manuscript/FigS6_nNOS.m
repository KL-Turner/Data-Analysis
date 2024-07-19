function [] = FigS6_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Transitions_Ephys';
load(resultsStruct);
cd(rootFolder)
groups = {'Blank_SAP','SSP_SAP'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
%% IOS mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Transitions_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(transitions)
            transition = transitions{1,cc};
            % pre-allocate necessary variable fields
            data.(group).(transition).dummCheck = 1;
            if isfield(data.(group).(transition),'EMG') == false
                data.(group).(transition).EMG = [];
                data.(group).(transition).Hip = [];
                data.(group).(transition).T = [];
                data.(group).(transition).F = [];
                data.(group).(transition).Cort = [];
                data.(group).(transition).HbT = [];
            end
            data.(group).(transition).EMG = cat(1,data.(group).(transition).EMG,Results_Transitions_Ephys.(animalID).(transition).EMG);
            data.(group).(transition).Hip = cat(3,data.(group).(transition).Hip,Results_Transitions_Ephys.(animalID).(transition).Hip);
            data.(group).(transition).T = cat(1,data.(group).(transition).T,Results_Transitions_Ephys.(animalID).(transition).T);
            data.(group).(transition).F = cat(1,data.(group).(transition).F,Results_Transitions_Ephys.(animalID).(transition).F);
            data.(group).(transition).Cort = cat(3,data.(group).(transition).Cort,Results_Transitions_Ephys.(animalID).(transition).Cort);
            data.(group).(transition).HbT = cat(1,data.(group).(transition).HbT,Results_Transitions_Ephys.(animalID).(transition).HbT);
        end
    end
end
% take average for each behavioral transition
for qq = 1:length(groups)
    group = groups{1,qq};
    for cc = 1:length(transitions)
        transition = transitions{1,cc};
        data.(group).(transition).meanEMG = mean(data.(group).(transition).EMG,1);
        data.(group).(transition).stdEMG = std(data.(group).(transition).EMG,0,1);
        data.(group).(transition).meanHbT = mean(data.(group).(transition).HbT,1);
        data.(group).(transition).stdHbT = std(data.(group).(transition).HbT,0,1);
        data.(group).(transition).meanHip = mean(data.(group).(transition).Hip,3)*100;
        data.(group).(transition).meanCort = mean(data.(group).(transition).Cort,3)*100;
        data.(group).(transition).T = mean(data.(group).(transition).T,1);
        data.(group).(transition).F = mean(data.(group).(transition).F,1);
    end
end
T1 = -30 + (1/30):(1/30):30;
T2 = -30 + (1/10):(1/10):30;
%% Arousal-state
for aa = 1:length(groups)
    group = groups{1,aa};
    figName = ['summaryFigure' num2str(aa)]; %#ok<NASGU>
    figName = figure;
    sgtitle([group ' arousal-state transitions'])
    %% [4b] Awake to NREM
    ax1 = subplot(6,2,1);
    % HbT and EMG
    p1 = plot(T1,data.(group).AWAKEtoNREM.meanHbT,'-','color',colors('dark candy apple red'),'LineWidth',2);
    hold on
    plot(T1,data.(group).AWAKEtoNREM.meanHbT + data.(group).AWAKEtoNREM.stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    plot(T1,data.(group).AWAKEtoNREM.meanHbT - data.(group).AWAKEtoNREM.stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    ylabel('\Delta[HbT] (\muM)')
    xlim([-30,30])
    %     ylim([-5,50])
    yyaxis right
    p2 = plot(T1,data.(group).AWAKEtoNREM.meanEMG,'-','color',colors('rich black'),'LineWidth',2);
    hold on
    plot(T1,data.(group).AWAKEtoNREM.meanEMG + data.(group).AWAKEtoNREM.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
    plot(T1,data.(group).AWAKEtoNREM.meanEMG - data.(group).AWAKEtoNREM.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
    title('[4b] Awake to NREM transition')
    xlabel('Time (s)')
    ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
    set(gca,'box','off')
    legend([p1,p2],'HbT','EMG')
    ax1.YAxis(1).Color = colors('dark candy apple red');
    ax1.YAxis(2).Color = colors('rich black');
    %     ylim([-1,0.5])
    ax1.TickLength = [0.03,0.03];
    % cort neural
    ax2 = subplot(6,2,3);
    Semilog_ImageSC(T2,data.(group).AWAKEtoNREM.F,data.(group).AWAKEtoNREM.meanCort,'y')
    axis xy
    c1 = colorbar;
    ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    %     caxis([-100,200])
    xlabel('Time (s)')
    ylabel({'Cortical LFP';'Frequency (Hz)'})
    set(gca,'Yticklabel','10^1')
    xlim([-30,30])
    set(gca,'box','off')
    ax2.TickLength = [0.03,0.03];
    % hippocampal neural
    ax3 = subplot(6,2,5);
    Semilog_ImageSC(T2,data.(group).AWAKEtoNREM.F,data.(group).AWAKEtoNREM.meanHip,'y')
    c2 = colorbar;
    ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    %     caxis([-100,200])
    xlabel('Time (s)')
    ylabel({'Hippocampal LFP';'Frequency (Hz)'})
    set(gca,'Yticklabel','10^1')
    xlim([-30,30])
    set(gca,'box','off')
    ax3.TickLength = [0.03,0.03];
    %% [4c] NREM to Awake
    ax4 = subplot(6,2,2);
    % HbT and EMG
    plot(T1,data.(group).NREMtoAWAKE.meanHbT,'-','color',colors('dark candy apple red'),'LineWidth',2);
    hold on
    plot(T1,data.(group).NREMtoAWAKE.meanHbT + data.(group).NREMtoAWAKE.stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    plot(T1,data.(group).NREMtoAWAKE.meanHbT - data.(group).NREMtoAWAKE.stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    ylabel('\Delta[HbT] (\muM)')
    xlim([-30,30])
    %     ylim([-5,50])
    yyaxis right
    plot(T1,data.(group).NREMtoAWAKE.meanEMG ,'-','color',colors('rich black'),'LineWidth',2);
    hold on
    plot(T1,data.(group).NREMtoAWAKE.meanEMG + data.(group).NREMtoAWAKE.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
    plot(T1,data.(group).NREMtoAWAKE.meanEMG - data.(group).NREMtoAWAKE.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
    title('[4c] NREM to Awake transition')
    xlabel('Time (s)')
    ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
    set(gca,'box','off')
    ax4.YAxis(1).Color = colors('dark candy apple red');
    ax4.YAxis(2).Color = colors('rich black');
    %     ylim([-1,0.5])
    ax4.TickLength = [0.03,0.03];
    % cort neural
    ax5 = subplot(6,2,4);
    Semilog_ImageSC(T2,data.(group).NREMtoAWAKE.F,data.(group).NREMtoAWAKE.meanCort,'y')
    axis xy
    c3 = colorbar;
    ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    %     caxis([-100,200])
    xlabel('Time (s)')
    ylabel({'Cortical LFP';'Frequency (Hz)'})
    set(gca,'Yticklabel','10^1')
    xlim([-30,30])
    set(gca,'box','off')
    ax5.TickLength = [0.03,0.03];
    % hippocampal neural
    ax6 = subplot(6,2,6);
    Semilog_ImageSC(T2,data.(group).NREMtoAWAKE.F,data.(group).NREMtoAWAKE.meanHip,'y')
    c4 = colorbar;
    ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    %     caxis([-100,200])
    xlabel('Time (s)')
    ylabel({'Hippocampal LFP';'Frequency (Hz)'})
    set(gca,'Yticklabel','10^1')
    xlim([-30,30])
    set(gca,'box','off')
    ax6.TickLength = [0.03,0.03];
    %% [4d] NREM to REM
    ax7 = subplot(6,2,7);
    % HbT and EMG
    plot(T1,data.(group).NREMtoREM.meanHbT,'-','color',colors('dark candy apple red'),'LineWidth',2);
    hold on
    plot(T1,data.(group).NREMtoREM.meanHbT + data.(group).NREMtoREM.stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    plot(T1,data.(group).NREMtoREM.meanHbT - data.(group).NREMtoREM.stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    ylabel('\Delta[HbT] (\muM)')
    xlim([-30,30])
    %     ylim([35,90])
    yyaxis right
    plot(T1,data.(group).NREMtoREM.meanEMG ,'-','color',colors('rich black'),'LineWidth',2);
    hold on
    plot(T1,data.(group).NREMtoREM.meanEMG + data.(group).NREMtoREM.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
    plot(T1,data.(group).NREMtoREM.meanEMG - data.(group).NREMtoREM.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
    title('[4d] NREM to REM transition')
    xlabel('Time (s)')
    ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
    set(gca,'box','off')
    ax7.YAxis(1).Color = colors('dark candy apple red');
    ax7.YAxis(2).Color = colors('rich black');
    %     ylim([-2,-0.5])
    ax7.TickLength = [0.03,0.03];
    % cort neural
    ax8 = subplot(6,2,9);
    Semilog_ImageSC(T2,data.(group).NREMtoREM.F,data.(group).NREMtoREM.meanCort,'y')
    axis xy
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    %     caxis([-100,300])
    xlabel('Time (s)')
    ylabel({'Cortical LFP';'Frequency (Hz)'})
    set(gca,'Yticklabel','10^1')
    xlim([-30,30])
    set(gca,'box','off')
    ax8.TickLength = [0.03,0.03];
    % hippocampal neural
    ax9 = subplot(6,2,11);
    Semilog_ImageSC(T2,data.(group).NREMtoREM.F,data.(group).NREMtoREM.meanHip,'y')
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    %     caxis([-100,300])
    xlabel('Time (s)')
    ylabel({'Hippocampal LFP';'Frequency (Hz)'})
    set(gca,'Yticklabel','10^1')
    xlim([-30,30])
    set(gca,'box','off')
    ax9.TickLength = [0.03,0.03];
    %% [4e] REM to Awake
    ax10 = subplot(6,2,8);
    plot(T1,data.(group).REMtoAWAKE.meanHbT,'-','color',colors('dark candy apple red'),'LineWidth',2);
    hold on
    plot(T1,data.(group).REMtoAWAKE.meanHbT + data.(group).REMtoAWAKE.stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    plot(T1,data.(group).REMtoAWAKE.meanHbT - data.(group).REMtoAWAKE.stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    ylabel('\Delta[HbT] (\muM)')
    xlim([-30,30])
    %     ylim([0,90])
    yyaxis right
    plot(T1,data.(group).REMtoAWAKE.meanEMG ,'-','color',colors('rich black'),'LineWidth',2);
    hold on
    plot(T1,data.(group).REMtoAWAKE.meanEMG + data.(group).REMtoAWAKE.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
    plot(T1,data.(group).REMtoAWAKE.meanEMG - data.(group).REMtoAWAKE.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
    title('[4e] REM to Awake transition')
    xlabel('Time (s)')
    ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
    set(gca,'box','off')
    ax10.YAxis(1).Color = colors('dark candy apple red');
    ax10.YAxis(2).Color = colors('rich black');
    %     ylim([-2,1])
    ax10.TickLength = [0.03,0.03];
    % cort neural
    ax11 = subplot(6,2,10);
    Semilog_ImageSC(T2,data.(group).REMtoAWAKE.F,data.(group).REMtoAWAKE.meanCort,'y')
    axis xy
    c7 = colorbar;
    ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    %     caxis([-100,300])
    xlabel('Time (s)')
    ylabel({'Cortical LFP';'Frequency (Hz)'})
    set(gca,'Yticklabel','10^1')
    xlim([-30,30])
    set(gca,'box','off')
    ax11.TickLength = [0.03,0.03];
    % hippocampal neural
    ax12 = subplot(6,2,12);
    Semilog_ImageSC(T2,data.(group).REMtoAWAKE.F,data.(group).REMtoAWAKE.meanHip,'y')
    c8 = colorbar;
    ylabel(c8,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    %     caxis([-100,300])
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
    %     %% save figure(s)
    %     if saveFigs == true
    %         dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Arousal Transitions - Bilateral IOS' delim];
    %         if ~exist(dirpath,'dir')
    %             mkdir(dirpath);
    %         end
    %         savefig(figName,[dirpath group '_ArousalTransitions_IOS']);
    %         set(figName,'PaperPositionMode','auto');
    %         print('-painters','-dpdf','-fillpage',[dirpath group '_ArousalTransitions_IOS'])
    %     end
end
%% comparison of each behavior
summaryFigure4 = figure;
sgtitle('arousal-state transitions')
%% [4b] Awake to NREM
ax1 = subplot(2,2,1);
% HbT and EMG
p1 = plot(T1,data.Naive.AWAKEtoNREM.meanHbT,'-','color',colors('sapphire'),'LineWidth',1);
hold on
% plot(T1,data.Naive.AWAKEtoNREM.meanHbT + data.Naive.AWAKEtoNREM.stdHbT,'-','color',colors('sapphire'),'LineWidth',0.5);
% plot(T1,data.Naive.AWAKEtoNREM.meanHbT - data.Naive.AWAKEtoNREM.stdHbT,'-','color',colors('sapphire'),'LineWidth',0.5);
p2 = plot(T1,data.SSP_SAP.AWAKEtoNREM.meanHbT,'-','color',colors('electric purple'),'LineWidth',1);
% plot(T1,data.SSP_SAP.AWAKEtoNREM.meanHbT + data.SSP_SAP.AWAKEtoNREM.stdHbT,'-','color',colors('electric purple'),'LineWidth',0.5);
% plot(T1,data.SSP_SAP.AWAKEtoNREM.meanHbT - data.SSP_SAP.AWAKEtoNREM.stdHbT,'-','color',colors('electric purple'),'LineWidth',0.5);
p3 = plot(T1,data.Blank_SAP.AWAKEtoNREM.meanHbT,'-','color',colors('north texas green'),'LineWidth',1);
% plot(T1,data.Blank_SAP.AWAKEtoNREM.meanHbT + data.Blank_SAP.AWAKEtoNREM.stdHbT,'-','color',colors('north texas green'),'LineWidth',0.5);
% plot(T1,data.Blank_SAP.AWAKEtoNREM.meanHbT - data.Blank_SAP.AWAKEtoNREM.stdHbT,'-','color',colors('north texas green'),'LineWidth',0.5);
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
%     ylim([-5,50])
yyaxis right
plot(T1,data.Naive.AWAKEtoNREM.meanEMG,'--','color',colors('sapphire'),'LineWidth',1);
% plot(T1,data.Naive.AWAKEtoNREM.meanEMG + data.Naive.AWAKEtoNREM.stdEMG,'-','color',colors('sapphire'),'LineWidth',0.5);
% plot(T1,data.Naive.AWAKEtoNREM.meanEMG - data.Naive.AWAKEtoNREM.stdEMG,'-','color',colors('sapphire'),'LineWidth',0.5);
hold on
plot(T1,data.SSP_SAP.AWAKEtoNREM.meanEMG,'--','color',colors('electric purple'),'LineWidth',1);
% plot(T1,data.SSP_SAP.AWAKEtoNREM.meanEMG + data.SSP_SAP.AWAKEtoNREM.stdEMG,'-','color',colors('electric purple'),'LineWidth',0.5);
% plot(T1,data.SSP_SAP.AWAKEtoNREM.meanEMG - data.SSP_SAP.AWAKEtoNREM.stdEMG,'-','color',colors('electric purple'),'LineWidth',0.5);
plot(T1,data.Blank_SAP.AWAKEtoNREM.meanEMG,'--','color',colors('north texas green'),'LineWidth',1);
% plot(T1,data.Blank_SAP.AWAKEtoNREM.meanEMG + data.Blank_SAP.AWAKEtoNREM.stdEMG,'-','color',colors('north texas green'),'LineWidth',0.5);
% plot(T1,data.Blank_SAP.AWAKEtoNREM.meanEMG - data.Blank_SAP.AWAKEtoNREM.stdEMG,'-','color',colors('north texas green'),'LineWidth',0.5);
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2,p3],'Naive','SSP-SAP','Blank-SAP')
ax1.YAxis(1).Color = colors('rich black');
ax1.YAxis(2).Color = colors('rich black');
%     ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];
%% [4c] NREM to Awake
ax4 = subplot(2,2,2);
% HbT and EMG
plot(T1,data.Naive.NREMtoAWAKE.meanHbT,'-','color',colors('sapphire'),'LineWidth',1);
hold on
% plot(T1,data.Naive.NREMtoAWAKE.meanHbT + data.Naive.NREMtoAWAKE.stdHbT,'-','color',colors('sapphire'),'LineWidth',0.5);
% plot(T1,data.Naive.NREMtoAWAKE.meanHbT - data.Naive.NREMtoAWAKE.stdHbT,'-','color',colors('sapphire'),'LineWidth',0.5);
plot(T1,data.SSP_SAP.NREMtoAWAKE.meanHbT,'-','color',colors('electric purple'),'LineWidth',1);
% plot(T1,data.SSP_SAP.NREMtoAWAKE.meanHbT + data.SSP_SAP.NREMtoAWAKE.stdHbT,'-','color',colors('electric purple'),'LineWidth',0.5);
% plot(T1,data.SSP_SAP.NREMtoAWAKE.meanHbT - data.SSP_SAP.NREMtoAWAKE.stdHbT,'-','color',colors('electric purple'),'LineWidth',0.5);
plot(T1,data.Blank_SAP.NREMtoAWAKE.meanHbT,'-','color',colors('north texas green'),'LineWidth',1);
% plot(T1,data.Blank_SAP.NREMtoAWAKE.meanHbT + data.Blank_SAP.NREMtoAWAKE.stdHbT,'-','color',colors('north texas green'),'LineWidth',0.5);
% plot(T1,data.Blank_SAP.NREMtoAWAKE.meanHbT - data.Blank_SAP.NREMtoAWAKE.stdHbT,'-','color',colors('north texas green'),'LineWidth',0.5);
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
%     ylim([-5,50])
yyaxis right
plot(T1,data.Naive.NREMtoAWAKE.meanEMG,'--','color',colors('sapphire'),'LineWidth',1);
hold on
% plot(T1,data.Naive.NREMtoAWAKE.meanEMG + data.Naive.NREMtoAWAKE.stdEMG,'-','color',colors('sapphire'),'LineWidth',0.5);
% plot(T1,data.Naive.NREMtoAWAKE.meanEMG - data.Naive.NREMtoAWAKE.stdEMG,'-','color',colors('sapphire'),'LineWidth',0.5);
plot(T1,data.SSP_SAP.NREMtoAWAKE.meanEMG,'--','color',colors('electric purple'),'LineWidth',1);
% plot(T1,data.SSP_SAP.NREMtoAWAKE.meanEMG + data.SSP_SAP.NREMtoAWAKE.stdEMG,'-','color',colors('electric purple'),'LineWidth',0.5);
% plot(T1,data.SSP_SAP.NREMtoAWAKE.meanEMG - data.SSP_SAP.NREMtoAWAKE.stdEMG,'-','color',colors('electric purple'),'LineWidth',0.5);
plot(T1,data.Blank_SAP.NREMtoAWAKE.meanEMG,'--','color',colors('north texas green'),'LineWidth',1);
% plot(T1,data.Blank_SAP.NREMtoAWAKE.meanEMG + data.Blank_SAP.NREMtoAWAKE.stdEMG,'-','color',colors('north texas green'),'LineWidth',0.5);
% plot(T1,data.Blank_SAP.NREMtoAWAKE.meanEMG - data.Blank_SAP.NREMtoAWAKE.stdEMG,'-','color',colors('north texas green'),'LineWidth',0.5);
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax4.YAxis(1).Color = colors('rich black');
ax4.YAxis(2).Color = colors('rich black');
%     ylim([-1,0.5])
ax4.TickLength = [0.03,0.03];
%% [4d] NREM to REM
ax7 = subplot(2,2,3);
% HbT and EMG
plot(T1,data.Naive.NREMtoREM.meanHbT,'-','color',colors('sapphire'),'LineWidth',1);
hold on
% plot(T1,data.Naive.NREMtoREM.meanHbT + data.Naive.NREMtoREM.stdHbT,'-','color',colors('sapphire'),'LineWidth',0.5);
% plot(T1,data.Naive.NREMtoREM.meanHbT - data.Naive.NREMtoREM.stdHbT,'-','color',colors('sapphire'),'LineWidth',0.5);
plot(T1,data.SSP_SAP.NREMtoREM.meanHbT,'-','color',colors('electric purple'),'LineWidth',1);
% plot(T1,data.SSP_SAP.NREMtoREM.meanHbT + data.SSP_SAP.NREMtoREM.stdHbT,'-','color',colors('electric purple'),'LineWidth',0.5);
% plot(T1,data.SSP_SAP.NREMtoREM.meanHbT - data.SSP_SAP.NREMtoREM.stdHbT,'-','color',colors('electric purple'),'LineWidth',0.5);
plot(T1,data.Blank_SAP.NREMtoREM.meanHbT,'-','color',colors('north texas green'),'LineWidth',1);
% plot(T1,data.Blank_SAP.NREMtoREM.meanHbT + data.Blank_SAP.NREMtoREM.stdHbT,'-','color',colors('north texas green'),'LineWidth',0.5);
% plot(T1,data.Blank_SAP.NREMtoREM.meanHbT - data.Blank_SAP.NREMtoREM.stdHbT,'-','color',colors('north texas green'),'LineWidth',0.5);
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
%     ylim([-5,50])
yyaxis right
plot(T1,data.Naive.NREMtoREM.meanEMG,'--','color',colors('sapphire'),'LineWidth',1);
hold on
% plot(T1,data.Naive.NREMtoREM.meanEMG + data.Naive.NREMtoREM.stdEMG,'-','color',colors('sapphire'),'LineWidth',0.5);
% plot(T1,data.Naive.NREMtoREM.meanEMG - data.Naive.NREMtoREM.stdEMG,'-','color',colors('sapphire'),'LineWidth',0.5);
plot(T1,data.SSP_SAP.NREMtoREM.meanEMG,'--','color',colors('electric purple'),'LineWidth',1);
% plot(T1,data.SSP_SAP.NREMtoREM.meanEMG + data.SSP_SAP.NREMtoREM.stdEMG,'-','color',colors('electric purple'),'LineWidth',0.5);
% plot(T1,data.SSP_SAP.NREMtoREM.meanEMG - data.SSP_SAP.NREMtoREM.stdEMG,'-','color',colors('electric purple'),'LineWidth',0.5);
plot(T1,data.Blank_SAP.NREMtoREM.meanEMG,'--','color',colors('north texas green'),'LineWidth',1);
% plot(T1,data.Blank_SAP.NREMtoREM.meanEMG + data.Blank_SAP.NREMtoREM.stdEMG,'-','color',colors('north texas green'),'LineWidth',0.5);
% plot(T1,data.Blank_SAP.NREMtoREM.meanEMG - data.Blank_SAP.NREMtoREM.stdEMG,'-','color',colors('north texas green'),'LineWidth',0.5);
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax7.YAxis(1).Color = colors('rich black');
ax7.YAxis(2).Color = colors('rich black');
%     ylim([-2,-0.5])
ax7.TickLength = [0.03,0.03];
%% [4e] REM to Awake
ax10 = subplot(2,2,4);
plot(T1,data.Naive.REMtoAWAKE.meanHbT,'-','color',colors('sapphire'),'LineWidth',1);
hold on
% plot(T1,data.Naive.REMtoAWAKE.meanHbT + data.Naive.REMtoAWAKE.stdHbT,'-','color',colors('sapphire'),'LineWidth',0.5);
% plot(T1,data.Naive.REMtoAWAKE.meanHbT - data.Naive.REMtoAWAKE.stdHbT,'-','color',colors('sapphire'),'LineWidth',0.5);
plot(T1,data.SSP_SAP.REMtoAWAKE.meanHbT,'-','color',colors('electric purple'),'LineWidth',1);
% plot(T1,data.SSP_SAP.REMtoAWAKE.meanHbT + data.SSP_SAP.REMtoAWAKE.stdHbT,'-','color',colors('electric purple'),'LineWidth',0.5);
% plot(T1,data.SSP_SAP.REMtoAWAKE.meanHbT - data.SSP_SAP.REMtoAWAKE.stdHbT,'-','color',colors('electric purple'),'LineWidth',0.5);
plot(T1,data.Blank_SAP.REMtoAWAKE.meanHbT,'-','color',colors('north texas green'),'LineWidth',1);
% plot(T1,data.Blank_SAP.REMtoAWAKE.meanHbT + data.Blank_SAP.REMtoAWAKE.stdHbT,'-','color',colors('north texas green'),'LineWidth',0.5);
% plot(T1,data.Blank_SAP.REMtoAWAKE.meanHbT - data.Blank_SAP.REMtoAWAKE.stdHbT,'-','color',colors('north texas green'),'LineWidth',0.5);
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
%     ylim([-5,50])
yyaxis right
plot(T1,data.Naive.REMtoAWAKE.meanEMG,'--','color',colors('sapphire'),'LineWidth',1);
hold on
% plot(T1,data.Naive.REMtoAWAKE.meanEMG + data.Naive.REMtoAWAKE.stdEMG,'-','color',colors('sapphire'),'LineWidth',0.5);
% plot(T1,data.Naive.REMtoAWAKE.meanEMG - data.Naive.REMtoAWAKE.stdEMG,'-','color',colors('sapphire'),'LineWidth',0.5);
plot(T1,data.SSP_SAP.REMtoAWAKE.meanEMG,'--','color',colors('electric purple'),'LineWidth',1);
% plot(T1,data.SSP_SAP.REMtoAWAKE.meanEMG + data.SSP_SAP.REMtoAWAKE.stdEMG,'-','color',colors('electric purple'),'LineWidth',0.5);
% plot(T1,data.SSP_SAP.REMtoAWAKE.meanEMG - data.SSP_SAP.REMtoAWAKE.stdEMG,'-','color',colors('electric purple'),'LineWidth',0.5);
plot(T1,data.Blank_SAP.REMtoAWAKE.meanEMG,'--','color',colors('north texas green'),'LineWidth',1);
% plot(T1,data.Blank_SAP.REMtoAWAKE.meanEMG + data.Blank_SAP.REMtoAWAKE.stdEMG,'-','color',colors('north texas green'),'LineWidth',0.5);
% plot(T1,data.Blank_SAP.REMtoAWAKE.meanEMG - data.Blank_SAP.REMtoAWAKE.stdEMG,'-','color',colors('north texas green'),'LineWidth',0.5);
title('REM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax10.YAxis(1).Color = colors('rich black');
ax10.YAxis(2).Color = colors('rich black');
%     ylim([-2,1])
ax10.TickLength = [0.03,0.03];
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Arousal Transitions - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure4,[dirpath 'ArousalTransitions_Comparison']);
    set(summaryFigure4,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'ArousalTransitions_Comparison'])
end

end