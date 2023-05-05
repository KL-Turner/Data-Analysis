function [] = ArousalTransitions_GCaMP_SI_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Transitions_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
dataTypes = {'HbT','HbO','HbR','GCaMP','EMG','T','F','Cort','Hip'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Transitions_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(transitions)
                transition = transitions{1,dd};
                data.(group).(hemisphere).(transition).dummCheck = 1;
                for ee = 1:length(dataTypes)
                    dataType = dataTypes{1,ee};
                    if isfield(data.(group).(hemisphere).(transition),dataType) == false
                        data.(group).(hemisphere).(transition).(dataType) = [];
                    end
                    % concatenate data across animals
                    if isempty(Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).(dataType)) == false
                        if any(strcmp(dataType,{'Cort','Hip'})) == true
                            data.(group).(hemisphere).(transition).(dataType) = cat(3,data.(group).(hemisphere).(transition).(dataType),Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).(dataType)*100);
                        else
                            data.(group).(hemisphere).(transition).(dataType) = cat(1,data.(group).(hemisphere).(transition).(dataType),Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).(dataType));
                        end
                    end
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(groups)
        group = groups{1,bb};
        for cc = 1:length(transitions)
            transition = transitions{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if any(strcmp(dataType,{'Cort','Hip'})) == true
                    data.(group).(hemisphere).(transition).(['mean_' dataType]) = mean(data.(group).(hemisphere).(transition).(dataType),3);
                else
                    data.(group).(hemisphere).(transition).(['mean_' dataType]) = mean(data.(group).(hemisphere).(transition).(dataType),1);
                    data.(group).(hemisphere).(transition).(['stdErr_' dataType]) = std(data.(group).(hemisphere).(transition).(dataType),1)./sqrt(size(data.(group).(hemisphere).(transition).(dataType),1));
                end
            end
        end
    end
end
% figure
T0 = -30 + (1/10):(1/10):30;
T1 = -30 + (1/30):(1/30):30;
T2 = -30 + (1/10):(1/10):30;
hemoDataTypes = {'HbT','HbO','HbR','GCaMP'};
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(hemoDataTypes)
            hemoDataType = hemoDataTypes{1,cc};
            summaryFigure = figure;
            sgtitle([strrep(group,'_',' ') ' ' hemisphere ' ' hemoDataType ' arousal-state transitions [GCaMP SI]'])
            %% Awake to NREM
            ax1 = subplot(6,2,1);
            % HbT and EMG
            p1 = plot(T0,data.(group).(hemisphere).AWAKEtoNREM.(['mean_' hemoDataType]),'-','color',colors('dark candy apple red'),'LineWidth',2);
            hold on
            plot(T0,data.(group).(hemisphere).AWAKEtoNREM.(['mean_' hemoDataType]) + data.(group).(hemisphere).AWAKEtoNREM.(['stdErr_' hemoDataType]),'-','color',colors('dark candy apple red'),'LineWidth',0.25);
            plot(T0,data.(group).(hemisphere).AWAKEtoNREM.(['mean_' hemoDataType]) - data.(group).(hemisphere).AWAKEtoNREM.(['stdErr_' hemoDataType]),'-','color',colors('dark candy apple red'),'LineWidth',0.25);
            if strcmp(hemoDataType,'GCaMP') == true
                ylabel('\DeltaF/F (%)')
            else
                ylabel(['\Delta[' hemoDataType '] (\muM)'])
            end
            xlim([-30,30])
            yyaxis right
            p2 = plot(T1,data.(group).(hemisphere).AWAKEtoNREM.mean_EMG,'-','color',colors('rich black'),'LineWidth',2);
            hold on
            plot(T1,data.(group).(hemisphere).AWAKEtoNREM.mean_EMG + data.(group).(hemisphere).AWAKEtoNREM.stdErr_EMG,'-','color',colors('rich black'),'LineWidth',0.25);
            plot(T1,data.(group).(hemisphere).AWAKEtoNREM.mean_EMG - data.(group).(hemisphere).AWAKEtoNREM.stdErr_EMG,'-','color',colors('rich black'),'LineWidth',0.25);
            title('Awake to NREM')
            xlabel('Time (s)')
            ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
            set(gca,'box','off')
            legend([p1,p2],'HbT','EMG')
            ax1.YAxis(1).Color = colors('dark candy apple red');
            ax1.YAxis(2).Color = colors('rich black');
            ylim([-1,1])
            ax1.TickLength = [0.03,0.03];
            % cort neural
            ax2 = subplot(6,2,3);
            Semilog_ImageSC(T2,data.(group).(hemisphere).AWAKEtoNREM.mean_F,data.(group).(hemisphere).AWAKEtoNREM.mean_Cort,'y')
            axis xy
            c1 = colorbar;
            ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
            caxis([-100,100])
            xlabel('Time (s)')
            ylabel({'Cortical LFP';'Frequency (Hz)'})
            set(gca,'Yticklabel','10^1')
            xlim([-30,30])
            set(gca,'box','off')
            ax2.TickLength = [0.03,0.03];
            % hippocampal neural
            ax3 = subplot(6,2,5);
            Semilog_ImageSC(T2,data.(group).(hemisphere).AWAKEtoNREM.mean_F,data.(group).(hemisphere).AWAKEtoNREM.mean_Hip,'y')
            c2 = colorbar;
            ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
            caxis([-100,100])
            xlabel('Time (s)')
            ylabel({'Hippocampal LFP';'Frequency (Hz)'})
            set(gca,'Yticklabel','10^1')
            xlim([-30,30])
            set(gca,'box','off')
            ax3.TickLength = [0.03,0.03];
            %% NREM to Awake
            ax4 = subplot(6,2,2);
            % HbT and EMG
            plot(T0,data.(group).(hemisphere).NREMtoAWAKE.(['mean_' hemoDataType]),'-','color',colors('dark candy apple red'),'LineWidth',2);
            hold on
            plot(T0,data.(group).(hemisphere).NREMtoAWAKE.(['mean_' hemoDataType]) + data.(group).(hemisphere).NREMtoAWAKE.(['stdErr_' hemoDataType]),'-','color',colors('dark candy apple red'),'LineWidth',0.25);
            plot(T0,data.(group).(hemisphere).NREMtoAWAKE.(['mean_' hemoDataType]) - data.(group).(hemisphere).NREMtoAWAKE.(['stdErr_' hemoDataType]),'-','color',colors('dark candy apple red'),'LineWidth',0.25);
            if strcmp(hemoDataType,'GCaMP') == true
                ylabel('\DeltaF/F (%)')
            else
                ylabel(['\Delta[' hemoDataType '] (\muM)'])
            end
            xlim([-30,30])
            yyaxis right
            plot(T1,data.(group).(hemisphere).NREMtoAWAKE.mean_EMG ,'-','color',colors('rich black'),'LineWidth',2);
            hold on
            plot(T1,data.(group).(hemisphere).NREMtoAWAKE.mean_EMG + data.(group).(hemisphere).NREMtoAWAKE.stdErr_EMG,'-','color',colors('rich black'),'LineWidth',0.25);
            plot(T1,data.(group).(hemisphere).NREMtoAWAKE.mean_EMG - data.(group).(hemisphere).NREMtoAWAKE.stdErr_EMG,'-','color',colors('rich black'),'LineWidth',0.25);
            title('NREM to Awake')
            xlabel('Time (s)')
            ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
            set(gca,'box','off')
            ax4.YAxis(1).Color = colors('dark candy apple red');
            ax4.YAxis(2).Color = colors('rich black');
            ylim([-1,1])
            ax4.TickLength = [0.03,0.03];
            % cort neural
            ax5 = subplot(6,2,4);
            Semilog_ImageSC(T2,data.(group).(hemisphere).NREMtoAWAKE.mean_F,data.(group).(hemisphere).NREMtoAWAKE.mean_Cort,'y')
            axis xy
            c3 = colorbar;
            ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
            caxis([-100,100])
            xlabel('Time (s)')
            ylabel({'Cortical LFP';'Frequency (Hz)'})
            set(gca,'Yticklabel','10^1')
            xlim([-30,30])
            set(gca,'box','off')
            ax5.TickLength = [0.03,0.03];
            % hippocampal neural
            ax6 = subplot(6,2,6);
            Semilog_ImageSC(T2,data.(group).(hemisphere).NREMtoAWAKE.mean_F,data.(group).(hemisphere).NREMtoAWAKE.mean_Hip,'y')
            c4 = colorbar;
            ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
            caxis([-100,100])
            xlabel('Time (s)')
            ylabel({'Hippocampal LFP';'Frequency (Hz)'})
            set(gca,'Yticklabel','10^1')
            xlim([-30,30])
            set(gca,'box','off')
            ax6.TickLength = [0.03,0.03];
            % NREM to REM
            ax7 = subplot(6,2,7);
            % HbT and EMG
            plot(T0,data.(group).(hemisphere).NREMtoREM.(['mean_' hemoDataType]),'-','color',colors('dark candy apple red'),'LineWidth',2);
            hold on
            plot(T0,data.(group).(hemisphere).NREMtoREM.(['mean_' hemoDataType]) + data.(group).(hemisphere).NREMtoREM.(['stdErr_' hemoDataType]),'-','color',colors('dark candy apple red'),'LineWidth',0.25);
            plot(T0,data.(group).(hemisphere).NREMtoREM.(['mean_' hemoDataType]) - data.(group).(hemisphere).NREMtoREM.(['stdErr_' hemoDataType]),'-','color',colors('dark candy apple red'),'LineWidth',0.25);
            if strcmp(hemoDataType,'GCaMP') == true
                ylabel('\DeltaF/F (%)')
            else
                ylabel(['\Delta[' hemoDataType '] (\muM)'])
            end
            xlim([-30,30])
            yyaxis right
            plot(T1,data.(group).(hemisphere).NREMtoREM.mean_EMG ,'-','color',colors('rich black'),'LineWidth',2);
            hold on
            plot(T1,data.(group).(hemisphere).NREMtoREM.mean_EMG + data.(group).(hemisphere).NREMtoREM.stdErr_EMG,'-','color',colors('rich black'),'LineWidth',0.25);
            plot(T1,data.(group).(hemisphere).NREMtoREM.mean_EMG - data.(group).(hemisphere).NREMtoREM.stdErr_EMG,'-','color',colors('rich black'),'LineWidth',0.25);
            title('NREM to REM')
            xlabel('Time (s)')
            ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
            set(gca,'box','off')
            ax7.YAxis(1).Color = colors('dark candy apple red');
            ax7.YAxis(2).Color = colors('rich black');
            ylim([-1.75,-0.5])
            ax7.TickLength = [0.03,0.03];
            % cort neural
            ax8 = subplot(6,2,9);
            Semilog_ImageSC(T2,data.(group).(hemisphere).NREMtoREM.mean_F,data.(group).(hemisphere).NREMtoREM.mean_Cort,'y')
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
            Semilog_ImageSC(T2,data.(group).(hemisphere).NREMtoREM.mean_F,data.(group).(hemisphere).NREMtoREM.mean_Hip,'y')
            c6 = colorbar;
            ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
            caxis([-100,300])
            xlabel('Time (s)')
            ylabel({'Hippocampal LFP';'Frequency (Hz)'})
            set(gca,'Yticklabel','10^1')
            xlim([-30,30])
            set(gca,'box','off')
            ax9.TickLength = [0.03,0.03];
            % REM to Awake
            ax10 = subplot(6,2,8);
            plot(T0,data.(group).(hemisphere).REMtoAWAKE.(['mean_' hemoDataType]),'-','color',colors('dark candy apple red'),'LineWidth',2);
            hold on
            plot(T0,data.(group).(hemisphere).REMtoAWAKE.(['mean_' hemoDataType]) + data.(group).(hemisphere).REMtoAWAKE.(['stdErr_' hemoDataType]),'-','color',colors('dark candy apple red'),'LineWidth',0.25);
            plot(T0,data.(group).(hemisphere).REMtoAWAKE.(['mean_' hemoDataType]) - data.(group).(hemisphere).REMtoAWAKE.(['stdErr_' hemoDataType]),'-','color',colors('dark candy apple red'),'LineWidth',0.25);
            if strcmp(hemoDataType,'GCaMP') == true
                ylabel('\DeltaF/F (%)')
            else
                ylabel(['\Delta[' hemoDataType '] (\muM)'])
            end
            xlim([-30,30])
            yyaxis right
            plot(T1,data.(group).(hemisphere).REMtoAWAKE.mean_EMG ,'-','color',colors('rich black'),'LineWidth',2);
            hold on
            plot(T1,data.(group).(hemisphere).REMtoAWAKE.mean_EMG + data.(group).(hemisphere).REMtoAWAKE.stdErr_EMG,'-','color',colors('rich black'),'LineWidth',0.25);
            plot(T1,data.(group).(hemisphere).REMtoAWAKE.mean_EMG - data.(group).(hemisphere).REMtoAWAKE.stdErr_EMG,'-','color',colors('rich black'),'LineWidth',0.25);
            title('REM to Awake')
            xlabel('Time (s)')
            ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
            set(gca,'box','off')
            ax10.YAxis(1).Color = colors('dark candy apple red');
            ax10.YAxis(2).Color = colors('rich black');
            ylim([-2,1.25])
            ax10.TickLength = [0.03,0.03];
            % cort neural
            ax11 = subplot(6,2,10);
            Semilog_ImageSC(T2,data.(group).(hemisphere).REMtoAWAKE.mean_F,data.(group).(hemisphere).REMtoAWAKE.mean_Cort,'y')
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
            Semilog_ImageSC(T2,data.(group).(hemisphere).REMtoAWAKE.mean_F,data.(group).(hemisphere).REMtoAWAKE.mean_Hip,'y')
            c8 = colorbar;
            ylabel(c8,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
            caxis([-100,300])
            xlabel('Time (s)')
            ylabel({'Hippocampal LFP';'Frequency (Hz)'})
            set(gca,'Yticklabel','10^1')
            xlim([-30,30])
            set(gca,'box','off')
            ax12.TickLength = [0.03,0.03];
            % axes positionns
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
            % save figure(s)
            if saveFigs == true
                dirpath = [rootFolder delim 'Summary Figures' delim 'Arousal Transitions' delim];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(summaryFigure,[dirpath 'ArousalTransitions_GCaMP_SI_' group '_' hemisphere '_' hemoDataType]);
            end
        end
    end
end
% figure
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(hemoDataTypes)
        hemoDataType = hemoDataTypes{1,bb};
        summaryFigure = figure;
        sgtitle([hemisphere ' ' hemoDataType ' arousal-state transitions [GCaMP SI]'])
        %% Awake to NREM
        ax1 = subplot(2,2,1);
        % HbT and EMG
        p2 = plot(T0,data.SSP_SAP.(hemisphere).AWAKEtoNREM.(['mean_' hemoDataType]),'-','color',colors('electric purple'),'LineWidth',2);
        hold on;
        plot(T0,data.SSP_SAP.(hemisphere).AWAKEtoNREM.(['mean_' hemoDataType]) + data.SSP_SAP.(hemisphere).AWAKEtoNREM.(['stdErr_' hemoDataType]),'-','color',colors('electric purple'),'LineWidth',0.25);
        plot(T0,data.SSP_SAP.(hemisphere).AWAKEtoNREM.(['mean_' hemoDataType]) - data.SSP_SAP.(hemisphere).AWAKEtoNREM.(['stdErr_' hemoDataType]),'-','color',colors('electric purple'),'LineWidth',0.25);
        p3 = plot(T0,data.Blank_SAP.(hemisphere).AWAKEtoNREM.(['mean_' hemoDataType]),'-','color',colors('north texas green'),'LineWidth',2);
        plot(T0,data.Blank_SAP.(hemisphere).AWAKEtoNREM.(['mean_' hemoDataType]) + data.Blank_SAP.(hemisphere).AWAKEtoNREM.(['stdErr_' hemoDataType]),'-','color',colors('north texas green'),'LineWidth',0.25);
        plot(T0,data.Blank_SAP.(hemisphere).AWAKEtoNREM.(['mean_' hemoDataType]) - data.Blank_SAP.(hemisphere).AWAKEtoNREM.(['stdErr_' hemoDataType]),'-','color',colors('north texas green'),'LineWidth',0.25);
        if strcmp(hemoDataType,'GCaMP') == true
            ylabel('\DeltaF/F (%)')
        else
            ylabel(['\Delta[' hemoDataType '] (\muM)'])
        end
        xlim([-30,30])
        title('Awake to NREM transition')
        xlabel('Time (s)')
        set(gca,'box','off')
        legend([p2,p3],'SSP-SAP','Blank-SAP')
        ax1.TickLength = [0.03,0.03];
        %% NREM to Awake
        ax4 = subplot(2,2,2);
        % HbT and EMG
        plot(T0,data.SSP_SAP.(hemisphere).NREMtoAWAKE.(['mean_' hemoDataType]),'-','color',colors('electric purple'),'LineWidth',2);
        hold on;
        plot(T0,data.SSP_SAP.(hemisphere).NREMtoAWAKE.(['mean_' hemoDataType]) + data.SSP_SAP.(hemisphere).NREMtoAWAKE.(['stdErr_' hemoDataType]),'-','color',colors('electric purple'),'LineWidth',0.25);
        plot(T0,data.SSP_SAP.(hemisphere).NREMtoAWAKE.(['mean_' hemoDataType]) - data.SSP_SAP.(hemisphere).NREMtoAWAKE.(['stdErr_' hemoDataType]),'-','color',colors('electric purple'),'LineWidth',0.25);
        plot(T0,data.Blank_SAP.(hemisphere).NREMtoAWAKE.(['mean_' hemoDataType]),'-','color',colors('north texas green'),'LineWidth',2);
        plot(T0,data.Blank_SAP.(hemisphere).NREMtoAWAKE.(['mean_' hemoDataType]) + data.Blank_SAP.(hemisphere).NREMtoAWAKE.(['stdErr_' hemoDataType]),'-','color',colors('north texas green'),'LineWidth',0.25);
        plot(T0,data.Blank_SAP.(hemisphere).NREMtoAWAKE.(['mean_' hemoDataType]) - data.Blank_SAP.(hemisphere).NREMtoAWAKE.(['stdErr_' hemoDataType]),'-','color',colors('north texas green'),'LineWidth',0.25);
        title('NREM to Awake')
        if strcmp(hemoDataType,'GCaMP') == true
            ylabel('\DeltaF/F (%)')
        else
            ylabel(['\Delta[' hemoDataType '] (\muM)'])
        end
        xlim([-30,30])
        xlabel('Time (s)')
        set(gca,'box','off')
        ax4.TickLength = [0.03,0.03];
        %% NREM to REM
        ax7 = subplot(2,2,3);
        % HbT and EMG
        plot(T0,data.SSP_SAP.(hemisphere).NREMtoREM.(['mean_' hemoDataType]),'-','color',colors('electric purple'),'LineWidth',2);
        hold on;
        plot(T0,data.SSP_SAP.(hemisphere).NREMtoREM.(['mean_' hemoDataType]) + data.SSP_SAP.(hemisphere).NREMtoREM.(['stdErr_' hemoDataType]),'-','color',colors('electric purple'),'LineWidth',0.25);
        plot(T0,data.SSP_SAP.(hemisphere).NREMtoREM.(['mean_' hemoDataType]) - data.SSP_SAP.(hemisphere).NREMtoREM.(['stdErr_' hemoDataType]),'-','color',colors('electric purple'),'LineWidth',0.25);
        plot(T0,data.Blank_SAP.(hemisphere).NREMtoREM.(['mean_' hemoDataType]),'-','color',colors('north texas green'),'LineWidth',2);
        plot(T0,data.Blank_SAP.(hemisphere).NREMtoREM.(['mean_' hemoDataType]) + data.Blank_SAP.(hemisphere).NREMtoREM.(['stdErr_' hemoDataType]),'-','color',colors('north texas green'),'LineWidth',0.25);
        plot(T0,data.Blank_SAP.(hemisphere).NREMtoREM.(['mean_' hemoDataType]) - data.Blank_SAP.(hemisphere).NREMtoREM.(['stdErr_' hemoDataType]),'-','color',colors('north texas green'),'LineWidth',0.25);
        if strcmp(hemoDataType,'GCaMP') == true
            ylabel('\DeltaF/F (%)')
        else
            ylabel(['\Delta[' hemoDataType '] (\muM)'])
        end
        xlim([-30,30])
        title('NREM to REM')
        xlabel('Time (s)')
        set(gca,'box','off')
        ax7.TickLength = [0.03,0.03];
        %% REM to Awake
        ax10 = subplot(2,2,4);
        plot(T0,data.SSP_SAP.(hemisphere).REMtoAWAKE.(['mean_' hemoDataType]),'-','color',colors('electric purple'),'LineWidth',2);
        hold on;
        plot(T0,data.SSP_SAP.(hemisphere).REMtoAWAKE.(['mean_' hemoDataType]) + data.SSP_SAP.(hemisphere).REMtoAWAKE.(['stdErr_' hemoDataType]),'-','color',colors('electric purple'),'LineWidth',0.25);
        plot(T0,data.SSP_SAP.(hemisphere).REMtoAWAKE.(['mean_' hemoDataType]) - data.SSP_SAP.(hemisphere).REMtoAWAKE.(['stdErr_' hemoDataType]),'-','color',colors('electric purple'),'LineWidth',0.25);
        plot(T0,data.Blank_SAP.(hemisphere).REMtoAWAKE.(['mean_' hemoDataType]),'-','color',colors('north texas green'),'LineWidth',2);
        plot(T0,data.Blank_SAP.(hemisphere).REMtoAWAKE.(['mean_' hemoDataType]) + data.Blank_SAP.(hemisphere).REMtoAWAKE.(['stdErr_' hemoDataType]),'-','color',colors('north texas green'),'LineWidth',0.25);
        plot(T0,data.Blank_SAP.(hemisphere).REMtoAWAKE.(['mean_' hemoDataType]) - data.Blank_SAP.(hemisphere).REMtoAWAKE.(['stdErr_' hemoDataType]),'-','color',colors('north texas green'),'LineWidth',0.25);
        if strcmp(hemoDataType,'GCaMP') == true
            ylabel('\DeltaF/F (%)')
        else
            ylabel(['\Delta[' hemoDataType '] (\muM)'])
        end
        xlim([-30,30])
        title('REM to Awake')
        xlabel('Time (s)')
        set(gca,'box','off')
        ax10.TickLength = [0.03,0.03];
        % save figure(s)
        if saveFigs == true
            savefig(summaryFigure,[dirpath 'ArousalTransitions_GCaMP_SI_ ' hemisphere '_' hemoDataType]);
        end
    end
end