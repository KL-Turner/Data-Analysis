function [AnalysisResults] = Fig1_S3_Test(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S3 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
animalIDs = {'T135','T141','T142','T144','T151','T155','T156','T157','T159'};
C57BL6J_IDs = {'T141','T155','T156','T157'};
SSP_SAP_IDs = {'T135','T142','T144','T151','T159'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
hemispheres = {'adjLH','adjRH'};
treatments = {'C57BL6J','SSP_SAP'};
% cd through each animal's directory and extract the appropriate analysis results
xx = 1;
zz = 1;
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    if ismember(animalID,C57BL6J_IDs) == true
        treatment = 'C57BL6J';
        for cc = 1:length(whiskDataTypes)
            whiskDataType = whiskDataTypes{1,cc};
            % left cortical
            data.(treatment).(whiskDataType).adjLH.HbT(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).CBV_HbT.HbT;
            data.(treatment).(whiskDataType).adjLH.CBV(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).CBV.CBV;
            data.(treatment).(whiskDataType).adjLH.cortMUA(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).MUA.corticalData;
            data.(treatment).(whiskDataType).adjLH.cortGam(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).Gam.corticalData;
            data.(treatment).(whiskDataType).adjLH.cortS(:,:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.corticalS;
            data.(treatment).(whiskDataType).adjLH.cortS_Gam(:,:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.corticalS(49:end,20:23);
            data.(treatment).(whiskDataType).adjLH.cortT(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.T;
            data.(treatment).(whiskDataType).adjLH.cortF(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.F;
            % right cortical
            data.(treatment).(whiskDataType).adjRH.HbT(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).CBV_HbT.HbT;
            data.(treatment).(whiskDataType).adjRH.CBV(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).CBV.CBV;
            data.(treatment).(whiskDataType).adjRH.cortMUA(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).MUA.corticalData;
            data.(treatment).(whiskDataType).adjRH.cortGam(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).Gam.corticalData;
            data.(treatment).(whiskDataType).adjRH.cortS(:,:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.corticalS;
            data.(treatment).(whiskDataType).adjRH.cortS_Gam(:,:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.corticalS(49:end,20:23);
            data.(treatment).(whiskDataType).adjRH.cortT(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.T;
            data.(treatment).(whiskDataType).adjRH.cortF(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.F;
            % hippocampal
            data.(treatment).(whiskDataType).Hip.hipMUA(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).MUA.hippocampalData;
            data.(treatment).(whiskDataType).Hip.hipGam(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).Gam.hippocampalData;
            data.(treatment).(whiskDataType).Hip.hipS(:,:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.hippocampalS;
            data.(treatment).(whiskDataType).Hip.hipS_Gam(:,:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.hippocampalS(49:end,20:23);
            data.(treatment).(whiskDataType).Hip.hipT(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.T;
            data.(treatment).(whiskDataType).Hip.hipF(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.F;
            % time vector
            data.(treatment).(whiskDataType).timeVector(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).timeVector;
        end
        xx = xx + 1;
    elseif ismember(animalID,SSP_SAP_IDs) == true
        treatment = 'SSP_SAP';
        for cc = 1:length(whiskDataTypes)
            whiskDataType = whiskDataTypes{1,cc};
            % left cortical
            data.(treatment).(whiskDataType).adjLH.HbT(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).CBV_HbT.HbT;
            data.(treatment).(whiskDataType).adjLH.CBV(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).CBV.CBV;
            data.(treatment).(whiskDataType).adjLH.cortMUA(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).MUA.corticalData;
            data.(treatment).(whiskDataType).adjLH.cortGam(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).Gam.corticalData;
            data.(treatment).(whiskDataType).adjLH.cortS(:,:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.corticalS;
            data.(treatment).(whiskDataType).adjLH.cortS_Gam(:,:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.corticalS(49:end,20:23);
            data.(treatment).(whiskDataType).adjLH.cortT(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.T;
            data.(treatment).(whiskDataType).adjLH.cortF(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.F;
            % right cortical
            data.(treatment).(whiskDataType).adjRH.HbT(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).CBV_HbT.HbT;
            data.(treatment).(whiskDataType).adjRH.CBV(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).CBV.CBV;
            data.(treatment).(whiskDataType).adjRH.cortMUA(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).MUA.corticalData;
            data.(treatment).(whiskDataType).adjRH.cortGam(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).Gam.corticalData;
            data.(treatment).(whiskDataType).adjRH.cortS(:,:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.corticalS;
            data.(treatment).(whiskDataType).adjRH.cortS_Gam(:,:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.corticalS(49:end,20:23);
            data.(treatment).(whiskDataType).adjRH.cortT(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.T;
            data.(treatment).(whiskDataType).adjRH.cortF(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.F;
            % hippocampal
            data.(treatment).(whiskDataType).Hip.hipMUA(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).MUA.hippocampalData;
            data.(treatment).(whiskDataType).Hip.hipGam(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).Gam.hippocampalData;
            data.(treatment).(whiskDataType).Hip.hipS(:,:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.hippocampalS;
            data.(treatment).(whiskDataType).Hip.hipS_Gam(:,:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.hippocampalS(49:end,20:23);
            data.(treatment).(whiskDataType).Hip.hipT(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.T;
            data.(treatment).(whiskDataType).Hip.hipF(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.F;
            % time vector
            data.(treatment).(whiskDataType).timeVector(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).timeVector;
        end
        zz = zz + 1;
    end
end
%% concatenate the data from the contra and ipsi data
% for ee = 1:length(whiskDataTypes)
%     whiskDataType = whiskDataTypes{1,ee};
%     data.(treatment).(whiskDataType).HbT = cat(2,data.(treatment).(whiskDataType).adjLH.HbT,data.(treatment).(whiskDataType).adjRH.HbT);
%     data.(treatment).(whiskDataType).CBV = cat(2,data.(treatment).(whiskDataType).adjLH.CBV,data.(treatment).(whiskDataType).adjRH.CBV);
%     data.(treatment).(whiskDataType).cortMUA = cat(2,data.(treatment).(whiskDataType).adjLH.cortMUA,data.(treatment).(whiskDataType).adjRH.cortMUA);
%     data.(treatment).(whiskDataType).cortGam = cat(2,data.(treatment).(whiskDataType).adjLH.cortGam,data.(treatment).(whiskDataType).adjRH.cortGam);
%     data.(treatment).(whiskDataType).cortS = cat(3,data.(treatment).(whiskDataType).adjLH.cortS,data.(treatment).(whiskDataType).adjRH.cortS);
%     data.(treatment).(whiskDataType).cortS_Gam = cat(3,data.(treatment).(whiskDataType).adjLH.cortS_Gam,data.(treatment).(whiskDataType).adjRH.cortS_Gam);
%     data.(treatment).(whiskDataType).cortT = cat(2,data.(treatment).(whiskDataType).adjLH.cortT,data.(treatment).(whiskDataType).adjRH.cortT);
%     data.(treatment).(whiskDataType).cortF = cat(2,data.(treatment).(whiskDataType).adjLH.cortF,data.(treatment).(whiskDataType).adjRH.cortF);
% end
%% concatenate the data from the contra and ipsi data
for ff = 1:length(treatments)
    treatment = treatments{1,ff};
    for ee = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,ee};
        for gg = 1:length(hemispheres)
            hemisphere = hemispheres{1,gg};
            data.(treatment).(whiskDataType).(hemisphere).meanHbT = mean(data.(treatment).(whiskDataType).(hemisphere).HbT,2);
            data.(treatment).(whiskDataType).(hemisphere).stdHbT = std(data.(treatment).(whiskDataType).(hemisphere).HbT,0,2);
            data.(treatment).(whiskDataType).(hemisphere).meanCBV = mean(data.(treatment).(whiskDataType).(hemisphere).CBV,2);
            data.(treatment).(whiskDataType).(hemisphere).stdCBV = std(data.(treatment).(whiskDataType).(hemisphere).CBV,0,2);
            data.(treatment).(whiskDataType).(hemisphere).meanCortMUA = mean(data.(treatment).(whiskDataType).(hemisphere).cortMUA,2);
            data.(treatment).(whiskDataType).(hemisphere).stdCortMUA = std(data.(treatment).(whiskDataType).(hemisphere).cortMUA,0,2);
            data.(treatment).(whiskDataType).(hemisphere).meanCortGam = mean(data.(treatment).(whiskDataType).(hemisphere).cortGam,2);
            data.(treatment).(whiskDataType).(hemisphere).stdCortGam = std(data.(treatment).(whiskDataType).(hemisphere).cortGam,0,2);
            data.(treatment).(whiskDataType).(hemisphere).meanCortS = mean(data.(treatment).(whiskDataType).(hemisphere).cortS,3).*100;
            data.(treatment).(whiskDataType).(hemisphere).mean_CortS_Gam = mean(mean(mean(data.(treatment).(whiskDataType).(hemisphere).cortS_Gam.*100,2),1),3);
            data.(treatment).(whiskDataType).(hemisphere).std_CortS_Gam = std(mean(mean(data.(treatment).(whiskDataType).(hemisphere).cortS_Gam.*100,2),1),0,3);
            data.(treatment).(whiskDataType).(hemisphere).meanCortT = mean(data.(treatment).(whiskDataType).(hemisphere).cortT,2);
            data.(treatment).(whiskDataType).(hemisphere).meanCortF = mean(data.(treatment).(whiskDataType).(hemisphere).cortF,2);
%             data.(treatment).(whiskDataType).(hemisphere).meanHipMUA = mean(data.(treatment).(whiskDataType).(hemisphere).Hip.hipMUA,2);
%             data.(treatment).(whiskDataType).(hemisphere).stdHipMUA = std(data.(treatment).(whiskDataType).(hemisphere).Hip.hipMUA,0,2);
%             data.(treatment).(whiskDataType).(hemisphere).meanHipGam = mean(data.(treatment).(whiskDataType).(hemisphere).Hip.hipGam,2);
%             data.(treatment).(whiskDataType).(hemisphere).stdHipGam = std(data.(treatment).(whiskDataType).(hemisphere).Hip.hipGam,0,2);
%             data.(treatment).(whiskDataType).(hemisphere).meanHipS = mean(data.(treatment).(whiskDataType).(hemisphere).Hip.hipS,3).*100;
%             data.(treatment).(whiskDataType).(hemisphere).mean_HipS_Gam = mean(mean(mean(data.(treatment).(whiskDataType).(hemisphere).Hip.hipS_Gam.*100,2),1),3);
%             data.(treatment).(whiskDataType).(hemisphere).std_HipS_Gam = std(mean(mean(data.(treatment).(whiskDataType).(hemisphere).Hip.hipS_Gam.*100,2),1),0,3);
%             data.(treatment).(whiskDataType).(hemisphere).meanHipT = mean(data.(treatment).(whiskDataType).(hemisphere).Hip.hipT,2);
%             data.(treatment).(whiskDataType).(hemisphere).meanHipF = mean(data.(treatment).(whiskDataType).(hemisphere).Hip.hipF,2);
            data.(treatment).(whiskDataType).meanTimeVector = mean(data.(treatment).(whiskDataType).timeVector,2);
        end
    end
end
%% Fig. 1-S3
summaryFigure = figure('Name','Fig1-S3 (a-r)');
sgtitle('Figure 1-S3 - Turner et al. 2020')


ax1 = subplot(2,2,1);
plot(data.C57BL6J.LongWhisks.meanTimeVector,data.C57BL6J.LongWhisks.adjLH.meanHbT,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.C57BL6J.LongWhisks.meanTimeVector,data.C57BL6J.LongWhisks.adjLH.meanHbT + data.C57BL6J.LongWhisks.adjLH.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.C57BL6J.LongWhisks.meanTimeVector,data.C57BL6J.LongWhisks.adjLH.meanHbT - data.C57BL6J.LongWhisks.adjLH.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
title('Control LH Extended whisk \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
xlim([-2,10])

ax2 = subplot(2,2,2);
plot(data.C57BL6J.LongWhisks.meanTimeVector,data.C57BL6J.LongWhisks.adjRH.meanHbT,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.C57BL6J.LongWhisks.meanTimeVector,data.C57BL6J.LongWhisks.adjRH.meanHbT + data.C57BL6J.LongWhisks.adjRH.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.C57BL6J.LongWhisks.meanTimeVector,data.C57BL6J.LongWhisks.adjRH.meanHbT - data.C57BL6J.LongWhisks.adjRH.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
title('Control RH Extended whisk \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
xlim([-2,10])

ax3 = subplot(2,2,3);
plot(data.SSP_SAP.LongWhisks.meanTimeVector,data.SSP_SAP.LongWhisks.adjLH.meanHbT,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.SSP_SAP.LongWhisks.meanTimeVector,data.SSP_SAP.LongWhisks.adjLH.meanHbT + data.SSP_SAP.LongWhisks.adjLH.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.SSP_SAP.LongWhisks.meanTimeVector,data.SSP_SAP.LongWhisks.adjLH.meanHbT - data.SSP_SAP.LongWhisks.adjLH.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
title('SSP-SAP UnRx LH Extended whisk \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
xlim([-2,10])

ax4 = subplot(2,2,4);
plot(data.SSP_SAP.LongWhisks.meanTimeVector,data.SSP_SAP.LongWhisks.adjRH.meanHbT,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.SSP_SAP.LongWhisks.meanTimeVector,data.SSP_SAP.LongWhisks.adjRH.meanHbT + data.SSP_SAP.LongWhisks.adjRH.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.SSP_SAP.LongWhisks.meanTimeVector,data.SSP_SAP.LongWhisks.adjRH.meanHbT - data.SSP_SAP.LongWhisks.adjRH.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
title('SSP-SAP treated RH Extended whisk \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
xlim([-2,10])

linkaxes([ax1,ax2,ax3,ax4],'xy')








% 
% 
% 
% %% [1-S3a] brief whisks cortical MUA
% ax1 = subplot(6,3,1);
% plot(data.(treatment).ShortWhisks.meanTimeVector,data.(treatment).ShortWhisks.meanCortMUA,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.(treatment).ShortWhisks.meanTimeVector,data.(treatment).ShortWhisks.meanCortMUA + data.(treatment).ShortWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.(treatment).ShortWhisks.meanTimeVector,data.(treatment).ShortWhisks.meanCortMUA - data.(treatment).ShortWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S3a] Brief whisk cortical MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax1.TickLength = [0.03,0.03];
% %% [1-S3b] moderate whisks cortical MUA
% ax2 = subplot(6,3,2);
% plot(data.(treatment).IntermediateWhisks.meanTimeVector,data.(treatment).IntermediateWhisks.meanCortMUA,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.(treatment).IntermediateWhisks.meanTimeVector,data.(treatment).IntermediateWhisks.meanCortMUA + data.(treatment).IntermediateWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.(treatment).IntermediateWhisks.meanTimeVector,data.(treatment).IntermediateWhisks.meanCortMUA - data.(treatment).IntermediateWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S3b] Moderate whisk cortical MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
% %% [1-S3c] extended whisks cortical MUA
% ax3 = subplot(6,3,3);
% plot(data.(treatment).LongWhisks.meanTimeVector,data.(treatment).LongWhisks.meanCortMUA,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.(treatment).LongWhisks.meanTimeVector,data.(treatment).LongWhisks.meanCortMUA + data.(treatment).LongWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.(treatment).LongWhisks.meanTimeVector,data.(treatment).LongWhisks.meanCortMUA - data.(treatment).LongWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S3c] Extended whisk cortical MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
% %% [1-S3d] brief whisks cortical LFP
% ax4 = subplot(6,3,4);
% imagesc(data.(treatment).ShortWhisks.meanCortT,data.(treatment).ShortWhisks.meanCortF,data.(treatment).ShortWhisks.meanCortS)
% title('[1-S3d] Brief whisk cortical LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-whisk time (s)')
% c4 = colorbar;
% ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-25,25])
% set(gca,'Ticklength',[0,0])
% axis square
% axis xy
% set(gca,'box','off')
% ax4.TickLength = [0.03,0.03];
% %% [1-S3e] moderate whisks cortical LFP
% ax5 = subplot(6,3,5);
% imagesc(data.(treatment).IntermediateWhisks.meanCortT,data.(treatment).IntermediateWhisks.meanCortF,data.(treatment).IntermediateWhisks.meanCortS)
% title('[1-S3e] Moderate whisk cortical LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-whisk time (s)')
% c5 = colorbar;
% ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-25,25])
% set(gca,'Ticklength',[0,0])
% axis square
% axis xy
% set(gca,'box','off')
% ax5.TickLength = [0.03,0.03];
% %% [1-S3f] extended whisks cortical LFP
% ax6 = subplot(6,3,6);
% imagesc(data.(treatment).LongWhisks.meanCortT,data.(treatment).LongWhisks.meanCortF,data.(treatment).LongWhisks.meanCortS)
% title('[1-S3f] Extended whisk cortical LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-whisk time (s)')
% c6 = colorbar;
% ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-25,25])
% set(gca,'Ticklength',[0,0])
% axis square
% axis xy
% set(gca,'box','off')
% ax6.TickLength = [0.03,0.03];
% %% [1-S3g] brief whisks hippocampal MUA
% ax7 = subplot(6,3,7);
% plot(data.(treatment).ShortWhisks.meanTimeVector,data.(treatment).ShortWhisks.meanHipMUA,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.(treatment).ShortWhisks.meanTimeVector,data.(treatment).ShortWhisks.meanHipMUA + data.(treatment).ShortWhisks.stdHipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.(treatment).ShortWhisks.meanTimeVector,data.(treatment).ShortWhisks.meanHipMUA - data.(treatment).ShortWhisks.stdHipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S3g] Brief whisk hippocampal MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax7.TickLength = [0.03,0.03];
% %% [1-S3h] moderate whisks hippocampal MUA
% ax8 = subplot(6,3,8);
% plot(data.(treatment).IntermediateWhisks.meanTimeVector,data.(treatment).IntermediateWhisks.meanHipMUA,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.(treatment).IntermediateWhisks.meanTimeVector,data.(treatment).IntermediateWhisks.meanHipMUA + data.(treatment).IntermediateWhisks.stdHipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.(treatment).IntermediateWhisks.meanTimeVector,data.(treatment).IntermediateWhisks.meanHipMUA - data.(treatment).IntermediateWhisks.stdHipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S3h] Moderate whisk hippocampal MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax8.TickLength = [0.03,0.03];
% %% [1-S3i] extended whisks hippocampal MUA
% ax9 = subplot(6,3,9);
% plot(data.(treatment).LongWhisks.meanTimeVector,data.(treatment).LongWhisks.meanHipMUA,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.(treatment).LongWhisks.meanTimeVector,data.(treatment).LongWhisks.meanHipMUA + data.(treatment).LongWhisks.stdHipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.(treatment).LongWhisks.meanTimeVector,data.(treatment).LongWhisks.meanHipMUA - data.(treatment).LongWhisks.stdHipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S3i] Extended whisk hippocampal MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax9.TickLength = [0.03,0.03];
% %% [1-S3j] brief whisks hippocampal LFP
% ax10 = subplot(6,3,10);
% imagesc(data.(treatment).ShortWhisks.meanHipT,data.(treatment).ShortWhisks.meanHipF,data.(treatment).ShortWhisks.meanHipS)
% title('[1-S3j] Brief whisk hippocampal LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-whisk time (s)')
% c10 = colorbar;
% ylabel(c10,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-25,25])
% set(gca,'Ticklength',[0,0])
% axis square
% axis xy
% set(gca,'box','off')
% ax10.TickLength = [0.03,0.03];
% %% [1-S3k] moderate whisks hippocampal LFP
% ax11 = subplot(6,3,11);
% imagesc(data.(treatment).IntermediateWhisks.meanHipT,data.(treatment).IntermediateWhisks.meanHipF,data.(treatment).IntermediateWhisks.meanHipS)
% title('[1-S3k] Moderate whisk hippocampal LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-whisk time (s)')
% c11 = colorbar;
% ylabel(c11,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-25,25])
% set(gca,'Ticklength',[0 0])
% axis square
% axis xy
% set(gca,'box','off')
% ax11.TickLength = [0.03,0.03];
% %% [1-S3l] extended whisks hippocampal LFP
% ax12 = subplot(6,3,12);
% imagesc(data.(treatment).LongWhisks.meanHipT,data.(treatment).LongWhisks.meanHipF,data.(treatment).LongWhisks.meanHipS)
% title('[1-S3l] Extended whisk hippocampal LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-whisk time (s)')
% c12 = colorbar;
% ylabel(c12,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-25,25])
% set(gca,'Ticklength',[0,0])
% axis square
% axis xy
% set(gca,'box','off')
% ax12.TickLength = [0.03,0.03];
% %% [1-S3m] brief whisks HbT
% ax13 = subplot(6,3,13);
% plot(data.(treatment).ShortWhisks.meanTimeVector,data.(treatment).ShortWhisks.meanHbT,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.(treatment).ShortWhisks.meanTimeVector,data.(treatment).ShortWhisks.meanHbT + data.(treatment).ShortWhisks.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.(treatment).ShortWhisks.meanTimeVector,data.(treatment).ShortWhisks.meanHbT - data.(treatment).ShortWhisks.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S3m] Brief whisk \Delta[HbT] (\muM)')
% ylabel('\Delta[HbT] (\muM)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax13.TickLength = [0.03,0.03];
% %% [1-S3n] moderate whisks HbT
% ax14 = subplot(6,3,14);
% plot(data.(treatment).IntermediateWhisks.meanTimeVector,data.(treatment).IntermediateWhisks.meanHbT,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.(treatment).IntermediateWhisks.meanTimeVector,data.(treatment).IntermediateWhisks.meanHbT + data.(treatment).IntermediateWhisks.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.(treatment).IntermediateWhisks.meanTimeVector,data.(treatment).IntermediateWhisks.meanHbT - data.(treatment).IntermediateWhisks.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S3n] Moderate whisk \Delta[HbT] (\muM)')
% ylabel('\Delta[HbT] (\muM)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax14.TickLength = [0.03,0.03];
% %% [1-S3o] extended whisks HbT
% ax15 = subplot(6,3,15);
% plot(data.(treatment).LongWhisks.meanTimeVector,data.(treatment).LongWhisks.meanHbT,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.(treatment).LongWhisks.meanTimeVector,data.(treatment).LongWhisks.meanHbT + data.(treatment).LongWhisks.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.(treatment).LongWhisks.meanTimeVector,data.(treatment).LongWhisks.meanHbT - data.(treatment).LongWhisks.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S3o] Extended whisk \Delta[HbT] (\muM)')
% ylabel('\Delta[HbT] (\muM)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax15.TickLength = [0.03,0.03];
% %% [1-S3p] brief whisks refl
% ax16 = subplot(6,3,16);
% plot(data.(treatment).ShortWhisks.meanTimeVector,data.(treatment).ShortWhisks.meanCBV,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.(treatment).ShortWhisks.meanTimeVector,data.(treatment).ShortWhisks.meanCBV + data.(treatment).ShortWhisks.stdCBV,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.(treatment).ShortWhisks.meanTimeVector,data.(treatment).ShortWhisks.meanCBV - data.(treatment).ShortWhisks.stdCBV,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S3p] Brief whisk reflectance')
% ylabel('\DeltaR/R (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax16.TickLength = [0.03,0.03];
% %% [1-S3q] moderate whisks refl
% ax17 = subplot(6,3,17);
% plot(data.(treatment).IntermediateWhisks.meanTimeVector,data.(treatment).IntermediateWhisks.meanCBV,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.(treatment).IntermediateWhisks.meanTimeVector,data.(treatment).IntermediateWhisks.meanCBV + data.(treatment).IntermediateWhisks.stdCBV,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.(treatment).IntermediateWhisks.meanTimeVector,data.(treatment).IntermediateWhisks.meanCBV - data.(treatment).IntermediateWhisks.stdCBV,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S3q] Moderate whisk reflectance')
% ylabel('\DeltaR/R (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax17.TickLength = [0.03,0.03];
% %% [1-S3r] extended whisks refl
% ax18 = subplot(6,3,18);
% plot(data.(treatment).LongWhisks.meanTimeVector,data.(treatment).LongWhisks.meanCBV,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.(treatment).LongWhisks.meanTimeVector,data.(treatment).LongWhisks.meanCBV + data.(treatment).LongWhisks.stdCBV,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.(treatment).LongWhisks.meanTimeVector,data.(treatment).LongWhisks.meanCBV - data.(treatment).LongWhisks.stdCBV,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S3r] Extended whisk reflectance')
% ylabel('\DeltaR/R (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax18.TickLength = [0.03,0.03];
% %% axes positions
% linkaxes([ax1,ax2,ax3,ax7,ax8,ax9],'xy')
% linkaxes([ax4,ax5,ax6,ax10,ax11,ax12],'xy')
% linkaxes([ax13,ax14,ax15],'xy')
% linkaxes([ax16,ax17,ax18],'xy')
% ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
% ax3Pos = get(ax3,'position');
% ax4Pos = get(ax4,'position');
% ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
% ax10Pos = get(ax10,'position');
% ax11Pos = get(ax11,'position');
% ax12Pos = get(ax12,'position');
% ax4Pos(3:4) = ax1Pos(3:4);
% ax5Pos(3:4) = ax2Pos(3:4);
% ax6Pos(3:4) = ax3Pos(3:4);
% ax10Pos(3:4) = ax1Pos(3:4);
% ax11Pos(3:4) = ax2Pos(3:4);
% ax12Pos(3:4) = ax3Pos(3:4);
% set(ax4,'position',ax4Pos);
% set(ax5,'position',ax5Pos);
% set(ax6,'position',ax6Pos);
% set(ax10,'position',ax10Pos);
% set(ax11,'position',ax11Pos);
% set(ax12,'position',ax12Pos);
% %% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig1-S3']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig1-S3'])
    %% Text diary
    diaryFile = [dirpath 'Fig1-S3_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
end
%     diary(diaryFile)
%     diary on
%     % text values
%     disp('======================================================================================================================')
%     disp('[1-S3] Text values for gamma/HbT changes')
%     disp('======================================================================================================================')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     % cortical MUA/LFP
%     [~,index] = max(data.(treatment).ShortWhisks.meanCortMUA);
%     disp(['Brief whisk Cort gamma MUA P/P (%): ' num2str(round(data.(treatment).ShortWhisks.meanCortMUA(index),1)) ' +/- ' num2str(round(data.(treatment).ShortWhisks.stdCortMUA(index),1))]); disp(' ')
%     [~,index] = max(data.(treatment).IntermediateWhisks.meanCortMUA);
%     disp(['Moderate whisk Cort gamma MUA P/P (%): ' num2str(round(data.(treatment).IntermediateWhisks.meanCortMUA(index),1)) ' +/- ' num2str(round(data.(treatment).IntermediateWhisks.stdCortMUA(index),1))]); disp(' ')
%     [~,index] = max(data.(treatment).LongWhisks.meanCortMUA);
%     disp(['Extended whisk Cort gamma MUA P/P (%): ' num2str(round(data.(treatment).LongWhisks.meanCortMUA(index),1)) ' +/- ' num2str(round(data.(treatment).LongWhisks.stdCortMUA(index),1))]); disp(' ')
%     % cortical LFP
%     disp(['Brief whisk Cort gamma LFP P/P (%): ' num2str(round(data.(treatment).ShortWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.(treatment).ShortWhisks.std_CortS_Gam,1))]); disp(' ')
%     disp(['Moderate whisk Cort gamma LFP P/P (%): ' num2str(round(data.(treatment).IntermediateWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.(treatment).IntermediateWhisks.std_CortS_Gam,1))]); disp(' ')
%     disp(['Extended whisk Cort gamma LFP P/P (%): ' num2str(round(data.(treatment).LongWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.(treatment).LongWhisks.std_CortS_Gam,1))]); disp(' ')
%     % hippocampal MUA
%     [~,index] = max(data.(treatment).ShortWhisks.meanHipMUA);
%     disp(['Brief whisk Hip gamma MUA P/P (%): ' num2str(round(data.(treatment).ShortWhisks.meanHipMUA(index),1)) ' +/- ' num2str(round(data.(treatment).ShortWhisks.stdHipMUA(index),1))]); disp(' ')
%     [~,index] = max(data.(treatment).IntermediateWhisks.meanHipMUA);
%     disp(['Moderate whisk Hip gamma MUA P/P (%): ' num2str(round(data.(treatment).IntermediateWhisks.meanHipMUA(index),1)) ' +/- ' num2str(round(data.(treatment).IntermediateWhisks.stdHipMUA(index),1))]); disp(' ')
%     [~,index] = max(data.(treatment).LongWhisks.meanHipMUA);
%     disp(['Extended whisk Hip gamma MUA P/P (%): ' num2str(round(data.(treatment).LongWhisks.meanHipMUA(index),1)) ' +/- ' num2str(round(data.(treatment).LongWhisks.stdHipMUA(index),1))]); disp(' ')
%     % hippocampal LFP
%     disp(['Brief whisk Hip gamma LFP P/P (%): ' num2str(round(data.(treatment).ShortWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.(treatment).ShortWhisks.std_HipS_Gam,1))]); disp(' ')
%     disp(['Moderate whisk Hip gamma LFP P/P (%): ' num2str(round(data.(treatment).IntermediateWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.(treatment).IntermediateWhisks.std_HipS_Gam,1))]); disp(' ')
%     disp(['Extended whisk Hip gamma LFP P/P (%): ' num2str(round(data.(treatment).LongWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.(treatment).LongWhisks.std_HipS_Gam,1))]); disp(' ')
%     % HbT
%     [~,index] = max(data.(treatment).ShortWhisks.meanHbT);
%     disp(['Brief whisk [HbT] (uM): ' num2str(round(data.(treatment).ShortWhisks.meanHbT(index),1)) ' +/- ' num2str(round(data.(treatment).ShortWhisks.stdHbT(index),1))]); disp(' ')
%     [~,index] = max(data.(treatment).IntermediateWhisks.meanHbT);
%     disp(['Moderate whisk [HbT] (uM): ' num2str(round(data.(treatment).IntermediateWhisks.meanHbT(index),1)) ' +/- ' num2str(round(data.(treatment).IntermediateWhisks.stdHbT(index),1))]); disp(' ')
%     [~,index] = max(data.(treatment).LongWhisks.meanHbT);
%     disp(['Extended whisk [HbT] (uM): ' num2str(round(data.(treatment).LongWhisks.meanHbT(index),1)) ' +/- ' num2str(round(data.(treatment).LongWhisks.stdHbT(index),1))]); disp(' ')
%     % R/R
%     [~,index] = min(data.(treatment).ShortWhisks.meanCBV);
%     disp(['Brief whisk refl R/R (%): ' num2str(round(data.(treatment).ShortWhisks.meanCBV(index),1)) ' +/- ' num2str(round(data.(treatment).ShortWhisks.stdCBV(index),1))]); disp(' ')
%     [~,index] = min(data.(treatment).IntermediateWhisks.meanCBV);
%     disp(['Moderate whisk refl R/R (%): ' num2str(round(data.(treatment).IntermediateWhisks.meanCBV(index),1)) ' +/- ' num2str(round(data.(treatment).IntermediateWhisks.stdCBV(index),1))]); disp(' ')
%     [~,index] = min(data.(treatment).LongWhisks.meanCBV);
%     disp(['Extended whisk refl R/R (%): ' num2str(round(data.(treatment).LongWhisks.meanCBV(index),1)) ' +/- ' num2str(round(data.(treatment).LongWhisks.stdCBV(index),1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     diary off
% end
% 
% end
