function [AnalysisResults] = Fig1_S2_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S2 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
IOSanimalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
dataTypes = {'adjLH','adjRH'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for dd = 1:length(solenoidNames)
            solenoidName = solenoidNames{1,dd};
            data.(dataType).(solenoidName).count(:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).count;
            data.(dataType).(solenoidName).HbT(:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).CBV_HbT.HbT;
            data.(dataType).(solenoidName).CBV(:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).CBV.CBV;
            data.(dataType).(solenoidName).cortMUA(:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).MUA.corticalData;
            data.(dataType).(solenoidName).hipMUA(:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).MUA.hippocampalData;
            data.(dataType).(solenoidName).cortGam(:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).Gam.corticalData;
            data.(dataType).(solenoidName).hipGam(:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).Gam.hippocampalData;
            data.(dataType).(solenoidName).timeVector(:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).timeVector;
            data.(dataType).(solenoidName).cortS(:,:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.corticalS;
            data.(dataType).(solenoidName).cortS_Gam(:,:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.corticalS(49:end,20:23);
            data.(dataType).(solenoidName).hipS(:,:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.hippocampalS;
            data.(dataType).(solenoidName).hipS_Gam(:,:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.hippocampalS(49:end,20:23);
            data.(dataType).(solenoidName).T(:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.T;
            data.(dataType).(solenoidName).F(:,aa) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.F;
        end
    end
end
% concatenate the data from the contra and ipsi data
data.Contra.count = cat(2,data.adjLH.RPadSol.count,data.adjRH.LPadSol.count);
data.Contra.HbT = cat(2,data.adjLH.RPadSol.HbT,data.adjRH.LPadSol.HbT);
data.Contra.CBV = cat(2,data.adjLH.RPadSol.CBV,data.adjRH.LPadSol.CBV);
data.Contra.cortMUA = cat(2,data.adjLH.RPadSol.cortMUA,data.adjRH.LPadSol.cortMUA);
data.Contra.hipMUA = data.adjRH.RPadSol.hipMUA;
data.Contra.cortGam = cat(2,data.adjLH.RPadSol.cortGam,data.adjRH.LPadSol.cortGam);
data.Contra.hipGam = data.adjRH.RPadSol.hipGam;
data.Contra.timeVector = cat(2,data.adjLH.RPadSol.timeVector,data.adjRH.LPadSol.timeVector);
data.Contra.cortS = cat(3,data.adjLH.RPadSol.cortS,data.adjRH.LPadSol.cortS);
data.Contra.cortS_Gam = cat(3,data.adjLH.RPadSol.cortS_Gam,data.adjRH.LPadSol.cortS_Gam);
data.Contra.hipS = data.adjRH.RPadSol.hipS;
data.Contra.hipS_Gam = data.adjRH.RPadSol.hipS_Gam;
data.Contra.T = cat(2,data.adjLH.RPadSol.T,data.adjRH.LPadSol.T);
data.Contra.F = cat(2,data.adjLH.RPadSol.F,data.adjRH.LPadSol.F);
data.Ipsi.count = cat(2,data.adjLH.LPadSol.count,data.adjRH.RPadSol.count);
data.Ipsi.HbT = cat(2,data.adjLH.LPadSol.HbT,data.adjRH.RPadSol.HbT);
data.Ipsi.CBV = cat(2,data.adjLH.LPadSol.CBV,data.adjRH.RPadSol.CBV);
data.Ipsi.cortMUA = cat(2,data.adjLH.LPadSol.cortMUA,data.adjRH.RPadSol.cortMUA);
data.Ipsi.hipMUA = data.adjRH.LPadSol.hipMUA;
data.Ipsi.cortGam = cat(2,data.adjLH.LPadSol.cortGam,data.adjRH.RPadSol.cortGam);
data.Ipsi.hipGam = data.adjRH.LPadSol.hipGam;
data.Ipsi.timeVector = cat(2,data.adjLH.LPadSol.timeVector,data.adjRH.RPadSol.timeVector);
data.Ipsi.cortS = cat(3,data.adjLH.LPadSol.cortS,data.adjRH.RPadSol.cortS);
data.Ipsi.cortS_Gam = cat(3,data.adjLH.LPadSol.cortS_Gam,data.adjRH.RPadSol.cortS_Gam);
data.Ipsi.hipS = data.adjRH.LPadSol.hipS;
data.Ipsi.hipS_Gam = data.adjRH.LPadSol.hipS_Gam;
data.Ipsi.T = cat(2,data.adjLH.LPadSol.T,data.adjRH.RPadSol.T);
data.Ipsi.F = cat(2,data.adjLH.LPadSol.F,data.adjRH.RPadSol.F);
data.Auditory.count = cat(2,data.adjLH.AudSol.count,data.adjRH.AudSol.count);
data.Auditory.HbT = cat(2,data.adjLH.AudSol.HbT,data.adjRH.AudSol.HbT);
data.Auditory.CBV = cat(2,data.adjLH.AudSol.CBV,data.adjRH.AudSol.CBV);
data.Auditory.cortMUA = cat(2,data.adjLH.AudSol.cortMUA,data.adjRH.AudSol.cortMUA);
data.Auditory.hipMUA = data.adjRH.AudSol.hipMUA;
data.Auditory.cortGam = cat(2,data.adjLH.AudSol.cortGam,data.adjRH.AudSol.cortGam);
data.Auditory.hipGam = data.adjRH.AudSol.hipGam;
data.Auditory.timeVector = cat(2,data.adjLH.AudSol.timeVector,data.adjRH.AudSol.timeVector);
data.Auditory.cortS = cat(3,data.adjLH.AudSol.cortS,data.adjRH.AudSol.cortS);
data.Auditory.cortS_Gam = cat(3,data.adjLH.AudSol.cortS_Gam,data.adjRH.AudSol.cortS_Gam);
data.Auditory.hipS = data.adjRH.AudSol.hipS;
data.Auditory.hipS_Gam = data.adjRH.AudSol.hipS_Gam;
data.Auditory.T = cat(2,data.adjLH.AudSol.T,data.adjRH.AudSol.T);
data.Auditory.F = cat(2,data.adjLH.AudSol.F,data.adjRH.AudSol.F);
% take the averages of each field through the proper dimension
for ff = 1:length(compDataTypes)
    compDataType = compDataTypes{1,ff};
    data.(compDataType).mean_Count = mean(data.(compDataType).count,2);
    data.(compDataType).std_Count = std(data.(compDataType).count,0,2);   
    data.(compDataType).mean_HbT = mean(data.(compDataType).HbT,2);
    data.(compDataType).std_HbT = std(data.(compDataType).HbT,0,2);
    data.(compDataType).mean_CBV = mean(data.(compDataType).CBV,2);
    data.(compDataType).std_CBV = std(data.(compDataType).CBV,0,2);
    data.(compDataType).mean_CortMUA = mean(data.(compDataType).cortMUA,2);
    data.(compDataType).std_CortMUA = std(data.(compDataType).cortMUA,0,2);
    data.(compDataType).mean_HipMUA = mean(data.(compDataType).hipMUA,2);
    data.(compDataType).std_HipMUA = std(data.(compDataType).hipMUA,0,2);
    data.(compDataType).mean_CortGam = mean(data.(compDataType).cortGam,2);
    data.(compDataType).std_CortGam = std(data.(compDataType).cortGam,0,2);
    data.(compDataType).mean_HipGam = mean(data.(compDataType).hipGam,2);
    data.(compDataType).std_HipGam = std(data.(compDataType).hipGam,0,2);
    data.(compDataType).mean_timeVector = mean(data.(compDataType).timeVector,2);
    data.(compDataType).mean_CortS = mean(data.(compDataType).cortS,3).*100;
    data.(compDataType).mean_CortS_Gam = mean(mean(mean(data.(compDataType).cortS_Gam.*100,2),1),3);
    data.(compDataType).std_CortS_Gam = std(mean(mean(data.(compDataType).cortS_Gam.*100,2),1),0,3);
    data.(compDataType).mean_HipS = mean(data.(compDataType).hipS,3).*100;
    data.(compDataType).mean_HipS_Gam = mean(mean(mean(data.(compDataType).hipS_Gam.*100,2),1),3);
    data.(compDataType).std_HipS_Gam = std(mean(mean(data.(compDataType).hipS_Gam.*100,2),1),0,3);
    data.(compDataType).mean_T = mean(data.(compDataType).T,2);
    data.(compDataType).mean_F = mean(data.(compDataType).F,2);
end
%% Fig. 1-S2
summaryFigure = figure('Name','Fig1-S2 (a-r)');
sgtitle('Figure 1-S2 (& 1d) - Turner et al. 2020')
%% [1-S2a] cortical MUA contra stim
ax1 = subplot(6,3,1);
plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA,'color',colors_eLife2020('rich black'),'LineWidth',1);
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA + data.Contra.std_CortMUA,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA - data.Contra.std_CortMUA,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[1d,1-S2a] Contra stim cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [1-S2b] cortical MUA ispi stim
ax2 = subplot(6,3,2);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA + data.Ipsi.std_CortMUA,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA - data.Ipsi.std_CortMUA,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[1-S2b] Ipsi stim cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [1-S2c] cortical MUA auditory stim
ax3 = subplot(6,3,3);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA + data.Auditory.std_CortMUA,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA - data.Auditory.std_CortMUA,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[1-S2c] Aud stim cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [1-S2d] cortical LFP contra stim
ax4 = subplot(6,3,4);
imagesc(data.Contra.mean_T,data.Contra.mean_F,data.Contra.mean_CortS)
title('[1d,1-S2d] Contra stim cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [1-S2e] cortical LFP ispi stim
ax5 = subplot(6,3,5);
imagesc(data.Ipsi.mean_T,data.Ipsi.mean_F,data.Ipsi.mean_CortS)
title('[1-S2e] Ipsi stim cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [1-S2f] cortical LFP auditory stim
ax6 = subplot(6,3,6);
imagesc(data.Auditory.mean_T,data.Auditory.mean_F,data.Auditory.mean_CortS)
title('[1-S2f] Aud stim cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [1-S2g] hippocampal MUA contra stim
ax7 = subplot(6,3,7);
plot(data.Contra.mean_timeVector,data.Contra.mean_HipMUA,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_HipMUA + data.Contra.std_HipMUA,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Contra.mean_timeVector,data.Contra.mean_HipMUA - data.Contra.std_HipMUA,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[1-S2g] Contra stim hippocampal MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% [1-S2h] hippocampal MUA ispi stim
ax8 = subplot(6,3,8);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipMUA,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipMUA + data.Ipsi.std_HipMUA,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipMUA - data.Ipsi.std_HipMUA,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[1-S2h] Ipsi stim hippocampal MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% [1-S2i] hippocampal MUA auditory stim
ax9 = subplot(6,3,9);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipMUA,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipMUA + data.Auditory.std_HipMUA,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipMUA - data.Auditory.std_HipMUA,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[1-S2i] Aud stim hippocampal MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [1-S2j] hippocampal LFP contra stim
ax10 = subplot(6,3,10);
imagesc(data.Contra.mean_T,data.Contra.mean_F,data.Contra.mean_HipS)
title('[1-S2j] Contra stim hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c10 = colorbar;
ylabel(c10,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax10.TickLength = [0.03,0.03];
%% [1-S2j] hippocampal LFP ispi stim
ax11 = subplot(6,3,11);
imagesc(data.Ipsi.mean_T,data.Ipsi.mean_F,data.Ipsi.mean_HipS)
title('[1-S2j] Ipsi stim hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c11 = colorbar;
ylabel(c11,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
%% [1-S2l] hippocampal LFP auditory stim
ax12 = subplot(6,3,12);
imagesc(data.Auditory.mean_T,data.Auditory.mean_F,data.Auditory.mean_HipS)
title('[1-S2l] Aud stim hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c12 = colorbar;
ylabel(c12,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];
%% [1-S2m] HbT contra stim
ax13 = subplot(6,3,13);
plot(data.Contra.mean_timeVector,data.Contra.mean_HbT,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_HbT + data.Contra.std_HbT,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Contra.mean_timeVector,data.Contra.mean_HbT - data.Contra.std_HbT,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[1d,1-S2m] Contra stim \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax13.TickLength = [0.03,0.03];
%% [1-S2n] HbT ispi stim
ax14 = subplot(6,3,14);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HbT,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HbT + data.Ipsi.std_HbT,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HbT - data.Ipsi.std_HbT,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[1-S2n] Ipsi stim \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];
%% [1-S2o] HbT auditory stim
ax15 = subplot(6,3,15);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HbT,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HbT + data.Auditory.std_HbT,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HbT - data.Auditory.std_HbT,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[1-S2o] Aud stim \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];
%% [1-S2p] refl contra stim
ax16 = subplot(6,3,16);
plot(data.Contra.mean_timeVector,data.Contra.mean_CBV,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_CBV + data.Contra.std_CBV,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Contra.mean_timeVector,data.Contra.mean_CBV - data.Contra.std_CBV,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[1-S2p] Contra stim reflectance')
ylabel('\DeltaR/R (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
%% [1-S2q] refl ispi stim
ax17 = subplot(6,3,17);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CBV,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CBV + data.Ipsi.std_CBV,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CBV - data.Ipsi.std_CBV,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[1-S2q] Ipsi stim reflectance')
ylabel('\DeltaR/R (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax17.TickLength = [0.03,0.03];
%% [1-S2r] refl auditory stim
ax18 = subplot(6,3,18);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CBV,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CBV + data.Auditory.std_CBV,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CBV - data.Auditory.std_CBV,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[1-S2r] Aud stim reflectance')
ylabel('\DeltaR/R (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax18.TickLength = [0.03,0.03];
%% adjust and link axes
linkaxes([ax1,ax2,ax3,ax7,ax8,ax9],'xy')
linkaxes([ax4,ax5,ax6,ax10,ax11,ax12],'xy')
linkaxes([ax13,ax14,ax15],'xy')
linkaxes([ax16,ax17,ax18],'xy')
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
ax12Pos = get(ax12,'position');
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax2Pos(3:4);
ax6Pos(3:4) = ax3Pos(3:4);
ax10Pos(3:4) = ax1Pos(3:4);
ax11Pos(3:4) = ax2Pos(3:4);
ax12Pos(3:4) = ax3Pos(3:4);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax10,'position',ax10Pos);
set(ax11,'position',ax11Pos);
set(ax12,'position',ax12Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig1-S2']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig1-S2'])
    %% text diary
    diaryFile = [dirpath 'Fig1-S2_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % text values
    disp('======================================================================================================================')
    disp('[1-S2] Text values for gamma/HbT/reflectance changes')
    disp('======================================================================================================================')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % cortical MUA/LFP
    [~,index] = max(data.Contra.mean_CortMUA);
    disp(['Contra stim Cort gamma MUA P/P (%): ' num2str(round(data.Contra.mean_CortMUA(index),1)) ' +/- ' num2str(round(data.Contra.std_CortMUA(index),1))]); disp(' ')
    [~,index] = max(data.Ipsi.mean_CortMUA);
    disp(['Ipsil stim Cort gamma MUA P/P (%): ' num2str(round(data.Ipsi.mean_CortMUA(index),1)) ' +/- ' num2str(round(data.Ipsi.std_CortMUA(index),1))]); disp(' ')
    [~,index] = max(data.Auditory.mean_CortMUA);
    disp(['Audit stim Cort gamma MUA P/P (%): ' num2str(round(data.Auditory.mean_CortMUA(index),1)) ' +/- ' num2str(round(data.Auditory.std_CortMUA(index),1))]); disp(' ')
    % cortical LFP
    disp(['Contra stim Cort gamma LFP P/P (%): ' num2str(round(data.Contra.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Contra.std_CortS_Gam,1))]); disp(' ')
    disp(['Ipsil stim Cort gamma LFP P/P (%): ' num2str(round(data.Ipsi.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Ipsi.std_CortS_Gam,1))]); disp(' ')
    disp(['Audit stim Cort gamma LFP P/P (%): ' num2str(round(data.Auditory.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Auditory.std_CortS_Gam,1))]); disp(' ')
    % hippocampal MUA
    [~,index] = max(data.Contra.mean_HipMUA);
    disp(['Contra stim Hip gamma MUA P/P (%): ' num2str(round(data.Contra.mean_HipMUA(index),1)) ' +/- ' num2str(round(data.Contra.std_HipMUA(index),1))]); disp(' ')
    [~,index] = max(data.Ipsi.mean_HipMUA);
    disp(['Ipsil stim Hip gamma MUA P/P (%): ' num2str(round(data.Ipsi.mean_HipMUA(index),1)) ' +/- ' num2str(round(data.Ipsi.std_HipMUA(index),1))]); disp(' ')
    [~,index] = max(data.Auditory.mean_HipMUA);
    disp(['Audit stim Hip gamma MUA P/P (%): ' num2str(round(data.Auditory.mean_HipMUA(index),1)) ' +/- ' num2str(round(data.Auditory.std_HipMUA(index),1))]); disp(' ')
    % hippocampal LFP
    disp(['Contra stim Hip gamma LFP P/P (%): ' num2str(round(data.Contra.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Contra.std_HipS_Gam,1))]); disp(' ')
    disp(['Ipsil stim Hip gamma LFP P/P (%): ' num2str(round(data.Ipsi.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Ipsi.std_HipS_Gam,1))]); disp(' ')
    disp(['Auditory stim Hip gamma LFP P/P (%): ' num2str(round(data.Auditory.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Auditory.std_HipS_Gam,1))]); disp(' ')
    % HbT
    [~,index] = max(data.Contra.mean_HbT);
    disp(['Contra stim [HbT] (uM): ' num2str(round(data.Contra.mean_HbT(index),1)) ' +/- ' num2str(round(data.Contra.std_HbT(index),1))]); disp(' ')
    [~,index] = max(data.Ipsi.mean_HbT);
    disp(['Ipsil stim [HbT] (uM): ' num2str(round(data.Ipsi.mean_HbT(index),1)) ' +/- ' num2str(round(data.Ipsi.std_HbT(index),1))]); disp(' ')
    [~,index] = max(data.Auditory.mean_HbT);
    disp(['Audit stim [HbT] (uM): ' num2str(round(data.Auditory.mean_HbT(index),1)) ' +/- ' num2str(round(data.Auditory.std_HbT(index),1))]); disp(' ')
    % R/R
    [~,index] = min(data.Contra.mean_CBV);
    disp(['Contra stim refl R/R (%): ' num2str(round(data.Contra.mean_CBV(index),1)) ' +/- ' num2str(round(data.Contra.std_CBV(index),1))]); disp(' ')
    [~,index] = min(data.Ipsi.mean_CBV);
    disp(['Ipsil stim refl R/R (%): ' num2str(round(data.Ipsi.mean_CBV(index),1)) ' +/- ' num2str(round(data.Ipsi.std_CBV(index),1))]); disp(' ')
    [~,index] = min(data.Auditory.mean_CBV);
    disp(['Audit stim refl R/R (%): ' num2str(round(data.Auditory.mean_CBV(index),1)) ' +/- ' num2str(round(data.Auditory.std_CBV(index),1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
end

end
