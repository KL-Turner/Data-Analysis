function [AnalysisResults] = Fig6_Test(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generate figure panel 6 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
animalIDs = {'T135','T141','T142','T144','T151','T155','T156','T157','T159'};
behavFields = {'Rest','NREM','REM'};
treatments = {'C57BL6J','SSP_SAP'};
%% cd through each animal's directory and extract the appropriate analysis results
xx = 1;
zz = 1;
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    if strcmp(animalID,'T141') == true || strcmp(animalID,'T155') == true || strcmp(animalID,'T156') == true || strcmp(animalID,'T157') == true
        for bb = 1:length(behavFields)
            behavField = behavFields{1,bb};
            data.C57BL6J.(behavField).adjLH.HbTvLFPxcVals(:,:,xx) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvLFPxcVals;
            data.C57BL6J.(behavField).adjLH.LFP_lags(:,:,xx) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.LFP_lags;
            data.C57BL6J.(behavField).adjLH.F(:,:,xx) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.F;
            data.C57BL6J.(behavField).adjRH.HbTvLFPxcVals(:,:,xx) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvLFPxcVals;
            data.C57BL6J.(behavField).adjRH.LFP_lags(:,:,xx) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.LFP_lags;
            data.C57BL6J.(behavField).adjRH.F(:,:,xx) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.F;
            data.C57BL6J.(behavField).adjLH.HbTvMUAxcVals(:,xx) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvMUAxcVals;
            data.C57BL6J.(behavField).adjLH.HbTvMUAxcVals_std(:,xx) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvMUAxcVals_std;
            data.C57BL6J.(behavField).adjLH.MUA_lags(:,xx) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.LFP_lags;
            data.C57BL6J.(behavField).adjRH.HbTvMUAxcVals(:,xx) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvMUAxcVals;
            data.C57BL6J.(behavField).adjRH.HbTvMUAxcVals_std(:,xx) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvMUAxcVals_std;
            data.C57BL6J.(behavField).adjRH.MUA_lags(:,xx) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.LFP_lags;
            data.C57BL6J.(behavField).animalID{xx,1} = animalID;
            data.C57BL6J.(behavField).behavior{xx,1} = behavField;
            data.C57BL6J.(behavField).LH{xx,1} = 'LH';
            data.C57BL6J.(behavField).RH{xx,1} = 'RH';
        end
        xx = xx + 1;
    elseif strcmp(animalID,'T135') == true || strcmp(animalID,'T142') == true || strcmp(animalID,'T144') == true || strcmp(animalID,'T151') == true || strcmp(animalID,'T159') == true
        for bb = 1:length(behavFields)
            behavField = behavFields{1,bb};
            data.SSP_SAP.(behavField).adjLH.HbTvLFPxcVals(:,:,zz) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvLFPxcVals;
            data.SSP_SAP.(behavField).adjLH.LFP_lags(:,:,zz) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.LFP_lags;
            data.SSP_SAP.(behavField).adjLH.F(:,:,zz) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.F;
            data.SSP_SAP.(behavField).adjRH.HbTvLFPxcVals(:,:,zz) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvLFPxcVals;
            data.SSP_SAP.(behavField).adjRH.LFP_lags(:,:,zz) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.LFP_lags;
            data.SSP_SAP.(behavField).adjRH.F(:,:,zz) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.F;
            data.SSP_SAP.(behavField).adjLH.HbTvMUAxcVals(:,zz) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvMUAxcVals;
            data.SSP_SAP.(behavField).adjLH.HbTvMUAxcVals_std(:,zz) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvMUAxcVals_std;
            data.SSP_SAP.(behavField).adjLH.MUA_lags(:,zz) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.LFP_lags;
            data.SSP_SAP.(behavField).adjRH.HbTvMUAxcVals(:,zz) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvMUAxcVals;
            data.SSP_SAP.(behavField).adjRH.HbTvMUAxcVals_std(:,zz) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvMUAxcVals_std;
            data.SSP_SAP.(behavField).adjRH.MUA_lags(:,zz) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.LFP_lags;
            data.SSP_SAP.(behavField).animalID{zz,1} = animalID;
            data.SSP_SAP.(behavField).behavior{zz,1} = behavField;
            data.SSP_SAP.(behavField).LH{zz,1} = 'LH';
            data.SSP_SAP.(behavField).RH{zz,1} = 'RH';
        end
        zz = zz + 1;
    end
end
%% concatenate the data from the left and right hemispheres
% for dd = 1:length(behavFields)
%     behavField = behavFields{1,dd};
%     data.(behavField).cat_HbTvLFPxcVals = cat(3,data.(behavField).adjLH.HbTvLFPxcVals, data.(behavField).adjRH.HbTvLFPxcVals);
%     data.(behavField).cat_LFP_lags = cat(3,data.(behavField).adjLH.LFP_lags, data.(behavField).adjRH.LFP_lags);
%     data.(behavField).cat_LFP_F = cat(3,data.(behavField).adjLH.F, data.(behavField).adjRH.F);
%     data.(behavField).cat_HbTvMUAxcVals = cat(2,data.(behavField).adjLH.HbTvMUAxcVals, data.(behavField).adjRH.HbTvMUAxcVals);
%     data.(behavField).cat_MUA_lags = cat(2,data.(behavField).adjLH.MUA_lags, data.(behavField).adjRH.MUA_lags);
% end
%% take the averages of each field through the proper dimension
for dd = 1:length(treatments)
    treatment = treatments{1,dd};
    for ff = 1:length(behavFields)
        behavField = behavFields{1,ff};
        % LH
        data.(treatment).(behavField).adjLH.meanHbTvLFPxcVals = mean(data.(treatment).(behavField).adjLH.HbTvLFPxcVals,3);
        data.(treatment).(behavField).adjLH.meanLFP_lags = mean(data.(treatment).(behavField).adjLH.LFP_lags,3);
        data.(treatment).(behavField).adjLH.meanLFP_F = mean(data.(treatment).(behavField).adjLH.F,3);
        data.(treatment).(behavField).adjLH.meanHbTvMUAxcVals = mean(data.(treatment).(behavField).adjLH.HbTvMUAxcVals,2);
        data.(treatment).(behavField).adjLH.stdHbTvMUAxcVals = std(data.(treatment).(behavField).adjLH.HbTvMUAxcVals,0,2);
        data.(treatment).(behavField).adjLH.meanMUA_lags = mean(data.(treatment).(behavField).adjLH.MUA_lags,2);
        % RH
        data.(treatment).(behavField).adjRH.meanHbTvLFPxcVals = mean(data.(treatment).(behavField).adjRH.HbTvLFPxcVals,3);
        data.(treatment).(behavField).adjRH.meanLFP_lags = mean(data.(treatment).(behavField).adjRH.LFP_lags,3);
        data.(treatment).(behavField).adjRH.meanLFP_F = mean(data.(treatment).(behavField).adjRH.F,3);
        data.(treatment).(behavField).adjRH.meanHbTvMUAxcVals = mean(data.(treatment).(behavField).adjRH.HbTvMUAxcVals,2);
        data.(treatment).(behavField).adjRH.stdHbTvMUAxcVals = std(data.(treatment).(behavField).adjRH.HbTvMUAxcVals,0,2);
        data.(treatment).(behavField).adjRH.meanMUA_lags = mean(data.(treatment).(behavField).adjRH.MUA_lags,2);
    end
end
%% find max/time to peak for MUA/Gamma-band power
% for gg = 1:length(behavFields)
%     behavField = behavFields{1,gg};
%     for hh = 1:size(data.(behavField).cat_HbTvLFPxcVals,3)
%         % gamma-band
%         LFPmat = data.(behavField).cat_HbTvLFPxcVals(:,:,hh);
%         gammaArray = mean(LFPmat(49:end,:),1);
%         [gammaMax,gammaIndex] = max(gammaArray);
%         data.(behavField).gammaPeak(hh,1) = gammaMax;
%         data.(behavField).gammaTTP(hh,1) = data.(behavField).meanLFP_lags(gammaIndex)/30;
%         % mua
%         muaArray = data.(behavField).cat_HbTvMUAxcVals(:,hh)';
%         [muaMax,muaIndex] = max(muaArray);
%         data.(behavField).muaPeak(hh,1) = muaMax;
%         data.(behavField).muaTTP(hh,1) = data.(behavField).meanMUA_lags(muaIndex)/30;
%     end
% end
%% statistics - generalized linear mixed effects model
% muaPeakTableSize = cat(1,data.Rest.muaPeak,data.NREM.muaPeak,data.REM.muaPeak);
% muaPeakTable = table('Size',[size(muaPeakTableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Peak','Behavior','Hemisphere'});
% muaPeakTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
% muaPeakTable.Peak = cat(1,data.Rest.muaPeak,data.NREM.muaPeak,data.REM.muaPeak);
% muaPeakTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
% muaPeakTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
% muaPeakFitFormula = 'Peak ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
% muaPeakStats = fitglme(muaPeakTable,muaPeakFitFormula); %#ok<*NASGU>
%% statistics - generalized linear mixed effects model
% muaTTPTableSize = cat(1,data.Rest.muaTTP,data.NREM.muaTTP,data.REM.muaTTP);
% muaTTPTable = table('Size',[size(muaTTPTableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Peak','Behavior','Hemisphere'});
% muaTTPTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
% muaTTPTable.TTP = cat(1,data.Rest.muaTTP,data.NREM.muaTTP,data.REM.muaTTP);
% muaTTPTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
% muaTTPTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
% muaTTPFitFormula = 'TTP ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
% muaTTPStats = fitglme(muaTTPTable,muaTTPFitFormula); %#ok<*NASGU>
%% statistics - generalized linear mixed effects model
% gammaPeakTableSize = cat(1,data.Rest.gammaPeak,data.NREM.gammaPeak,data.REM.gammaPeak);
% gammaPeakTable = table('Size',[size(gammaPeakTableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Peak','Behavior','Hemisphere'});
% gammaPeakTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
% gammaPeakTable.Peak = cat(1,data.Rest.gammaPeak,data.NREM.gammaPeak,data.REM.gammaPeak);
% gammaPeakTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
% gammaPeakTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
% gammaPeakFitFormula = 'Peak ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
% gammaPeakStats = fitglme(gammaPeakTable,gammaPeakFitFormula); %#ok<*NASGU>
%% statistics - generalized linear mixed effects model
% gammaTTPTableSize = cat(1,data.Rest.gammaTTP,data.NREM.gammaTTP,data.REM.gammaTTP);
% gammaTTPTable = table('Size',[size(gammaTTPTableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Peak','Behavior','Hemisphere'});
% gammaTTPTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
% gammaTTPTable.TTP = cat(1,data.Rest.gammaTTP,data.NREM.gammaTTP,data.REM.gammaTTP);
% gammaTTPTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
% gammaTTPTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
% gammaTTPFitFormula = 'TTP ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
% gammaTTPStats = fitglme(gammaTTPTable,gammaTTPFitFormula); %#ok<*NASGU>
%% Fig. 6
% summaryFigure = figure('Name','Fig6 (a-c)');
% sgtitle('Figure 6 - Turner et al. 2020')
%% [6a] rest MUA-HbT XCorr
freq = 30;
restLag = 5;
sleepLag = 5;
% ax1 = subplot(2,3,1);
% plot(data.Rest.meanMUA_lags,data.Rest.meanHbTvMUAxcVals,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.Rest.meanMUA_lags,data.Rest.meanHbTvMUAxcVals + data.Rest.stdHbTvMUAxcVals,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Rest.meanMUA_lags,data.Rest.meanHbTvMUAxcVals - data.Rest.stdHbTvMUAxcVals,'color',colors('battleship grey'),'LineWidth',0.5)
% title({'[6a] Awake Rest','MUA-[HbT] XCorr'})
% xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
% xticklabels({'-5','-2.5','0','2.5','5'})
% xlim([-restLag*freq,restLag*freq])
% xlabel('Lags (s)')
% ylabel({'Corr. coefficient';'MUA vs. \Delta[HbT] (\muM)'})
% axis square
% set(gca,'box','off')
% ax1.TickLength = [0.03,0.03];
% ylim([-0.1,0.5])
% %% [6b] NREM MUA-HbT XCorr
% ax2 = subplot(2,3,2);
% plot(data.NREM.meanMUA_lags,data.NREM.meanHbTvMUAxcVals,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.NREM.meanMUA_lags,data.NREM.meanHbTvMUAxcVals + data.NREM.stdHbTvMUAxcVals,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.NREM.meanMUA_lags,data.NREM.meanHbTvMUAxcVals - data.NREM.stdHbTvMUAxcVals,'color',colors('battleship grey'),'LineWidth',0.5)
% title({'[6b] NREM','MUA-[HbT] XCorr'})
% xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
% xticklabels({'-5','-2.5','0','2.5','5'})
% xlim([-sleepLag*freq,sleepLag*freq])
% xlabel('Lags (s)')
% ylabel({'Corr. coefficient';'MUA vs. \Delta[HbT] (\muM)'})
% axis square
% set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
% ylim([-0.1,0.5])
% %% [6c] REM MUA-HbT XCorr
% ax3 = subplot(2,3,3);
% plot(data.REM.meanMUA_lags,data.REM.meanHbTvMUAxcVals,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.REM.meanMUA_lags,data.REM.meanHbTvMUAxcVals + data.REM.stdHbTvMUAxcVals,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.REM.meanMUA_lags,data.REM.meanHbTvMUAxcVals - data.REM.stdHbTvMUAxcVals,'color',colors('battleship grey'),'LineWidth',0.5)
% title({'[6c] REM','MUA-[HbT] XCorr'})
% xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
% xticklabels({'-5','-2.5','0','2.5','5'})
% xlim([-sleepLag*freq,sleepLag*freq])
% xlabel('Lags (s)')
% ylabel({'Corr. coefficient';'MUA vs. \Delta[HbT] (\muM)'})
% axis square
% set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
% ylim([-0.1,0.5])
fig1 = figure;   %('Name','Fig6 (a-c)');
% sgtitle('Figure 6 - Turner et al. 2020')
%% [6a bottom] rest LFP-HbT XCorr
ax1 = subplot(2,2,1);
imagesc(data.C57BL6J.Rest.adjLH.meanLFP_lags,data.C57BL6J.Rest.adjLH.meanLFP_F,data.C57BL6J.Rest.adjLH.meanHbTvLFPxcVals)
title({'Control LH Awake Rest','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c1 = colorbar;
caxis([-0.1,0.1])
ylabel(c1,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [6a bottom] rest LFP-HbT XCorr
ax2 = subplot(2,2,2);
imagesc(data.C57BL6J.Rest.adjRH.meanLFP_lags,data.C57BL6J.Rest.adjRH.meanLFP_F,data.C57BL6J.Rest.adjRH.meanHbTvLFPxcVals)
title({'Control RH Awake Rest','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c2 = colorbar;
caxis([-0.1,0.1])
ylabel(c2,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [6a bottom] rest LFP-HbT XCorr
ax7 = subplot(2,2,3);
imagesc(data.SSP_SAP.Rest.adjLH.meanLFP_lags,data.SSP_SAP.Rest.adjLH.meanLFP_F,data.SSP_SAP.Rest.adjLH.meanHbTvLFPxcVals)
title({'SSP-SAP UnRx LH Awake Rest','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c7 = colorbar;
caxis([-0.1,0.1])
ylabel(c7,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% [6a bottom] rest LFP-HbT XCorr
ax8 = subplot(2,2,4);
imagesc(data.SSP_SAP.Rest.adjRH.meanLFP_lags,data.SSP_SAP.Rest.adjRH.meanLFP_F,data.SSP_SAP.Rest.adjRH.meanHbTvLFPxcVals)
title({'SSP-SAP treated RH Awake Rest','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c8 = colorbar;
caxis([-0.1,0.1])
ylabel(c8,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%
fig2 = figure;
%% [6a bottom] rest LFP-HbT XCorr
ax3 = subplot(2,2,1);
imagesc(data.C57BL6J.NREM.adjLH.meanLFP_lags,data.C57BL6J.NREM.adjLH.meanLFP_F,data.C57BL6J.NREM.adjLH.meanHbTvLFPxcVals)
title({'Control LH NREM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c3 = colorbar;
caxis([-0.4,0.4])
ylabel(c3,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [6a bottom] rest LFP-HbT XCorr
ax4 = subplot(2,2,2);
imagesc(data.C57BL6J.NREM.adjRH.meanLFP_lags,data.C57BL6J.NREM.adjRH.meanLFP_F,data.C57BL6J.NREM.adjRH.meanHbTvLFPxcVals)
title({'Control RH NREM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c4 = colorbar;
caxis([-0.4,0.4])
ylabel(c4,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [6a bottom] rest LFP-HbT XCorr
ax9 = subplot(2,2,3);
imagesc(data.SSP_SAP.NREM.adjLH.meanLFP_lags,data.SSP_SAP.NREM.adjLH.meanLFP_F,data.SSP_SAP.NREM.adjLH.meanHbTvLFPxcVals)
title({'SSP-SAP UnRx LH NREM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c9 = colorbar;
caxis([-0.4,0.4])
ylabel(c9,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [6a bottom] rest LFP-HbT XCorr
ax10 = subplot(2,2,4);
imagesc(data.SSP_SAP.NREM.adjRH.meanLFP_lags,data.SSP_SAP.NREM.adjRH.meanLFP_F,data.SSP_SAP.NREM.adjRH.meanHbTvLFPxcVals)
title({'SSP-SAP treated RH NREM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c10 = colorbar;
caxis([-0.4,0.4])
ylabel(c10,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax10.TickLength = [0.03,0.03];
%
fig3 = figure;
%% [6a bottom] rest LFP-HbT XCorr
ax5 = subplot(2,2,1);
imagesc(data.C57BL6J.REM.adjLH.meanLFP_lags,data.C57BL6J.REM.adjLH.meanLFP_F,data.C57BL6J.REM.adjLH.meanHbTvLFPxcVals)
title({'Control LH REM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c5 = colorbar;
caxis([-0.2,0.2])
ylabel(c5,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [6a bottom] rest LFP-HbT XCorr
ax6 = subplot(2,2,2);
imagesc(data.C57BL6J.REM.adjRH.meanLFP_lags,data.C57BL6J.REM.adjRH.meanLFP_F,data.C57BL6J.REM.adjRH.meanHbTvLFPxcVals)
title({'Control RH REM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c6 = colorbar;
caxis([-0.2,0.2])
ylabel(c6,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [6a bottom] rest LFP-HbT XCorr
ax11 = subplot(2,2,3);
imagesc(data.SSP_SAP.REM.adjLH.meanLFP_lags,data.SSP_SAP.REM.adjLH.meanLFP_F,data.SSP_SAP.REM.adjLH.meanHbTvLFPxcVals)
title({'SSP-SAP UnRx LH REM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c11 = colorbar;
caxis([-0.2,0.2])
ylabel(c11,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
%% [6a bottom] rest LFP-HbT XCorr
ax12 = subplot(2,2,4);
imagesc(data.SSP_SAP.REM.adjRH.meanLFP_lags,data.SSP_SAP.REM.adjRH.meanLFP_F,data.SSP_SAP.REM.adjRH.meanHbTvLFPxcVals)
title({'SSP-SAP treated RH REM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c12 = colorbar;
caxis([-0.2,0.2])
ylabel(c12,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];
%% link axes and adjust positions
% linkaxes([ax1,ax2,ax3],'y')
% ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
% ax3Pos = get(ax3,'position');
% ax4Pos = get(ax1,'position');
% ax5Pos = get(ax1,'position');
% ax6Pos = get(ax6,'position');
% ax4Pos(3:4) = ax1Pos(3:4);
% ax5Pos(3:4) = ax2Pos(3:4);
% ax6Pos(3:4) = ax3Pos(3:4);
% set(ax1,'position',ax4Pos);
% set(ax1,'position',ax5Pos);
% set(ax6,'position',ax6Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
%     savefig(summaryFigure, [dirpath 'Fig6_Test']);
    savefig(fig1, [dirpath 'Fig6_Test']);
    savefig(fig2, [dirpath 'Fig6_Test']);
    savefig(fig3, [dirpath 'Fig6_Test']);
% 
%     set(summaryFigure,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-bestfit',[dirpath 'Fig6'])
end

end
