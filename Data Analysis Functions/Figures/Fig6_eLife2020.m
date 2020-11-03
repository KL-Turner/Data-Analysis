function [AnalysisResults] = Fig6_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generate figure panel 6 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Rest','NREM','REM'};
%% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        data.(behavField).adjLH.HbTvLFPxcVals(:,:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvLFPxcVals;
        data.(behavField).adjLH.LFP_lags(:,:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.LFP_lags;
        data.(behavField).adjLH.F(:,:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.F;
        data.(behavField).adjRH.HbTvLFPxcVals(:,:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvLFPxcVals;
        data.(behavField).adjRH.LFP_lags(:,:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.LFP_lags;
        data.(behavField).adjRH.F(:,:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.F;
        data.(behavField).adjLH.HbTvMUAxcVals(:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvMUAxcVals;
        data.(behavField).adjLH.HbTvMUAxcVals_std(:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvMUAxcVals_std;
        data.(behavField).adjLH.MUA_lags(:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.LFP_lags;
        data.(behavField).adjRH.HbTvMUAxcVals(:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvMUAxcVals;
        data.(behavField).adjRH.HbTvMUAxcVals_std(:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvMUAxcVals_std;
        data.(behavField).adjRH.MUA_lags(:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.LFP_lags;
        data.(behavField).animalID{aa,1} = animalID;
        data.(behavField).behavior{aa,1} = behavField;
        data.(behavField).LH{aa,1} = 'LH';
        data.(behavField).RH{aa,1} = 'RH';
    end
end
% concatenate the data from the left and right hemispheres
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    data.(behavField).cat_HbTvLFPxcVals = cat(3,data.(behavField).adjLH.HbTvLFPxcVals, data.(behavField).adjRH.HbTvLFPxcVals);
    data.(behavField).cat_LFP_lags = cat(3,data.(behavField).adjLH.LFP_lags, data.(behavField).adjRH.LFP_lags);
    data.(behavField).cat_LFP_F = cat(3,data.(behavField).adjLH.F, data.(behavField).adjRH.F);
    data.(behavField).cat_HbTvMUAxcVals = cat(2,data.(behavField).adjLH.HbTvMUAxcVals, data.(behavField).adjRH.HbTvMUAxcVals);
    data.(behavField).cat_MUA_lags = cat(2,data.(behavField).adjLH.MUA_lags, data.(behavField).adjRH.MUA_lags);
end
% take the averages of each field through the proper dimension
for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.(behavField).meanHbTvLFPxcVals = mean(data.(behavField).cat_HbTvLFPxcVals,3);
    data.(behavField).meanLFP_lags = mean(data.(behavField).cat_LFP_lags,3);
    data.(behavField).meanLFP_F = mean(data.(behavField).cat_LFP_F,3);
    data.(behavField).meanHbTvMUAxcVals = mean(data.(behavField).cat_HbTvMUAxcVals,2);
    data.(behavField).stdHbTvMUAxcVals = std(data.(behavField).cat_HbTvMUAxcVals,0,2);
    data.(behavField).meanMUA_lags = mean(data.(behavField).cat_MUA_lags,2);
end
%% find max/time to peak for MUA/Gamma-band power
for gg = 1:length(behavFields)
    behavField = behavFields{1,gg};
    for hh = 1:size(data.(behavField).cat_HbTvLFPxcVals,3)
        % gamma-band
        LFPmat = data.(behavField).cat_HbTvLFPxcVals(:,:,hh);
        gammaArray = mean(LFPmat(49:end,:),1);
        [gammaMax,gammaIndex] = max(gammaArray);
        data.(behavField).gammaPeak(hh,1) = gammaMax;
        data.(behavField).gammaTTP(hh,1) = data.(behavField).meanLFP_lags(gammaIndex)/30;
        % mua
        muaArray = data.(behavField).cat_HbTvMUAxcVals(:,hh)';
        [muaMax,muaIndex] = max(muaArray);
        data.(behavField).muaPeak(hh,1) = muaMax;
        data.(behavField).muaTTP(hh,1) = data.(behavField).meanMUA_lags(muaIndex)/30;
    end
end
%% statistics - generalized linear mixed effects model
muaPeakTableSize = cat(1,data.Rest.muaPeak,data.NREM.muaPeak,data.REM.muaPeak);
muaPeakTable = table('Size',[size(muaPeakTableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Peak','Behavior','Hemisphere'});
muaPeakTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
muaPeakTable.Peak = cat(1,data.Rest.muaPeak,data.NREM.muaPeak,data.REM.muaPeak);
muaPeakTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
muaPeakTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
muaPeakFitFormula = 'Peak ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
muaPeakStats = fitglme(muaPeakTable,muaPeakFitFormula); %#ok<*NASGU>
%% statistics - generalized linear mixed effects model
muaTTPTableSize = cat(1,data.Rest.muaTTP,data.NREM.muaTTP,data.REM.muaTTP);
muaTTPTable = table('Size',[size(muaTTPTableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Peak','Behavior','Hemisphere'});
muaTTPTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
muaTTPTable.TTP = cat(1,data.Rest.muaTTP,data.NREM.muaTTP,data.REM.muaTTP);
muaTTPTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
muaTTPTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
muaTTPFitFormula = 'TTP ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
muaTTPStats = fitglme(muaTTPTable,muaTTPFitFormula); %#ok<*NASGU>
%% statistics - generalized linear mixed effects model
gammaPeakTableSize = cat(1,data.Rest.gammaPeak,data.NREM.gammaPeak,data.REM.gammaPeak);
gammaPeakTable = table('Size',[size(gammaPeakTableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Peak','Behavior','Hemisphere'});
gammaPeakTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
gammaPeakTable.Peak = cat(1,data.Rest.gammaPeak,data.NREM.gammaPeak,data.REM.gammaPeak);
gammaPeakTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
gammaPeakTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
gammaPeakFitFormula = 'Peak ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
gammaPeakStats = fitglme(gammaPeakTable,gammaPeakFitFormula); %#ok<*NASGU>
%% statistics - generalized linear mixed effects model
gammaTTPTableSize = cat(1,data.Rest.gammaTTP,data.NREM.gammaTTP,data.REM.gammaTTP);
gammaTTPTable = table('Size',[size(gammaTTPTableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Peak','Behavior','Hemisphere'});
gammaTTPTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
gammaTTPTable.TTP = cat(1,data.Rest.gammaTTP,data.NREM.gammaTTP,data.REM.gammaTTP);
gammaTTPTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
gammaTTPTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
gammaTTPFitFormula = 'TTP ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
gammaTTPStats = fitglme(gammaTTPTable,gammaTTPFitFormula); %#ok<*NASGU>
%% Fig. 6
summaryFigure = figure('Name','Fig6 (a-c)');
sgtitle('Figure 6 - Turner et al. 2020')
%% [6a] rest MUA-HbT XCorr
freq = 30;
restLag = 5;
sleepLag = 5;
ax1 = subplot(2,3,1);
plot(data.Rest.meanMUA_lags,data.Rest.meanHbTvMUAxcVals,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.Rest.meanMUA_lags,data.Rest.meanHbTvMUAxcVals + data.Rest.stdHbTvMUAxcVals,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Rest.meanMUA_lags,data.Rest.meanHbTvMUAxcVals - data.Rest.stdHbTvMUAxcVals,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title({'[6a] Awake Rest','MUA-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel({'Corr. coefficient';'MUA vs. \Delta[HbT] (\muM)'})
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ylim([-0.1,0.5])
%% [6b] NREM MUA-HbT XCorr
ax2 = subplot(2,3,2);
plot(data.NREM.meanMUA_lags,data.NREM.meanHbTvMUAxcVals,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.NREM.meanMUA_lags,data.NREM.meanHbTvMUAxcVals + data.NREM.stdHbTvMUAxcVals,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.NREM.meanMUA_lags,data.NREM.meanHbTvMUAxcVals - data.NREM.stdHbTvMUAxcVals,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title({'[6b] NREM','MUA-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel({'Corr. coefficient';'MUA vs. \Delta[HbT] (\muM)'})
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
ylim([-0.1,0.5])
%% [6c] REM MUA-HbT XCorr
ax3 = subplot(2,3,3);
plot(data.REM.meanMUA_lags,data.REM.meanHbTvMUAxcVals,'color',colors_eLife2020('rich black'),'LineWidth',1)
hold on
plot(data.REM.meanMUA_lags,data.REM.meanHbTvMUAxcVals + data.REM.stdHbTvMUAxcVals,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.REM.meanMUA_lags,data.REM.meanHbTvMUAxcVals - data.REM.stdHbTvMUAxcVals,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title({'[6c] REM','MUA-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel({'Corr. coefficient';'MUA vs. \Delta[HbT] (\muM)'})
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
ylim([-0.1,0.5])
%% [6a bottom] rest LFP-HbT XCorr
ax4 = subplot(2,3,4);
imagesc(data.Rest.meanLFP_lags,data.Rest.meanLFP_F,data.Rest.meanHbTvLFPxcVals)
title({'Awake Rest','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c4 = colorbar;
caxis([-0.2,0.2])
ylabel(c4,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [6b bottom] NREM LFP-HbT XCorr
ax5 = subplot(2,3,5);
imagesc(data.NREM.meanLFP_lags,data.NREM.meanLFP_F,data.NREM.meanHbTvLFPxcVals)
title({'NREM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c5 = colorbar;
caxis([-0.4,0.4])
ylabel(c5,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [6c bottom] REM LFP-HbT XCorr
ax6 = subplot(2,3,6);
imagesc(data.REM.meanLFP_lags,data.REM.meanLFP_F,data.REM.meanHbTvLFPxcVals)
title({'REM','LFP-[HbT] XCorr'})
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
%% link axes and adjust positions
linkaxes([ax1,ax2,ax3],'y')
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax2Pos(3:4);
ax6Pos(3:4) = ax3Pos(3:4);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure, [dirpath 'Fig6']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig6'])
end

end
