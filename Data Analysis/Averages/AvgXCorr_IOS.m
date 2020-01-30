%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Across-animal average of the cross correlation between neural activity and lagged hemodynamics between
%            different arousal states.
%________________________________________________________________________________________________________________________

clear
clc

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111'};
driveLetters = {'M','M','M','M','M','M','M','M','M'};
behavFields = {'Rest','NREM','REM'};

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\Turner_Manuscript_Summer2020\' animalID '\Bilateral Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        data.(behavField).adjLH.HbTvLFPxcVals(:,:,a) = AnalysisResults.XCorr.(behavField).adjLH.HbTvLFPxcVals;
        data.(behavField).adjLH.LFP_lags(:,:,a) = AnalysisResults.XCorr.(behavField).adjLH.LFP_lags;
        data.(behavField).adjLH.F(:,:,a) = AnalysisResults.XCorr.(behavField).adjLH.F;
        data.(behavField).adjRH.HbTvLFPxcVals(:,:,a) = AnalysisResults.XCorr.(behavField).adjRH.HbTvLFPxcVals;
        data.(behavField).adjRH.LFP_lags(:,:,a) = AnalysisResults.XCorr.(behavField).adjRH.LFP_lags;
        data.(behavField).adjRH.F(:,:,a) = AnalysisResults.XCorr.(behavField).adjRH.F;
        data.(behavField).adjLH.HbTvMUAxcVals(:,a) = AnalysisResults.XCorr.(behavField).adjLH.HbTvMUAxcVals;
        data.(behavField).adjLH.HbTvMUAxcVals_std(:,a) = AnalysisResults.XCorr.(behavField).adjLH.HbTvMUAxcVals_std;
        data.(behavField).adjLH.MUA_lags(:,a) = AnalysisResults.XCorr.(behavField).adjLH.LFP_lags;
        data.(behavField).adjRH.HbTvMUAxcVals(:,a) = AnalysisResults.XCorr.(behavField).adjRH.HbTvMUAxcVals;
        data.(behavField).adjRH.HbTvMUAxcVals_std(:,a) = AnalysisResults.XCorr.(behavField).adjRH.HbTvMUAxcVals_std;
        data.(behavField).adjRH.MUA_lags(:,a) = AnalysisResults.XCorr.(behavField).adjRH.LFP_lags;
    end
end

% concatenate the data from the left and right hemispheres
for d = 1:length(behavFields)
    behavField = behavFields{1,d};
    data.(behavField).cat_HbTvLFPxcVals = cat(3,data.(behavField).adjLH.HbTvLFPxcVals, data.(behavField).adjRH.HbTvLFPxcVals);
    data.(behavField).cat_LFP_lags = cat(3,data.(behavField).adjLH.LFP_lags, data.(behavField).adjRH.LFP_lags);
    data.(behavField).cat_LFP_F = cat(3,data.(behavField).adjLH.F, data.(behavField).adjRH.F);
    data.(behavField).cat_HbTvMUAxcVals = cat(2,data.(behavField).adjLH.HbTvMUAxcVals, data.(behavField).adjRH.HbTvMUAxcVals);
    data.(behavField).cat_MUA_lags = cat(2,data.(behavField).adjLH.MUA_lags, data.(behavField).adjRH.MUA_lags);
end

% take the averages of each field through the proper dimension
for f = 1:length(behavFields)
    behavField = behavFields{1,f};
    data.(behavField).meanHbTvLFPxcVals = mean(data.(behavField).cat_HbTvLFPxcVals,3);
    data.(behavField).meanLFP_lags = mean(data.(behavField).cat_LFP_lags,3);
    data.(behavField).meanLFP_F = mean(data.(behavField).cat_LFP_F,3);
    data.(behavField).meanHbTvMUAxcVals = mean(data.(behavField).cat_HbTvMUAxcVals,2);
    data.(behavField).stdHbTvMUAxcVals = std(data.(behavField).cat_HbTvMUAxcVals,0,2);
    data.(behavField).meanMUA_lags = mean(data.(behavField).cat_MUA_lags,2);
end

%% summary figure(s)
% lagTime2 = 15;
summaryFigure = figure;
sgtitle('Neural-Hemo Cross-Correlations')

%% Rest MUA
freq = 10;
restLag = 5;
sleepLag = 5;
ax1 = subplot(2,3,1);
plot(data.Rest.meanMUA_lags,data.Rest.meanHbTvMUAxcVals,'k')
hold on
plot(data.Rest.meanMUA_lags,data.Rest.meanHbTvMUAxcVals + data.Rest.stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
plot(data.Rest.meanMUA_lags,data.Rest.meanHbTvMUAxcVals - data.Rest.stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
title('Awake Rest MUA')
xticks([-restLag*freq -restLag*freq/2 0 restLag*freq/2 restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq restLag*freq])
xlabel('Lags (s)')
ylabel({'Corr. Coefficient';'MUA vs. \DeltaHbT (\muM)'})
axis square
set(gca,'box','off')

%% NREM MUA
ax2 = subplot(2,3,2);
plot(data.NREM.meanMUA_lags,data.NREM.meanHbTvMUAxcVals,'k')
hold on
plot(data.NREM.meanMUA_lags,data.NREM.meanHbTvMUAxcVals + data.NREM.stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
plot(data.NREM.meanMUA_lags,data.NREM.meanHbTvMUAxcVals - data.NREM.stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
title('NREM MUA')
xticks([-sleepLag*freq -sleepLag*freq/2 0 sleepLag*freq/2 sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq sleepLag*freq])
xlabel('Lags (s)')
ylabel({'Corr. Coefficient';'MUA vs. \DeltaHbT (\muM)'})
axis square
set(gca,'box','off')

%% REM MUA
ax3 = subplot(2,3,3);
plot(data.REM.meanMUA_lags,data.REM.meanHbTvMUAxcVals,'k')
hold on
plot(data.REM.meanMUA_lags,data.REM.meanHbTvMUAxcVals + data.REM.stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
plot(data.REM.meanMUA_lags,data.REM.meanHbTvMUAxcVals - data.REM.stdHbTvMUAxcVals,'color',colors_IOS('battleship grey'))
title('REM MUA')
xticks([-sleepLag*freq -sleepLag*freq/2 0 sleepLag*freq/2 sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq sleepLag*freq])
xlabel('Lags (s)')
ylabel({'Corr. Coefficient';'MUA vs. \DeltaHbT (\muM)'})
axis square
set(gca,'box','off')

%% Rest LFP
ax4 = subplot(2,3,4);
imagesc(data.Rest.meanLFP_lags,data.Rest.meanLFP_F,data.Rest.meanHbTvLFPxcVals)
title('Awake Rest LFP')
xticks([-restLag*freq -restLag*freq/2 0 restLag*freq/2 restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq restLag*freq])
xlabel('Lags (s)')
ylabel('Frequency (Hz)')
ylim([1 100])
c4 = colorbar;
ylabel(c4,{'Corr. Coefficient';'LFP vs. \DeltaHbT (\muM)'})
axis xy
axis square
set(gca,'box','off')

%% NREM LFP
ax5 = subplot(2,3,5);
imagesc(data.NREM.meanLFP_lags,data.NREM.meanLFP_F,data.NREM.meanHbTvLFPxcVals)
title('NREM LFP')
xticks([-sleepLag*freq -sleepLag*freq/2 0 sleepLag*freq/2 sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq sleepLag*freq])
xlabel('Lags (s)')
ylabel('Frequency (Hz)')
ylim([1 100])
c5 = colorbar;
ylabel(c5,{'Corr. Coefficient';'LFP vs. \DeltaHbT (\muM)'})
axis xy
axis square
set(gca,'box','off')

%% REM LFP
ax6 = subplot(2,3,6);
imagesc(data.REM.meanLFP_lags,data.REM.meanLFP_F,data.REM.meanHbTvLFPxcVals)
title('REM LFP')
xticks([-sleepLag*freq -sleepLag*freq/2 0 sleepLag*freq/2 sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq sleepLag*freq])
xlabel('Lags (s)')
ylabel('Frequency (Hz)')
ylim([1 100])
c6 = colorbar;
ylabel(c6,{'Corr. Coefficient';'LFP vs. \DeltaHbT (\muM)'})
axis xy
axis square
set(gca,'box','off')

linkaxes([ax1 ax2 ax3],'y')

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

% save figure(s)
dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\';
if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end
savefig(summaryFigure, [dirpath 'Summary Figure - Cross Correlation']);

