%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

clear
clc

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111'};
driveLetters = {'E','E','E','F','F','F','D','D','D'};
behavFields = {'Rest','NREM','REM','Unstim','All'};
baselineTypes = {'manualSelection','setDuration','entireDuration'};

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        if strcmp(behavField,'Rest') == true
            for c = 1:length(baselineTypes)
                baselineType = baselineTypes{1,c};
                data.(behavField).(baselineType).LH.CBVvLFPxcVals(:,:,a) = AnalysisResults.XCorr.(behavField).LH.(baselineType).CBVvLFPxcVals;
                data.(behavField).(baselineType).LH.HbTvLFPxcVals(:,:,a) = AnalysisResults.XCorr.(behavField).LH.(baselineType).HbTvLFPxcVals;
                data.(behavField).(baselineType).LH.LFP_lags(:,:,a) = AnalysisResults.XCorr.(behavField).LH.(baselineType).LFP_lags;
                data.(behavField).(baselineType).LH.F(:,:,a) = AnalysisResults.XCorr.(behavField).LH.(baselineType).F;
                
                data.(behavField).(baselineType).LH.CBVvMUAxcVals(:,a) = AnalysisResults.XCorr.(behavField).LH.(baselineType).CBVvMUAxcVals;
                data.(behavField).(baselineType).LH.CBVvMUAxcVals_std(:,a) = AnalysisResults.XCorr.(behavField).LH.(baselineType).CBVvMUAxcVals_std;
                data.(behavField).(baselineType).LH.HbTvMUAxcVals(:,a) = AnalysisResults.XCorr.(behavField).LH.(baselineType).HbTvMUAxcVals;
                data.(behavField).(baselineType).LH.HbTvMUAxcVals_std(:,a) = AnalysisResults.XCorr.(behavField).LH.(baselineType).HbTvMUAxcVals_std;
                data.(behavField).(baselineType).LH.MUA_lags(:,a) = AnalysisResults.XCorr.(behavField).LH.(baselineType).LFP_lags;
                
                data.(behavField).(baselineType).RH.CBVvLFPxcVals(:,:,a) = AnalysisResults.XCorr.(behavField).RH.(baselineType).CBVvLFPxcVals;
                data.(behavField).(baselineType).RH.HbTvLFPxcVals(:,:,a) = AnalysisResults.XCorr.(behavField).RH.(baselineType).HbTvLFPxcVals;
                data.(behavField).(baselineType).RH.LFP_lags(:,:,a) = AnalysisResults.XCorr.(behavField).RH.(baselineType).LFP_lags;
                data.(behavField).(baselineType).RH.F(:,:,a) = AnalysisResults.XCorr.(behavField).RH.(baselineType).F;
                
                data.(behavField).(baselineType).RH.CBVvMUAxcVals(:,a) = AnalysisResults.XCorr.(behavField).RH.(baselineType).CBVvMUAxcVals;
                data.(behavField).(baselineType).RH.CBVvMUAxcVals_std(:,a) = AnalysisResults.XCorr.(behavField).RH.(baselineType).CBVvMUAxcVals_std;
                data.(behavField).(baselineType).RH.HbTvMUAxcVals(:,a) = AnalysisResults.XCorr.(behavField).RH.(baselineType).HbTvMUAxcVals;
                data.(behavField).(baselineType).RH.HbTvMUAxcVals_std(:,a) = AnalysisResults.XCorr.(behavField).RH.(baselineType).HbTvMUAxcVals_std;
                data.(behavField).(baselineType).RH.MUA_lags(:,a) = AnalysisResults.XCorr.(behavField).RH.(baselineType).LFP_lags;
            end
        else
                data.(behavField).LH.CBVvLFPxcVals(:,:,a) = AnalysisResults.XCorr.(behavField).LH.CBVvLFPxcVals;
                data.(behavField).LH.HbTvLFPxcVals(:,:,a) = AnalysisResults.XCorr.(behavField).LH.HbTvLFPxcVals;
                data.(behavField).LH.LFP_lags(:,:,a) = AnalysisResults.XCorr.(behavField).LH.LFP_lags;
                data.(behavField).LH.F(:,:,a) = AnalysisResults.XCorr.(behavField).LH.F;
                
                data.(behavField).LH.CBVvMUAxcVals(:,a) = AnalysisResults.XCorr.(behavField).LH.CBVvMUAxcVals;
                data.(behavField).LH.CBVvMUAxcVals_std(:,a) = AnalysisResults.XCorr.(behavField).LH.CBVvMUAxcVals_std;
                data.(behavField).LH.HbTvMUAxcVals(:,a) = AnalysisResults.XCorr.(behavField).LH.HbTvMUAxcVals;
                data.(behavField).LH.HbTvMUAxcVals_std(:,a) = AnalysisResults.XCorr.(behavField).LH.HbTvMUAxcVals_std;
                data.(behavField).LH.MUA_lags(:,a) = AnalysisResults.XCorr.(behavField).LH.LFP_lags;
                
                data.(behavField).RH.CBVvLFPxcVals(:,:,a) = AnalysisResults.XCorr.(behavField).RH.CBVvLFPxcVals;
                data.(behavField).RH.HbTvLFPxcVals(:,:,a) = AnalysisResults.XCorr.(behavField).RH.HbTvLFPxcVals;
                data.(behavField).RH.LFP_lags(:,:,a) = AnalysisResults.XCorr.(behavField).RH.LFP_lags;
                data.(behavField).RH.F(:,:,a) = AnalysisResults.XCorr.(behavField).RH.F;
                
                data.(behavField).RH.CBVvMUAxcVals(:,a) = AnalysisResults.XCorr.(behavField).RH.CBVvMUAxcVals;
                data.(behavField).RH.CBVvMUAxcVals_std(:,a) = AnalysisResults.XCorr.(behavField).RH.CBVvMUAxcVals_std;
                data.(behavField).RH.HbTvMUAxcVals(:,a) = AnalysisResults.XCorr.(behavField).RH.HbTvMUAxcVals;
                data.(behavField).RH.HbTvMUAxcVals_std(:,a) = AnalysisResults.XCorr.(behavField).RH.HbTvMUAxcVals_std;
                data.(behavField).RH.MUA_lags(:,a) = AnalysisResults.XCorr.(behavField).RH.LFP_lags;
        end
    end
end

% concatenate the data from the left and right hemispheres
for d = 1:length(behavFields)
    behavField = behavFields{1,d};
    if strcmp(behavField,'Rest') == true
        for e = 1:length(baselineTypes)
            baselineType = baselineTypes{1,e};
            data.(behavField).(baselineType).cat_CBVvLFPxcVals = cat(3,data.(behavField).(baselineType).LH.CBVvLFPxcVals, data.(behavField).(baselineType).RH.CBVvLFPxcVals);
            data.(behavField).(baselineType).cat_HbTvLFPxcVals = cat(3,data.(behavField).(baselineType).LH.HbTvLFPxcVals, data.(behavField).(baselineType).RH.HbTvLFPxcVals);
            data.(behavField).(baselineType).cat_LFP_lags = cat(3,data.(behavField).(baselineType).LH.LFP_lags, data.(behavField).(baselineType).RH.LFP_lags);
            data.(behavField).(baselineType).cat_LFP_F = cat(3,data.(behavField).(baselineType).LH.F, data.(behavField).(baselineType).RH.F);
            
            data.(behavField).(baselineType).cat_CBVvMUAxcVals = cat(2,data.(behavField).(baselineType).LH.CBVvMUAxcVals, data.(behavField).(baselineType).RH.CBVvMUAxcVals);
            data.(behavField).(baselineType).cat_HbTvMUAxcVals = cat(2,data.(behavField).(baselineType).LH.HbTvMUAxcVals, data.(behavField).(baselineType).RH.HbTvMUAxcVals);
            data.(behavField).(baselineType).cat_MUA_lags = cat(2,data.(behavField).(baselineType).LH.MUA_lags, data.(behavField).(baselineType).RH.MUA_lags); 
        end
    else
        data.(behavField).cat_CBVvLFPxcVals = cat(3,data.(behavField).LH.CBVvLFPxcVals, data.(behavField).RH.CBVvLFPxcVals);
        data.(behavField).cat_HbTvLFPxcVals = cat(3,data.(behavField).LH.HbTvLFPxcVals, data.(behavField).RH.HbTvLFPxcVals);
        data.(behavField).cat_LFP_lags = cat(3,data.(behavField).LH.LFP_lags, data.(behavField).RH.LFP_lags);
        data.(behavField).cat_LFP_F = cat(3,data.(behavField).LH.F, data.(behavField).RH.F);
        
        data.(behavField).cat_CBVvMUAxcVals = cat(2,data.(behavField).LH.CBVvMUAxcVals, data.(behavField).RH.CBVvMUAxcVals);
        data.(behavField).cat_HbTvMUAxcVals = cat(2,data.(behavField).LH.HbTvMUAxcVals, data.(behavField).RH.HbTvMUAxcVals);
        data.(behavField).cat_MUA_lags = cat(2,data.(behavField).LH.MUA_lags, data.(behavField).RH.MUA_lags);
    end
end

% take the averages of each field through the proper dimension
for f = 1:length(behavFields)
    behavField = behavFields{1,f};
    if strcmp(behavField,'Rest') == true
        for g = 1:length(baselineTypes)
            baselineType = baselineTypes{1,g};
            data.(behavField).(baselineType).meanCBVvLFPxcVals = mean(data.(behavField).(baselineType).cat_CBVvLFPxcVals,3);
            data.(behavField).(baselineType).meanHbTvLFPxcVals = mean(data.(behavField).(baselineType).cat_HbTvLFPxcVals,3);
            data.(behavField).(baselineType).meanLFP_lags = mean(data.(behavField).(baselineType).cat_LFP_lags,3);
            data.(behavField).(baselineType).meanLFP_F = mean(data.(behavField).(baselineType).cat_LFP_F,3);
            
            data.(behavField).(baselineType).meanCBVvMUAxcVals = mean(data.(behavField).(baselineType).cat_CBVvMUAxcVals,2);
            data.(behavField).(baselineType).stdCBVvMUAxcVals = std(data.(behavField).(baselineType).cat_CBVvMUAxcVals,0,2);
            data.(behavField).(baselineType).meanHbTvMUAxcVals = mean(data.(behavField).(baselineType).cat_HbTvMUAxcVals,2);
            data.(behavField).(baselineType).stdHbTvMUAxcVals = std(data.(behavField).(baselineType).cat_HbTvMUAxcVals,0,2);
            data.(behavField).(baselineType).meanMUA_lags = mean(data.(behavField).(baselineType).cat_MUA_lags,2);
        end
    else
        data.(behavField).meanCBVvLFPxcVals = mean(data.(behavField).cat_CBVvLFPxcVals,3);
        data.(behavField).meanHbTvLFPxcVals = mean(data.(behavField).cat_HbTvLFPxcVals,3);
        data.(behavField).meanLFP_lags = mean(data.(behavField).cat_LFP_lags,3);
        data.(behavField).meanLFP_F = mean(data.(behavField).cat_LFP_F,3);
        
        data.(behavField).meanCBVvMUAxcVals = mean(data.(behavField).cat_CBVvMUAxcVals,2);
        data.(behavField).stdCBVvMUAxcVals = std(data.(behavField).cat_CBVvMUAxcVals,0,2);
        data.(behavField).meanHbTvMUAxcVals = mean(data.(behavField).cat_HbTvMUAxcVals,2);
        data.(behavField).stdHbTvMUAxcVals = std(data.(behavField).cat_HbTvMUAxcVals,0,2);
        data.(behavField).meanMUA_lags = mean(data.(behavField).cat_MUA_lags,2);
    end
end

%% summary figure(s)
hemoDataTypes = {'CBV','HbT'};
freq = 5;
lagTime1 = 5;
% lagTime2 = 15;
for h = 1:length(baselineTypes)
    baselineType = baselineTypes{1,h};
    for j = 1:length(hemoDataTypes)
        hemoDataType = hemoDataTypes{1,j};
        summaryFigure = figure;
        sgtitle({[hemoDataType '-neural cross-correlations - ' baselineType ' for resting data'],' '})
        
        %% Rest MUA
        ax1 = subplot(2,5,1);
        plot(data.Rest.(baselineType).meanMUA_lags,data.Rest.(baselineType).(['mean' hemoDataType 'vMUAxcVals']),'k')
        hold on
        plot(data.Rest.(baselineType).meanMUA_lags,data.Rest.(baselineType).(['mean' hemoDataType 'vMUAxcVals']) + data.Rest.(baselineType).(['std' hemoDataType 'vMUAxcVals']),'color',colors_IOS('battleship grey'))
        plot(data.Rest.(baselineType).meanMUA_lags,data.Rest.(baselineType).(['mean' hemoDataType 'vMUAxcVals']) - data.Rest.(baselineType).(['std' hemoDataType 'vMUAxcVals']),'color',colors_IOS('battleship grey'))
        title('Rest')
        xticks([-lagTime1*freq -lagTime1*freq/2 0 lagTime1*freq/2 lagTime1*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lagTime1*freq lagTime1*freq])
        xlabel('lags (sec)')
        ylabel('cross-correlation')
        axis square
             
        %% NREM MUA
        ax2 = subplot(2,5,2);
        plot(data.NREM.meanMUA_lags,data.NREM.(['mean' hemoDataType 'vMUAxcVals']),'k')
        hold on
        plot(data.NREM.meanMUA_lags,data.NREM.(['mean' hemoDataType 'vMUAxcVals']) + data.NREM.(['std' hemoDataType 'vMUAxcVals']),'color',colors_IOS('battleship grey'))
        plot(data.NREM.meanMUA_lags,data.NREM.(['mean' hemoDataType 'vMUAxcVals']) - data.NREM.(['std' hemoDataType 'vMUAxcVals']),'color',colors_IOS('battleship grey'))
        title('NREM')
        xticks([-lagTime1*freq -lagTime1*freq/2 0 lagTime1*freq/2 lagTime1*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lagTime1*freq lagTime1*freq])
        xlabel('lags (sec)')
        ylabel('cross-correlation')
        axis square
        
        %% REM MUA
        ax3 = subplot(2,5,3);
        plot(data.REM.meanMUA_lags,data.REM.(['mean' hemoDataType 'vMUAxcVals']),'k')
        hold on
        plot(data.REM.meanMUA_lags,data.REM.(['mean' hemoDataType 'vMUAxcVals']) + data.REM.(['std' hemoDataType 'vMUAxcVals']),'color',colors_IOS('battleship grey'))
        plot(data.REM.meanMUA_lags,data.REM.(['mean' hemoDataType 'vMUAxcVals']) - data.REM.(['std' hemoDataType 'vMUAxcVals']),'color',colors_IOS('battleship grey'))
        title('REM')
        xticks([-lagTime1*freq -lagTime1*freq/2 0 lagTime1*freq/2 lagTime1*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lagTime1*freq lagTime1*freq])
        xlabel('lags (sec)')
        ylabel('cross-correlation')
        axis square

        %% Unstim Data MUA
        ax4 = subplot(2,5,4);
        plot(data.Unstim.meanMUA_lags,data.Unstim.(['mean' hemoDataType 'vMUAxcVals']),'k')
        hold on
        plot(data.Unstim.meanMUA_lags,data.Unstim.(['mean' hemoDataType 'vMUAxcVals']) + data.Unstim.(['std' hemoDataType 'vMUAxcVals']),'color',colors_IOS('battleship grey'))
        plot(data.Unstim.meanMUA_lags,data.Unstim.(['mean' hemoDataType 'vMUAxcVals']) - data.Unstim.(['std' hemoDataType 'vMUAxcVals']),'color',colors_IOS('battleship grey'))
        title('Unstim data')
        xticks([-lagTime1*freq -lagTime1*freq/2 0 lagTime1*freq/2 lagTime1*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lagTime1*freq lagTime1*freq])
        xlabel('lags (sec)')
        ylabel('cross-correlation')
        axis square
        
        %% All Data MUA
        ax5 = subplot(2,5,5);
        plot(data.All.meanMUA_lags,data.All.(['mean' hemoDataType 'vMUAxcVals']),'k')
        hold on
        plot(data.All.meanMUA_lags,data.All.(['mean' hemoDataType 'vMUAxcVals']) + data.All.(['std' hemoDataType 'vMUAxcVals']),'color',colors_IOS('battleship grey'))
        plot(data.All.meanMUA_lags,data.All.(['mean' hemoDataType 'vMUAxcVals']) - data.All.(['std' hemoDataType 'vMUAxcVals']),'color',colors_IOS('battleship grey'))
        title('All data')
        xticks([-lagTime1*freq -lagTime1*freq/2 0 lagTime1*freq/2 lagTime1*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lagTime1*freq lagTime1*freq])
        xlabel('lags (sec)')
        ylabel('cross-correlation')
        axis square
        
        linkaxes([ax1 ax2 ax3 ax4 ax5],'y')

        %% Rest LFP
        ax6 = subplot(2,5,6);
        imagesc(data.Rest.(baselineType).meanLFP_lags,data.Rest.(baselineType).meanLFP_F,data.Rest.(baselineType).(['mean' hemoDataType 'vLFPxcVals']))
        xticks([-lagTime1*freq -lagTime1*freq/2 0 lagTime1*freq/2 lagTime1*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lagTime1*freq lagTime1*freq])
        xlabel('lags (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        colorbar
        axis xy
        axis square
        
        %% NREM LFP
        ax7 = subplot(2,5,7);
        imagesc(data.NREM.meanLFP_lags,data.NREM.meanLFP_F,data.NREM.(['mean' hemoDataType 'vLFPxcVals']))
        xticks([-lagTime1*freq -lagTime1*freq/2 0 lagTime1*freq/2 lagTime1*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lagTime1*freq lagTime1*freq])
        xlabel('lags (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        colorbar
        axis xy
        axis square
        
        %% REM LFP
        ax8 = subplot(2,5,8);
        imagesc(data.REM.meanLFP_lags,data.REM.meanLFP_F,data.REM.(['mean' hemoDataType 'vLFPxcVals']))
        xticks([-lagTime1*freq -lagTime1*freq/2 0 lagTime1*freq/2 lagTime1*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lagTime1*freq lagTime1*freq])
        xlabel('lags (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        colorbar
        axis xy
        axis square
        
        %% Unstim data LFP
        ax9 = subplot(2,5,9);
        imagesc(data.Unstim.meanLFP_lags,data.Unstim.meanLFP_F,data.Unstim.(['mean' hemoDataType 'vLFPxcVals']))
        xticks([-lagTime1*freq -lagTime1*freq/2 0 lagTime1*freq/2 lagTime1*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lagTime1*freq lagTime1*freq])
        xlabel('lags (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        colorbar
        axis xy
        axis square
        
        %% All data LFP
        ax10 = subplot(2,5,10);
        imagesc(data.All.meanLFP_lags,data.All.meanLFP_F,data.All.(['mean' hemoDataType 'vLFPxcVals']))
        xticks([-lagTime1*freq -lagTime1*freq/2 0 lagTime1*freq/2 lagTime1*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lagTime1*freq lagTime1*freq])
        xlabel('lags (sec)')
        ylabel('Freq (Hz)')
        ylim([1 100])
        colorbar
        axis xy
        axis square
        
        % save figure(s)
        dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Cross Correlation\';
        if ~exist(dirpath, 'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigure, [dirpath hemoDataType '_' baselineType '_AverageXCorr']);
    end
end
