%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Calculate the average coherence of different behavioral states
%________________________________________________________________________________________________________________________
%
%   Inputs: none
%
%   Outputs: Generates summary figures saved to C: drive Documents folder
%
%   Last Revised: Oct 1st, 2019
%________________________________________________________________________________________________________________________
clear
clc

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110'};
driveLetters = {'E','E','E','F','F','F','D','D'};
behavFields = {'Rest','NREM','REM','Unstim','All'};
baselineTypes = {'manualSelection','setDuration','entireDuration'};
coherr_dataTypes = {'CBV','CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','muaPower'};
colorbrewer_setA_colorA = [0.520000 0.520000 0.510000];
colorbrewer_setA_colorB = [(31/256) (120/256) (180/256)];
colorbrewer_setA_colorC = [(255/256) (0/256) (115/256)];
colorbrewer_setA_colorD = [(51/256) (160/256) (44/256)];
colorbrewer_setA_colorE = [(255/256) (140/256) (0/256)];

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        for c = 1:length(coherr_dataTypes)
            coherr_dataType = coherr_dataTypes{1,c};
            if strcmp(behavField,'Rest') == true
                for d = 1:length(baselineTypes)
                    baselineType = baselineTypes{1,d};
                    data.(behavField).(coherr_dataType).(baselineType).C(:,a) = AnalysisResults.Coherence.(behavField).(coherr_dataType).(baselineType).C;
                    data.(behavField).(coherr_dataType).(baselineType).f(:,a) = AnalysisResults.Coherence.(behavField).(coherr_dataType).(baselineType).f;
                    data.(behavField).(coherr_dataType).(baselineType).confC(:,a) = AnalysisResults.Coherence.(behavField).(coherr_dataType).(baselineType).confC;
                end
            else
                data.(behavField).(coherr_dataType).C(:,a) = AnalysisResults.Coherence.(behavField).(coherr_dataType).C;
                data.(behavField).(coherr_dataType).f(:,a) = AnalysisResults.Coherence.(behavField).(coherr_dataType).f;
                data.(behavField).(coherr_dataType).confC(:,a) = AnalysisResults.Coherence.(behavField).(coherr_dataType).confC;
            end
        end
    end
end

% take the mean and standard deviation of each set of signals
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(coherr_dataTypes)
        coherr_dataType = coherr_dataTypes{1,f};
        if strcmp(behavField,'Rest')
            for g = 1:length(baselineTypes)
                baselineType = baselineTypes{1,g};
                data.(behavField).(coherr_dataType).(baselineType).meanC = mean(data.(behavField).(coherr_dataType).(baselineType).C,2);
                data.(behavField).(coherr_dataType).(baselineType).stdC = std(data.(behavField).(coherr_dataType).(baselineType).C,0,2);
                data.(behavField).(coherr_dataType).(baselineType).meanf = mean(data.(behavField).(coherr_dataType).(baselineType).f,2);
                data.(behavField).(coherr_dataType).(baselineType).maxConfC = max(data.(behavField).(coherr_dataType).(baselineType).confC);
                data.(behavField).(coherr_dataType).(baselineType).maxConfC_Y = ones(length(data.(behavField).(coherr_dataType).(baselineType).meanf),1)*data.(behavField).(coherr_dataType).(baselineType).maxConfC;
            end
        else
            data.(behavField).(coherr_dataType).meanC = mean(data.(behavField).(coherr_dataType).C,2);
            data.(behavField).(coherr_dataType).stdC = std(data.(behavField).(coherr_dataType).C,0,2);
            data.(behavField).(coherr_dataType).meanf = mean(data.(behavField).(coherr_dataType).f,2);
            data.(behavField).(coherr_dataType).maxConfC = max(data.(behavField).(coherr_dataType).confC);
            data.(behavField).(coherr_dataType).maxConfC_Y = ones(length(data.(behavField).(coherr_dataType).meanf),1)*data.(behavField).(coherr_dataType).maxConfC;
        end
    end
end

%% summary figure(s)
 for h = 1:length(baselineTypes)
     baselineType = baselineTypes{1,h};  
     summaryFigure = figure;
     sgtitle({['L/R Coherence - ' baselineType ' for resting data'],' '})
     %% CBV
     subplot(2,4,1);
     plot(data.Rest.CBV.(baselineType).meanf,data.Rest.CBV.(baselineType).meanC,'color',colorbrewer_setA_colorA,'LineWidth',2)
     hold on
     plot(data.NREM.CBV.meanf,data.NREM.CBV.meanC,'color',colorbrewer_setA_colorB,'LineWidth',2)
     plot(data.REM.CBV.meanf,data.REM.CBV.meanC,'color',colorbrewer_setA_colorC,'LineWidth',2)
     plot(data.Unstim.CBV.meanf,data.Unstim.CBV.meanC,'color',colorbrewer_setA_colorD,'LineWidth',2)
     plot(data.All.CBV.meanf,data.All.CBV.meanC,'color',colorbrewer_setA_colorE,'LineWidth',2)
     plot(data.Rest.CBV.(baselineType).meanf,data.Rest.CBV.(baselineType).maxConfC_Y,'color',colorbrewer_setA_colorA,'LineWidth',1)
     plot(data.NREM.CBV.meanf,data.NREM.CBV.maxConfC_Y,'color',colorbrewer_setA_colorB,'LineWidth',1)
     plot(data.REM.CBV.meanf,data.REM.CBV.maxConfC_Y,'color',colorbrewer_setA_colorC,'LineWidth',1)
     plot(data.Unstim.CBV.meanf,data.Unstim.CBV.maxConfC_Y,'color',colorbrewer_setA_colorD','LineWidth',1)
     plot(data.All.CBV.meanf,data.All.CBV.maxConfC_Y,'color',colorbrewer_setA_colorE','LineWidth',1)
     title('CBV Reflectance')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     legend('Rest','NREM','REM','Unstim','All','Rest 95% conf','NREM 95% conf','REM 95% conf','Unstim 95% conf','All data 95% conf')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% CBV HbT
     subplot(2,4,2);
     plot(data.Rest.CBV_HbT.(baselineType).meanf,data.Rest.CBV_HbT.(baselineType).meanC,'color',colorbrewer_setA_colorA,'LineWidth',2)
     hold on
     plot(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.meanC,'color',colorbrewer_setA_colorB,'LineWidth',2)
     plot(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.meanC,'color',colorbrewer_setA_colorC,'LineWidth',2)
     plot(data.Unstim.CBV_HbT.meanf,data.Unstim.CBV_HbT.meanC,'color',colorbrewer_setA_colorD,'LineWidth',2)
     plot(data.All.CBV_HbT.meanf,data.All.CBV_HbT.meanC,'color',colorbrewer_setA_colorE,'LineWidth',2)
     plot(data.Rest.CBV_HbT.(baselineType).meanf,data.Rest.CBV_HbT.(baselineType).maxConfC_Y,'color',colorbrewer_setA_colorA,'LineWidth',1)
     plot(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.maxConfC_Y,'color',colorbrewer_setA_colorB,'LineWidth',1)
     plot(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.maxConfC_Y,'color',colorbrewer_setA_colorC,'LineWidth',1)
     plot(data.Unstim.CBV_HbT.meanf,data.Unstim.CBV_HbT.maxConfC_Y,'color',colorbrewer_setA_colorD','LineWidth',1)
     plot(data.All.CBV_HbT.meanf,data.All.CBV_HbT.maxConfC_Y,'color',colorbrewer_setA_colorE','LineWidth',1)
     title('CBV HbT')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% Delta-band power
     subplot(2,4,3);
     plot(data.Rest.deltaBandPower.(baselineType).meanf,data.Rest.deltaBandPower.(baselineType).meanC,'color',colorbrewer_setA_colorA,'LineWidth',2)
     hold on
     plot(data.NREM.deltaBandPower.meanf,data.NREM.deltaBandPower.meanC,'color',colorbrewer_setA_colorB,'LineWidth',2)
     plot(data.REM.deltaBandPower.meanf,data.REM.deltaBandPower.meanC,'color',colorbrewer_setA_colorC,'LineWidth',2)
     plot(data.Unstim.deltaBandPower.meanf,data.Unstim.deltaBandPower.meanC,'color',colorbrewer_setA_colorD,'LineWidth',2)
     plot(data.All.deltaBandPower.meanf,data.All.deltaBandPower.meanC,'color',colorbrewer_setA_colorE,'LineWidth',2)
     plot(data.Rest.deltaBandPower.(baselineType).meanf,data.Rest.deltaBandPower.(baselineType).maxConfC_Y,'color',colorbrewer_setA_colorA,'LineWidth',1)
     plot(data.NREM.deltaBandPower.meanf,data.NREM.deltaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorB,'LineWidth',1)
     plot(data.REM.deltaBandPower.meanf,data.REM.deltaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorC,'LineWidth',1)
     plot(data.Unstim.deltaBandPower.meanf,data.Unstim.deltaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorD','LineWidth',1)
     plot(data.All.deltaBandPower.meanf,data.All.deltaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorE','LineWidth',1)
     title('Delta-band power')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% Theta-band power
     subplot(2,4,4);
     plot(data.Rest.thetaBandPower.(baselineType).meanf,data.Rest.thetaBandPower.(baselineType).meanC,'color',colorbrewer_setA_colorA,'LineWidth',2)
     hold on
     plot(data.NREM.thetaBandPower.meanf,data.NREM.thetaBandPower.meanC,'color',colorbrewer_setA_colorB,'LineWidth',2)
     plot(data.REM.thetaBandPower.meanf,data.REM.thetaBandPower.meanC,'color',colorbrewer_setA_colorC,'LineWidth',2)
     plot(data.Unstim.thetaBandPower.meanf,data.Unstim.thetaBandPower.meanC,'color',colorbrewer_setA_colorD,'LineWidth',2)
     plot(data.All.thetaBandPower.meanf,data.All.thetaBandPower.meanC,'color',colorbrewer_setA_colorE,'LineWidth',2)
     plot(data.Rest.thetaBandPower.(baselineType).meanf,data.Rest.thetaBandPower.(baselineType).maxConfC_Y,'color',colorbrewer_setA_colorA,'LineWidth',1)
     plot(data.NREM.thetaBandPower.meanf,data.NREM.thetaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorB,'LineWidth',1)
     plot(data.REM.thetaBandPower.meanf,data.REM.thetaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorC,'LineWidth',1)
     plot(data.Unstim.thetaBandPower.meanf,data.Unstim.thetaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorD','LineWidth',1)
     plot(data.All.thetaBandPower.meanf,data.All.thetaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorE','LineWidth',1)
     title('Theta-band power')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% Alpha-band power
     subplot(2,4,5);
     plot(data.alphaBandPower.CBV.(baselineType).meanf,data.Rest.alphaBandPower.(baselineType).meanC,'color',colorbrewer_setA_colorA,'LineWidth',2)
     hold on
     plot(data.NREM.alphaBandPower.meanf,data.NREM.alphaBandPower.meanC,'color',colorbrewer_setA_colorB,'LineWidth',2)
     plot(data.REM.alphaBandPower.meanf,data.REM.alphaBandPower.meanC,'color',colorbrewer_setA_colorC,'LineWidth',2)
     plot(data.Unstim.alphaBandPower.meanf,data.Unstim.alphaBandPower.meanC,'color',colorbrewer_setA_colorD,'LineWidth',2)
     plot(data.All.alphaBandPower.meanf,data.All.alphaBandPower.meanC,'color',colorbrewer_setA_colorE,'LineWidth',2)
     plot(data.Rest.alphaBandPower.(baselineType).meanf,data.Rest.alphaBandPower.(baselineType).maxConfC_Y,'color',colorbrewer_setA_colorA,'LineWidth',1)
     plot(data.NREM.alphaBandPower.meanf,data.NREM.alphaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorB,'LineWidth',1)
     plot(data.REM.alphaBandPower.meanf,data.REM.alphaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorC,'LineWidth',1)
     plot(data.Unstim.alphaBandPower.meanf,data.Unstim.alphaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorD','LineWidth',1)
     plot(data.All.alphaBandPower.meanf,data.All.alphaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorE','LineWidth',1)
     title('Alpha-band power')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% Beta-band power
     subplot(2,4,6);
     plot(data.Rest.betaBandPower.(baselineType).meanf,data.Rest.betaBandPower.(baselineType).meanC,'color',colorbrewer_setA_colorA,'LineWidth',2)
     hold on
     plot(data.NREM.betaBandPower.meanf,data.NREM.betaBandPower.meanC,'color',colorbrewer_setA_colorB,'LineWidth',2)
     plot(data.REM.betaBandPower.meanf,data.REM.betaBandPower.meanC,'color',colorbrewer_setA_colorC,'LineWidth',2)
     plot(data.Unstim.betaBandPower.meanf,data.Unstim.betaBandPower.meanC,'color',colorbrewer_setA_colorD,'LineWidth',2)
     plot(data.All.betaBandPower.meanf,data.All.betaBandPower.meanC,'color',colorbrewer_setA_colorE,'LineWidth',2)
     plot(data.Rest.betaBandPower.(baselineType).meanf,data.Rest.betaBandPower.(baselineType).maxConfC_Y,'color',colorbrewer_setA_colorA,'LineWidth',1)
     plot(data.NREM.betaBandPower.meanf,data.NREM.betaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorB,'LineWidth',1)
     plot(data.REM.betaBandPower.meanf,data.REM.betaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorC,'LineWidth',1)
     plot(data.Unstim.betaBandPower.meanf,data.Unstim.betaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorD','LineWidth',1)
     plot(data.All.betaBandPower.meanf,data.All.betaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorE','LineWidth',1)
     title('Beta-band power')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% Gamma-band power
     subplot(2,4,7);
     plot(data.Rest.gammaBandPower.(baselineType).meanf,data.Rest.gammaBandPower.(baselineType).meanC,'color',colorbrewer_setA_colorA,'LineWidth',2)
     hold on
     plot(data.NREM.gammaBandPower.meanf,data.NREM.gammaBandPower.meanC,'color',colorbrewer_setA_colorB,'LineWidth',2)
     plot(data.REM.gammaBandPower.meanf,data.REM.gammaBandPower.meanC,'color',colorbrewer_setA_colorC,'LineWidth',2)
     plot(data.Unstim.gammaBandPower.meanf,data.Unstim.gammaBandPower.meanC,'color',colorbrewer_setA_colorD,'LineWidth',2)
     plot(data.All.gammaBandPower.meanf,data.All.gammaBandPower.meanC,'color',colorbrewer_setA_colorE,'LineWidth',2)
     plot(data.Rest.gammaBandPower.(baselineType).meanf,data.Rest.gammaBandPower.(baselineType).maxConfC_Y,'color',colorbrewer_setA_colorA,'LineWidth',1)
     plot(data.NREM.gammaBandPower.meanf,data.NREM.gammaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorB,'LineWidth',1)
     plot(data.REM.gammaBandPower.meanf,data.REM.gammaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorC,'LineWidth',1)
     plot(data.Unstim.gammaBandPower.meanf,data.Unstim.gammaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorD','LineWidth',1)
     plot(data.All.gammaBandPower.meanf,data.All.gammaBandPower.maxConfC_Y,'color',colorbrewer_setA_colorE','LineWidth',1)
     title('Gamma-band power')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% MUA power
     subplot(2,4,8);
     plot(data.Rest.muaPower.(baselineType).meanf,data.Rest.muaPower.(baselineType).meanC,'color',colorbrewer_setA_colorA,'LineWidth',2)
     hold on
     plot(data.NREM.muaPower.meanf,data.NREM.muaPower.meanC,'color',colorbrewer_setA_colorB,'LineWidth',2)
     plot(data.REM.muaPower.meanf,data.REM.muaPower.meanC,'color',colorbrewer_setA_colorC,'LineWidth',2)
     plot(data.Unstim.muaPower.meanf,data.Unstim.muaPower.meanC,'color',colorbrewer_setA_colorD,'LineWidth',2)
     plot(data.All.muaPower.meanf,data.All.muaPower.meanC,'color',colorbrewer_setA_colorE,'LineWidth',2)
     plot(data.Rest.muaPower.(baselineType).meanf,data.Rest.muaPower.(baselineType).maxConfC_Y,'color',colorbrewer_setA_colorA,'LineWidth',1)
     plot(data.NREM.muaPower.meanf,data.NREM.muaPower.maxConfC_Y,'color',colorbrewer_setA_colorB,'LineWidth',1)
     plot(data.REM.muaPower.meanf,data.REM.muaPower.maxConfC_Y,'color',colorbrewer_setA_colorC,'LineWidth',1)
     plot(data.Unstim.muaPower.meanf,data.Unstim.muaPower.maxConfC_Y,'color',colorbrewer_setA_colorD','LineWidth',1)
     plot(data.All.muaPower.meanf,data.All.muaPower.maxConfC_Y,'color',colorbrewer_setA_colorE','LineWidth',1)
     title('MUA power')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     % save figure(s)
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Coherence\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end 
    savefig(summaryFigure, [dirpath baselineType '_AverageCoherence']);
 end
