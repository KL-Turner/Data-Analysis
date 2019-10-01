%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised: Oct 1st, 2019
%________________________________________________________________________________________________________________________
clear
clc

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110'};
driveLetters = {'E','E','E','F','F','F','D','D'};
behavFields = {'Rest','AllData','NREM','REM'};
baselineTypes = {'manualSelection','setDuration','entireDuration'};
coherr_dataTypes = {'CBV','CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','muaPower'};
graphColors = {'sapphire','electric purple','coral red','vegas gold','jungle green','deep carrot orange','battleship grey'}; 

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
     plot(data.Rest.CBV.(baselineType).meanf,data.Rest.CBV.(baselineType).meanC,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.CBV.meanf,data.AllData.CBV.meanC,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.CBV.meanf,data.NREM.CBV.meanC,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.CBV.meanf,data.REM.CBV.meanC,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     plot(data.Rest.CBV.(baselineType).meanf,data.Rest.CBV.(baselineType).maxConfC_Y,'--','color',colors_IOS(graphColors{1,1}),'LineWidth',1)
     plot(data.AllData.CBV.meanf,data.AllData.CBV.maxConfC_Y,'--','color',colors_IOS(graphColors{1,2}),'LineWidth',1)
     plot(data.NREM.CBV.meanf,data.NREM.CBV.maxConfC_Y,'--','color',colors_IOS(graphColors{1,3}),'LineWidth',1)
     plot(data.REM.CBV.meanf,data.REM.CBV.maxConfC_Y,'--','color',colors_IOS(graphColors{1,4}),'LineWidth',1)
     title('CBV Reflectance')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     legend('Rest','AllData','NREM','REM','Rest 95% conf','All Data 95% conf','NREM 95% conf','REM 95% conf')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% CBV HbT
     subplot(2,4,2);
     plot(data.Rest.CBV_HbT.(baselineType).meanf,data.Rest.CBV_HbT.(baselineType).meanC,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.CBV_HbT.meanf,data.AllData.CBV_HbT.meanC,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.meanC,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.meanC,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     plot(data.Rest.CBV_HbT.(baselineType).meanf,data.Rest.CBV_HbT.(baselineType).maxConfC_Y,'--','color',colors_IOS(graphColors{1,1}),'LineWidth',1)
     plot(data.AllData.CBV_HbT.meanf,data.AllData.CBV_HbT.maxConfC_Y,'--','color',colors_IOS(graphColors{1,2}),'LineWidth',1)
     plot(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.maxConfC_Y,'--','color',colors_IOS(graphColors{1,3}),'LineWidth',1)
     plot(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.maxConfC_Y,'--','color',colors_IOS(graphColors{1,4}),'LineWidth',1)
     title('CBV HbT')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% Delta-band power
     subplot(2,4,3);
     plot(data.Rest.deltaBandPower.(baselineType).meanf,data.Rest.deltaBandPower.(baselineType).meanC,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.deltaBandPower.meanf,data.AllData.deltaBandPower.meanC,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.deltaBandPower.meanf,data.NREM.deltaBandPower.meanC,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.deltaBandPower.meanf,data.REM.deltaBandPower.meanC,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     plot(data.Rest.deltaBandPower.(baselineType).meanf,data.Rest.deltaBandPower.(baselineType).maxConfC_Y,'--','color',colors_IOS(graphColors{1,1}),'LineWidth',1)
     plot(data.AllData.deltaBandPower.meanf,data.AllData.deltaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,2}),'LineWidth',1)
     plot(data.NREM.deltaBandPower.meanf,data.NREM.deltaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,3}),'LineWidth',1)
     plot(data.REM.deltaBandPower.meanf,data.REM.deltaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,4}),'LineWidth',1)
     title('Delta-band power')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% Theta-band power
     subplot(2,4,4);
     plot(data.Rest.thetaBandPower.(baselineType).meanf,data.Rest.thetaBandPower.(baselineType).meanC,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.thetaBandPower.meanf,data.AllData.thetaBandPower.meanC,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.thetaBandPower.meanf,data.NREM.thetaBandPower.meanC,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.thetaBandPower.meanf,data.REM.thetaBandPower.meanC,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     plot(data.Rest.thetaBandPower.(baselineType).meanf,data.Rest.thetaBandPower.(baselineType).maxConfC_Y,'--','color',colors_IOS(graphColors{1,1}),'LineWidth',1)
     plot(data.AllData.thetaBandPower.meanf,data.AllData.thetaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,2}),'LineWidth',1)
     plot(data.NREM.thetaBandPower.meanf,data.NREM.thetaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,3}),'LineWidth',1)
     plot(data.REM.thetaBandPower.meanf,data.REM.thetaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,4}),'LineWidth',1)
     title('Theta-band power')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% Alpha-band power
     subplot(2,4,5);
     plot(data.Rest.alphaBandPower.(baselineType).meanf,data.Rest.alphaBandPower.(baselineType).meanC,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.alphaBandPower.meanf,data.AllData.alphaBandPower.meanC,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.alphaBandPower.meanf,data.NREM.alphaBandPower.meanC,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.alphaBandPower.meanf,data.REM.alphaBandPower.meanC,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     plot(data.Rest.alphaBandPower.(baselineType).meanf,data.Rest.alphaBandPower.(baselineType).maxConfC_Y,'--','color',colors_IOS(graphColors{1,1}),'LineWidth',1)
     plot(data.AllData.alphaBandPower.meanf,data.AllData.alphaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,2}),'LineWidth',1)
     plot(data.NREM.alphaBandPower.meanf,data.NREM.alphaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,3}),'LineWidth',1)
     plot(data.REM.alphaBandPower.meanf,data.REM.alphaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,4}),'LineWidth',1)
     title('Alpha-band power')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% Beta-band power
     subplot(2,4,6);
     plot(data.Rest.betaBandPower.(baselineType).meanf,data.Rest.betaBandPower.(baselineType).meanC,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.betaBandPower.meanf,data.AllData.betaBandPower.meanC,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.betaBandPower.meanf,data.NREM.betaBandPower.meanC,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.betaBandPower.meanf,data.REM.betaBandPower.meanC,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     plot(data.Rest.betaBandPower.(baselineType).meanf,data.Rest.betaBandPower.(baselineType).maxConfC_Y,'--','color',colors_IOS(graphColors{1,1}),'LineWidth',1)
     plot(data.AllData.betaBandPower.meanf,data.AllData.betaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,2}),'LineWidth',1)
     plot(data.NREM.betaBandPower.meanf,data.NREM.betaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,3}),'LineWidth',1)
     plot(data.REM.betaBandPower.meanf,data.REM.betaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,4}),'LineWidth',1)
     title('Beta-band power')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% Gamma-band power
     subplot(2,4,7);
     plot(data.Rest.gammaBandPower.(baselineType).meanf,data.Rest.gammaBandPower.(baselineType).meanC,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     plot(data.AllData.gammaBandPower.meanf,data.AllData.gammaBandPower.meanC,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     hold on
     plot(data.NREM.gammaBandPower.meanf,data.NREM.gammaBandPower.meanC,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.gammaBandPower.meanf,data.REM.gammaBandPower.meanC,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     plot(data.Rest.gammaBandPower.(baselineType).meanf,data.Rest.gammaBandPower.(baselineType).maxConfC_Y,'--','color',colors_IOS(graphColors{1,1}),'LineWidth',1)
     plot(data.AllData.gammaBandPower.meanf,data.AllData.gammaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,2}),'LineWidth',1)
     plot(data.NREM.gammaBandPower.meanf,data.NREM.gammaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,3}),'LineWidth',1)
     plot(data.REM.gammaBandPower.meanf,data.REM.gammaBandPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,4}),'LineWidth',1)
     title('Gamma-band power')
     ylabel('Coherence')
     xlabel('Freq (Hz)')
     axis square
     ylim([0 1])
     xlim([0 1])
     
     %% MUA power
     subplot(2,4,8);
     plot(data.Rest.muaPower.(baselineType).meanf,data.Rest.muaPower.(baselineType).meanC,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.muaPower.meanf,data.AllData.muaPower.meanC,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.muaPower.meanf,data.NREM.muaPower.meanC,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.muaPower.meanf,data.REM.muaPower.meanC,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     plot(data.Rest.muaPower.(baselineType).meanf,data.Rest.muaPower.(baselineType).maxConfC_Y,'--','color',colors_IOS(graphColors{1,1}),'LineWidth',1)
     plot(data.AllData.muaPower.meanf,data.AllData.muaPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,2}),'LineWidth',1)
     plot(data.NREM.muaPower.meanf,data.NREM.muaPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,3}),'LineWidth',1)
     plot(data.REM.muaPower.meanf,data.REM.muaPower.maxConfC_Y,'--','color',colors_IOS(graphColors{1,4}),'LineWidth',1)
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
