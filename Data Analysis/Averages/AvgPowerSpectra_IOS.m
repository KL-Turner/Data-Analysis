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
powerSpec_dataTypes = {'CBV','CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','muaPower'};

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        for c = 1:length(powerSpec_dataTypes)
            powerSpec_dataType = powerSpec_dataTypes{1,c};
            if strcmp(behavField,'Rest') == true
                for d = 1:length(baselineTypes)
                    baselineType = baselineTypes{1,d};
                    data.(behavField).(powerSpec_dataType).(baselineType).LH.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).(baselineType).LH.S;
                    data.(behavField).(powerSpec_dataType).(baselineType).LH.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).(baselineType).LH.f;
                    data.(behavField).(powerSpec_dataType).(baselineType).RH.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).(baselineType).RH.S;
                    data.(behavField).(powerSpec_dataType).(baselineType).RH.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).(baselineType).RH.f;
                    if strmcp(powerSpec_dataType,'CBV') == false && strcmp(powerSpec_dataType,'CBV_HbT') == false
                        data.(behavField).(powerSpec_dataType).(baselineType).Hip.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).(baselineType).Hip.S;
                        data.(behavField).(powerSpec_dataType).(baselineType).Hip.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).(baselineType).Hip.f;
                    end
                end
            else
                data.(behavField).(powerSpec_dataType).LH.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).LH.S;
                data.(behavField).(powerSpec_dataType).LH.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).LH.f;
                data.(behavField).(powerSpec_dataType).RH.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).RH.S;
                data.(behavField).(powerSpec_dataType).RH.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).RH.f;
                if strmcp(powerSpec_dataType,'CBV') == false && strcmp(powerSpec_dataType,'CBV_HbT') == false
                    data.(behavField).(powerSpec_dataType).Hip.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).Hip.S;
                    data.(behavField).(powerSpec_dataType).Hip.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).Hip.f;
                end
            end
        end
    end
end

% concatenate the data from the left and right hemispheres
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(powerSpec_dataTypes)
        powerSpec_dataType = powerSpec_dataTypes{1,f};
        if strcmp(behavField,'Rest') == true
            for g = 1:length(baselineTypes)
                baselineType = baselineTypes{1,g;
                data.(behavField).(powerSpec_dataType).(baselineType).cat_S = cat(2,data.(behavField).(powerSpec_dataType).(baselineType).LH.S,data.(behavField).(powerSpec_dataType).(baselineType).RH.S);
                data.(behavField).(powerSpec_dataType).(baselineType).cat_f = cat(2,data.(behavField).(powerSpec_dataType).(baselineType).LH.f,data.(behavField).(powerSpec_dataType).(baselineType).RH.f);
            end
        else
            data.(behavField).(powerSpec_dataType).cat_S = cat(2,data.(behavField).(powerSpec_dataType).LH.S,data.(behavField).(powerSpec_dataType).RH.S);
            data.(behavField).(powerSpec_dataType).cat_f = cat(2,data.(behavField).(powerSpec_dataType).LH.f,data.(behavField).(powerSpec_dataType).RH.f);
        end
    end
end

% take the mean and standard deviation of each set of signals
for h = 1:length(behavFields)
    behavField = behavFields{1,h};
    for j = 1:length(powerSpec_dataTypes)
        powerSpec_dataType = powerSpec_dataTypes{1,j};
        if strcmp(behavField,'Rest')
            for k = 1:length(baselineTypes)
                baselineType = baselineTypes{1,k};
                data.(behavField).(powerSpec_dataType).(baselineType).meanCortS = mean(data.(behavField).(powerSpec_dataType).(baselineType).cat_S,2);
                data.(behavField).(powerSpec_dataType).(baselineType).stdCortS = std(data.(behavField).(powerSpec_dataType).(baselineType).cat_S,0,2);
                data.(behavField).(powerSpec_dataType).(baselineType).meanCortf = mean(data.(behavField).(powerSpec_dataType).(baselineType).cat_f,2);
                if strmcp(powerSpec_dataType,'CBV') == false && strcmp(powerSpec_dataType,'CBV_HbT') == false
                    data.(behavField).(powerSpec_dataType).(baselineType).meanHipS = mean(data.(behavField).(powerSpec_dataType).(baselineType).Hip.S,2);
                    data.(behavField).(powerSpec_dataType).(baselineType).stdHipS = std(data.(behavField).(powerSpec_dataType).(baselineType).Hip.S,0,2);
                    data.(behavField).(powerSpec_dataType).(baselineType).meanHipf = mean(data.(behavField).(powerSpec_dataType).(baselineType).Hip.f,2);
                end
            end
        else
            data.(behavField).(powerSpec_dataType).meanCortS = mean(data.(behavField).(powerSpec_dataType).cat_S,2);
            data.(behavField).(powerSpec_dataType).stdCortS = std(data.(behavField).(powerSpec_dataType).cat_S,0,2);
            data.(behavField).(powerSpec_dataType).meanCortf = mean(data.(behavField).(powerSpec_dataType).cat_f,2);
            if strmcp(powerSpec_dataType,'CBV') == false && strcmp(powerSpec_dataType,'CBV_HbT') == false
                data.(behavField).(powerSpec_dataType).meanHipS = mean(data.(behavField).(powerSpec_dataType).Hip.S,2);
                data.(behavField).(powerSpec_dataType).stdHipS = std(data.(behavField).(powerSpec_dataType).Hip.S,0,2);
                data.(behavField).(powerSpec_dataType).meanHipf = mean(data.(behavField).(powerSpec_dataType).Hip.f,2);
            end
        end
    end
end

%% summary figure(s) - cortex
 for m = 1:length(baselineTypes)
     baselineType = baselineTypes{1,m};
     
     summaryFigure = figure;
     sgtitle({['CBV and cortical power spectra - ' baselineType ' for resting data'],' '})
     %% CBV
     ax1 = subplot(2,4,1);
     plot(data.Rest.CBV.(baselineType).meanCortf,data.Rest.CBV.(baselineType).meanCortS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.CBV.meanCortf,data.AllData.CBV.meanCortS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.CBV.meanCortf,data.NREM.CBV.meanCortS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.CBV.meanCortf,data.REM.CBV.meanCortS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('CBV Reflectance')
     ylabel('Power')
     xlabel('Freq (Hz)')
     legend('Rest','AllData','NREM','REM')
     axis square
     xlim([0 1])
     
     %% CBV HbT
     ax2 = subplot(2,4,2);
     plot(data.Rest.CBV_HbT.(baselineType).meanCortf,data.Rest.CBV_HbT.(baselineType).meanCortS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.CBV_HbT.meanCortf,data.AllData.CBV_HbT.meanCortS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.CBV_HbT.meanCortf,data.NREM.CBV_HbT.meanCortS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.CBV_HbT.meanCortf,data.REM.CBV_HbT.meanCortS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('CBV HbT')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     %% Delta-band power
     ax3 = subplot(2,4,3);
     plot(data.Rest.deltaBandPower.(baselineType).meanCortf,data.Rest.deltaBandPower.(baselineType).meanCortS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.deltaBandPower.meanCortf,data.AllData.deltaBandPower.meanCortS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.deltaBandPower.meanCortf,data.NREM.deltaBandPower.meanCortS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.deltaBandPower.meanCortf,data.REM.deltaBandPower.meanCortS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('Delta-band power')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     %% Theta-band power
     ax4 = subplot(2,4,4);
     plot(data.Rest.thetaBandPower.(baselineType).meanCortf,data.Rest.thetaBandPower.(baselineType).meanCortS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.thetaBandPower.meanCortf,data.AllData.thetaBandPower.meanCortS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.thetaBandPower.meanCortf,data.NREM.thetaBandPower.meanCortS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.thetaBandPower.meanCortf,data.REM.thetaBandPower.meanCortS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('Theta-band power')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     %% Alpha-band power
     ax5 = subplot(2,4,5);
     plot(data.Rest.alphaBandPower.(baselineType).meanCortf,data.Rest.alphaBandPower.(baselineType).meanCortS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.alphaBandPower.meanCortf,data.AllData.alphaBandPower.meanCortS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.alphaBandPower.meanCortf,data.NREM.alphaBandPower.meanCortS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.alphaBandPower.meanCortf,data.REM.alphaBandPower.meanCortS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('Alpha-band power')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     %% Beta-band power
     ax6 = subplot(2,4,6);
     plot(data.Rest.betaBandPower.(baselineType).meanCortf,data.Rest.betaBandPower.(baselineType).meanCortS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.betaBandPower.meanCortf,data.AllData.betaBandPower.meanCortS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.betaBandPower.meanCortf,data.NREM.betaBandPower.meanCortS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.betaBandPower.meanCortf,data.REM.betaBandPower.meanCortS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('Beta-band power')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     %% Gamma-band power
     ax7 = subplot(2,4,7);
     plot(data.Rest.gammaBandPower.(baselineType).meanCortf,data.Rest.gammaBandPower.(baselineType).meanCortS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     plot(data.AllData.gammaBandPower.meanCortf,data.AllData.gammaBandPower.meanCortS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     hold on
     plot(data.NREM.gammaBandPower.meanCortf,data.NREM.gammaBandPower.meanCortS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.gammaBandPower.meanCortf,data.REM.gammaBandPower.meanCortS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('Gamma-band power')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     %% MUA power
     ax8 = subplot(2,4,8);
     plot(data.Rest.muaPower.(baselineType).meanCortf,data.Rest.muaPower.(baselineType).meanCortS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.muaPower.meanCortf,data.AllData.muaPower.meanCortS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.muaPower.meanCortf,data.NREM.muaPower.meanCortS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.muaPower.meanCortf,data.REM.muaPower.meanCortS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('MUA power')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],'y')
     
     % save figure(s)
    dirpath = 'S:\Users\klt8\Documents\Analysis Average Figures\Power Spectra\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end 
    savefig(summaryFigure, [dirpath baselineType '_CorticalPowerSpectra']);
 end
 
 %% summary figure(s) - hippocampus
 for n = 1:length(baselineTypes)
     baselineType = baselineTypes{1,n};
     
     summaryFigure = figure;
     sgtitle({['CBV and hippocampal power spectra - ' baselineType ' for resting data'],' '})
     %% CBV
     ax1 = subplot(2,4,1);
     plot(data.Rest.CBV.(baselineType).meanHipf,data.Rest.CBV.(baselineType).meanHipS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.CBV.meanHipf,data.AllData.CBV.meanHipS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.CBV.meanHipf,data.NREM.CBV.meanHipS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.CBV.meanHipf,data.REM.CBV.meanHipS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('CBV Reflectance')
     ylabel('Power')
     xlabel('Freq (Hz)')
     legend('Rest','AllData','NREM','REM')
     axis square
     xlim([0 1])
     
     %% CBV HbT
     ax2 = subplot(2,4,2);
     plot(data.Rest.CBV_HbT.(baselineType).meanHipf,data.Rest.CBV_HbT.(baselineType).meanHipS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.CBV_HbT.meanHipf,data.AllData.CBV_HbT.meanCortS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.CBV_HbT.meanHipf,data.NREM.CBV_HbT.meanHipS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.CBV_HbT.meanHipf,data.REM.CBV_HbT.meanHipS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('CBV HbT')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     %% Delta-band power
     ax3 = subplot(2,4,3);
     plot(data.Rest.deltaBandPower.(baselineType).meanHipf,data.Rest.deltaBandPower.(baselineType).meanHipS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.deltaBandPower.meanHipf,data.AllData.deltaBandPower.meanHipS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.deltaBandPower.meanHipf,data.NREM.deltaBandPower.meanHipS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.deltaBandPower.meanHipf,data.REM.deltaBandPower.meanHipS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('Delta-band power')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     %% Theta-band power
     ax4 = subplot(2,4,4);
     plot(data.Rest.thetaBandPower.(baselineType).meanHipf,data.Rest.thetaBandPower.(baselineType).meanHipS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.thetaBandPower.meanHipf,data.AllData.thetaBandPower.meanHipS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.thetaBandPower.meanHipf,data.NREM.thetaBandPower.meanHipS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.thetaBandPower.meanHipf,data.REM.thetaBandPower.meanHipS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('Theta-band power')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     %% Alpha-band power
     ax5 = subplot(2,4,5);
     plot(data.Rest.alphaBandPower.(baselineType).meanHipf,data.Rest.alphaBandPower.(baselineType).meanHipS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.alphaBandPower.meanHipf,data.AllData.alphaBandPower.meanHipS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.alphaBandPower.meanHipf,data.NREM.alphaBandPower.meanHipS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.alphaBandPower.meanHipf,data.REM.alphaBandPower.meanHipS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('Alpha-band power')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     %% Beta-band power
     ax6 = subplot(2,4,6);
     plot(data.Rest.betaBandPower.(baselineType).meanHipf,data.Rest.betaBandPower.(baselineType).meanHipS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.betaBandPower.meanHipf,data.AllData.betaBandPower.meanHipS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.betaBandPower.meanHipf,data.NREM.betaBandPower.meanHipS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.betaBandPower.meanHipf,data.REM.betaBandPower.meanHipS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('Beta-band power')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     %% Gamma-band power
     ax7 = subplot(2,4,7);
     plot(data.Rest.gammaBandPower.(baselineType).meanHipf,data.Rest.gammaBandPower.(baselineType).meanHipS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     plot(data.AllData.gammaBandPower.meanHipf,data.AllData.gammaBandPower.meanHipS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     hold on
     plot(data.NREM.gammaBandPower.meanHipf,data.NREM.gammaBandPower.meanHipS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.gammaBandPower.meanHipf,data.REM.gammaBandPower.meanHipS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('Gamma-band power')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     %% MUA power
     ax8 = subplot(2,4,8);
     plot(data.Rest.muaPower.(baselineType).meanHipf,data.Rest.muaPower.(baselineType).meanHipS,'color',colors_IOS(graphColors{1,1}),'LineWidth',2)
     hold on
     plot(data.AllData.muaPower.meanHipf,data.AllData.muaPower.meanHipS,'color',colors_IOS(graphColors{1,2}),'LineWidth',2)
     plot(data.NREM.muaPower.meanHipf,data.NREM.muaPower.meanHipS,'color',colors_IOS(graphColors{1,3}),'LineWidth',2)
     plot(data.REM.muaPower.meanHipf,data.REM.muaPower.meanHipS,'color',colors_IOS(graphColors{1,4}),'LineWidth',2)
     title('MUA power')
     ylabel('Power')
     xlabel('Freq (Hz)')
     axis square
     xlim([0 1])
     
     linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],'y')

     % save figure(s)
    dirpath = 'S:\Users\klt8\Documents\Analysis Average Figures\Power Spectra\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end 
    savefig(summaryFigure, [dirpath baselineType '_HippocampalPowerSpectra']);
 end