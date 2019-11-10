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

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111'};
driveLetters = {'E','E','E','F','F','F','D','D','D'};
behavFields = {'Rest','NREM','REM','Unstim','All'};
baselineTypes = {'manualSelection','setDuration','entireDuration'};
powerSpec_dataTypes = {'CBV','CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','muaPower'};
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
        for c = 1:length(powerSpec_dataTypes)
            powerSpec_dataType = powerSpec_dataTypes{1,c};
            if strcmp(behavField,'Rest') == true
                for d = 1:length(baselineTypes)
                    baselineType = baselineTypes{1,d};
                    data.(behavField).(powerSpec_dataType).(baselineType).adjLH.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).(baselineType).adjLH.S;
                    data.(behavField).(powerSpec_dataType).(baselineType).adjLH.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).(baselineType).adjLH.f;
                    data.(behavField).(powerSpec_dataType).(baselineType).adjRH.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).(baselineType).adjRH.S;
                    data.(behavField).(powerSpec_dataType).(baselineType).adjRH.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).(baselineType).adjRH.f;
                    if strcmp(powerSpec_dataType,'CBV') == false && strcmp(powerSpec_dataType,'CBV_HbT') == false
                        data.(behavField).(powerSpec_dataType).(baselineType).Hip.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).(baselineType).Hip.S;
                        data.(behavField).(powerSpec_dataType).(baselineType).Hip.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).(baselineType).Hip.f;
                    end
                end
            else
                data.(behavField).(powerSpec_dataType).adjLH.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).adjLH.S;
                data.(behavField).(powerSpec_dataType).adjLH.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).adjLH.f;
                data.(behavField).(powerSpec_dataType).adjRH.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).adjRH.S;
                data.(behavField).(powerSpec_dataType).adjRH.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).adjRH.f;
                if strcmp(powerSpec_dataType,'CBV') == false && strcmp(powerSpec_dataType,'CBV_HbT') == false
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
                baselineType = baselineTypes{1,g};
                data.(behavField).(powerSpec_dataType).(baselineType).cat_S = cat(2,data.(behavField).(powerSpec_dataType).(baselineType).adjLH.S,data.(behavField).(powerSpec_dataType).(baselineType).adjRH.S);
                data.(behavField).(powerSpec_dataType).(baselineType).cat_f = cat(2,data.(behavField).(powerSpec_dataType).(baselineType).adjLH.f,data.(behavField).(powerSpec_dataType).(baselineType).adjRH.f);
            end
        else
            data.(behavField).(powerSpec_dataType).cat_S = cat(2,data.(behavField).(powerSpec_dataType).adjLH.S,data.(behavField).(powerSpec_dataType).adjRH.S);
            data.(behavField).(powerSpec_dataType).cat_f = cat(2,data.(behavField).(powerSpec_dataType).adjLH.f,data.(behavField).(powerSpec_dataType).adjRH.f);
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
                if strcmp(powerSpec_dataType,'CBV') == false && strcmp(powerSpec_dataType,'CBV_HbT') == false
                    data.(behavField).(powerSpec_dataType).(baselineType).meanHipS = mean(data.(behavField).(powerSpec_dataType).(baselineType).Hip.S,2);
                    data.(behavField).(powerSpec_dataType).(baselineType).stdHipS = std(data.(behavField).(powerSpec_dataType).(baselineType).Hip.S,0,2);
                    data.(behavField).(powerSpec_dataType).(baselineType).meanHipf = mean(data.(behavField).(powerSpec_dataType).(baselineType).Hip.f,2);
                end
            end
        else
            data.(behavField).(powerSpec_dataType).meanCortS = mean(data.(behavField).(powerSpec_dataType).cat_S,2);
            data.(behavField).(powerSpec_dataType).stdCortS = std(data.(behavField).(powerSpec_dataType).cat_S,0,2);
            data.(behavField).(powerSpec_dataType).meanCortf = mean(data.(behavField).(powerSpec_dataType).cat_f,2);
            if strcmp(powerSpec_dataType,'CBV') == false && strcmp(powerSpec_dataType,'CBV_HbT') == false
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
    L1 = loglog(data.Rest.CBV.(baselineType).meanCortf,data.Rest.CBV.(baselineType).meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',2);
    hold on
    loglog(data.Rest.CBV.(baselineType).meanCortf,data.Rest.CBV.(baselineType).meanCortS + data.Rest.CBV.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.CBV.(baselineType).meanCortf,data.Rest.CBV.(baselineType).meanCortS - data.Rest.CBV.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    L2 = loglog(data.NREM.CBV.meanCortf,data.NREM.CBV.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.CBV.meanCortf,data.NREM.CBV.meanCortS + data.NREM.CBV.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.CBV.meanCortf,data.NREM.CBV.meanCortS - data.NREM.CBV.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    L3 = loglog(data.REM.CBV.meanCortf,data.REM.CBV.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.CBV.meanCortf,data.REM.CBV.meanCortS + data.REM.CBV.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.CBV.meanCortf,data.REM.CBV.meanCortS - data.REM.CBV.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    L4 = loglog(data.Unstim.CBV.meanCortf,data.Unstim.CBV.meanCortS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.CBV.meanCortf,data.Unstim.CBV.meanCortS + data.Unstim.CBV.stdCortS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.CBV.meanCortf,data.Unstim.CBV.meanCortS - data.Unstim.CBV.stdCortS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    L5 = loglog(data.All.CBV.meanCortf,data.All.CBV.meanCortS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.CBV.meanCortf,data.All.CBV.meanCortS + data.All.CBV.stdCortS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.CBV.meanCortf,data.All.CBV.meanCortS - data.All.CBV.stdCortS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('CBV Reflectance')
    ylabel('Power')
    xlabel('Freq (Hz)')
    legend([L1,L2,L3,L4,L5],'Rest','NREM','REM','Unstim','All')
    axis square
    xlim([0.05 1])
    
    %% CBV HbT
    ax2 = subplot(2,4,2);
    loglog(data.Rest.CBV_HbT.(baselineType).meanCortf,data.Rest.CBV_HbT.(baselineType).meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.CBV_HbT.(baselineType).meanCortf,data.Rest.CBV_HbT.(baselineType).meanCortS + data.Rest.CBV_HbT.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.CBV_HbT.(baselineType).meanCortf,data.Rest.CBV_HbT.(baselineType).meanCortS - data.Rest.CBV_HbT.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.CBV_HbT.meanCortf,data.NREM.CBV_HbT.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.CBV_HbT.meanCortf,data.NREM.CBV_HbT.meanCortS + data.NREM.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.CBV_HbT.meanCortf,data.NREM.CBV_HbT.meanCortS - data.NREM.CBV_HbT.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.CBV_HbT.meanCortf,data.REM.CBV_HbT.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.CBV_HbT.meanCortf,data.REM.CBV_HbT.meanCortS + data.REM.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.CBV_HbT.meanCortf,data.REM.CBV_HbT.meanCortS - data.REM.CBV_HbT.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.CBV_HbT.meanCortf,data.Unstim.CBV_HbT.meanCortS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.CBV_HbT.meanCortf,data.Unstim.CBV_HbT.meanCortS + data.Unstim.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.CBV_HbT.meanCortf,data.Unstim.CBV_HbT.meanCortS - data.Unstim.CBV_HbT.stdCortS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.CBV_HbT.meanCortf,data.All.CBV_HbT.meanCortS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.CBV_HbT.meanCortf,data.All.CBV_HbT.meanCortS + data.All.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.CBV_HbT.meanCortf,data.All.CBV_HbT.meanCortS - data.All.CBV_HbT.stdCortS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('CBV HbT')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    %% Delta-band power
    ax3 = subplot(2,4,3);
    loglog(data.Rest.deltaBandPower.(baselineType).meanCortf,data.Rest.deltaBandPower.(baselineType).meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.deltaBandPower.(baselineType).meanCortf,data.Rest.deltaBandPower.(baselineType).meanCortS + data.Rest.deltaBandPower.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.deltaBandPower.(baselineType).meanCortf,data.Rest.deltaBandPower.(baselineType).meanCortS - data.Rest.deltaBandPower.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.deltaBandPower.meanCortf,data.NREM.deltaBandPower.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.deltaBandPower.meanCortf,data.NREM.deltaBandPower.meanCortS + data.NREM.deltaBandPower.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.deltaBandPower.meanCortf,data.NREM.deltaBandPower.meanCortS - data.NREM.deltaBandPower.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.deltaBandPower.meanCortf,data.REM.deltaBandPower.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.deltaBandPower.meanCortf,data.REM.deltaBandPower.meanCortS + data.REM.deltaBandPower.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.deltaBandPower.meanCortf,data.REM.deltaBandPower.meanCortS - data.REM.deltaBandPower.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.deltaBandPower.meanCortf,data.Unstim.deltaBandPower.meanCortS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.deltaBandPower.meanCortf,data.Unstim.deltaBandPower.meanCortS + data.Unstim.deltaBandPower.stdCortS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.deltaBandPower.meanCortf,data.Unstim.deltaBandPower.meanCortS - data.Unstim.deltaBandPower.stdCortS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.deltaBandPower.meanCortf,data.All.deltaBandPower.meanCortS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.deltaBandPower.meanCortf,data.All.deltaBandPower.meanCortS + data.All.deltaBandPower.stdCortS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.deltaBandPower.meanCortf,data.All.deltaBandPower.meanCortS - data.All.deltaBandPower.stdCortS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('Delta-band power')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    %% Theta-band power
    ax4 = subplot(2,4,4);
    loglog(data.Rest.thetaBandPower.(baselineType).meanCortf,data.Rest.thetaBandPower.(baselineType).meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.thetaBandPower.(baselineType).meanCortf,data.Rest.thetaBandPower.(baselineType).meanCortS + data.Rest.thetaBandPower.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.thetaBandPower.(baselineType).meanCortf,data.Rest.thetaBandPower.(baselineType).meanCortS - data.Rest.thetaBandPower.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.thetaBandPower.meanCortf,data.NREM.thetaBandPower.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.thetaBandPower.meanCortf,data.NREM.thetaBandPower.meanCortS + data.NREM.thetaBandPower.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.thetaBandPower.meanCortf,data.NREM.thetaBandPower.meanCortS - data.NREM.thetaBandPower.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.thetaBandPower.meanCortf,data.REM.thetaBandPower.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.thetaBandPower.meanCortf,data.REM.thetaBandPower.meanCortS + data.REM.thetaBandPower.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.thetaBandPower.meanCortf,data.REM.thetaBandPower.meanCortS - data.REM.thetaBandPower.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.thetaBandPower.meanCortf,data.Unstim.thetaBandPower.meanCortS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.thetaBandPower.meanCortf,data.Unstim.thetaBandPower.meanCortS + data.Unstim.thetaBandPower.stdCortS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.thetaBandPower.meanCortf,data.Unstim.thetaBandPower.meanCortS - data.Unstim.thetaBandPower.stdCortS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.thetaBandPower.meanCortf,data.All.thetaBandPower.meanCortS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.thetaBandPower.meanCortf,data.All.thetaBandPower.meanCortS + data.All.thetaBandPower.stdCortS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.thetaBandPower.meanCortf,data.All.thetaBandPower.meanCortS - data.All.thetaBandPower.stdCortS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('Theta-band power')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    %% Alpha-band power
    ax5 = subplot(2,4,5);
    loglog(data.Rest.alphaBandPower.(baselineType).meanCortf,data.Rest.alphaBandPower.(baselineType).meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.alphaBandPower.(baselineType).meanCortf,data.Rest.alphaBandPower.(baselineType).meanCortS + data.Rest.alphaBandPower.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.alphaBandPower.(baselineType).meanCortf,data.Rest.alphaBandPower.(baselineType).meanCortS - data.Rest.alphaBandPower.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.alphaBandPower.meanCortf,data.NREM.alphaBandPower.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.alphaBandPower.meanCortf,data.NREM.alphaBandPower.meanCortS + data.NREM.alphaBandPower.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.alphaBandPower.meanCortf,data.NREM.alphaBandPower.meanCortS - data.NREM.alphaBandPower.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.alphaBandPower.meanCortf,data.REM.alphaBandPower.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.alphaBandPower.meanCortf,data.REM.alphaBandPower.meanCortS + data.REM.alphaBandPower.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.alphaBandPower.meanCortf,data.REM.alphaBandPower.meanCortS - data.REM.alphaBandPower.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.alphaBandPower.meanCortf,data.Unstim.alphaBandPower.meanCortS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.alphaBandPower.meanCortf,data.Unstim.alphaBandPower.meanCortS + data.Unstim.alphaBandPower.stdCortS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.alphaBandPower.meanCortf,data.Unstim.alphaBandPower.meanCortS - data.Unstim.alphaBandPower.stdCortS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.alphaBandPower.meanCortf,data.All.alphaBandPower.meanCortS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.alphaBandPower.meanCortf,data.All.alphaBandPower.meanCortS + data.All.alphaBandPower.stdCortS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.alphaBandPower.meanCortf,data.All.alphaBandPower.meanCortS - data.All.alphaBandPower.stdCortS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('Alpha-band power')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    %% Beta-band power
    ax6 = subplot(2,4,6);
    loglog(data.Rest.betaBandPower.(baselineType).meanCortf,data.Rest.betaBandPower.(baselineType).meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.betaBandPower.(baselineType).meanCortf,data.Rest.betaBandPower.(baselineType).meanCortS + data.Rest.betaBandPower.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.betaBandPower.(baselineType).meanCortf,data.Rest.betaBandPower.(baselineType).meanCortS - data.Rest.betaBandPower.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.betaBandPower.meanCortf,data.NREM.betaBandPower.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.betaBandPower.meanCortf,data.NREM.betaBandPower.meanCortS + data.NREM.betaBandPower.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.betaBandPower.meanCortf,data.NREM.betaBandPower.meanCortS - data.NREM.betaBandPower.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.betaBandPower.meanCortf,data.REM.betaBandPower.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.betaBandPower.meanCortf,data.REM.betaBandPower.meanCortS + data.REM.betaBandPower.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.betaBandPower.meanCortf,data.REM.betaBandPower.meanCortS - data.REM.betaBandPower.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.betaBandPower.meanCortf,data.Unstim.betaBandPower.meanCortS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.betaBandPower.meanCortf,data.Unstim.betaBandPower.meanCortS + data.Unstim.betaBandPower.stdCortS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.betaBandPower.meanCortf,data.Unstim.betaBandPower.meanCortS - data.Unstim.betaBandPower.stdCortS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.betaBandPower.meanCortf,data.All.betaBandPower.meanCortS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.betaBandPower.meanCortf,data.All.betaBandPower.meanCortS + data.All.betaBandPower.stdCortS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.betaBandPower.meanCortf,data.All.betaBandPower.meanCortS - data.All.betaBandPower.stdCortS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('Beta-band power')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    %% Gamma-band power
    ax7 = subplot(2,4,7);
    loglog(data.Rest.gammaBandPower.(baselineType).meanCortf,data.Rest.gammaBandPower.(baselineType).meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.gammaBandPower.(baselineType).meanCortf,data.Rest.gammaBandPower.(baselineType).meanCortS + data.Rest.gammaBandPower.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.gammaBandPower.(baselineType).meanCortf,data.Rest.gammaBandPower.(baselineType).meanCortS - data.Rest.gammaBandPower.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.gammaBandPower.meanCortf,data.NREM.gammaBandPower.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.gammaBandPower.meanCortf,data.NREM.gammaBandPower.meanCortS + data.NREM.gammaBandPower.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.gammaBandPower.meanCortf,data.NREM.gammaBandPower.meanCortS - data.NREM.gammaBandPower.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.gammaBandPower.meanCortf,data.REM.gammaBandPower.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.gammaBandPower.meanCortf,data.REM.gammaBandPower.meanCortS + data.REM.gammaBandPower.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.gammaBandPower.meanCortf,data.REM.gammaBandPower.meanCortS - data.REM.gammaBandPower.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.gammaBandPower.meanCortf,data.Unstim.gammaBandPower.meanCortS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.gammaBandPower.meanCortf,data.Unstim.gammaBandPower.meanCortS + data.Unstim.gammaBandPower.stdCortS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.gammaBandPower.meanCortf,data.Unstim.gammaBandPower.meanCortS - data.Unstim.gammaBandPower.stdCortS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.gammaBandPower.meanCortf,data.All.gammaBandPower.meanCortS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.gammaBandPower.meanCortf,data.All.gammaBandPower.meanCortS + data.All.gammaBandPower.stdCortS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.gammaBandPower.meanCortf,data.All.gammaBandPower.meanCortS - data.All.gammaBandPower.stdCortS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('Gamma-band power')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    %% MUA power
    ax8 = subplot(2,4,8);
    loglog(data.Rest.muaPower.(baselineType).meanCortf,data.Rest.muaPower.(baselineType).meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.muaPower.(baselineType).meanCortf,data.Rest.muaPower.(baselineType).meanCortS + data.Rest.muaPower.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.muaPower.(baselineType).meanCortf,data.Rest.muaPower.(baselineType).meanCortS - data.Rest.muaPower.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.muaPower.meanCortf,data.NREM.muaPower.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.muaPower.meanCortf,data.NREM.muaPower.meanCortS + data.NREM.muaPower.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.muaPower.meanCortf,data.NREM.muaPower.meanCortS - data.NREM.muaPower.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.muaPower.meanCortf,data.REM.muaPower.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.muaPower.meanCortf,data.REM.muaPower.meanCortS + data.REM.muaPower.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.muaPower.meanCortf,data.REM.muaPower.meanCortS - data.REM.muaPower.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.muaPower.meanCortf,data.Unstim.muaPower.meanCortS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.muaPower.meanCortf,data.Unstim.muaPower.meanCortS + data.Unstim.muaPower.stdCortS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.muaPower.meanCortf,data.Unstim.muaPower.meanCortS - data.Unstim.muaPower.stdCortS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.muaPower.meanCortf,data.All.muaPower.meanCortS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.muaPower.meanCortf,data.All.muaPower.meanCortS + data.All.muaPower.stdCortS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.muaPower.meanCortf,data.All.muaPower.meanCortS - data.All.muaPower.stdCortS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('MUA power')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    % save figure(s)
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Power Spectra\';
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
    L1 = loglog(data.Rest.CBV.(baselineType).meanCortf,data.Rest.CBV.(baselineType).meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',2);
    hold on
    loglog(data.Rest.CBV.(baselineType).meanCortf,data.Rest.CBV.(baselineType).meanCortS + data.Rest.CBV.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.CBV.(baselineType).meanCortf,data.Rest.CBV.(baselineType).meanCortS - data.Rest.CBV.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    L2 = loglog(data.NREM.CBV.meanCortf,data.NREM.CBV.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.CBV.meanCortf,data.NREM.CBV.meanCortS + data.NREM.CBV.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.CBV.meanCortf,data.NREM.CBV.meanCortS - data.NREM.CBV.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    L3 = loglog(data.REM.CBV.meanCortf,data.REM.CBV.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.CBV.meanCortf,data.REM.CBV.meanCortS + data.REM.CBV.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.CBV.meanCortf,data.REM.CBV.meanCortS - data.REM.CBV.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    L4 = loglog(data.Unstim.CBV.meanCortf,data.Unstim.CBV.meanCortS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.CBV.meanCortf,data.Unstim.CBV.meanCortS + data.Unstim.CBV.stdCortS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.CBV.meanCortf,data.Unstim.CBV.meanCortS - data.Unstim.CBV.stdCortS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    L5 = loglog(data.All.CBV.meanCortf,data.All.CBV.meanCortS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.CBV.meanCortf,data.All.CBV.meanCortS + data.All.CBV.stdCortS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.CBV.meanCortf,data.All.CBV.meanCortS - data.All.CBV.stdCortS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('CBV Reflectance')
    ylabel('Power')
    xlabel('Freq (Hz)')
    legend([L1,L2,L3,L4,L5],'Rest','NREM','REM','Unstim','All')
    axis square
    xlim([0.05 1])
    
    %% CBV HbT
    ax2 = subplot(2,4,2);
    loglog(data.Rest.CBV_HbT.(baselineType).meanCortf,data.Rest.CBV_HbT.(baselineType).meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.CBV_HbT.(baselineType).meanCortf,data.Rest.CBV_HbT.(baselineType).meanCortS + data.Rest.CBV_HbT.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.CBV_HbT.(baselineType).meanCortf,data.Rest.CBV_HbT.(baselineType).meanCortS - data.Rest.CBV_HbT.(baselineType).stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.CBV_HbT.meanCortf,data.NREM.CBV_HbT.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.CBV_HbT.meanCortf,data.NREM.CBV_HbT.meanCortS + data.NREM.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.CBV_HbT.meanCortf,data.NREM.CBV_HbT.meanCortS - data.NREM.CBV_HbT.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.CBV_HbT.meanCortf,data.REM.CBV_HbT.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.CBV_HbT.meanCortf,data.REM.CBV_HbT.meanCortS + data.REM.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.CBV_HbT.meanCortf,data.REM.CBV_HbT.meanCortS - data.REM.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.CBV_HbT.meanCortf,data.Unstim.CBV_HbT.meanCortS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.CBV_HbT.meanCortf,data.Unstim.CBV_HbT.meanCortS + data.Unstim.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.CBV_HbT.meanCortf,data.Unstim.CBV_HbT.meanCortS - data.Unstim.CBV_HbT.stdCortS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.CBV_HbT.meanCortf,data.All.CBV_HbT.meanCortS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.CBV_HbT.meanCortf,data.All.CBV_HbT.meanCortS + data.All.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.CBV_HbT.meanCortf,data.All.CBV_HbT.meanCortS - data.All.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('CBV HbT')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    %% Delta-band power
    ax3 = subplot(2,4,3);
    loglog(data.Rest.deltaBandPower.(baselineType).meanHipf,data.Rest.deltaBandPower.(baselineType).meanHipS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.deltaBandPower.(baselineType).meanHipf,data.Rest.deltaBandPower.(baselineType).meanHipS + data.Rest.deltaBandPower.(baselineType).stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.deltaBandPower.(baselineType).meanHipf,data.Rest.deltaBandPower.(baselineType).meanHipS - data.Rest.deltaBandPower.(baselineType).stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.deltaBandPower.meanHipf,data.NREM.deltaBandPower.meanHipS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.deltaBandPower.meanHipf,data.NREM.deltaBandPower.meanHipS + data.NREM.deltaBandPower.stdHipS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.deltaBandPower.meanHipf,data.NREM.deltaBandPower.meanHipS - data.NREM.deltaBandPower.stdHipS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.deltaBandPower.meanHipf,data.REM.deltaBandPower.meanHipS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.deltaBandPower.meanHipf,data.REM.deltaBandPower.meanHipS + data.REM.deltaBandPower.stdHipS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.deltaBandPower.meanHipf,data.REM.deltaBandPower.meanHipS - data.REM.deltaBandPower.stdHipS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.deltaBandPower.meanHipf,data.Unstim.deltaBandPower.meanHipS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.deltaBandPower.meanHipf,data.Unstim.deltaBandPower.meanHipS + data.Unstim.deltaBandPower.stdCortS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.deltaBandPower.meanHipf,data.Unstim.deltaBandPower.meanHipS - data.Unstim.deltaBandPower.stdHipS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.deltaBandPower.meanHipf,data.All.deltaBandPower.meanHipS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.deltaBandPower.meanHipf,data.All.deltaBandPower.meanHipS + data.All.deltaBandPower.stdHipS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.deltaBandPower.meanHipf,data.All.deltaBandPower.meanHipS - data.All.deltaBandPower.stdHipS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('Delta-band power')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    %% Theta-band power
    ax4 = subplot(2,4,4);
    loglog(data.Rest.thetaBandPower.(baselineType).meanHipf,data.Rest.thetaBandPower.(baselineType).meanHipS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.thetaBandPower.(baselineType).meanHipf,data.Rest.thetaBandPower.(baselineType).meanHipS + data.Rest.thetaBandPower.(baselineType).stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.thetaBandPower.(baselineType).meanHipf,data.Rest.thetaBandPower.(baselineType).meanHipS - data.Rest.thetaBandPower.(baselineType).stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.thetaBandPower.meanHipf,data.NREM.thetaBandPower.meanHipS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.thetaBandPower.meanHipf,data.NREM.thetaBandPower.meanHipS + data.NREM.thetaBandPower.stdHipS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.thetaBandPower.meanHipf,data.NREM.thetaBandPower.meanHipS - data.NREM.thetaBandPower.stdHipS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.thetaBandPower.meanHipf,data.REM.thetaBandPower.meanHipS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.thetaBandPower.meanHipf,data.REM.thetaBandPower.meanHipS + data.REM.thetaBandPower.stdHipS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.thetaBandPower.meanHipf,data.REM.thetaBandPower.meanHipS - data.REM.thetaBandPower.stdHipS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.thetaBandPower.meanHipf,data.Unstim.thetaBandPower.meanHipS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.thetaBandPower.meanHipf,data.Unstim.thetaBandPower.meanHipS + data.Unstim.thetaBandPower.stdHipS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.thetaBandPower.meanHipf,data.Unstim.thetaBandPower.meanHipS - data.Unstim.thetaBandPower.stdHipS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.thetaBandPower.meanHipf,data.All.thetaBandPower.meanHipS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.thetaBandPower.meanHipf,data.All.thetaBandPower.meanHipS + data.All.thetaBandPower.stdHipS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.thetaBandPower.meanHipf,data.All.thetaBandPower.meanHipS - data.All.thetaBandPower.stdHipS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('Theta-band power')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    %% Alpha-band power
    ax5 = subplot(2,4,5);
    loglog(data.Rest.alphaBandPower.(baselineType).meanHipf,data.Rest.alphaBandPower.(baselineType).meanHipS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.alphaBandPower.(baselineType).meanHipf,data.Rest.alphaBandPower.(baselineType).meanHipS + data.Rest.alphaBandPower.(baselineType).stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.alphaBandPower.(baselineType).meanHipf,data.Rest.alphaBandPower.(baselineType).meanHipS - data.Rest.alphaBandPower.(baselineType).stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.alphaBandPower.meanHipf,data.NREM.alphaBandPower.meanHipS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.alphaBandPower.meanHipf,data.NREM.alphaBandPower.meanHipS + data.NREM.alphaBandPower.stdHipS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.alphaBandPower.meanHipf,data.NREM.alphaBandPower.meanHipS - data.NREM.alphaBandPower.stdHipS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.alphaBandPower.meanHipf,data.REM.alphaBandPower.meanHipS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.alphaBandPower.meanHipf,data.REM.alphaBandPower.meanHipS + data.REM.alphaBandPower.stdHipS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.alphaBandPower.meanHipf,data.REM.alphaBandPower.meanHipS - data.REM.alphaBandPower.stdHipS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.alphaBandPower.meanHipf,data.Unstim.alphaBandPower.meanHipS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.alphaBandPower.meanHipf,data.Unstim.alphaBandPower.meanHipS + data.Unstim.alphaBandPower.stdHipS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.alphaBandPower.meanHipf,data.Unstim.alphaBandPower.meanHipS - data.Unstim.alphaBandPower.stdHipS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.alphaBandPower.meanHipf,data.All.alphaBandPower.meanHipS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.alphaBandPower.meanHipf,data.All.alphaBandPower.meanHipS + data.All.alphaBandPower.stdHipS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.alphaBandPower.meanHipf,data.All.alphaBandPower.meanHipS - data.All.alphaBandPower.stdHipS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('Alpha-band power')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    %% Beta-band power
    ax6 = subplot(2,4,6);
    loglog(data.Rest.betaBandPower.(baselineType).meanHipf,data.Rest.betaBandPower.(baselineType).meanHipS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.betaBandPower.(baselineType).meanHipf,data.Rest.betaBandPower.(baselineType).meanHipS + data.Rest.betaBandPower.(baselineType).stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.betaBandPower.(baselineType).meanHipf,data.Rest.betaBandPower.(baselineType).meanHipS - data.Rest.betaBandPower.(baselineType).stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.betaBandPower.meanHipf,data.NREM.betaBandPower.meanHipS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.betaBandPower.meanHipf,data.NREM.betaBandPower.meanHipS + data.NREM.betaBandPower.stdHipS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.betaBandPower.meanHipf,data.NREM.betaBandPower.meanHipS - data.NREM.betaBandPower.stdHipS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.betaBandPower.meanHipf,data.REM.betaBandPower.meanHipS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.betaBandPower.meanHipf,data.REM.betaBandPower.meanHipS + data.REM.betaBandPower.stdHipS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.betaBandPower.meanHipf,data.REM.betaBandPower.meanHipS - data.REM.betaBandPower.stdHipS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.betaBandPower.meanHipf,data.Unstim.betaBandPower.meanHipS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.betaBandPower.meanHipf,data.Unstim.betaBandPower.meanHipS + data.Unstim.betaBandPower.stdHipS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.betaBandPower.meanHipf,data.Unstim.betaBandPower.meanHipS - data.Unstim.betaBandPower.stdHipS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.betaBandPower.meanHipf,data.All.betaBandPower.meanHipS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.betaBandPower.meanHipf,data.All.betaBandPower.meanHipS + data.All.betaBandPower.stdHipS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.betaBandPower.meanHipf,data.All.betaBandPower.meanHipS - data.All.betaBandPower.stdHipS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('Beta-band power')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    %% Gamma-band power
    ax7 = subplot(2,4,7);
    loglog(data.Rest.gammaBandPower.(baselineType).meanHipf,data.Rest.gammaBandPower.(baselineType).meanHipS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.gammaBandPower.(baselineType).meanHipf,data.Rest.gammaBandPower.(baselineType).meanHipS + data.Rest.gammaBandPower.(baselineType).stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.gammaBandPower.(baselineType).meanHipf,data.Rest.gammaBandPower.(baselineType).meanHipS - data.Rest.gammaBandPower.(baselineType).stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.gammaBandPower.meanHipf,data.NREM.gammaBandPower.meanHipS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.gammaBandPower.meanHipf,data.NREM.gammaBandPower.meanHipS + data.NREM.gammaBandPower.stdHipS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.gammaBandPower.meanHipf,data.NREM.gammaBandPower.meanHipS - data.NREM.gammaBandPower.stdHipS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.gammaBandPower.meanHipf,data.REM.gammaBandPower.meanHipS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.gammaBandPower.meanHipf,data.REM.gammaBandPower.meanHipS + data.REM.gammaBandPower.stdHipS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.gammaBandPower.meanHipf,data.REM.gammaBandPower.meanHipS - data.REM.gammaBandPower.stdHipS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.gammaBandPower.meanHipf,data.Unstim.gammaBandPower.meanHipS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.gammaBandPower.meanHipf,data.Unstim.gammaBandPower.meanHipS + data.Unstim.gammaBandPower.stdHipS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.gammaBandPower.meanHipf,data.Unstim.gammaBandPower.meanHipS - data.Unstim.gammaBandPower.stdHipS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.gammaBandPower.meanHipf,data.All.gammaBandPower.meanHipS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.gammaBandPower.meanHipf,data.All.gammaBandPower.meanHipS + data.All.gammaBandPower.stdHipS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.gammaBandPower.meanHipf,data.All.gammaBandPower.meanHipS - data.All.gammaBandPower.stdHipS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('Gamma-band power')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
    %% MUA power
    ax8 = subplot(2,4,8);
    loglog(data.Rest.muaPower.(baselineType).meanHipf,data.Rest.muaPower.(baselineType).meanHipS,'color',colorbrewer_setA_colorA,'LineWidth',2)
    hold on
    loglog(data.Rest.muaPower.(baselineType).meanHipf,data.Rest.muaPower.(baselineType).meanHipS + data.Rest.muaPower.(baselineType).stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.Rest.muaPower.(baselineType).meanHipf,data.Rest.muaPower.(baselineType).meanHipS - data.Rest.muaPower.(baselineType).stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
    loglog(data.NREM.muaPower.meanHipf,data.NREM.muaPower.meanHipS,'color',colorbrewer_setA_colorB,'LineWidth',2);
    loglog(data.NREM.muaPower.meanHipf,data.NREM.muaPower.meanHipS + data.NREM.muaPower.stdHipS,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.NREM.muaPower.meanHipf,data.NREM.muaPower.meanHipS - data.NREM.muaPower.stdHipS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
    loglog(data.REM.muaPower.meanHipf,data.REM.muaPower.meanHipS,'color',colorbrewer_setA_colorC,'LineWidth',2);
    loglog(data.REM.muaPower.meanHipf,data.REM.muaPower.meanHipS + data.REM.muaPower.stdHipS,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.REM.muaPower.meanHipf,data.REM.muaPower.meanHipS - data.REM.muaPower.stdHipS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
    loglog(data.Unstim.muaPower.meanHipf,data.Unstim.muaPower.meanHipS,'color',colorbrewer_setA_colorD,'LineWidth',2);
    loglog(data.Unstim.muaPower.meanHipf,data.Unstim.muaPower.meanHipS + data.Unstim.muaPower.stdHipS,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.Unstim.muaPower.meanHipf,data.Unstim.muaPower.meanHipS - data.Unstim.muaPower.stdHipS ,'color',colorbrewer_setA_colorD,'LineWidth',1)
    loglog(data.All.muaPower.meanHipf,data.All.muaPower.meanHipS,'color',colorbrewer_setA_colorE,'LineWidth',2);
    loglog(data.All.muaPower.meanHipf,data.All.muaPower.meanHipS + data.All.muaPower.stdHipS,'color',colorbrewer_setA_colorE,'LineWidth',1)
    loglog(data.All.muaPower.meanHipf,data.All.muaPower.meanHipS - data.All.muaPower.stdHipS ,'color',colorbrewer_setA_colorE,'LineWidth',1)
    title('MUA power')
    ylabel('Power')
    xlabel('Freq (Hz)')
    axis square
    xlim([0.05 1])
    
     % save figure(s)
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Power Spectra\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end 
    savefig(summaryFigure, [dirpath baselineType '_HippocampalPowerSpectra']);
end
 