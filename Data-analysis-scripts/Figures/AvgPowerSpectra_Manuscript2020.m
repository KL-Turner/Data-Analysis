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
driveLetters = {'M','M','M','M','M','M','M','M','M'};
behavFields = {'Rest','NREM','REM'};
powerSpec_dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
colorbrewer_setA_colorB = [(31/256) (120/256) (180/256)];
colorbrewer_setA_colorA = [(51/256) (160/256) (44/256)];
colorbrewer_setA_colorC = [(255/256) (140/256) (0/256)];

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\Turner_Manuscript_Summer2020\' animalID '\Bilateral Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        for c = 1:length(powerSpec_dataTypes)
            powerSpec_dataType = powerSpec_dataTypes{1,c};
            data.(behavField).(powerSpec_dataType).adjLH.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).adjLH.S;
            data.(behavField).(powerSpec_dataType).adjLH.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).adjLH.f;
            data.(behavField).(powerSpec_dataType).adjRH.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).adjRH.S;
            data.(behavField).(powerSpec_dataType).adjRH.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).adjRH.f;
            if strcmp(powerSpec_dataType,'CBV_HbT') == false
                data.(behavField).(powerSpec_dataType).Hip.S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).Hip.S;
                data.(behavField).(powerSpec_dataType).Hip.f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerSpec_dataType).Hip.f;
            end
        end
    end
end

% concatenate the data from the left and right hemispheres
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(powerSpec_dataTypes)
        powerSpec_dataType = powerSpec_dataTypes{1,f};
        data.(behavField).(powerSpec_dataType).cat_S = cat(2,data.(behavField).(powerSpec_dataType).adjLH.S,data.(behavField).(powerSpec_dataType).adjRH.S);
        data.(behavField).(powerSpec_dataType).cat_f = cat(2,data.(behavField).(powerSpec_dataType).adjLH.f,data.(behavField).(powerSpec_dataType).adjRH.f);
    end
end

% take the mean and standard deviation of each set of signals
for h = 1:length(behavFields)
    behavField = behavFields{1,h};
    for j = 1:length(powerSpec_dataTypes)
        powerSpec_dataType = powerSpec_dataTypes{1,j};
        data.(behavField).(powerSpec_dataType).meanCortS = mean(data.(behavField).(powerSpec_dataType).cat_S,2);
        data.(behavField).(powerSpec_dataType).stdCortS = std(data.(behavField).(powerSpec_dataType).cat_S,0,2);
        data.(behavField).(powerSpec_dataType).meanCortf = mean(data.(behavField).(powerSpec_dataType).cat_f,2);
        if strcmp(powerSpec_dataType,'CBV_HbT') == false
            data.(behavField).(powerSpec_dataType).meanHipS = mean(data.(behavField).(powerSpec_dataType).Hip.S,2);
            data.(behavField).(powerSpec_dataType).stdHipS = std(data.(behavField).(powerSpec_dataType).Hip.S,0,2);
            data.(behavField).(powerSpec_dataType).meanHipf = mean(data.(behavField).(powerSpec_dataType).Hip.f,2);
        end
    end
end

%% summary figure(s)
summaryFigure = figure;
sgtitle('Power Spectral Density')
%% CBV HbT
ax1 = subplot(2,6,[1,7]);
L1 = loglog(data.Rest.CBV_HbT.meanCortf,data.Rest.CBV_HbT.meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',3);
hold on
% loglog(data.Rest.CBV_HbT.meanCortf,data.Rest.CBV_HbT.meanCortS + data.Rest.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
% loglog(data.Rest.CBV_HbT.meanCortf,data.Rest.CBV_HbT.meanCortS - data.Rest.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
L2 = loglog(data.NREM.CBV_HbT.meanCortf,data.NREM.CBV_HbT.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',3);
% loglog(data.NREM.CBV_HbT.meanCortf,data.NREM.CBV_HbT.meanCortS + data.NREM.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
% loglog(data.NREM.CBV_HbT.meanCortf,data.NREM.CBV_HbT.meanCortS - data.NREM.CBV_HbT.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
L3 = loglog(data.REM.CBV_HbT.meanCortf,data.REM.CBV_HbT.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',3);
% loglog(data.REM.CBV_HbT.meanCortf,data.REM.CBV_HbT.meanCortS + data.REM.CBV_HbT.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
% loglog(data.REM.CBV_HbT.meanCortf,data.REM.CBV_HbT.meanCortS - data.REM.CBV_HbT.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
title('\DeltaHbT (\muM)')
ylabel('Power')
xlabel('Frequency (Hz)')
legend([L1,L2,L3],'Awake Rest','NREM','REM','Location','SouthWest')
axis square
xlim([0.05 1])
set(gca,'box','off')

%% Delta-band power
ax2 = subplot(2,6,2);
loglog(data.Rest.deltaBandPower.meanCortf,data.Rest.deltaBandPower.meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',3)
hold on
% loglog(data.Rest.deltaBandPower.meanCortf,data.Rest.deltaBandPower.meanCortS + data.Rest.deltaBandPower.stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
% loglog(data.Rest.deltaBandPower.meanCortf,data.Rest.deltaBandPower.meanCortS - data.Rest.deltaBandPower.stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
loglog(data.NREM.deltaBandPower.meanCortf,data.NREM.deltaBandPower.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',3);
% loglog(data.NREM.deltaBandPower.meanCortf,data.NREM.deltaBandPower.meanCortS + data.NREM.deltaBandPower.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
% loglog(data.NREM.deltaBandPower.meanCortf,data.NREM.deltaBandPower.meanCortS - data.NREM.deltaBandPower.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
loglog(data.REM.deltaBandPower.meanCortf,data.REM.deltaBandPower.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',3);
% loglog(data.REM.deltaBandPower.meanCortf,data.REM.deltaBandPower.meanCortS + data.REM.deltaBandPower.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
% loglog(data.REM.deltaBandPower.meanCortf,data.REM.deltaBandPower.meanCortS - data.REM.deltaBandPower.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
title({'Cortical';'Delta-band [1-4 Hz]'})
ylabel('Power')
xlabel('Frequency (Hz)')
axis square
xlim([0.05 1])
set(gca,'box','off')

%% Theta-band power
ax3 = subplot(2,6,3);
loglog(data.Rest.thetaBandPower.meanCortf,data.Rest.thetaBandPower.meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',3)
hold on
% loglog(data.Rest.thetaBandPower.meanCortf,data.Rest.thetaBandPower.meanCortS + data.Rest.thetaBandPower.stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
% loglog(data.Rest.thetaBandPower.meanCortf,data.Rest.thetaBandPower.meanCortS - data.Rest.thetaBandPower.stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
loglog(data.NREM.thetaBandPower.meanCortf,data.NREM.thetaBandPower.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',3);
% loglog(data.NREM.thetaBandPower.meanCortf,data.NREM.thetaBandPower.meanCortS + data.NREM.thetaBandPower.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
% loglog(data.NREM.thetaBandPower.meanCortf,data.NREM.thetaBandPower.meanCortS - data.NREM.thetaBandPower.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
loglog(data.REM.thetaBandPower.meanCortf,data.REM.thetaBandPower.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',3);
% loglog(data.REM.thetaBandPower.meanCortf,data.REM.thetaBandPower.meanCortS + data.REM.thetaBandPower.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
% loglog(data.REM.thetaBandPower.meanCortf,data.REM.thetaBandPower.meanCortS - data.REM.thetaBandPower.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
title({'Cortical';'Theta-band [4-10 Hz]'})
ylabel('Power')
xlabel('Frequency (Hz)')
axis square
xlim([0.05 1])
set(gca,'box','off')

%% Alpha-band power
ax4 = subplot(2,6,4);
loglog(data.Rest.alphaBandPower.meanCortf,data.Rest.alphaBandPower.meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',3)
hold on
% loglog(data.Rest.alphaBandPower.meanCortf,data.Rest.alphaBandPower.meanCortS + data.Rest.alphaBandPower.stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
% loglog(data.Rest.alphaBandPower.meanCortf,data.Rest.alphaBandPower.meanCortS - data.Rest.alphaBandPower.stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
loglog(data.NREM.alphaBandPower.meanCortf,data.NREM.alphaBandPower.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',3);
% loglog(data.NREM.alphaBandPower.meanCortf,data.NREM.alphaBandPower.meanCortS + data.NREM.alphaBandPower.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
% loglog(data.NREM.alphaBandPower.meanCortf,data.NREM.alphaBandPower.meanCortS - data.NREM.alphaBandPower.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
loglog(data.REM.alphaBandPower.meanCortf,data.REM.alphaBandPower.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',3);
% loglog(data.REM.alphaBandPower.meanCortf,data.REM.alphaBandPower.meanCortS + data.REM.alphaBandPower.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
% loglog(data.REM.alphaBandPower.meanCortf,data.REM.alphaBandPower.meanCortS - data.REM.alphaBandPower.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
title({'Cortical';'Alpha-band [10-13 Hz]'})
ylabel('Power')
xlabel('Frequency (Hz)')
axis square
xlim([0.05 1])
set(gca,'box','off')

%% Beta-band power
ax5 = subplot(2,6,5);
loglog(data.Rest.betaBandPower.meanCortf,data.Rest.betaBandPower.meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',3)
hold on
% loglog(data.Rest.betaBandPower.meanCortf,data.Rest.betaBandPower.meanCortS + data.Rest.betaBandPower.stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
% loglog(data.Rest.betaBandPower.meanCortf,data.Rest.betaBandPower.meanCortS - data.Rest.betaBandPower.stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
loglog(data.NREM.betaBandPower.meanCortf,data.NREM.betaBandPower.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',3);
% loglog(data.NREM.betaBandPower.meanCortf,data.NREM.betaBandPower.meanCortS + data.NREM.betaBandPower.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
% loglog(data.NREM.betaBandPower.meanCortf,data.NREM.betaBandPower.meanCortS - data.NREM.betaBandPower.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
loglog(data.REM.betaBandPower.meanCortf,data.REM.betaBandPower.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',3);
% loglog(data.REM.betaBandPower.meanCortf,data.REM.betaBandPower.meanCortS + data.REM.betaBandPower.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
% loglog(data.REM.betaBandPower.meanCortf,data.REM.betaBandPower.meanCortS - data.REM.betaBandPower.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
title({'Cortical';'Beta-band [13-30 Hz]'})
ylabel('Power')
xlabel('Frequency (Hz)')
axis square
xlim([0.05 1])
set(gca,'box','off')

%% Gamma-band power
ax6 = subplot(2,6,6);
loglog(data.Rest.gammaBandPower.meanCortf,data.Rest.gammaBandPower.meanCortS,'color',colorbrewer_setA_colorA,'LineWidth',3)
hold on
% loglog(data.Rest.gammaBandPower.meanCortf,data.Rest.gammaBandPower.meanCortS + data.Rest.gammaBandPower.stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
% loglog(data.Rest.gammaBandPower.meanCortf,data.Rest.gammaBandPower.meanCortS - data.Rest.gammaBandPower.stdCortS,'color',colorbrewer_setA_colorA,'LineWidth',1)
loglog(data.NREM.gammaBandPower.meanCortf,data.NREM.gammaBandPower.meanCortS,'color',colorbrewer_setA_colorB,'LineWidth',3);
% loglog(data.NREM.gammaBandPower.meanCortf,data.NREM.gammaBandPower.meanCortS + data.NREM.gammaBandPower.stdCortS,'color',colorbrewer_setA_colorB,'LineWidth',1)
% loglog(data.NREM.gammaBandPower.meanCortf,data.NREM.gammaBandPower.meanCortS - data.NREM.gammaBandPower.stdCortS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
loglog(data.REM.gammaBandPower.meanCortf,data.REM.gammaBandPower.meanCortS,'color',colorbrewer_setA_colorC,'LineWidth',3);
% loglog(data.REM.gammaBandPower.meanCortf,data.REM.gammaBandPower.meanCortS + data.REM.gammaBandPower.stdCortS,'color',colorbrewer_setA_colorC,'LineWidth',1)
% loglog(data.REM.gammaBandPower.meanCortf,data.REM.gammaBandPower.meanCortS - data.REM.gammaBandPower.stdCortS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
title({'Cortical';'Gamma-band [30-100 Hz]'})
ylabel('Power')
xlabel('Frequency (Hz)')
axis square
xlim([0.05 1])
set(gca,'box','off')

%% Delta-band power
ax8 = subplot(2,6,8);
loglog(data.Rest.deltaBandPower.meanHipf,data.Rest.deltaBandPower.meanHipS,'color',colorbrewer_setA_colorA,'LineWidth',3)
hold on
% loglog(data.Rest.deltaBandPower.meanHipf,data.Rest.deltaBandPower.meanHipS + data.Rest.deltaBandPower.stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
% loglog(data.Rest.deltaBandPower.meanHipf,data.Rest.deltaBandPower.meanHipS - data.Rest.deltaBandPower.stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
loglog(data.NREM.deltaBandPower.meanHipf,data.NREM.deltaBandPower.meanHipS,'color',colorbrewer_setA_colorB,'LineWidth',3);
% loglog(data.NREM.deltaBandPower.meanHipf,data.NREM.deltaBandPower.meanHipS + data.NREM.deltaBandPower.stdHipS,'color',colorbrewer_setA_colorB,'LineWidth',1)
% loglog(data.NREM.deltaBandPower.meanHipf,data.NREM.deltaBandPower.meanHipS - data.NREM.deltaBandPower.stdHipS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
loglog(data.REM.deltaBandPower.meanHipf,data.REM.deltaBandPower.meanHipS,'color',colorbrewer_setA_colorC,'LineWidth',3);
% loglog(data.REM.deltaBandPower.meanHipf,data.REM.deltaBandPower.meanHipS + data.REM.deltaBandPower.stdHipS,'color',colorbrewer_setA_colorC,'LineWidth',1)
% loglog(data.REM.deltaBandPower.meanHipf,data.REM.deltaBandPower.meanHipS - data.REM.deltaBandPower.stdHipS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
title({'Hippocampal';'Delta-band [1-4 Hz]'})
ylabel('Power')
xlabel('Frequency (Hz)')
axis square
xlim([0.05 1])
set(gca,'box','off')

%% Theta-band power
ax9 = subplot(2,6,9);
loglog(data.Rest.thetaBandPower.meanHipf,data.Rest.thetaBandPower.meanHipS,'color',colorbrewer_setA_colorA,'LineWidth',3)
hold on
% loglog(data.Rest.thetaBandPower.meanHipf,data.Rest.thetaBandPower.meanHipS + data.Rest.thetaBandPower.stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
% loglog(data.Rest.thetaBandPower.meanHipf,data.Rest.thetaBandPower.meanHipS - data.Rest.thetaBandPower.stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
loglog(data.NREM.thetaBandPower.meanHipf,data.NREM.thetaBandPower.meanHipS,'color',colorbrewer_setA_colorB,'LineWidth',3);
% loglog(data.NREM.thetaBandPower.meanHipf,data.NREM.thetaBandPower.meanHipS + data.NREM.thetaBandPower.stdHipS,'color',colorbrewer_setA_colorB,'LineWidth',1)
% loglog(data.NREM.thetaBandPower.meanHipf,data.NREM.thetaBandPower.meanHipS - data.NREM.thetaBandPower.stdHipS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
loglog(data.REM.thetaBandPower.meanHipf,data.REM.thetaBandPower.meanHipS,'color',colorbrewer_setA_colorC,'LineWidth',3);
% loglog(data.REM.thetaBandPower.meanHipf,data.REM.thetaBandPower.meanHipS + data.REM.thetaBandPower.stdHipS,'color',colorbrewer_setA_colorC,'LineWidth',1)
% loglog(data.REM.thetaBandPower.meanHipf,data.REM.thetaBandPower.meanHipS - data.REM.thetaBandPower.stdHipS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
title({'Hippocampal';'Theta-band [4-10 Hz]'})
ylabel('Power')
xlabel('Frequency (Hz)')
axis square
xlim([0.05 1])
set(gca,'box','off')

%% Alpha-band power
ax10 = subplot(2,6,10);
loglog(data.Rest.alphaBandPower.meanHipf,data.Rest.alphaBandPower.meanHipS,'color',colorbrewer_setA_colorA,'LineWidth',3)
hold on
% loglog(data.Rest.alphaBandPower.meanHipf,data.Rest.alphaBandPower.meanHipS + data.Rest.alphaBandPower.stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
% loglog(data.Rest.alphaBandPower.meanHipf,data.Rest.alphaBandPower.meanHipS - data.Rest.alphaBandPower.stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
loglog(data.NREM.alphaBandPower.meanHipf,data.NREM.alphaBandPower.meanHipS,'color',colorbrewer_setA_colorB,'LineWidth',3);
% loglog(data.NREM.alphaBandPower.meanHipf,data.NREM.alphaBandPower.meanHipS + data.NREM.alphaBandPower.stdHipS,'color',colorbrewer_setA_colorB,'LineWidth',1)
% loglog(data.NREM.alphaBandPower.meanHipf,data.NREM.alphaBandPower.meanHipS - data.NREM.alphaBandPower.stdHipS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
loglog(data.REM.alphaBandPower.meanHipf,data.REM.alphaBandPower.meanHipS,'color',colorbrewer_setA_colorC,'LineWidth',3);
% loglog(data.REM.alphaBandPower.meanHipf,data.REM.alphaBandPower.meanHipS + data.REM.alphaBandPower.stdHipS,'color',colorbrewer_setA_colorC,'LineWidth',1)
% loglog(data.REM.alphaBandPower.meanHipf,data.REM.alphaBandPower.meanHipS - data.REM.alphaBandPower.stdHipS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
title({'Hippocampal';'Alpha-band [10-13 Hz]'})
ylabel('Power')
xlabel('Frequency (Hz)')
axis square
xlim([0.05 1])
set(gca,'box','off')

%% Beta-band power
ax11 = subplot(2,6,11);
loglog(data.Rest.betaBandPower.meanHipf,data.Rest.betaBandPower.meanHipS,'color',colorbrewer_setA_colorA,'LineWidth',3)
hold on
% loglog(data.Rest.betaBandPower.meanHipf,data.Rest.betaBandPower.meanHipS + data.Rest.betaBandPower.stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
% loglog(data.Rest.betaBandPower.meanHipf,data.Rest.betaBandPower.meanHipS - data.Rest.betaBandPower.stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
loglog(data.NREM.betaBandPower.meanHipf,data.NREM.betaBandPower.meanHipS,'color',colorbrewer_setA_colorB,'LineWidth',3);
% loglog(data.NREM.betaBandPower.meanHipf,data.NREM.betaBandPower.meanHipS + data.NREM.betaBandPower.stdHipS,'color',colorbrewer_setA_colorB,'LineWidth',1)
% loglog(data.NREM.betaBandPower.meanHipf,data.NREM.betaBandPower.meanHipS - data.NREM.betaBandPower.stdHipS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
loglog(data.REM.betaBandPower.meanHipf,data.REM.betaBandPower.meanHipS,'color',colorbrewer_setA_colorC,'LineWidth',3);
% loglog(data.REM.betaBandPower.meanHipf,data.REM.betaBandPower.meanHipS + data.REM.betaBandPower.stdHipS,'color',colorbrewer_setA_colorC,'LineWidth',1)
% loglog(data.REM.betaBandPower.meanHipf,data.REM.betaBandPower.meanHipS - data.REM.betaBandPower.stdHipS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
title({'Hippocampal';'Beta-band [13-30 Hz]'})
ylabel('Power')
xlabel('Frequency (Hz)')
axis square
xlim([0.05 1])
set(gca,'box','off')

%% Gamma-band power
ax12 = subplot(2,6,12);
loglog(data.Rest.gammaBandPower.meanHipf,data.Rest.gammaBandPower.meanHipS,'color',colorbrewer_setA_colorA,'LineWidth',3)
hold on
% loglog(data.Rest.gammaBandPower.meanHipf,data.Rest.gammaBandPower.meanHipS + data.Rest.gammaBandPower.stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
% loglog(data.Rest.gammaBandPower.meanHipf,data.Rest.gammaBandPower.meanHipS - data.Rest.gammaBandPower.stdHipS,'color',colorbrewer_setA_colorA,'LineWidth',1)
loglog(data.NREM.gammaBandPower.meanHipf,data.NREM.gammaBandPower.meanHipS,'color',colorbrewer_setA_colorB,'LineWidth',3);
% loglog(data.NREM.gammaBandPower.meanHipf,data.NREM.gammaBandPower.meanHipS + data.NREM.gammaBandPower.stdHipS,'color',colorbrewer_setA_colorB,'LineWidth',1)
% loglog(data.NREM.gammaBandPower.meanHipf,data.NREM.gammaBandPower.meanHipS - data.NREM.gammaBandPower.stdHipS ,'color',colorbrewer_setA_colorB,'LineWidth',1)
loglog(data.REM.gammaBandPower.meanHipf,data.REM.gammaBandPower.meanHipS,'color',colorbrewer_setA_colorC,'LineWidth',3);
% loglog(data.REM.gammaBandPower.meanHipf,data.REM.gammaBandPower.meanHipS + data.REM.gammaBandPower.stdHipS,'color',colorbrewer_setA_colorC,'LineWidth',1)
% loglog(data.REM.gammaBandPower.meanHipf,data.REM.gammaBandPower.meanHipS - data.REM.gammaBandPower.stdHipS ,'color',colorbrewer_setA_colorC,'LineWidth',1)
title({'Hippocampal';'Gamma-band [30-100 Hz]'})
ylabel('Power')
xlabel('Frequency (Hz)')
axis square
xlim([0.05 1])
set(gca,'box','off')

linkaxes([ax2,ax3,ax4,ax5,ax6],'xy')
linkaxes([ax8,ax9,ax10,ax11,ax12],'xy')

% save figure(s)
dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\';
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath  'Summary Figure - Power Spectra']);
