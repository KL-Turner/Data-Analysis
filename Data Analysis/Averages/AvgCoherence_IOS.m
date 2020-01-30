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

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111'};
driveLetters = {'M','M','M','M','M','M','M','M','M'};
behavFields = {'Rest','NREM','REM'};
coherr_dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
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
        for c = 1:length(coherr_dataTypes)
            coherr_dataType = coherr_dataTypes{1,c};
            data.(behavField).(coherr_dataType).C(:,a) = AnalysisResults.Coherence.(behavField).(coherr_dataType).C;
            data.(behavField).(coherr_dataType).f(:,a) = AnalysisResults.Coherence.(behavField).(coherr_dataType).f;
            data.(behavField).(coherr_dataType).confC(:,a) = AnalysisResults.Coherence.(behavField).(coherr_dataType).confC;
        end
    end
end

% take the mean and standard deviation of each set of signals
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(coherr_dataTypes)
        coherr_dataType = coherr_dataTypes{1,f};
        data.(behavField).(coherr_dataType).meanC = mean(data.(behavField).(coherr_dataType).C,2);
        data.(behavField).(coherr_dataType).stdC = std(data.(behavField).(coherr_dataType).C,0,2);
        data.(behavField).(coherr_dataType).meanf = mean(data.(behavField).(coherr_dataType).f,2);
        data.(behavField).(coherr_dataType).maxConfC = max(data.(behavField).(coherr_dataType).confC);
        data.(behavField).(coherr_dataType).maxConfC_Y = ones(length(data.(behavField).(coherr_dataType).meanf),1)*data.(behavField).(coherr_dataType).maxConfC;
    end
end

%% summary figure(s)
summaryFigure = figure;
sgtitle('Spectral Coherence Between Bilateral Signals')
%% CBV_HbT
subplot(2,3,1);
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.meanC,'color',colorbrewer_setA_colorA,'LineWidth',3);
hold on
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.meanC + data.Rest.CBV_HbT.stdC,'color',colorbrewer_setA_colorA,'LineWidth',1)
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.meanC - data.Rest.CBV_HbT.stdC,'color',colorbrewer_setA_colorA,'LineWidth',1)
semilogx(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.meanC,'color',colorbrewer_setA_colorB,'LineWidth',3);
semilogx(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.meanC + data.NREM.CBV_HbT.stdC,'color',colorbrewer_setA_colorB,'LineWidth',1)
semilogx(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.meanC - data.NREM.CBV_HbT.stdC,'color',colorbrewer_setA_colorB,'LineWidth',1)
semilogx(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.meanC,'color',colorbrewer_setA_colorC,'LineWidth',3);
semilogx(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.meanC + data.REM.CBV_HbT.stdC,'color',colorbrewer_setA_colorC,'LineWidth',1)
semilogx(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.meanC - data.REM.CBV_HbT.stdC,'color',colorbrewer_setA_colorC,'LineWidth',1)
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorA,'LineWidth',1);
semilogx(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorB,'LineWidth',1);
semilogx(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorC,'LineWidth',1);
title('\DeltaHbT (\muM)')
ylabel('Coherence')
xlabel('Frequency (Hz)')
axis square
ylim([0 1])
xlim([0.05 1])
set(gca,'box','off')

%% Delta-band power
subplot(2,3,2);
semilogx(data.Rest.deltaBandPower.meanf,data.Rest.deltaBandPower.meanC,'color',colorbrewer_setA_colorA,'LineWidth',3);
hold on
semilogx(data.Rest.deltaBandPower.meanf,data.Rest.deltaBandPower.meanC + data.Rest.deltaBandPower.stdC,'color',colorbrewer_setA_colorA,'LineWidth',1)
semilogx(data.Rest.deltaBandPower.meanf,data.Rest.deltaBandPower.meanC - data.Rest.deltaBandPower.stdC,'color',colorbrewer_setA_colorA,'LineWidth',1)
semilogx(data.NREM.deltaBandPower.meanf,data.NREM.deltaBandPower.meanC,'color',colorbrewer_setA_colorB,'LineWidth',3);
semilogx(data.NREM.deltaBandPower.meanf,data.NREM.deltaBandPower.meanC + data.NREM.deltaBandPower.stdC,'color',colorbrewer_setA_colorB,'LineWidth',1)
semilogx(data.NREM.deltaBandPower.meanf,data.NREM.deltaBandPower.meanC - data.NREM.deltaBandPower.stdC,'color',colorbrewer_setA_colorB,'LineWidth',1)
semilogx(data.REM.deltaBandPower.meanf,data.REM.deltaBandPower.meanC,'color',colorbrewer_setA_colorC,'LineWidth',3);
semilogx(data.REM.deltaBandPower.meanf,data.REM.deltaBandPower.meanC + data.REM.deltaBandPower.stdC,'color',colorbrewer_setA_colorC,'LineWidth',1)
semilogx(data.REM.deltaBandPower.meanf,data.REM.deltaBandPower.meanC - data.REM.deltaBandPower.stdC,'color',colorbrewer_setA_colorC,'LineWidth',1)
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorA,'LineWidth',1);
semilogx(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorB,'LineWidth',1);
semilogx(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorC,'LineWidth',1);
title('Delta-band [1-4 Hz]')
ylabel('Coherence')
xlabel('Frequency (Hz)')
axis square
ylim([0 1])
xlim([0.05 1])
set(gca,'box','off')

%% Theta-band power
subplot(2,3,3)
L1 = semilogx(data.Rest.thetaBandPower.meanf,data.Rest.thetaBandPower.meanC,'color',colorbrewer_setA_colorA,'LineWidth',3);
hold on
semilogx(data.Rest.thetaBandPower.meanf,data.Rest.thetaBandPower.meanC + data.Rest.thetaBandPower.stdC,'color',colorbrewer_setA_colorA,'LineWidth',1)
semilogx(data.Rest.thetaBandPower.meanf,data.Rest.thetaBandPower.meanC - data.Rest.thetaBandPower.stdC,'color',colorbrewer_setA_colorA,'LineWidth',1)
L2 = semilogx(data.NREM.thetaBandPower.meanf,data.NREM.thetaBandPower.meanC,'color',colorbrewer_setA_colorB,'LineWidth',3);
semilogx(data.NREM.thetaBandPower.meanf,data.NREM.thetaBandPower.meanC + data.NREM.thetaBandPower.stdC,'color',colorbrewer_setA_colorB,'LineWidth',1)
semilogx(data.NREM.thetaBandPower.meanf,data.NREM.thetaBandPower.meanC - data.NREM.thetaBandPower.stdC,'color',colorbrewer_setA_colorB,'LineWidth',1)
L3 = semilogx(data.REM.thetaBandPower.meanf,data.REM.thetaBandPower.meanC,'color',colorbrewer_setA_colorC,'LineWidth',3);
semilogx(data.REM.thetaBandPower.meanf,data.REM.thetaBandPower.meanC + data.REM.thetaBandPower.stdC,'color',colorbrewer_setA_colorC,'LineWidth',1)
semilogx(data.REM.thetaBandPower.meanf,data.REM.thetaBandPower.meanC - data.REM.thetaBandPower.stdC,'color',colorbrewer_setA_colorC,'LineWidth',1)
L4 = semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorA,'LineWidth',1);
L5 = semilogx(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorB,'LineWidth',1);
L6 = semilogx(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorC,'LineWidth',1);
title('Theta-band [4-10 Hz]')
ylabel('Coherence')
xlabel('Frequency (Hz)')
legend([L1,L2,L3,L4,L5,L6],'Rest','NREM','REM','Rest 95% conf','NREM 95% conf','REM 95% conf')
axis square
ylim([0 1])
xlim([0.05 1])
set(gca,'box','off')

%% Alpha-band power
subplot(2,3,4);
semilogx(data.Rest.alphaBandPower.meanf,data.Rest.alphaBandPower.meanC,'color',colorbrewer_setA_colorA,'LineWidth',3);
hold on
semilogx(data.Rest.alphaBandPower.meanf,data.Rest.alphaBandPower.meanC + data.Rest.alphaBandPower.stdC,'color',colorbrewer_setA_colorA,'LineWidth',1)
semilogx(data.Rest.alphaBandPower.meanf,data.Rest.alphaBandPower.meanC - data.Rest.alphaBandPower.stdC,'color',colorbrewer_setA_colorA,'LineWidth',1)
semilogx(data.NREM.alphaBandPower.meanf,data.NREM.alphaBandPower.meanC,'color',colorbrewer_setA_colorB,'LineWidth',3);
semilogx(data.NREM.alphaBandPower.meanf,data.NREM.alphaBandPower.meanC + data.NREM.alphaBandPower.stdC,'color',colorbrewer_setA_colorB,'LineWidth',1)
semilogx(data.NREM.alphaBandPower.meanf,data.NREM.alphaBandPower.meanC - data.NREM.alphaBandPower.stdC,'color',colorbrewer_setA_colorB,'LineWidth',1)
semilogx(data.REM.alphaBandPower.meanf,data.REM.alphaBandPower.meanC,'color',colorbrewer_setA_colorC,'LineWidth',3);
semilogx(data.REM.alphaBandPower.meanf,data.REM.alphaBandPower.meanC + data.REM.alphaBandPower.stdC,'color',colorbrewer_setA_colorC,'LineWidth',1)
semilogx(data.REM.alphaBandPower.meanf,data.REM.alphaBandPower.meanC - data.REM.alphaBandPower.stdC,'color',colorbrewer_setA_colorC,'LineWidth',1)
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorA,'LineWidth',1);
semilogx(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorB,'LineWidth',1);
semilogx(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorC,'LineWidth',1);
title('Alpha-band [10-13 Hz]')
ylabel('Coherence')
xlabel('Frequency (Hz)')
axis square
ylim([0 1])
xlim([0.05 1])
set(gca,'box','off')

%% Beta-band power
subplot(2,3,5);
semilogx(data.Rest.betaBandPower.meanf,data.Rest.betaBandPower.meanC,'color',colorbrewer_setA_colorA,'LineWidth',3);
hold on
semilogx(data.Rest.betaBandPower.meanf,data.Rest.betaBandPower.meanC + data.Rest.betaBandPower.stdC,'color',colorbrewer_setA_colorA,'LineWidth',1)
semilogx(data.Rest.betaBandPower.meanf,data.Rest.betaBandPower.meanC - data.Rest.betaBandPower.stdC,'color',colorbrewer_setA_colorA,'LineWidth',1)
semilogx(data.NREM.betaBandPower.meanf,data.NREM.betaBandPower.meanC,'color',colorbrewer_setA_colorB,'LineWidth',3);
semilogx(data.NREM.betaBandPower.meanf,data.NREM.betaBandPower.meanC + data.NREM.betaBandPower.stdC,'color',colorbrewer_setA_colorB,'LineWidth',1)
semilogx(data.NREM.betaBandPower.meanf,data.NREM.betaBandPower.meanC - data.NREM.betaBandPower.stdC,'color',colorbrewer_setA_colorB,'LineWidth',1)
semilogx(data.REM.betaBandPower.meanf,data.REM.betaBandPower.meanC,'color',colorbrewer_setA_colorC,'LineWidth',3);
semilogx(data.REM.betaBandPower.meanf,data.REM.betaBandPower.meanC + data.REM.betaBandPower.stdC,'color',colorbrewer_setA_colorC,'LineWidth',1)
semilogx(data.REM.betaBandPower.meanf,data.REM.betaBandPower.meanC - data.REM.betaBandPower.stdC,'color',colorbrewer_setA_colorC,'LineWidth',1)
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorA,'LineWidth',1);
semilogx(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorB,'LineWidth',1);
semilogx(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorC,'LineWidth',1);
title('Beta-band [13-30 Hz]')
ylabel('Coherence')
xlabel('Frequency (Hz)')
axis square
ylim([0 1])
xlim([0.05 1])
set(gca,'box','off')

%% Gamma-band power
subplot(2,3,6);
semilogx(data.Rest.gammaBandPower.meanf,data.Rest.gammaBandPower.meanC,'color',colorbrewer_setA_colorA,'LineWidth',3);
hold on
semilogx(data.Rest.gammaBandPower.meanf,data.Rest.gammaBandPower.meanC + data.Rest.gammaBandPower.stdC,'color',colorbrewer_setA_colorA,'LineWidth',1)
semilogx(data.Rest.gammaBandPower.meanf,data.Rest.gammaBandPower.meanC - data.Rest.gammaBandPower.stdC,'color',colorbrewer_setA_colorA,'LineWidth',1)
semilogx(data.NREM.gammaBandPower.meanf,data.NREM.gammaBandPower.meanC,'color',colorbrewer_setA_colorB,'LineWidth',3);
semilogx(data.NREM.gammaBandPower.meanf,data.NREM.gammaBandPower.meanC + data.NREM.gammaBandPower.stdC,'color',colorbrewer_setA_colorB,'LineWidth',1)
semilogx(data.NREM.gammaBandPower.meanf,data.NREM.gammaBandPower.meanC - data.NREM.gammaBandPower.stdC,'color',colorbrewer_setA_colorB,'LineWidth',1)
semilogx(data.REM.gammaBandPower.meanf,data.REM.gammaBandPower.meanC,'color',colorbrewer_setA_colorC,'LineWidth',3);
semilogx(data.REM.gammaBandPower.meanf,data.REM.gammaBandPower.meanC + data.REM.gammaBandPower.stdC,'color',colorbrewer_setA_colorC,'LineWidth',1)
semilogx(data.REM.gammaBandPower.meanf,data.REM.gammaBandPower.meanC - data.REM.gammaBandPower.stdC,'color',colorbrewer_setA_colorC,'LineWidth',1)
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorA,'LineWidth',1);
semilogx(data.NREM.CBV_HbT.meanf,data.NREM.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorB,'LineWidth',1);
semilogx(data.REM.CBV_HbT.meanf,data.REM.CBV_HbT.maxConfC_Y,'--','color',colorbrewer_setA_colorC,'LineWidth',1);
title('Gamma-band [30-100 Hz]')
ylabel('Coherence')
xlabel('Frequency (Hz)')
axis square
ylim([0 1])
xlim([0.05 1])
set(gca,'box','off')

% save figure(s)
dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\';
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Coherence']);
