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
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
dataTypes = {'adjLH','adjRH'};

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\Turner_Manuscript_Summer2020\' animalID '\Bilateral Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(dataTypes)
        dataType = dataTypes{1,b};
        for d = 1:length(solenoidNames)
            solenoidName = solenoidNames{1,d};
            data.(dataType).(solenoidName).HbT(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoidName).CBV_HbT.HbT;
            data.(dataType).(solenoidName).CBV(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoidName).CBV.CBV;
            data.(dataType).(solenoidName).cortMUA(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoidName).MUA.corticalData;
            data.(dataType).(solenoidName).hipMUA(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoidName).MUA.hippocampalData;
            data.(dataType).(solenoidName).timeVector(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoidName).timeVector;
            data.(dataType).(solenoidName).cortS(:,:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoidName).LFP.corticalS;
            data.(dataType).(solenoidName).hipS(:,:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoidName).LFP.hippocampalS;
            data.(dataType).(solenoidName).T(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoidName).LFP.T;
            data.(dataType).(solenoidName).F(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(solenoidName).LFP.F;
        end
    end
end

% concatenate the data from the contra and ipsi data
data.Contra.HbT = cat(2,data.adjLH.RPadSol.HbT,data.adjRH.LPadSol.HbT);
data.Contra.CBV = cat(2,data.adjLH.RPadSol.CBV,data.adjRH.LPadSol.CBV);
data.Contra.cortMUA = cat(2,data.adjLH.RPadSol.cortMUA,data.adjRH.LPadSol.cortMUA);
data.Contra.hipMUA = data.adjRH.RPadSol.hipMUA;
data.Contra.timeVector = cat(2,data.adjLH.RPadSol.timeVector,data.adjRH.LPadSol.timeVector);
data.Contra.cortS = cat(3,data.adjLH.RPadSol.cortS,data.adjRH.LPadSol.cortS);
data.Contra.hipS = data.adjRH.RPadSol.hipS;
data.Contra.T = cat(2,data.adjLH.RPadSol.T,data.adjRH.LPadSol.T);
data.Contra.F = cat(2,data.adjLH.RPadSol.F,data.adjRH.LPadSol.F);

data.Ipsi.HbT = cat(2,data.adjLH.LPadSol.HbT,data.adjRH.RPadSol.HbT);
data.Ipsi.CBV = cat(2,data.adjLH.LPadSol.CBV,data.adjRH.RPadSol.CBV);
data.Ipsi.cortMUA = cat(2,data.adjLH.LPadSol.cortMUA,data.adjRH.RPadSol.cortMUA);
data.Ipsi.hipMUA = data.adjRH.LPadSol.hipMUA;
data.Ipsi.timeVector = cat(2,data.adjLH.LPadSol.timeVector,data.adjRH.RPadSol.timeVector);
data.Ipsi.cortS = cat(3,data.adjLH.LPadSol.cortS,data.adjRH.RPadSol.cortS);
data.Ipsi.hipS = data.adjRH.LPadSol.hipS;
data.Ipsi.T = cat(2,data.adjLH.LPadSol.T,data.adjRH.RPadSol.T);
data.Ipsi.F = cat(2,data.adjLH.LPadSol.F,data.adjRH.RPadSol.F);

data.Auditory.HbT = cat(2,data.adjLH.AudSol.HbT,data.adjRH.AudSol.HbT);
data.Auditory.CBV = cat(2,data.adjLH.AudSol.CBV,data.adjRH.AudSol.CBV);
data.Auditory.cortMUA = cat(2,data.adjLH.AudSol.cortMUA,data.adjRH.AudSol.cortMUA);
data.Auditory.hipMUA = data.adjRH.AudSol.hipMUA;
data.Auditory.timeVector = cat(2,data.adjLH.AudSol.timeVector,data.adjRH.AudSol.timeVector);
data.Auditory.cortS = cat(3,data.adjLH.AudSol.cortS,data.adjRH.AudSol.cortS);
data.Auditory.hipS = data.adjRH.AudSol.hipS;
data.Auditory.T = cat(2,data.adjLH.AudSol.T,data.adjRH.AudSol.T);
data.Auditory.F = cat(2,data.adjLH.AudSol.F,data.adjRH.AudSol.F);

% take the averages of each field through the proper dimension
for f = 1:length(compDataTypes)
    compDataType = compDataTypes{1,f};
    data.(compDataType).mean_HbT = mean(data.(compDataType).HbT,2);
    data.(compDataType).std_HbT = std(data.(compDataType).HbT,0,2);
    data.(compDataType).mean_CBV = mean(data.(compDataType).CBV,2);
    data.(compDataType).std_CBV = std(data.(compDataType).CBV,0,2);
    data.(compDataType).mean_CortMUA = mean(data.(compDataType).cortMUA,2);
    data.(compDataType).std_CortMUA = std(data.(compDataType).cortMUA,0,2);
    data.(compDataType).mean_HipMUA = mean(data.(compDataType).hipMUA,2);
    data.(compDataType).std_HipMUA = std(data.(compDataType).hipMUA,0,2);
    data.(compDataType).mean_timeVector = mean(data.(compDataType).timeVector,2);
    data.(compDataType).mean_CortS = mean(data.(compDataType).cortS,3).*100;
    data.(compDataType).mean_HipS = mean(data.(compDataType).hipS,3).*100;
    data.(compDataType).mean_T = mean(data.(compDataType).T,2);
    data.(compDataType).mean_F = mean(data.(compDataType).F,2);
end

%% summary figure(s)
summaryFigure = figure;
sgtitle('Stimulus-evoked averages')

%% Cortical MUA Contra Stim
ax1 = subplot(5,3,1);
plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA,'k')
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA + data.Contra.std_CortMUA,'color',colors_IOS('battleship grey'))
plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA - data.Contra.std_CortMUA,'color',colors_IOS('battleship grey'))
title('Contra stim cortical MUA')
ylabel('\DeltaP/P (%)')
axis square

%% Cortical MUA Ispi Stim
ax2 = subplot(5,3,2);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA,'k')
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA + data.Ipsi.std_CortMUA,'color',colors_IOS('battleship grey'))
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA - data.Ipsi.std_CortMUA,'color',colors_IOS('battleship grey'))
title('Ipsi stim cortical MUA')
ylabel('\DeltaP/P (%)')
axis square

%% Cortical MUA Auditory Stim
ax3 = subplot(5,3,3);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA,'k')
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA + data.Auditory.std_CortMUA,'color',colors_IOS('battleship grey'))
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA - data.Auditory.std_CortMUA,'color',colors_IOS('battleship grey'))
title('Aud stim cortical MUA')
ylabel('\DeltaP/P (%)')
axis square

%% Cortical LFP Contra Stim
ax4 = subplot(5,3,4);
imagesc(data.Contra.mean_T,data.Contra.mean_F,data.Contra.mean_CortS)
title('Contra stim cortical LFP')
ylabel('Frequency (Hz)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)')
caxis([-50 100])
set(gca,'Ticklength',[0 0])
axis square
axis xy

%% Cortical LFP Ispi Stim
ax5 = subplot(5,3,5);
imagesc(data.Ipsi.mean_T,data.Ipsi.mean_F,data.Ipsi.mean_CortS)
title('Ipsi stim cortical LFP')
ylabel('Frequency (Hz)')
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)')
caxis([-50 100]) 
set(gca,'Ticklength',[0 0])
axis square
axis xy

%% Cortical LFP Auditory Stim
ax6 = subplot(5,3,6);
imagesc(data.Auditory.mean_T,data.Auditory.mean_F,data.Auditory.mean_CortS)
title('Aud stim cortical LFP')
ylabel('Frequency (Hz)')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)')
caxis([-50 100]) 
set(gca,'Ticklength',[0 0])
axis square
axis xy

%% Hippocampal MUA Contra Stim
ax7 = subplot(5,3,7);
plot(data.Contra.mean_timeVector,data.Contra.mean_HipMUA,'k')
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_HipMUA + data.Contra.std_HipMUA,'color',colors_IOS('battleship grey'))
plot(data.Contra.mean_timeVector,data.Contra.mean_HipMUA - data.Contra.std_HipMUA,'color',colors_IOS('battleship grey'))
title('Contra stim hippocampal MUA')
ylabel('\DeltaP/P (%)')
axis square

%% Hippocampal MUA Ispi Stim
ax8 = subplot(5,3,8);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipMUA,'k')
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipMUA + data.Ipsi.std_HipMUA,'color',colors_IOS('battleship grey'))
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipMUA - data.Ipsi.std_HipMUA,'color',colors_IOS('battleship grey'))
title('Ipsi stim hippocampal MUA')
ylabel('\DeltaP/P (%)')
axis square

%% Hippocampal MUA Auditory Stim
ax9 = subplot(5,3,9);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipMUA,'k')
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipMUA + data.Auditory.std_HipMUA,'color',colors_IOS('battleship grey'))
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipMUA - data.Auditory.std_HipMUA,'color',colors_IOS('battleship grey'))
title('Aud stim hippocampal MUA')
ylabel('\DeltaP/P (%)')
axis square

%% Hippocampal LFP Contra Stim
ax10 = subplot(5,3,10);
imagesc(data.Contra.mean_T,data.Contra.mean_F,data.Contra.mean_HipS)
title('Contra stim hippocampal LFP')
ylabel('Frequency (Hz)')
c10 = colorbar;
ylabel(c10,'\DeltaP/P (%)')
caxis([-50 100]) 
axis square
axis xy

%% Hippocampal LFP Ispi Stim
ax11 = subplot(5,3,11);
imagesc(data.Ipsi.mean_T,data.Ipsi.mean_F,data.Ipsi.mean_HipS)
title('Ipsi stim hippocampal LFP')
ylabel('Frequency (Hz)')
c11 = colorbar;
ylabel(c11,'\DeltaP/P (%)')
caxis([-50 100]) 
set(gca,'Ticklength',[0 0])
axis square
axis xy

%% Hippocampal LFP Auditory Stim
ax12 = subplot(5,3,12);
imagesc(data.Auditory.mean_T,data.Auditory.mean_F,data.Auditory.mean_HipS)
title('Aud stim hippocampal LFP')
ylabel('Frequency (Hz)')
c12 = colorbar;
ylabel(c12,'\DeltaP/P (%)')
caxis([-50 100]) 
set(gca,'Ticklength',[0 0])
axis square
axis xy

%% CBV HbT Contra Stim
ax13 = subplot(5,6,25);
plot(data.Contra.mean_timeVector,data.Contra.mean_HbT,'k')
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_HbT + data.Contra.std_HbT,'color',colors_IOS('battleship grey'))
plot(data.Contra.mean_timeVector,data.Contra.mean_HbT - data.Contra.std_HbT,'color',colors_IOS('battleship grey'))
title('Contra stim \DeltaHbT (\muM)')
ylabel('\DeltaHbT (\muM)')
xlabel('Peristimuls time (s)') 
axis square

%% CBV Refl Contra Stim
ax14 = subplot(5,6,26);
plot(data.Contra.mean_timeVector,data.Contra.mean_CBV,'k')
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_CBV + data.Contra.std_CBV,'color',colors_IOS('battleship grey'))
plot(data.Contra.mean_timeVector,data.Contra.mean_CBV - data.Contra.std_CBV,'color',colors_IOS('battleship grey'))
title('Contra stim reflectance')
ylabel('\DeltaR/R (%)')
xlabel('Peristimuls time (s)') 
axis square

%% CBV HbT Ispi Stim
ax15 = subplot(5,6,27);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HbT,'k')
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HbT + data.Ipsi.std_HbT,'color',colors_IOS('battleship grey'))
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HbT - data.Ipsi.std_HbT,'color',colors_IOS('battleship grey'))
title('Ipsi stim \DeltaHbT (\muM)')
ylabel('\DeltaHbT (\muM)')
xlabel('Peristimuls time (s)') 
axis square

%% CBV Refl Ispi Stim
ax16 = subplot(5,6,28);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CBV,'k')
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CBV + data.Ipsi.std_CBV,'color',colors_IOS('battleship grey'))
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CBV - data.Ipsi.std_CBV,'color',colors_IOS('battleship grey'))
title('Ipsi stim reflectance')
ylabel('\DeltaR/R (%)')
xlabel('Peristimuls time (s)') 
axis square

%% CBV HbT Auditory Stim
ax17 = subplot(5,6,29);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HbT,'k')
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HbT + data.Auditory.std_HbT,'color',colors_IOS('battleship grey'))
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HbT - data.Auditory.std_HbT,'color',colors_IOS('battleship grey'))
title('Aud stim \DeltaHbT (\muM)')
ylabel('\DeltaHbT (\muM)')
xlabel('Peristimuls time (s)') 
axis square

%% CBV Refl Auditory Stim
ax18 = subplot(5,6,30);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CBV,'k')
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CBV + data.Auditory.std_CBV,'color',colors_IOS('battleship grey'))
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CBV - data.Auditory.std_CBV,'color',colors_IOS('battleship grey'))
title('Aud stim reflectance')
ylabel('\DeltaR/R (%)')
xlabel('Peristimuls time (s)') 
axis square

linkaxes([ax1 ax2 ax3 ax7 ax8 ax9],'xy')
linkaxes([ax4 ax5 ax6 ax10 ax11 ax12],'xy')
linkaxes([ax13 ax15 ax17],'xy')
linkaxes([ax14 ax16 ax18],'xy')

ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
ax12Pos = get(ax12,'position');
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax2Pos(3:4);
ax6Pos(3:4) = ax3Pos(3:4);
ax10Pos(3:4) = ax1Pos(3:4);
ax11Pos(3:4) = ax2Pos(3:4);
ax12Pos(3:4) = ax3Pos(3:4);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax10,'position',ax10Pos);
set(ax11,'position',ax11Pos);
set(ax12,'position',ax12Pos);

% save figure(s)
dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\';
if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Stim Responses']);
