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
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
baselineTypes = {'manualSelection','setDuration','entireDuration'};
dataTypes = {'LH','RH'};

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(dataTypes)
        dataType = dataTypes{1,b};
        for c = 1:length(baselineTypes)
            baselineType = baselineTypes{1,c};
            for d = 1:length(solenoidNames)
                solenoidName = solenoidNames{1,d};
                data.(dataType).(baselineType).(solenoidName).CBV(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).CBV.Refl;
                data.(dataType).(baselineType).(solenoidName).HbT(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).CBV.HbT;
                data.(dataType).(baselineType).(solenoidName).cortMUA(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).MUA.corticalData;
                data.(dataType).(baselineType).(solenoidName).hipMUA(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).MUA.hippocampalData;
                data.(dataType).(baselineType).(solenoidName).timeVector(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).timeVector;
                data.(dataType).(baselineType).(solenoidName).cortS(:,:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).LFP.corticalS;
                data.(dataType).(baselineType).(solenoidName).hipS(:,:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).LFP.hippocampalS;
                data.(dataType).(baselineType).(solenoidName).T(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).LFP.T;
                data.(dataType).(baselineType).(solenoidName).F(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).LFP.F;
            end
        end
    end
end

% concatenate the data from the contra and ipsi data
for e = 1:length(baselineTypes)
    baselineType = baselineTypes{1,e};
    data.Contra.(baselineType).CBV = cat(2,data.LH.(baselineType).RPadSol.CBV,data.RH.(baselineType).LPadSol.CBV);
    data.Contra.(baselineType).HbT = cat(2,data.LH.(baselineType).RPadSol.HbT,data.RH.(baselineType).LPadSol.HbT);
    data.Contra.(baselineType).cortMUA = cat(2,data.LH.(baselineType).RPadSol.cortMUA,data.RH.(baselineType).LPadSol.cortMUA);
    data.Contra.(baselineType).hipMUA = data.RH.(baselineType).RPadSol.hipMUA;
    data.Contra.(baselineType).timeVector = cat(2,data.LH.(baselineType).RPadSol.timeVector,data.RH.(baselineType).LPadSol.timeVector);
    data.Contra.(baselineType).cortS = cat(3,data.LH.(baselineType).RPadSol.cortS,data.RH.(baselineType).LPadSol.cortS);
    data.Contra.(baselineType).hipS = data.RH.(baselineType).RPadSol.hipS;
    data.Contra.(baselineType).T = cat(2,data.LH.(baselineType).RPadSol.T,data.RH.(baselineType).LPadSol.T);
    data.Contra.(baselineType).F = cat(2,data.LH.(baselineType).RPadSol.F,data.RH.(baselineType).LPadSol.F);

    data.Ipsi.(baselineType).CBV = cat(2,data.LH.(baselineType).LPadSol.CBV,data.RH.(baselineType).RPadSol.CBV);
    data.Ipsi.(baselineType).HbT = cat(2,data.LH.(baselineType).LPadSol.HbT,data.RH.(baselineType).RPadSol.HbT);
    data.Ipsi.(baselineType).cortMUA = cat(2,data.LH.(baselineType).LPadSol.cortMUA,data.RH.(baselineType).RPadSol.cortMUA);
    data.Ipsi.(baselineType).hipMUA = data.RH.(baselineType).LPadSol.hipMUA;
    data.Ipsi.(baselineType).timeVector = cat(2,data.LH.(baselineType).LPadSol.timeVector,data.RH.(baselineType).RPadSol.timeVector);
    data.Ipsi.(baselineType).cortS = cat(3,data.LH.(baselineType).LPadSol.cortS,data.RH.(baselineType).RPadSol.cortS);
    data.Ipsi.(baselineType).hipS = data.RH.(baselineType).LPadSol.hipS;
    data.Ipsi.(baselineType).T = cat(2,data.LH.(baselineType).LPadSol.T,data.RH.(baselineType).RPadSol.T);
    data.Ipsi.(baselineType).F = cat(2,data.LH.(baselineType).LPadSol.F,data.RH.(baselineType).RPadSol.F);
    
    data.Auditory.(baselineType).CBV = cat(2,data.LH.(baselineType).AudSol.CBV,data.RH.(baselineType).AudSol.CBV);
    data.Auditory.(baselineType).HbT = cat(2,data.LH.(baselineType).AudSol.HbT,data.RH.(baselineType).AudSol.HbT);
    data.Auditory.(baselineType).cortMUA = cat(2,data.LH.(baselineType).AudSol.cortMUA,data.RH.(baselineType).AudSol.cortMUA);
    data.Auditory.(baselineType).hipMUA = data.RH.(baselineType).AudSol.hipMUA;
    data.Auditory.(baselineType).timeVector = cat(2,data.LH.(baselineType).AudSol.timeVector,data.RH.(baselineType).AudSol.timeVector);
    data.Auditory.(baselineType).cortS = cat(3,data.LH.(baselineType).AudSol.cortS,data.RH.(baselineType).AudSol.cortS);
    data.Auditory.(baselineType).hipS = data.RH.(baselineType).AudSol.hipS;
    data.Auditory.(baselineType).T = cat(2,data.LH.(baselineType).AudSol.T,data.RH.(baselineType).AudSol.T);
    data.Auditory.(baselineType).F = cat(2,data.LH.(baselineType).AudSol.F,data.RH.(baselineType).AudSol.F);
end

% take the averages of each field through the proper dimension
for f = 1:length(compDataTypes)
    compDataType = compDataTypes{1,f};
    for g = 1:length(baselineTypes)
        baselineType = baselineTypes{1,g};
        data.(compDataType).(baselineType).mean_CBV = mean(data.(compDataType).(baselineType).CBV,2);
        data.(compDataType).(baselineType).std_CBV = std(data.(compDataType).(baselineType).CBV,0,2);
        data.(compDataType).(baselineType).mean_HbT = mean(data.(compDataType).(baselineType).HbT,2);
        data.(compDataType).(baselineType).std_HbT = std(data.(compDataType).(baselineType).HbT,0,2);
        data.(compDataType).(baselineType).mean_CortMUA = mean(data.(compDataType).(baselineType).cortMUA,2);
        data.(compDataType).(baselineType).std_CortMUA = std(data.(compDataType).(baselineType).cortMUA,0,2);
        data.(compDataType).(baselineType).mean_HipMUA = mean(data.(compDataType).(baselineType).hipMUA,2);
        data.(compDataType).(baselineType).std_HipMUA = std(data.(compDataType).(baselineType).hipMUA,0,2);
        data.(compDataType).(baselineType).mean_timeVector = mean(data.(compDataType).(baselineType).timeVector,2);
        data.(compDataType).(baselineType).mean_CortS = mean(data.(compDataType).(baselineType).cortS,3);
        data.(compDataType).(baselineType).mean_HipS = mean(data.(compDataType).(baselineType).hipS,3);
        data.(compDataType).(baselineType).mean_T = mean(data.(compDataType).(baselineType).T,2);
        data.(compDataType).(baselineType).mean_F = mean(data.(compDataType).(baselineType).F,2);
    end
end

%% summary figure(s)
for h = 1:length(baselineTypes)
    baselineType = baselineTypes{1,h};
    summaryFigure = figure;
    sgtitle({['Cortical stimulus-evoked averages - ' baselineType],' '})
    
    %% Cortical MUA Contra Stim
    ax1 = subplot(4,3,1);
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_CortMUA,'k')
    hold on
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_CortMUA + data.Contra.(baselineType).std_CortMUA,'color',colors_IOS('battleship grey'))
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_CortMUA - data.Contra.(baselineType).std_CortMUA,'color',colors_IOS('battleship grey'))
    title('Contralateral stim cortical MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% Cortical MUA Ispi Stim
    ax2 = subplot(4,3,2);
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_CortMUA,'k')
    hold on
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_CortMUA + data.Ipsi.(baselineType).std_CortMUA,'color',colors_IOS('battleship grey'))
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_CortMUA - data.Ipsi.(baselineType).std_CortMUA,'color',colors_IOS('battleship grey'))
    title('Ipsilateral stim cortical MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% Cortical MUA Auditory Stim
    ax3 = subplot(4,3,3);
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_CortMUA,'k')
    hold on
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_CortMUA + data.Auditory.(baselineType).std_CortMUA,'color',colors_IOS('battleship grey'))
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_CortMUA - data.Auditory.(baselineType).std_CortMUA,'color',colors_IOS('battleship grey'))
    title('Auditory stim cortical MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% Cortical LFP Contra Stim
    ax4 = subplot(4,3,4);
    imagesc(data.Contra.(baselineType).mean_T,data.Contra.(baselineType).mean_F,data.Contra.(baselineType).mean_CortS)
    title('Contralateral stim cortical LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    colorbar
    caxis([-0.5 1])
    set(gca,'Ticklength',[0 0])
    axis square
    axis xy
    
    %% Cortical LFP Ispi Stim
    ax5 = subplot(4,3,5);
    imagesc(data.Ipsi.(baselineType).mean_T,data.Ipsi.(baselineType).mean_F,data.Ipsi.(baselineType).mean_CortS)
    title('Ipsilateral stim cortical LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    colorbar
    caxis([-0.5 1])
    set(gca,'Ticklength',[0 0])
    axis square
    axis xy
    
    %% Cortical LFP Auditory Stim
    ax6 = subplot(4,3,6);
    imagesc(data.Auditory.(baselineType).mean_T,data.Auditory.(baselineType).mean_F,data.Auditory.(baselineType).mean_CortS)
    title('Auditory stim cortical LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    colorbar
    caxis([-0.5 1])
    set(gca,'Ticklength',[0 0])
    axis square
    axis xy
    
    %% CBV HbT Contra Stim
    ax7 = subplot(4,3,7);
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_HbT,'k')
    hold on
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_HbT + data.Contra.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_HbT - data.Contra.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    title('Contralateral stim HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% CBV HbT Ispi Stim
    ax8 = subplot(4,3,8);
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_HbT,'k')
    hold on
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_HbT + data.Ipsi.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_HbT - data.Ipsi.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    title('Ipsilateral stim HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% CBV HbT Auditory Stim
    ax9 = subplot(4,3,9);
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_HbT,'k')
    hold on
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_HbT + data.Auditory.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_HbT - data.Auditory.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    title('Auditory stim HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% CBV Refl Contra Stim
    ax10 = subplot(4,3,10);
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_CBV,'k')
    hold on
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_CBV + data.Contra.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_CBV - data.Contra.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    title('Contralateral stim reflectance')
    ylabel('\DeltaR/R (%)')
    xlabel('Time (sec)')
    axis square
    
    %% CBV Refl Ispi Stim
    ax11 = subplot(4,3,11);
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_CBV,'k')
    hold on
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_CBV + data.Ipsi.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_CBV - data.Ipsi.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    title('Ipsilateral stim reflectance')
    ylabel('\DeltaR/R (%)')
    xlabel('Time (sec)')
    axis square
    
    %% CBV Refl Auditory Stim
    ax12 = subplot(4,3,12);
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_CBV,'k')
    hold on
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_CBV + data.Auditory.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_CBV - data.Auditory.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    title('Auditory stim feflectance')
    ylabel('\DeltaR/R (%)')
    xlabel('Time (sec)')
    axis square
    
    linkaxes([ax1 ax2 ax3],'xy')
    linkaxes([ax7 ax8 ax9],'xy')
    linkaxes([ax10 ax11 ax12],'xy')
    
    % save figure(s)
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Stimulus-evoked Responses\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure, [dirpath baselineType '_AverageCorticalStimResponses']);
end

%% summary figure(s)
for h = 1:length(baselineTypes)
    baselineType = baselineTypes{1,h};
    summaryFigure = figure;
    sgtitle({['Hippocampal stimulus-evoked averages - ' baselineType],' '})
    
    %% Hippocampal MUA Contra Stim
    ax1 = subplot(4,3,1);
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_HipMUA,'k')
    hold on
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_HipMUA + data.Contra.(baselineType).std_HipMUA,'color',colors_IOS('battleship grey'))
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_HipMUA - data.Contra.(baselineType).std_HipMUA,'color',colors_IOS('battleship grey'))
    title('Contralateral stim hippocampal MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% Hippocampal MUA Ispi Stim
    ax2 = subplot(4,3,2);
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_HipMUA,'k')
    hold on
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_HipMUA + data.Ipsi.(baselineType).std_HipMUA,'color',colors_IOS('battleship grey'))
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_HipMUA - data.Ipsi.(baselineType).std_HipMUA,'color',colors_IOS('battleship grey'))
    title('Ipsilateral stim hippocampal MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% Hippocampal MUA Auditory Stim
    ax3 = subplot(4,3,3);
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_HipMUA,'k')
    hold on
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_HipMUA + data.Auditory.(baselineType).std_HipMUA,'color',colors_IOS('battleship grey'))
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_HipMUA - data.Auditory.(baselineType).std_HipMUA,'color',colors_IOS('battleship grey'))
    title('Auditory stim hippocampal MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% Hippocampal LFP Contra Stim
    ax4 = subplot(4,3,4);
    imagesc(data.Contra.(baselineType).mean_T,data.Contra.(baselineType).mean_F,data.Contra.(baselineType).mean_HipS)
    title('Contralateral stim hippocampal LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    colorbar
    caxis([-0.5 1])
    set(gca,'Ticklength',[0 0])
    axis square
    axis xy
    
    %% Hippocampal LFP Ispi Stim
    ax5 = subplot(4,3,5);
    imagesc(data.Ipsi.(baselineType).mean_T,data.Ipsi.(baselineType).mean_F,data.Ipsi.(baselineType).mean_HipS)
    title('Ipsilateral stim hippocampal LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    colorbar
    caxis([-0.5 1])
    set(gca,'Ticklength',[0 0])
    axis square
    axis xy
    
    %% Hippocampal LFP Auditory Stim
    ax6 = subplot(4,3,6);
    imagesc(data.Auditory.(baselineType).mean_T,data.Auditory.(baselineType).mean_F,data.Auditory.(baselineType).mean_HipS)
    title('Auditory stim hippocampal LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    colorbar
    caxis([-0.5 1])
    set(gca,'Ticklength',[0 0])
    axis square
    axis xy
    
    %% CBV HbT Contra Stim
    ax7 = subplot(4,3,7);
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_HbT,'k')
    hold on
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_HbT + data.Contra.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_HbT - data.Contra.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    title('Contralateral stim HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% CBV HbT Ispi Stim
    ax8 = subplot(4,3,8);
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_HbT,'k')
    hold on
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_HbT + data.Ipsi.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_HbT - data.Ipsi.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    title('Ipsilateral stim HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% CBV HbT Auditory Stim
    ax9 = subplot(4,3,9);
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_HbT,'k')
    hold on
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_HbT + data.Auditory.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_HbT - data.Auditory.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    title('Auditory stim HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% CBV Refl Contra Stim
    ax10 = subplot(4,3,10);
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_CBV,'k')
    hold on
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_CBV + data.Contra.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_CBV - data.Contra.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    title('Contralateral stim reflectance')
    ylabel('\DeltaR/R (%)')
    xlabel('Time (sec)')
    axis square
    
    %% CBV Refl Ispi Stim
    ax11 = subplot(4,3,11);
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_CBV,'k')
    hold on
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_CBV + data.Ipsi.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_CBV - data.Ipsi.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    title('Ipsilateral stim reflectance')
    ylabel('\DeltaR/R (%)')
    xlabel('Time (sec)')
    axis square
    
    %% CBV Refl Auditory Stim
    ax12 = subplot(4,3,12);
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_CBV,'k')
    hold on
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_CBV + data.Auditory.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_CBV - data.Auditory.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    title('Auditory stim reflectance')
    ylabel('\DeltaR/R (%)')
    xlabel('Time (sec)')
    axis square
    
    linkaxes([ax1 ax2 ax3],'xy')
    linkaxes([ax7 ax8 ax9],'xy')
    linkaxes([ax10 ax11 ax12],'xy')
    
    % save figure(s)
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Stimulus-evoked Responses\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure, [dirpath baselineType '_AverageHippocampalStimResponses']);
end
