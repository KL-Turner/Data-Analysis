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
            for d = 1:length(solenoids)
                solenoidName = solenoids{1,d};
                data.(dataType).(baselineType).(solenoidName).CBV(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).CBV.Refl;
                data.(dataType).(baselineType).(solenoidName).HbT(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).CBV.HbT;
                data.(dataType).(baselineType).(solenoidName).MUA(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).MUA.data;
                data.(dataType).(baselineType).(solenoidName).timeVector = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).timeVector;
                data.(dataType).(baselineType).(solenoidName).S(:,:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).LFP.S;
                data.(dataType).(baselineType).(solenoidName).T(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).LFP.T;
                data.(dataType).(baselineType).(solenoidName).F(:,a) = AnalysisResults.EvokedAvgs.Stim.(dataType).(baselineType).(solenoidName).LFP.F;
            end
        end
    end
end

% concatenate the data from the contra and ipsi data
for e = 1:length(baselineTypes)
    baselineType = baselineTypes{1,e};
    data.Ispi.(baselineType).CBV = cat(2,data.LH.(baselineType).RPadSol.CBV,data.RH.(baselineType).LPadSol.CBV);
    data.Ispi.(baselineType).HbT = cat(2,data.LH.(baselineType).RPadSol.HbT,data.RH.(baselineType).LPadSol.HbT);
    data.Ispi.(baselineType).MUA = cat(2,data.LH.(baselineType).RPadSol.MUA,data.RH.(baselineType).LPadSol.MUA);
    data.Ispi.(baselineType).timeVector = cat(2,data.LH.(baselineType).RPadSol.timeVector,data.RH.(baselineType).LPadSol.timeVector);
    data.Ispi.(baselineType).S = cat(3,data.LH.(baselineType).RPadSol.S,data.RH.(baselineType).LPadSol.S);
    data.Ispi.(baselineType).T = cat(2,data.LH.(baselineType).RPadSol.T,data.RH.(baselineType).LPadSol.T);
    data.Ispi.(baselineType).F = cat(2,data.LH.(baselineType).RPadSol.F,data.RH.(baselineType).LPadSol.F);

    data.Contra.(baselineType).CBV = cat(2,data.LH.(baselineType).LPadSol.CBV,data.RH.(baselineType).RPadSol.CBV);
    data.Contra.(baselineType).HbT = cat(2,data.LH.(baselineType).LPadSol.HbT,data.RH.(baselineType).RPadSol.HbT);
    data.Contra.(baselineType).MUA = cat(2,data.LH.(baselineType).LPadSol.MUA,data.RH.(baselineType).RPadSol.MUA);
    data.Contra.(baselineType).timeVector = cat(2,data.LH.(baselineType).LPadSol.timeVector,data.RH.(baselineType).RPadSol.timeVector);
    data.Contra.(baselineType).S = cat(3,data.LH.(baselineType).LPadSol.S,data.RH.(baselineType).RPadSol.S);
    data.Contra.(baselineType).T = cat(2,data.LH.(baselineType).LPadSol.T,data.RH.(baselineType).RPadSol.T);
    data.Contra.(baselineType).F = cat(2,data.LH.(baselineType).LPadSol.F,data.RH.(baselineType).RPadSol.F);
    
    data.Auditory.(baselineType).CBV = cat(2,data.LH.(baselineType).AudSol.CBV,data.RH.(baselineType).AudSol.CBV);
    data.Auditory.(baselineType).HbT = cat(2,data.LH.(baselineType).AudSol.HbT,data.RH.(baselineType).AudSol.HbT);
    data.Auditory.(baselineType).MUA = cat(2,data.LH.(baselineType).AudSol.MUA,data.RH.(baselineType).AudSol.MUA);
    data.Auditory.(baselineType).timeVector = cat(2,data.LH.(baselineType).AudSol.timeVector,data.RH.(baselineType).AudSol.timeVector);
    data.Auditory.(baselineType).S = cat(3,data.LH.(baselineType).AudSol.S,data.RH.(baselineType).AudSol.S);
    data.Auditory.(baselineType).T = cat(2,data.LH.(baselineType).AudSol.T,data.RH.(baselineType).AudSol.T);
    data.Auditory.(baselineType).F = cat(2,data.LH.(baselineType).AudSol.F,data.RH.(baselineType).AudSol.F);
end

% take the averages of each field through the proper dimension
for f = 1:length(compDataTypes)
    compDataType = compDataTypes{1,f};
    for g = 1:length(baselineDataTypes)
        baselineType = baselineTypes{1,g};
        data.(compDataType).(baselineType).mean_CBV = mean(data.(compDataType).(baselineType).CBV,2);
        data.(compDataType).(baselineType).std_CBV = std(data.(compDataType).(baselineType).CBV,0,2);
        data.(compDataType).(baselineType).mean_HbT = mean(data.(compDataType).(baselineType).HbT,2);
        data.(compDataType).(baselineType).std_HbT = std(data.(compDataType).(baselineType).HbT,0,2);
        data.(compDataType).(baselineType).mean_MUA = mean(data.(compDataType).(baselineType).MUA,2);
        data.(compDataType).(baselineType).std_MUA = std(data.(compDataType).(baselineType).MUA,0,2);
        data.(compDataType).(baselineType).mean_timeVector = mean(data.(compDataType).(baselineType).timeVector,2);
        data.(compDataType).(baselineType).mean_S = mean(data.(compDataType).(baselineType).S,3);
        data.(compDataType).(baselineType).mean_T = mean(data.(compDataType).(baselineType).T,2);
        data.(compDataType).(baselineType).mean_F = mean(data.(compDataType).(baselineType).F,2);
    end
end

%% summary figure(s)
for h = 1:length(baselineTypes)
    baselineType = baselineTypes{1,h};
    summaryFigure = figure;
    sgtitle({['Stimulus-evoked averages - ' baselineType],' '})
    
    %% MUA Contra Stim
    ax1 = subplot(4,3,1);
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_MUA,'k')
    hold on
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_MUA + data.Contra.(baselineType).std_MUA,'color',colors_IOS('battleship grey'))
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_MUA - data.Contra.(baselineType).std_MUA,'color',colors_IOS('battleship grey'))
    title('Contralateral Stim MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% MUA Ispi Stim
    ax2 = subplot(4,3,3);
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_MUA,'k')
    hold on
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_MUA + data.Ipsi.(baselineType).std_MUA,'color',colors_IOS('battleship grey'))
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_MUA - data.Ipsi.(baselineType).std_MUA,'color',colors_IOS('battleship grey'))
    title('Ipsilateral Stim MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% MUA Auditory Stim
    ax3 = subplot(4,3,3);
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_MUA,'k')
    hold on
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_MUA + data.Auditory.(baselineType).std_MUA,'color',colors_IOS('battleship grey'))
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_MUA - data.Auditory.(baselineType).std_MUA,'color',colors_IOS('battleship grey'))
    title('Auditory Stim MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% LFP Contra Stim
    ax4 = subplot(4,3,4);
    imagesc(data.Contra.(baselineType).mean_T,data.Contra.(baselineType).mean_F,data.Contra.(baselineType).mean_S)
    title('Contralateral Stim LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    caxis([-0.5 1])
    set(gca,'Ticklength',[0 0])
    axis square
    
    %% LFP Ispi Stim
    ax5 = subplot(4,3,5);
    imagesc(data.Ipsi.(baselineType).mean_T,data.Ipsi.(baselineType).mean_F,data.Ipsi.(baselineType).mean_S)
    title('Ipsilateral Stim LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    caxis([-0.5 1])
    set(gca,'Ticklength',[0 0])
    axis square
    
    %% LFP Auditory Stim
    ax6 = subplot(4,3,6);
    imagesc(data.Auditory.(baselineType).mean_T,data.Auditory.(baselineType).mean_F,data.Auditory.(baselineType).mean_S)
    title('Auditory Stim LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    caxis([-0.5 1])
    set(gca,'Ticklength',[0 0])
    axis square
    
    %% CBV HbT Contra Stim
    ax7 = subplot(4,3,7);
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_HbT,'k')
    hold on
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_HbT + data.Contra.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_HbT - data.Contra.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    title('Contralateral Stim HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% CBV HbT Ispi Stim
    ax8 = subplot(4,3,8);
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_HbT,'k')
    hold on
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_HbT + data.Ipsi.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_HbT - data.Ipsi.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    title('Ipsilateral Stim HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% CBV HbT Auditory Stim
    ax9 = subplot(4,3,9);
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_HbT,'k')
    hold on
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_HbT + data.Auditory.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_HbT - data.Auditory.(baselineType).std_HbT,'color',colors_IOS('battleship grey'))
    title('Auditory Stim HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% CBV Refl Contra Stim
    ax10 = subplot(4,3,10);
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_CBV,'k')
    hold on
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_CBV + data.Contra.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    plot(data.Contra.(baselineType).mean_timeVector,data.Contra.(baselineType).mean_CBV - data.Contra.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    title('Contralateral Stim Reflectance')
    ylabel('\DeltaR/R (%)')
    xlabel('Time (sec)')
    axis square
    
    %% CBV Refl Ispi Stim
    ax11 = subplot(4,3,11);
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_CBV,'k')
    hold on
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_CBV + data.Ipsi.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    plot(data.Ipsi.(baselineType).mean_timeVector,data.Ipsi.(baselineType).mean_CBV - data.Ipsi.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    title('Ipsilateral Stim Reflectance')
    ylabel('\DeltaR/R (%)')
    xlabel('Time (sec)')
    axis square
    
    %% CBV Refl Auditory Stim
    ax12 = subplot(4,3,12);
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_CBV,'k')
    hold on
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_CBV + data.Auditory.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    plot(data.Auditory.(baselineType).mean_timeVector,data.Auditory.(baselineType).mean_CBV - data.Auditory.(baselineType).std_CBV,'color',colors_IOS('battleship grey'))
    title('Auditory Stim Reflectance')
    ylabel('\DeltaR/R (%)')
    xlabel('Time (sec)')
    axis square
    
    linkaxes([ax1 ax2 ax3],'y')
    linkaxes([ax7 ax8 ax9],'y')
    linkaxes([ax10 ax11 ax12],'y')
    
    % save figure(s)
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Stimulus Responses\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure, [dirpath baselineType '_AverageStimResponse']);
end
