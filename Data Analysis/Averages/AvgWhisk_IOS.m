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
whiskDataTypes = {'Short','Intermediate','Long'};
criteriaDataTypes = {'whiskCriteriaA','whiskCriteriaB','whiskCriteriaC'};
baselineTypes = {'manualSelection','setDuration','entireDuration'};
dataTypes = {'LH','RH'};

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(baselineTypes)
        baselineType = baselineTypes{1,b};
        for c = 1:length(whiskDataTypes)
            whiskDataType = whiskDataTypes{1,c};
            criteriaDataType = criteriaDataTypes{1,c};
            data.(baselineType).(whiskDataType).LH.CBV(:,a) = AnalysisResults.EvokedAvgs.Whisk.LH.(baselineType).(criteriaDataType).CBV.Refl;
            data.(baselineType).(whiskDataType).LH.HbT(:,a) = AnalysisResults.EvokedAvgs.Whisk.LH.(baselineType).(criteriaDataType).CBV.HbT;
            data.(baselineType).(whiskDataType).LH.cortMUA(:,a) = AnalysisResults.EvokedAvgs.Whisk.LH.(baselineType).(criteriaDataType).MUA.corticalData;
            data.(baselineType).(whiskDataType).LH.cortS(:,:,a) = AnalysisResults.EvokedAvgs.Whisk.LH.(baselineType).(criteriaDataType).LFP.corticalS;
            data.(baselineType).(whiskDataType).LH.cortT(:,a) = AnalysisResults.EvokedAvgs.Whisk.LH.(baselineType).(criteriaDataType).LFP.T;
            data.(baselineType).(whiskDataType).LH.cortF(:,a) = AnalysisResults.EvokedAvgs.Whisk.LH.(baselineType).(criteriaDataType).LFP.F;
            
            data.(baselineType).(whiskDataType).RH.CBV(:,a) = AnalysisResults.EvokedAvgs.Whisk.RH.(baselineType).(criteriaDataType).CBV.Refl;
            data.(baselineType).(whiskDataType).RH.HbT(:,a) = AnalysisResults.EvokedAvgs.Whisk.RH.(baselineType).(criteriaDataType).CBV.HbT;
            data.(baselineType).(whiskDataType).RH.cortMUA(:,a) = AnalysisResults.EvokedAvgs.Whisk.RH.(baselineType).(criteriaDataType).MUA.corticalData;
            data.(baselineType).(whiskDataType).RH.cortS(:,:,a) = AnalysisResults.EvokedAvgs.Whisk.RH.(baselineType).(criteriaDataType).LFP.corticalS;
            data.(baselineType).(whiskDataType).RH.cortT(:,a) = AnalysisResults.EvokedAvgs.Whisk.RH.(baselineType).(criteriaDataType).LFP.T;
            data.(baselineType).(whiskDataType).RH.cortF(:,a) = AnalysisResults.EvokedAvgs.Whisk.RH.(baselineType).(criteriaDataType).LFP.F;
            
            data.(baselineType).(whiskDataType).Hip.hipMUA(:,a) = AnalysisResults.EvokedAvgs.Whisk.LH.(baselineType).(criteriaDataType).MUA.hippocampalData;
            data.(baselineType).(whiskDataType).Hip.hipS(:,:,a) = AnalysisResults.EvokedAvgs.Whisk.LH.(baselineType).(criteriaDataType).LFP.hippocampalS;
            data.(baselineType).(whiskDataType).Hip.hipT(:,a) = AnalysisResults.EvokedAvgs.Whisk.LH.(baselineType).(criteriaDataType).LFP.T;
            data.(baselineType).(whiskDataType).Hip.hipF(:,a) = AnalysisResults.EvokedAvgs.Whisk.LH.(baselineType).(criteriaDataType).LFP.F;
            
            data.(baselineType).(whiskDataType).timeVector(:,a) = AnalysisResults.EvokedAvgs.Whisk.LH.(baselineType).(criteriaDataType).timeVector;
        end
    end
end

% concatenate the data from the contra and ipsi data
for d = 1:length(baselineTypes)
    baselineType = baselineTypes{1,d};
    for e = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,e};
        data.(baselineType).(whiskDataType).CBV = cat(2,data.(baselineType).(whiskDataType).LH.CBV,data.(baselineType).(whiskDataType).RH.CBV);
        data.(baselineType).(whiskDataType).HbT = cat(2,data.(baselineType).(whiskDataType).LH.HbT,data.(baselineType).(whiskDataType).RH.HbT);
        data.(baselineType).(whiskDataType).cortMUA = cat(2,data.(baselineType).(whiskDataType).LH.cortMUA,data.(baselineType).(whiskDataType).RH.cortMUA);
        data.(baselineType).(whiskDataType).cortS = cat(3,data.(baselineType).(whiskDataType).LH.cortS,data.(baselineType).(whiskDataType).RH.cortS);
        data.(baselineType).(whiskDataType).cortT = cat(2,data.(baselineType).(whiskDataType).LH.cortT,data.(baselineType).(whiskDataType).RH.cortT);
        data.(baselineType).(whiskDataType).cortF = cat(2,data.(baselineType).(whiskDataType).LH.cortF,data.(baselineType).(whiskDataType).RH.cortF);
    end
end

% concatenate the data from the contra and ipsi data
for d = 1:length(baselineTypes)
    baselineType = baselineTypes{1,d};
    for e = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,e};
        data.(baselineType).(whiskDataType).meanCBV = mean(data.(baselineType).(whiskDataType).CBV,2);
        data.(baselineType).(whiskDataType).stdCBV = std(data.(baselineType).(whiskDataType).CBV,0,2);
        data.(baselineType).(whiskDataType).meanHbT = mean(data.(baselineType).(whiskDataType).HbT,2);
        data.(baselineType).(whiskDataType).stdHbT = std(data.(baselineType).(whiskDataType).HbT,0,2);
        data.(baselineType).(whiskDataType).meanCortMUA = mean(data.(baselineType).(whiskDataType).cortMUA,2);
        data.(baselineType).(whiskDataType).stdCortMUA = std(data.(baselineType).(whiskDataType).cortMUA,0,2);
        data.(baselineType).(whiskDataType).meanCortS = mean(data.(baselineType).(whiskDataType).cortS,3);
        data.(baselineType).(whiskDataType).meanCortT = mean(data.(baselineType).(whiskDataType).cortT,2);
        data.(baselineType).(whiskDataType).meanCortF = mean(data.(baselineType).(whiskDataType).cortF,2);
        
        data.(baselineType).(whiskDataType).meanHipMUA = mean(data.(baselineType).(whiskDataType).Hip.hipMUA,2);
        data.(baselineType).(whiskDataType).stdHipMUA = std(data.(baselineType).(whiskDataType).Hip.hipMUA,0,2);
        data.(baselineType).(whiskDataType).meanHipS = mean(data.(baselineType).(whiskDataType).Hip.hipS,3);
        data.(baselineType).(whiskDataType).meanHipT = mean(data.(baselineType).(whiskDataType).Hip.hipT,2);
        data.(baselineType).(whiskDataType).meanHipF = mean(data.(baselineType).(whiskDataType).Hip.hipF,2);
        
        data.(baselineType).(whiskDataType).meanTimeVector = mean(data.(baselineType).(whiskDataType).timeVector(:,a),2);
    end
end

%% summary figure(s)
for f = 1:length(baselineTypes)
    baselineType = baselineTypes{1,f};
    summaryFigure = figure;
    sgtitle({['Cortical whisking-evoked averages - ' baselineType],' '})
    
    %% Short whisks cortical MUA
    ax1 = subplot(4,3,1);
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCortMUA,'k');
    hold on
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCortMUA + data.(baselineType).Short.stdCortMUA,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCortMUA - data.(baselineType).Short.stdCortMUA,'color',colors_IOS('battleship grey'))
    title('Short whisking cortical MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
     %% Intermediate whisks cortical MUA
    ax2 = subplot(4,3,2);
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCortMUA,'k');
    hold on
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCortMUA + data.(baselineType).Intermediate.stdCortMUA,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCortMUA - data.(baselineType).Intermediate.stdCortMUA,'color',colors_IOS('battleship grey'))
    title('Short whisking cortical MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% Long whisks cortical MUA
    ax3 = subplot(4,3,3);
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCortMUA,'k');
    hold on
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCortMUA + data.(baselineType).Long.stdCortMUA,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCortMUA - data.(baselineType).Long.stdCortMUA,'color',colors_IOS('battleship grey'))
    title('Short whisking cortical MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% Short whisks cortical LFP 
    ax4 = subplot(4,3,4);
    imagesc(data.(baselineType).Short.meanCortT,data.(baselineType).Short.meanCortF,data.(baselineType).Short.meanCortS)
    title('Short whisking cortical LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    colorbar
    caxis([-0.25 0.5])
    set(gca,'Ticklength',[0 0])
    axis square
    axis xy
    
    %% Intermediate whisks cortical LFP 
    ax5 = subplot(4,3,5);
    imagesc(data.(baselineType).Intermediate.meanCortT,data.(baselineType).Intermediate.meanCortF,data.(baselineType).Intermediate.meanCortS)
    title('Intermediate whisking cortical LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    colorbar
    caxis([-0.25 0.5])
    set(gca,'Ticklength',[0 0])
    axis square
    axis xy
    
    %% Long whisks cortical LFP 
    ax6 = subplot(4,3,6);
    imagesc(data.(baselineType).Long.meanCortT,data.(baselineType).Long.meanCortF,data.(baselineType).Long.meanCortS)
    title('Long whisking cortical LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    colorbar
    caxis([-0.25 0.5])
    set(gca,'Ticklength',[0 0])
    axis square
    axis xy
    
    %% Short whisks HbT
    ax7 = subplot(4,3,7);
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCBV,'k');
    hold on
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCBV + data.(baselineType).Short.stdCBV,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCBV - data.(baselineType).Short.stdCBV,'color',colors_IOS('battleship grey'))
    title('Short whisking HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% Intermediate whisks HbT
    ax8 = subplot(4,3,8);
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCBV,'k');
    hold on
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCBV + data.(baselineType).Intermediate.stdCBV,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCBV - data.(baselineType).Intermediate.stdCBV,'color',colors_IOS('battleship grey'))
    title('Intermediate whisking HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% Long whisks HbT
    ax9 = subplot(4,3,9);
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCBV,'k');
    hold on
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCBV + data.(baselineType).Long.stdCBV,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCBV - data.(baselineType).Long.stdCBV,'color',colors_IOS('battleship grey'))
    title('Long whisking HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% Short whisks reflectance
    ax10 = subplot(4,3,10);
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCBV,'k');
    hold on
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCBV + data.(baselineType).Short.stdCBV,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCBV - data.(baselineType).Short.stdCBV,'color',colors_IOS('battleship grey'))
    title('Short whisking reflectance')
    ylabel('\DeltaR/R')
    xlabel('Time (sec)')
    axis square
    
    %% Intermediate whisks reflectance
    ax11 = subplot(4,3,11);
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCBV,'k');
    hold on
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCBV + data.(baselineType).Intermediate.stdCBV,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCBV - data.(baselineType).Intermediate.stdCBV,'color',colors_IOS('battleship grey'))
    title('Intermediate whisking reflectance')
    ylabel('\DeltaR/R')
    xlabel('Time (sec)')
    axis square
    
    %% Long whisks reflectance
    ax12 = subplot(4,3,12);
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCBV,'k');
    hold on
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCBV + data.(baselineType).Long.stdCBV,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCBV - data.(baselineType).Long.stdCBV,'color',colors_IOS('battleship grey'))
    title('Long whisking reflectance')
    ylabel('\DeltaR/R')
    xlabel('Time (sec)')
    axis square
    
    linkaxes([ax1 ax2 ax3],'xy')
    linkaxes([ax7 ax8 ax9],'xy')
    linkaxes([ax10 ax11 ax12],'xy')
    
    % save figure(s)
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Stimulus Responses\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure, [dirpath baselineType '_AverageCorticalWhiskResponses']);
end

%% summary figure(s)
for f = 1:length(baselineTypes)
    baselineType = baselineTypes{1,f};
    summaryFigure = figure;
    sgtitle({['Hippocampal whisking-evoked averages - ' baselineType],' '})
    
    %% Short whisks hippocampal MUA
    ax1 = subplot(4,3,1);
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanHipMUA,'k');
    hold on
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanHipMUA + data.(baselineType).Short.stdHipMUA,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanHipMUA - data.(baselineType).Short.stdHipMUA,'color',colors_IOS('battleship grey'))
    title('Short whisking hippocampal MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
     %% Intermediate whisks hippocampal MUA
    ax2 = subplot(4,3,2);
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanHipMUA,'k');
    hold on
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanHipMUA + data.(baselineType).Intermediate.stdHipMUA,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanHipMUA - data.(baselineType).Intermediate.stdHipMUA,'color',colors_IOS('battleship grey'))
    title('Short whisking hippocampal MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% Long whisks hippocampal MUA
    ax3 = subplot(4,3,3);
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanHipMUA,'k');
    hold on
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanHipMUA + data.(baselineType).Long.stdHipMUA,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanHipMUA - data.(baselineType).Long.stdHipMUA,'color',colors_IOS('battleship grey'))
    title('Short whisking hippocampal MUA')
    ylabel('Fold-change')
    xlabel('Time (sec)')
    axis square
    
    %% Short whisks hippocampal LFP 
    ax4 = subplot(4,3,4);
    imagesc(data.(baselineType).Short.meanHipT,data.(baselineType).Short.meanHipF,data.(baselineType).Short.meanHipS)
    title('Short whisking hippocampal LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    colorbar
    caxis([-0.25 0.5])
    set(gca,'Ticklength',[0 0])
    axis square
    axis xy
    
    %% Intermediate whisks hippocampal LFP 
    ax5 = subplot(4,3,5);
    imagesc(data.(baselineType).Intermediate.meanHipT,data.(baselineType).Intermediate.meanHipF,data.(baselineType).Intermediate.meanHipS)
    title('Intermediate whisking hippocampal LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    colorbar
    caxis([-0.25 0.5])
    set(gca,'Ticklength',[0 0])
    axis square
    axis xy
    
    %% Long whisks hippocampal LFP 
    ax6 = subplot(4,3,6);
    imagesc(data.(baselineType).Long.meanHipT,data.(baselineType).Long.meanHipF,data.(baselineType).Long.meanHipS)
    title('Long whisking hippocampal LFP')
    ylabel('Freq (Hz)')
    xlabel('Time (sec)')
    colorbar
    caxis([-0.25 0.5])
    set(gca,'Ticklength',[0 0])
    axis square
    axis xy
    
    %% Short whisks HbT
    ax7 = subplot(4,3,7);
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCBV,'k');
    hold on
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCBV + data.(baselineType).Short.stdCBV,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCBV - data.(baselineType).Short.stdCBV,'color',colors_IOS('battleship grey'))
    title('Short whisking HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% Intermediate whisks HbT
    ax8 = subplot(4,3,8);
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCBV,'k');
    hold on
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCBV + data.(baselineType).Intermediate.stdCBV,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCBV - data.(baselineType).Intermediate.stdCBV,'color',colors_IOS('battleship grey'))
    title('Intermediate whisking HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% Long whisks HbT
    ax9 = subplot(4,3,9);
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCBV,'k');
    hold on
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCBV + data.(baselineType).Long.stdCBV,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCBV - data.(baselineType).Long.stdCBV,'color',colors_IOS('battleship grey'))
    title('Long whisking HbT')
    ylabel('\DeltaHbT')
    xlabel('Time (sec)')
    axis square
    
    %% Short whisks reflectance
    ax10 = subplot(4,3,10);
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCBV,'k');
    hold on
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCBV + data.(baselineType).Short.stdCBV,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCBV - data.(baselineType).Short.stdCBV,'color',colors_IOS('battleship grey'))
    title('Short whisking reflectance')
    ylabel('\DeltaR/R')
    xlabel('Time (sec)')
    axis square
    
    %% Intermediate whisks reflectance
    ax11 = subplot(4,3,11);
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCBV,'k');
    hold on
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCBV + data.(baselineType).Intermediate.stdCBV,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCBV - data.(baselineType).Intermediate.stdCBV,'color',colors_IOS('battleship grey'))
    title('Intermediate whisking reflectance')
    ylabel('\DeltaR/R')
    xlabel('Time (sec)')
    axis square
    
    %% Long whisks reflectance
    ax12 = subplot(4,3,12);
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCBV,'k');
    hold on
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCBV + data.(baselineType).Long.stdCBV,'color',colors_IOS('battleship grey'))
    plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCBV - data.(baselineType).Long.stdCBV,'color',colors_IOS('battleship grey'))
    title('Long whisking reflectance')
    ylabel('\DeltaR/R')
    xlabel('Time (sec)')
    axis square
    
    linkaxes([ax1 ax2 ax3],'xy')
    linkaxes([ax7 ax8 ax9],'xy')
    linkaxes([ax10 ax11 ax12],'xy')
    
    % save figure(s)
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Stimulus Responses\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure, [dirpath baselineType '_AverageHippocampalWhiskResponses']);
end
