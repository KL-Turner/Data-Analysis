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
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
dataTypes = {'adjLH','adjRH'};

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\Turner_Manuscript_Summer2020\' animalID '\Bilateral Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for c = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,c};
        data.(whiskDataType).adjLH.HbT(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).CBV_HbT.HbT;
        data.(whiskDataType).adjLH.CBV(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).CBV.CBV;
        data.(whiskDataType).adjLH.cortMUA(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).MUA.corticalData;
        data.(whiskDataType).adjLH.cortS(:,:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.corticalS;
        data.(whiskDataType).adjLH.cortT(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.T;
        data.(whiskDataType).adjLH.cortF(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.F;
        
        data.(whiskDataType).adjRH.HbT(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjRH.(whiskDataType).CBV_HbT.HbT;
        data.(whiskDataType).adjRH.CBV(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjRH.(whiskDataType).CBV.CBV;
        data.(whiskDataType).adjRH.cortMUA(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjRH.(whiskDataType).MUA.corticalData;
        data.(whiskDataType).adjRH.cortS(:,:,a) = AnalysisResults.EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.corticalS;
        data.(whiskDataType).adjRH.cortT(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.T;
        data.(whiskDataType).adjRH.cortF(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.F;
        
        data.(whiskDataType).Hip.hipMUA(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).MUA.hippocampalData;
        data.(whiskDataType).Hip.hipS(:,:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.hippocampalS;
        data.(whiskDataType).Hip.hipT(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.T;
        data.(whiskDataType).Hip.hipF(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.F;
        
        data.(whiskDataType).timeVector(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).timeVector;
    end
end

% concatenate the data from the contra and ipsi data
for e = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,e};
    data.(whiskDataType).HbT = cat(2,data.(whiskDataType).adjLH.HbT,data.(whiskDataType).adjRH.HbT);
    data.(whiskDataType).CBV = cat(2,data.(whiskDataType).adjLH.CBV,data.(whiskDataType).adjRH.CBV);
    data.(whiskDataType).cortMUA = cat(2,data.(whiskDataType).adjLH.cortMUA,data.(whiskDataType).adjRH.cortMUA);
    data.(whiskDataType).cortS = cat(3,data.(whiskDataType).adjLH.cortS,data.(whiskDataType).adjRH.cortS);
    data.(whiskDataType).cortT = cat(2,data.(whiskDataType).adjLH.cortT,data.(whiskDataType).adjRH.cortT);
    data.(whiskDataType).cortF = cat(2,data.(whiskDataType).adjLH.cortF,data.(whiskDataType).adjRH.cortF);
end

% concatenate the data from the contra and ipsi data
for e = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,e};
    data.(whiskDataType).meanHbT = mean(data.(whiskDataType).HbT,2);
    data.(whiskDataType).stdHbT = std(data.(whiskDataType).HbT,0,2);
    data.(whiskDataType).meanCBV = mean(data.(whiskDataType).CBV,2);
    data.(whiskDataType).stdCBV = std(data.(whiskDataType).CBV,0,2);
    data.(whiskDataType).meanCortMUA = mean(data.(whiskDataType).cortMUA,2);
    data.(whiskDataType).stdCortMUA = std(data.(whiskDataType).cortMUA,0,2);
    data.(whiskDataType).meanCortS = mean(data.(whiskDataType).cortS,3).*100;
    data.(whiskDataType).meanCortT = mean(data.(whiskDataType).cortT,2);
    data.(whiskDataType).meanCortF = mean(data.(whiskDataType).cortF,2);
    data.(whiskDataType).meanHipMUA = mean(data.(whiskDataType).Hip.hipMUA,2);
    data.(whiskDataType).stdHipMUA = std(data.(whiskDataType).Hip.hipMUA,0,2);
    data.(whiskDataType).meanHipS = mean(data.(whiskDataType).Hip.hipS,3).*100;
    data.(whiskDataType).meanHipT = mean(data.(whiskDataType).Hip.hipT,2);
    data.(whiskDataType).meanHipF = mean(data.(whiskDataType).Hip.hipF,2);
    data.(whiskDataType).meanTimeVector = mean(data.(whiskDataType).timeVector(:,a),2);
end

%% summary figure(s)
summaryFigure = figure;
sgtitle('Whisking-evoked averages')

%% ShortWhisks whisks cortical MUA
ax1 = subplot(5,3,1);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA,'k');
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA + data.ShortWhisks.stdCortMUA,'color',colors_IOS('battleship grey'))
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA - data.ShortWhisks.stdCortMUA,'color',colors_IOS('battleship grey'))
title('Short whisking cortical MUA')
ylabel('\DeltaP/P (%)')
axis square
set(gca,'box','off')

%% IntermediateWhisks whisks cortical MUA
ax2 = subplot(5,3,2);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA,'k');
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA + data.IntermediateWhisks.stdCortMUA,'color',colors_IOS('battleship grey'))
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA - data.IntermediateWhisks.stdCortMUA,'color',colors_IOS('battleship grey'))
title('Intermed. whisking cortical MUA')
ylabel('\DeltaP/P (%)')
axis square
set(gca,'box','off')

%% LongWhisks whisks cortical MUA
ax3 = subplot(5,3,3);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA,'k');
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA + data.LongWhisks.stdCortMUA,'color',colors_IOS('battleship grey'))
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA - data.LongWhisks.stdCortMUA,'color',colors_IOS('battleship grey'))
title('Long whisking cortical MUA')
ylabel('\DeltaP/P (%)')
axis square
set(gca,'box','off')

%% ShortWhisks whisks cortical LFP
ax4 = subplot(5,3,4);
imagesc(data.ShortWhisks.meanCortT,data.ShortWhisks.meanCortF,data.ShortWhisks.meanCortS)
title('Short whisking cortical LFP')
ylabel('Frequency (Hz)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)')
caxis([-25 25])
set(gca,'Ticklength',[0 0])
axis square
axis xy
set(gca,'box','off')

%% IntermediateWhisks whisks cortical LFP
ax5 = subplot(5,3,5);
imagesc(data.IntermediateWhisks.meanCortT,data.IntermediateWhisks.meanCortF,data.IntermediateWhisks.meanCortS)
title('Intermed. whisking cortical LFP')
ylabel('Frequency (Hz)')
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)')
caxis([-25 25])
set(gca,'Ticklength',[0 0])
axis square
axis xy
set(gca,'box','off')

%% LongWhisks whisks cortical LFP
ax6 = subplot(5,3,6);
imagesc(data.LongWhisks.meanCortT,data.LongWhisks.meanCortF,data.LongWhisks.meanCortS)
title('Long whisking cortical LFP')
ylabel('Frequency (Hz)')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)')
caxis([-25 25])
set(gca,'Ticklength',[0 0])
axis square
axis xy
set(gca,'box','off')

%% ShortWhisks whisks hippocampal MUA
ax7 = subplot(5,3,7);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHipMUA,'k');
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHipMUA + data.ShortWhisks.stdHipMUA,'color',colors_IOS('battleship grey'))
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHipMUA - data.ShortWhisks.stdHipMUA,'color',colors_IOS('battleship grey'))
title('Short whisking hippocampal MUA')
ylabel('\DeltaP/P (%)')
axis square
set(gca,'box','off')

%% IntermediateWhisks whisks hippocampal MUA
ax8 = subplot(5,3,8);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHipMUA,'k');
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHipMUA + data.IntermediateWhisks.stdHipMUA,'color',colors_IOS('battleship grey'))
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHipMUA - data.IntermediateWhisks.stdHipMUA,'color',colors_IOS('battleship grey'))
title('Intermed. whisking hippocampal MUA')
ylabel('\DeltaP/P (%)')
axis square
set(gca,'box','off')

%% LongWhisks whisks hippocampal MUA
ax9 = subplot(5,3,9);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHipMUA,'k');
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHipMUA + data.LongWhisks.stdHipMUA,'color',colors_IOS('battleship grey'))
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHipMUA - data.LongWhisks.stdHipMUA,'color',colors_IOS('battleship grey'))
title('Long whisking hippocampal MUA')
ylabel('\DeltaP/P (%)')
axis square
set(gca,'box','off')

%% ShortWhisks whisks hippocampal LFP
ax10 = subplot(5,3,10);
imagesc(data.ShortWhisks.meanHipT,data.ShortWhisks.meanHipF,data.ShortWhisks.meanHipS)
title('Short whisking hippocampal LFP')
ylabel('Frequency (Hz)')
c10 = colorbar;
ylabel(c10,'\DeltaP/P (%)')
caxis([-25 25])
set(gca,'Ticklength',[0 0])
axis square
axis xy
set(gca,'box','off')

%% IntermediateWhisks whisks hippocampal LFP
ax11 = subplot(5,3,11);
imagesc(data.IntermediateWhisks.meanHipT,data.IntermediateWhisks.meanHipF,data.IntermediateWhisks.meanHipS)
title('Intermed. whisking hippocampal LFP')
ylabel('Frequency (Hz)')
c11 = colorbar;
ylabel(c11,'\DeltaP/P (%)')
caxis([-25 25])
set(gca,'Ticklength',[0 0])
axis square
axis xy
set(gca,'box','off')

%% Long whisks hippocampal LFP
ax12 = subplot(5,3,12);
imagesc(data.LongWhisks.meanHipT,data.LongWhisks.meanHipF,data.LongWhisks.meanHipS)
title('Long whisking hippocampal LFP')
ylabel('Frequency (Hz)')
c12 = colorbar;
ylabel(c12,'\DeltaP/P (%)')
caxis([-25 25])
set(gca,'Ticklength',[0 0])
axis square
axis xy
set(gca,'box','off')

%% Short whisks HbT
ax13 = subplot(5,6,25);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHbT,'k');
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHbT + data.ShortWhisks.stdHbT,'color',colors_IOS('battleship grey'))
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHbT - data.ShortWhisks.stdHbT,'color',colors_IOS('battleship grey'))
title('Short whisking \DeltaHbT (\muM)')
ylabel('\DeltaHbT (\muM)')
xlabel('Peristimuls time (s)')
axis square
set(gca,'box','off')

%% Short whisks refl
ax14 = subplot(5,6,26);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCBV,'k');
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCBV + data.ShortWhisks.stdCBV,'color',colors_IOS('battleship grey'))
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCBV - data.ShortWhisks.stdCBV,'color',colors_IOS('battleship grey'))
title('Short whisking reflectance')
ylabel('\DeltaR/R (%)')
xlabel('Peristimuls time (s)')
axis square
set(gca,'box','off')

%% Intermediate whisks HbT
ax15 = subplot(5,6,27);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHbT,'k');
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHbT + data.IntermediateWhisks.stdHbT,'color',colors_IOS('battleship grey'))
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHbT - data.IntermediateWhisks.stdHbT,'color',colors_IOS('battleship grey'))
title('Intermed. whisking \DeltaHbT (\muM)')
ylabel('\DeltaHbT (\muM)')
xlabel('Peristimuls time (s)')
axis square
set(gca,'box','off')

%% Intermediate whisks refl
ax16 = subplot(5,6,28);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCBV,'k');
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCBV + data.IntermediateWhisks.stdCBV,'color',colors_IOS('battleship grey'))
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCBV - data.IntermediateWhisks.stdCBV,'color',colors_IOS('battleship grey'))
title('Intermed. whisking reflectance')
ylabel('\DeltaR/R (%)')
xlabel('Peristimuls time (s)')
axis square
set(gca,'box','off')

%% Long whisks HbT
ax17 = subplot(5,6,29);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHbT,'k');
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHbT + data.LongWhisks.stdHbT,'color',colors_IOS('battleship grey'))
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHbT - data.LongWhisks.stdHbT,'color',colors_IOS('battleship grey'))
title('Long whisking \DeltaHbT (\muM)')
ylabel('\DeltaHbT (\muM)')
xlabel('Peristimuls time (s)')
axis square
set(gca,'box','off')

%% Long whisks refl
ax18 = subplot(5,6,30);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCBV,'k');
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCBV + data.LongWhisks.stdCBV,'color',colors_IOS('battleship grey'))
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCBV - data.LongWhisks.stdCBV,'color',colors_IOS('battleship grey'))
title('Long whisking reflectance')
ylabel('\DeltaR/R (%)')
xlabel('Peristimuls time (s)')
axis square
set(gca,'box','off')

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
savefig(summaryFigure,[dirpath 'Summary Figure - Whisk Responses']);
