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
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
dataTypes = {'adjLH','adjRH'};

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for c = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,c};
        data.(whiskDataType).adjLH.HbT(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).HbT;
        data.(whiskDataType).adjLH.cortMUA(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).MUA.corticalData;
        data.(whiskDataType).adjLH.cortS(:,:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.corticalS;
        data.(whiskDataType).adjLH.cortT(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.T;
        data.(whiskDataType).adjLH.cortF(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.F;
        
        data.(whiskDataType).adjRH.HbT(:,a) = AnalysisResults.EvokedAvgs.Whisk.adjRH.(whiskDataType).HbT;
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
    data.(whiskDataType).meanCortMUA = mean(data.(whiskDataType).cortMUA,2);
    data.(whiskDataType).stdCortMUA = std(data.(whiskDataType).cortMUA,0,2);
    data.(whiskDataType).meanCortS = mean(data.(whiskDataType).cortS,3);
    data.(whiskDataType).meanCortT = mean(data.(whiskDataType).cortT,2);
    data.(whiskDataType).meanCortF = mean(data.(whiskDataType).cortF,2);
    data.(whiskDataType).meanHipMUA = mean(data.(whiskDataType).Hip.hipMUA,2);
    data.(whiskDataType).stdHipMUA = std(data.(whiskDataType).Hip.hipMUA,0,2);
    data.(whiskDataType).meanHipS = mean(data.(whiskDataType).Hip.hipS,3);
    data.(whiskDataType).meanHipT = mean(data.(whiskDataType).Hip.hipT,2);
    data.(whiskDataType).meanHipF = mean(data.(whiskDataType).Hip.hipF,2);
    data.(whiskDataType).meanTimeVector = mean(data.(whiskDataType).timeVector(:,a),2);
end

%% summary figure(s)
summaryFigure = figure;
sgtitle('Whisking-evoked averages')

%% Short whisks cortical MUA
ax1 = subplot(5,3,1);
plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCortMUA,'k');
hold on
plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCortMUA + data.(baselineType).Short.stdCortMUA,'color',colors_IOS('battleship grey'))
plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanCortMUA - data.(baselineType).Short.stdCortMUA,'color',colors_IOS('battleship grey'))
title('Short whisking cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peristimuls time (s)')
axis square

%% Intermediate whisks cortical MUA
ax2 = subplot(5,3,2);
plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCortMUA,'k');
hold on
plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCortMUA + data.(baselineType).Intermediate.stdCortMUA,'color',colors_IOS('battleship grey'))
plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanCortMUA - data.(baselineType).Intermediate.stdCortMUA,'color',colors_IOS('battleship grey'))
title('Short whisking cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peristimuls time (s)')
axis square

%% Long whisks cortical MUA
ax3 = subplot(5,3,3);
plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCortMUA,'k');
hold on
plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCortMUA + data.(baselineType).Long.stdCortMUA,'color',colors_IOS('battleship grey'))
plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanCortMUA - data.(baselineType).Long.stdCortMUA,'color',colors_IOS('battleship grey'))
title('Short whisking cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peristimuls time (s)')
axis square

%% Short whisks cortical LFP
ax4 = subplot(5,3,4);
imagesc(data.(baselineType).Short.meanCortT,data.(baselineType).Short.meanCortF,data.(baselineType).Short.meanCortS)
title('Short whisking cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peristimuls time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)')
caxis([-0.25 0.5]);
set(gca,'Ticklength',[0 0])
axis square
axis xy

%% Intermediate whisks cortical LFP
ax5 = subplot(5,3,5);
imagesc(data.(baselineType).Intermediate.meanCortT,data.(baselineType).Intermediate.meanCortF,data.(baselineType).Intermediate.meanCortS)
title('Intermediate whisking cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peristimuls time (s)')
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)')
caxis([-0.25 0.5])
set(gca,'Ticklength',[0 0])
axis square
axis xy

%% Long whisks cortical LFP
ax6 = subplot(5,3,6);
imagesc(data.(baselineType).Long.meanCortT,data.(baselineType).Long.meanCortF,data.(baselineType).Long.meanCortS)
title('Long whisking cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peristimuls time (s)')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)')
caxis([-0.25 0.5])
set(gca,'Ticklength',[0 0])
axis square
axis xy

%% Short whisks hippocampal MUA
ax7 = subplot(5,3,7);
plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanHipMUA,'k');
hold on
plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanHipMUA + data.(baselineType).Short.stdHipMUA,'color',colors_IOS('battleship grey'))
plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanHipMUA - data.(baselineType).Short.stdHipMUA,'color',colors_IOS('battleship grey'))
title('Short whisking hippocampal MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peristimuls time (s)')
axis square

%% Intermediate whisks hippocampal MUA
ax8 = subplot(5,3,8);
plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanHipMUA,'k');
hold on
plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanHipMUA + data.(baselineType).Intermediate.stdHipMUA,'color',colors_IOS('battleship grey'))
plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanHipMUA - data.(baselineType).Intermediate.stdHipMUA,'color',colors_IOS('battleship grey'))
title('Short whisking hippocampal MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peristimuls time (s)')
axis square

%% Long whisks hippocampal MUA
ax9 = subplot(5,3,9);
plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanHipMUA,'k');
hold on
plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanHipMUA + data.(baselineType).Long.stdHipMUA,'color',colors_IOS('battleship grey'))
plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanHipMUA - data.(baselineType).Long.stdHipMUA,'color',colors_IOS('battleship grey'))
title('Short whisking hippocampal MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peristimuls time (s)')
axis square

%% Short whisks hippocampal LFP
ax10 = subplot(5,3,10);
imagesc(data.(baselineType).Short.meanHipT,data.(baselineType).Short.meanHipF,data.(baselineType).Short.meanHipS)
title('Short whisking hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peristimuls time (s)')
c10 = colorbar;
ylabel(c10,'\DeltaP/P (%)')
caxis([-0.25 0.5])
set(gca,'Ticklength',[0 0])
axis square
axis xy

%% Intermediate whisks hippocampal LFP
ax11 = subplot(5,3,11);
imagesc(data.(baselineType).Intermediate.meanHipT,data.(baselineType).Intermediate.meanHipF,data.(baselineType).Intermediate.meanHipS)
title('Intermediate whisking hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peristimuls time (s)')
c11 = colorbar;
ylabel(c11,'\DeltaP/P (%)')
caxis([-0.25 0.5])
set(gca,'Ticklength',[0 0])
axis square
axis xy

%% Long whisks hippocampal LFP
ax12 = subplot(5,3,12);
imagesc(data.(baselineType).Long.meanHipT,data.(baselineType).Long.meanHipF,data.(baselineType).Long.meanHipS)
title('Long whisking hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peristimuls time (s)')
c12 = colorbar;
ylabel(c12,'\DeltaP/P (%)')
caxis([-0.25 0.5])
set(gca,'Ticklength',[0 0])
axis square
axis xy

%% Short whisks HbT
ax13 = subplot(5,3,13);
plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanHbT,'k');
hold on
plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanHbT + data.(baselineType).Short.stdHbT,'color',colors_IOS('battleship grey'))
plot(data.(baselineType).Short.meanTimeVector,data.(baselineType).Short.meanHbT - data.(baselineType).Short.stdHbT,'color',colors_IOS('battleship grey'))
title('Short whisking HbT')
ylabel('\DeltaHbT')
xlabel('Peristimuls time (s)')
axis square

%% Intermediate whisks HbT
ax14 = subplot(5,3,14);
plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanHbT,'k');
hold on
plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanHbT + data.(baselineType).Intermediate.stdHbT,'color',colors_IOS('battleship grey'))
plot(data.(baselineType).Intermediate.meanTimeVector,data.(baselineType).Intermediate.meanHbT - data.(baselineType).Intermediate.stdHbT,'color',colors_IOS('battleship grey'))
title('Intermediate whisking HbT')
ylabel('\DeltaHbT')
xlabel('Peristimuls time (s)')
axis square

%% Long whisks HbT
ax15 = subplot(5,3,15);
plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanHbT,'k');
hold on
plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanHbT + data.(baselineType).Long.stdHbT,'color',colors_IOS('battleship grey'))
plot(data.(baselineType).Long.meanTimeVector,data.(baselineType).Long.meanHbT - data.(baselineType).Long.stdHbT,'color',colors_IOS('battleship grey'))
title('Long whisking HbT')
ylabel('\DeltaHbT')
xlabel('Peristimuls time (s)')
axis square

linkaxes([ax1 ax2 ax3 ax7 ax8 ax9],'xy')
linkaxes([ax4 ax5 ax6 ax10 ax11 ax12],'xy')
linkaxes([ax13 ax14 ax15],'xy')

% save figure(s)
dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Whisking-evoked Responses\';
if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Average Whisk Responses']);
