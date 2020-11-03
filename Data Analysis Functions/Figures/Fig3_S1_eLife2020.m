function [AnalysisResults] = Fig3_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 3-S1 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
animalIDs = {'T115','T116','T117','T118','T125','T126'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
% cd through each animal's directory and extract the appropriate analysis results
data.EvokedAvgs.ShortWhisks.means = [];
data.EvokedAvgs.IntermediateWhisks.means = [];
data.EvokedAvgs.LongWhisks.means = [];
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for c = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,c};
        vesselIDs = fieldnames(AnalysisResults.(animalID).EvokedAvgs.(whiskDataType));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            if strcmp(vesselID(1),'V') == false
                % LH cortical
                data.EvokedAvgs.(whiskDataType).means = vertcat(data.EvokedAvgs.(whiskDataType).means,AnalysisResults.(animalID).EvokedAvgs.(whiskDataType).(vesselID).mean);
                data.EvokedAvgs.(whiskDataType).timeVector = AnalysisResults.(animalID).EvokedAvgs.(whiskDataType).(vesselID).timeVector;
            end
        end
    end
end
% mean/std of the data
for e = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,e};
    data.EvokedAvgs.(whiskDataType).mean = mean(data.EvokedAvgs.(whiskDataType).means,1);
    data.EvokedAvgs.(whiskDataType).StD = std(data.EvokedAvgs.(whiskDataType).means,0,1);
end
%% Fig. 3-S1
summaryFigure = figure('Name','Fig3-S1 (a-c)');
sgtitle('Figure 3-S1 - Turner et al. 2020')
%% [3-S1a] brief whisks
ax1 = subplot(1,3,1);
plot(data.EvokedAvgs.ShortWhisks.timeVector,data.EvokedAvgs.ShortWhisks.mean,'color',colors_eLife2020('rich black'),'LineWidth',1);
hold on
plot(data.EvokedAvgs.ShortWhisks.timeVector,data.EvokedAvgs.ShortWhisks.mean + data.EvokedAvgs.ShortWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.EvokedAvgs.ShortWhisks.timeVector,data.EvokedAvgs.ShortWhisks.mean - data.EvokedAvgs.ShortWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[3-S1a] Brief whisk response')
ylabel('\DeltaD/D (%)')
xlabel('Peri-whisk time (s)')
axis square
xlim([-2,10])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [3-S1b] moderate whisks
ax2 = subplot(1,3,2);
plot(data.EvokedAvgs.IntermediateWhisks.timeVector,data.EvokedAvgs.IntermediateWhisks.mean,'color',colors_eLife2020('rich black'),'LineWidth',1);
hold on
plot(data.EvokedAvgs.IntermediateWhisks.timeVector,data.EvokedAvgs.IntermediateWhisks.mean + data.EvokedAvgs.IntermediateWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.EvokedAvgs.IntermediateWhisks.timeVector,data.EvokedAvgs.IntermediateWhisks.mean - data.EvokedAvgs.IntermediateWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[3d,3-S1b] Moderate whisk response')
ylabel('\DeltaD/D (%)')
xlabel('Peri-whisk time (s)')
axis square
xlim([-2,10])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [3-S1c] extended whisks
ax3 = subplot(1,3,3);
plot(data.EvokedAvgs.LongWhisks.timeVector,data.EvokedAvgs.LongWhisks.mean,'color',colors_eLife2020('rich black'),'LineWidth',1);
hold on
plot(data.EvokedAvgs.LongWhisks.timeVector,data.EvokedAvgs.LongWhisks.mean + data.EvokedAvgs.LongWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.EvokedAvgs.LongWhisks.timeVector,data.EvokedAvgs.LongWhisks.mean - data.EvokedAvgs.LongWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
title('[3-S1c] Extended whisk response')
ylabel('\DeltaD/D (%)')
xlabel('Peri-whisk time (s)')
axis square
xlim([-2,10])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
linkaxes([ax1,ax2,ax3],'xy')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig3-S1']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig3-S1'])
    %% Text diary
    diaryFile = [dirpath 'Fig3-S1_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % text values
    disp('======================================================================================================================')
    disp('[3-S1] Text values for arteriole D/D changes')
    disp('======================================================================================================================')
    disp('----------------------------------------------------------------------------------------------------------------------')
    [~,index] = max(data.EvokedAvgs.ShortWhisks.mean);
    disp(['Brief whisk D/D (%): ' num2str(round(data.EvokedAvgs.ShortWhisks.mean(index),1)) ' +/- ' num2str(round(data.EvokedAvgs.ShortWhisks.StD(index),1))]); disp(' ')
    [~,index] = max(data.EvokedAvgs.IntermediateWhisks.mean);
    disp(['Moderate whisk D/D (%): ' num2str(round(data.EvokedAvgs.IntermediateWhisks.mean(index),1)) ' +/- ' num2str(round(data.EvokedAvgs.IntermediateWhisks.StD(index),1))]); disp(' ')
    [~,index] = max(data.EvokedAvgs.LongWhisks.mean);
    disp(['Extended whisk D/D (%): ' num2str(round(data.EvokedAvgs.LongWhisks.mean(index),1)) ' +/- ' num2str(round(data.EvokedAvgs.LongWhisks.StD(index),1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
end

end
