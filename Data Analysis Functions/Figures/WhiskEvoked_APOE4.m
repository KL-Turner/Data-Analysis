function [AnalysisResults] = WhiskEvoked_APOE4(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 3-S1 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
animalIDs = {'T115','T116','T117','T118','T125','T126','T201','T202'};
C57BL6J_IDs = {'T115','T116','T117','T118','T125'};
APOE4_IDs = {'T201','T202'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
treatments = {'C57BL6J','APOE4'};
% cd through each animal's directory and extract the appropriate analysis results
data.C57BL6J.EvokedAvgs.ShortWhisks.means = []; data.APOE4.EvokedAvgs.ShortWhisks.means = [];
data.C57BL6J.EvokedAvgs.IntermediateWhisks.means = []; data.APOE4.EvokedAvgs.IntermediateWhisks.means = [];
data.C57BL6J.EvokedAvgs.LongWhisks.means = []; data.APOE4.EvokedAvgs.LongWhisks.means = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    % recognize treatment based on animal group
    if ismember(animalIDs{1,aa},C57BL6J_IDs) == true
        treatment = 'C57BL6J';
    elseif ismember(animalIDs{1,aa},APOE4_IDs) == true
        treatment = 'APOE4';
    end
    for bb = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,bb};
        vesselIDs = fieldnames(AnalysisResults.(animalID).EvokedAvgs.(whiskDataType));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            if strcmp(vesselID(1),'V') == false
                data.(treatment).EvokedAvgs.(whiskDataType).means = vertcat(data.(treatment).EvokedAvgs.(whiskDataType).means,AnalysisResults.(animalID).EvokedAvgs.(whiskDataType).(vesselID).mean);
                data.(treatment).EvokedAvgs.(whiskDataType).timeVector = AnalysisResults.(animalID).EvokedAvgs.(whiskDataType).(vesselID).timeVector;
            end
        end
    end
end
% mean/std of the data
for dd = 1:length(treatments)
    treatment = treatments{1,dd};
    for ee = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,ee};
        data.(treatment).EvokedAvgs.(whiskDataType).mean = mean(data.(treatment).EvokedAvgs.(whiskDataType).means,1);
        data.(treatment).EvokedAvgs.(whiskDataType).StD = std(data.(treatment).EvokedAvgs.(whiskDataType).means,0,1);
    end
end
%% Fig. 3-S1
summaryFigure = figure;
sgtitle('Arteriole diameter whisking-evoked responses')
%% [3-S1a] brief whisks
ax1 = subplot(1,3,1);
p1 = plot(data.C57BL6J.EvokedAvgs.ShortWhisks.timeVector,data.C57BL6J.EvokedAvgs.ShortWhisks.mean,'color',colors_eLife2020('rich black'),'LineWidth',1);
hold on
plot(data.C57BL6J.EvokedAvgs.ShortWhisks.timeVector,data.C57BL6J.EvokedAvgs.ShortWhisks.mean + data.C57BL6J.EvokedAvgs.ShortWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.C57BL6J.EvokedAvgs.ShortWhisks.timeVector,data.C57BL6J.EvokedAvgs.ShortWhisks.mean - data.C57BL6J.EvokedAvgs.ShortWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
p2 = plot(data.APOE4.EvokedAvgs.ShortWhisks.timeVector,data.APOE4.EvokedAvgs.ShortWhisks.mean,'color',colors_eLife2020('dark candy apple red'),'LineWidth',1);
plot(data.APOE4.EvokedAvgs.ShortWhisks.timeVector,data.APOE4.EvokedAvgs.ShortWhisks.mean + data.APOE4.EvokedAvgs.ShortWhisks.StD,'color',colors_eLife2020('candy apple red'),'LineWidth',0.5)
plot(data.APOE4.EvokedAvgs.ShortWhisks.timeVector,data.APOE4.EvokedAvgs.ShortWhisks.mean - data.APOE4.EvokedAvgs.ShortWhisks.StD,'color',colors_eLife2020('candy apple red'),'LineWidth',0.5)
title('Brief whisk response')
ylabel('\DeltaD/D (%)')
xlabel('Peri-whisk time (s)')
legend([p1,p2],'C57BL6J','APOE4')
axis square
xlim([-2,10])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [3-S1b] moderate whisks
ax2 = subplot(1,3,2);
plot(data.C57BL6J.EvokedAvgs.IntermediateWhisks.timeVector,data.C57BL6J.EvokedAvgs.IntermediateWhisks.mean,'color',colors_eLife2020('rich black'),'LineWidth',1);
hold on
plot(data.C57BL6J.EvokedAvgs.IntermediateWhisks.timeVector,data.C57BL6J.EvokedAvgs.IntermediateWhisks.mean + data.C57BL6J.EvokedAvgs.IntermediateWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.C57BL6J.EvokedAvgs.IntermediateWhisks.timeVector,data.C57BL6J.EvokedAvgs.IntermediateWhisks.mean - data.C57BL6J.EvokedAvgs.IntermediateWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.APOE4.EvokedAvgs.IntermediateWhisks.timeVector,data.APOE4.EvokedAvgs.IntermediateWhisks.mean,'color',colors_eLife2020('dark candy apple red'),'LineWidth',1);
plot(data.APOE4.EvokedAvgs.IntermediateWhisks.timeVector,data.APOE4.EvokedAvgs.IntermediateWhisks.mean + data.APOE4.EvokedAvgs.IntermediateWhisks.StD,'color',colors_eLife2020('candy apple red'),'LineWidth',0.5)
plot(data.APOE4.EvokedAvgs.IntermediateWhisks.timeVector,data.APOE4.EvokedAvgs.IntermediateWhisks.mean - data.APOE4.EvokedAvgs.IntermediateWhisks.StD,'color',colors_eLife2020('candy apple red'),'LineWidth',0.5)
title('Moderate whisk response')
ylabel('\DeltaD/D (%)')
xlabel('Peri-whisk time (s)')
axis square
xlim([-2,10])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [3-S1c] extended whisks
ax3 = subplot(1,3,3);
plot(data.C57BL6J.EvokedAvgs.LongWhisks.timeVector,data.C57BL6J.EvokedAvgs.LongWhisks.mean,'color',colors_eLife2020('rich black'),'LineWidth',1);
hold on
plot(data.C57BL6J.EvokedAvgs.LongWhisks.timeVector,data.C57BL6J.EvokedAvgs.LongWhisks.mean + data.C57BL6J.EvokedAvgs.LongWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.C57BL6J.EvokedAvgs.LongWhisks.timeVector,data.C57BL6J.EvokedAvgs.LongWhisks.mean - data.C57BL6J.EvokedAvgs.LongWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.APOE4.EvokedAvgs.LongWhisks.timeVector,data.APOE4.EvokedAvgs.LongWhisks.mean,'color',colors_eLife2020('dark candy apple red'),'LineWidth',1);
plot(data.APOE4.EvokedAvgs.LongWhisks.timeVector,data.APOE4.EvokedAvgs.LongWhisks.mean + data.APOE4.EvokedAvgs.LongWhisks.StD,'color',colors_eLife2020('candy apple red'),'LineWidth',0.5)
plot(data.APOE4.EvokedAvgs.LongWhisks.timeVector,data.APOE4.EvokedAvgs.LongWhisks.mean - data.APOE4.EvokedAvgs.LongWhisks.StD,'color',colors_eLife2020('candy apple red'),'LineWidth',0.5)
title('Extended whisk response')
ylabel('\DeltaD/D (%)')
xlabel('Peri-whisk time (s)')
axis square
xlim([-2,10])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
linkaxes([ax1,ax2,ax3],'xy')
% %% [3-S1a] brief whisks
% ax4 = subplot(2,3,4);
% for aa = 1:size(data.C57BL6J.EvokedAvgs.ShortWhisks.means,1)
%     plot(data.C57BL6J.EvokedAvgs.ShortWhisks.timeVector,data.C57BL6J.EvokedAvgs.ShortWhisks.means(aa,:),'color',colors_eLife2020('rich black'),'LineWidth',0.25);
%     hold on
% end
% for bb = 1:size(data.APOE4.EvokedAvgs.ShortWhisks.means,1)
%     plot(data.APOE4.EvokedAvgs.ShortWhisks.timeVector,data.APOE4.EvokedAvgs.ShortWhisks.means(bb,:),'color',colors_eLife2020('dark candy apple red'),'LineWidth',0.25);
%     hold on
% end
% title('Brief whisk response')
% ylabel('\DeltaD/D (%)')
% xlabel('Peri-whisk time (s)')
% legend([p1,p2],'C57BL6J','APOE4')
% axis square
% xlim([-2,10])
% set(gca,'box','off')
% ax4.TickLength = [0.03,0.03];
% %% [3-S1b] moderate whisks
% ax5 = subplot(2,3,5);
% for aa = 1:size(data.C57BL6J.EvokedAvgs.IntermediateWhisks.means,1)
%     plot(data.C57BL6J.EvokedAvgs.IntermediateWhisks.timeVector,data.C57BL6J.EvokedAvgs.IntermediateWhisks.means(aa,:),'color',colors_eLife2020('rich black'),'LineWidth',0.25);
%     hold on
% end
% for bb = 1:size(data.APOE4.EvokedAvgs.IntermediateWhisks.means,1)
%     plot(data.APOE4.EvokedAvgs.IntermediateWhisks.timeVector,data.APOE4.EvokedAvgs.IntermediateWhisks.means(bb,:),'color',colors_eLife2020('dark candy apple red'),'LineWidth',0.25);
%     hold on
% end
% title('Moderate whisk response')
% ylabel('\DeltaD/D (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% xlim([-2,10])
% set(gca,'box','off')
% ax5.TickLength = [0.03,0.03];
% %% [3-S1c] extended whisks
% ax6 = subplot(2,3,6);
% for aa = 1:size(data.C57BL6J.EvokedAvgs.LongWhisks.means,1)
%     plot(data.C57BL6J.EvokedAvgs.LongWhisks.timeVector,data.C57BL6J.EvokedAvgs.LongWhisks.means(aa,:),'color',colors_eLife2020('rich black'),'LineWidth',0.25);
%     hold on
% end
% for bb = 1:size(data.APOE4.EvokedAvgs.LongWhisks.means,1)
%     plot(data.APOE4.EvokedAvgs.LongWhisks.timeVector,data.APOE4.EvokedAvgs.LongWhisks.means(bb,:),'color',colors_eLife2020('dark candy apple red'),'LineWidth',0.25);
%     hold on
% end
% title('Extended whisk response')
% ylabel('\DeltaD/D (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% xlim([-2,10])
% set(gca,'box','off')
% ax6.TickLength = [0.03,0.03];
% linkaxes([ax4,ax5,ax6],'xy')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'APOE_WhiskEvoked-S1']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'APOE_WhiskEvoked'])
end

end
