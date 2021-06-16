function [AnalysisResults] = WhiskEvoked_Wenke(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 3-S1 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
animalIDs = {'MK01_D','MK02_D','MK01_C','MK02_C'};
Dural_IDs = {'MK01_D','MK02_D'};
Capillary_IDs = {'MK01_C','MK02_C'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
treatments = {'Dural','Capillary'};
% cd through each animal's directory and extract the appropriate analysis results
data.Dural.EvokedAvgs.ShortWhisks.means = []; data.Capillary.EvokedAvgs.ShortWhisks.means = [];
data.Dural.EvokedAvgs.IntermediateWhisks.means = []; data.Capillary.EvokedAvgs.IntermediateWhisks.means = [];
data.Dural.EvokedAvgs.LongWhisks.means = []; data.Capillary.EvokedAvgs.LongWhisks.means = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    % recognize treatment based on animal group
    if ismember(animalIDs{1,aa},Dural_IDs) == true
        treatment = 'Dural';
    elseif ismember(animalIDs{1,aa},Capillary_IDs) == true
        treatment = 'Capillary';
    end
    for bb = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,bb};
        vesselIDs = fieldnames(AnalysisResults.(treatment).(animalID).EvokedAvgs.(whiskDataType));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            if strcmp(vesselID(1),'V') == false
                data.(treatment).EvokedAvgs.(whiskDataType).means = vertcat(data.(treatment).EvokedAvgs.(whiskDataType).means,AnalysisResults.(treatment).(animalID).EvokedAvgs.(whiskDataType).(vesselID).mean);
                data.(treatment).EvokedAvgs.(whiskDataType).timeVector = AnalysisResults.(treatment).(animalID).EvokedAvgs.(whiskDataType).(vesselID).timeVector;
            end
        end
    end
end
% mean/std of the data
for dd = 1:length(treatments)
    treatment = treatments{1,dd};
    for ee = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,ee};
        data.(treatment).EvokedAvgs.(whiskDataType).mean = mean(data.(treatment).EvokedAvgs.(whiskDataType).means,1)/1000;
        data.(treatment).EvokedAvgs.(whiskDataType).StD = std(data.(treatment).EvokedAvgs.(whiskDataType).means,0,1)/1000;
    end
end
%% Fig. 3-S1
summaryFigure = figure;
sgtitle('Arteriole diameter whisking-evoked responses')
%% [3-S1a] brief whisks
ax1 = subplot(1,3,1);
p1 = plot(data.Dural.EvokedAvgs.ShortWhisks.timeVector,data.Dural.EvokedAvgs.ShortWhisks.mean,'color',colors_eLife2020('rich black'),'LineWidth',3);
hold on
plot(data.Dural.EvokedAvgs.ShortWhisks.timeVector,data.Dural.EvokedAvgs.ShortWhisks.mean + data.Dural.EvokedAvgs.ShortWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Dural.EvokedAvgs.ShortWhisks.timeVector,data.Dural.EvokedAvgs.ShortWhisks.mean - data.Dural.EvokedAvgs.ShortWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
p2 = plot(data.Capillary.EvokedAvgs.ShortWhisks.timeVector,data.Capillary.EvokedAvgs.ShortWhisks.mean,'color',colors_eLife2020('dark candy apple red'),'LineWidth',3);
plot(data.Capillary.EvokedAvgs.ShortWhisks.timeVector,data.Capillary.EvokedAvgs.ShortWhisks.mean + data.Capillary.EvokedAvgs.ShortWhisks.StD,'color',colors_eLife2020('candy apple red'),'LineWidth',0.5)
plot(data.Capillary.EvokedAvgs.ShortWhisks.timeVector,data.Capillary.EvokedAvgs.ShortWhisks.mean - data.Capillary.EvokedAvgs.ShortWhisks.StD,'color',colors_eLife2020('candy apple red'),'LineWidth',0.5)
title('Brief whisk response')
ylabel('\DeltaV (mm/sec)')
xlabel('Peri-whisk time (s)')
legend([p1,p2],'Dural','Capillary')
axis square
xlim([-2,10])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [3-S1b] moderate whisks
ax2 = subplot(1,3,2);
plot(data.Dural.EvokedAvgs.IntermediateWhisks.timeVector,data.Dural.EvokedAvgs.IntermediateWhisks.mean,'color',colors_eLife2020('rich black'),'LineWidth',3);
hold on
plot(data.Dural.EvokedAvgs.IntermediateWhisks.timeVector,data.Dural.EvokedAvgs.IntermediateWhisks.mean + data.Dural.EvokedAvgs.IntermediateWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Dural.EvokedAvgs.IntermediateWhisks.timeVector,data.Dural.EvokedAvgs.IntermediateWhisks.mean - data.Dural.EvokedAvgs.IntermediateWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Capillary.EvokedAvgs.IntermediateWhisks.timeVector,data.Capillary.EvokedAvgs.IntermediateWhisks.mean,'color',colors_eLife2020('dark candy apple red'),'LineWidth',3);
plot(data.Capillary.EvokedAvgs.IntermediateWhisks.timeVector,data.Capillary.EvokedAvgs.IntermediateWhisks.mean + data.Capillary.EvokedAvgs.IntermediateWhisks.StD,'color',colors_eLife2020('candy apple red'),'LineWidth',0.5)
plot(data.Capillary.EvokedAvgs.IntermediateWhisks.timeVector,data.Capillary.EvokedAvgs.IntermediateWhisks.mean - data.Capillary.EvokedAvgs.IntermediateWhisks.StD,'color',colors_eLife2020('candy apple red'),'LineWidth',0.5)
title('Moderate whisk response')
ylabel('\DeltaV (mm/sec)')
xlabel('Peri-whisk time (s)')
axis square
xlim([-2,10])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [3-S1c] extended whisks
ax3 = subplot(1,3,3);
plot(data.Dural.EvokedAvgs.LongWhisks.timeVector,data.Dural.EvokedAvgs.LongWhisks.mean,'color',colors_eLife2020('rich black'),'LineWidth',3);
hold on
plot(data.Dural.EvokedAvgs.LongWhisks.timeVector,data.Dural.EvokedAvgs.LongWhisks.mean + data.Dural.EvokedAvgs.LongWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Dural.EvokedAvgs.LongWhisks.timeVector,data.Dural.EvokedAvgs.LongWhisks.mean - data.Dural.EvokedAvgs.LongWhisks.StD,'color',colors_eLife2020('battleship grey'),'LineWidth',0.5)
plot(data.Capillary.EvokedAvgs.LongWhisks.timeVector,data.Capillary.EvokedAvgs.LongWhisks.mean,'color',colors_eLife2020('dark candy apple red'),'LineWidth',3);
plot(data.Capillary.EvokedAvgs.LongWhisks.timeVector,data.Capillary.EvokedAvgs.LongWhisks.mean + data.Capillary.EvokedAvgs.LongWhisks.StD,'color',colors_eLife2020('candy apple red'),'LineWidth',0.5)
plot(data.Capillary.EvokedAvgs.LongWhisks.timeVector,data.Capillary.EvokedAvgs.LongWhisks.mean - data.Capillary.EvokedAvgs.LongWhisks.StD,'color',colors_eLife2020('candy apple red'),'LineWidth',0.5)
title('Extended whisk response')
ylabel('\DeltaV (mm/sec)')
xlabel('Peri-whisk time (s)')
axis square
xlim([-2,10])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
linkaxes([ax1,ax2,ax3],'xy')
% %% [3-S1a] brief whisks
% ax4 = subplot(2,3,4);
% for aa = 1:size(data.Dural.EvokedAvgs.ShortWhisks.means,1)
%     plot(data.Dural.EvokedAvgs.ShortWhisks.timeVector,data.Dural.EvokedAvgs.ShortWhisks.means(aa,:),'color',colors_eLife2020('rich black'),'LineWidth',0.25);
%     hold on
% end
% for bb = 1:size(data.Capillary.EvokedAvgs.ShortWhisks.means,1)
%     plot(data.Capillary.EvokedAvgs.ShortWhisks.timeVector,data.Capillary.EvokedAvgs.ShortWhisks.means(bb,:),'color',colors_eLife2020('dark candy apple red'),'LineWidth',0.25);
%     hold on
% end
% title('Brief whisk response')
% ylabel('\DeltaV (mm/sec)')
% xlabel('Peri-whisk time (s)')
% legend([p1,p2],'Dural','Capillary')
% axis square
% xlim([-2,10])
% set(gca,'box','off')
% ax4.TickLength = [0.03,0.03];
% %% [3-S1b] moderate whisks
% ax5 = subplot(2,3,5);
% for aa = 1:size(data.Dural.EvokedAvgs.IntermediateWhisks.means,1)
%     plot(data.Dural.EvokedAvgs.IntermediateWhisks.timeVector,data.Dural.EvokedAvgs.IntermediateWhisks.means(aa,:),'color',colors_eLife2020('rich black'),'LineWidth',0.25);
%     hold on
% end
% for bb = 1:size(data.Capillary.EvokedAvgs.IntermediateWhisks.means,1)
%     plot(data.Capillary.EvokedAvgs.IntermediateWhisks.timeVector,data.Capillary.EvokedAvgs.IntermediateWhisks.means(bb,:),'color',colors_eLife2020('dark candy apple red'),'LineWidth',0.25);
%     hold on
% end
% title('Moderate whisk response')
% ylabel('\DeltaV (mm/sec)')
% xlabel('Peri-whisk time (s)')
% axis square
% xlim([-2,10])
% set(gca,'box','off')
% ax5.TickLength = [0.03,0.03];
% %% [3-S1c] extended whisks
% ax6 = subplot(2,3,6);
% for aa = 1:size(data.Dural.EvokedAvgs.LongWhisks.means,1)
%     plot(data.Dural.EvokedAvgs.LongWhisks.timeVector,data.Dural.EvokedAvgs.LongWhisks.means(aa,:),'color',colors_eLife2020('rich black'),'LineWidth',0.25);
%     hold on
% end
% for bb = 1:size(data.Capillary.EvokedAvgs.LongWhisks.means,1)
%     plot(data.Capillary.EvokedAvgs.LongWhisks.timeVector,data.Capillary.EvokedAvgs.LongWhisks.means(bb,:),'color',colors_eLife2020('dark candy apple red'),'LineWidth',0.25);
%     hold on
% end
% title('Extended whisk response')
% ylabel('\DeltaV (mm/sec)')
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
