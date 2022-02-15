function [] = ArousalStateDiameter_Pupil(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

% behavior colors
colorRest = [(0/256),(166/256),(81/256)];
colorWhisk = [(31/256),(120/256),(179/256)];
colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
% set-up and process data
resultsStruct = 'Results_BehavData.mat';
load(resultsStruct);
animalIDs = fieldnames(Results_BehavData);
behavFields = {'Rest','Whisk','Stim','NREM','REM'};
% mean HbT comparison between behaviors
% pre-allocate the date for each day
% pre-allocate
for cc = 1:length(behavFields)
    behavField = behavFields{1,cc};
    data.(behavField).indMeanArea = [];
    data.(behavField).indArea = [];
    data.(behavField).indMeanDiameter = [];
    data.(behavField).indDiameter = [];
    data.(behavField).indMeanzArea = [];
    data.(behavField).indzArea = [];
    data.(behavField).indMeanzDiameter = [];
    data.(behavField).indzDiameter = [];
end
% concatenate
for cc = 1:length(animalIDs)
    animalID = animalIDs{cc,1};
    for dd = 1:length(behavFields)
        behavField = behavFields{1,dd};
        data.(behavField).indMeanArea = cat(1,data.(behavField).indMeanArea,mean(Results_BehavData.(animalID).(behavField).mmArea.eventMeans,'omitnan'));
        data.(behavField).indMeanDiameter = cat(1,data.(behavField).indMeanDiameter,mean(Results_BehavData.(animalID).(behavField).mmDiameter.eventMeans,'omitnan'));
        data.(behavField).indMeanzArea = cat(1,data.(behavField).indMeanzArea,mean(Results_BehavData.(animalID).(behavField).zArea.eventMeans,'omitnan'));
        data.(behavField).indMeanzDiameter = cat(1,data.(behavField).indMeanzDiameter,mean(Results_BehavData.(animalID).(behavField).zDiameter.eventMeans,'omitnan'));
        indArea = []; indDiameter = []; indzArea = []; indzDiameter = [];
        for ee = 1:length(Results_BehavData.(animalID).(behavField).mmArea.indData)
            indArea = cat(2,indArea,Results_BehavData.(animalID).(behavField).mmArea.indData{ee,1});
            indDiameter = cat(2,indDiameter,Results_BehavData.(animalID).(behavField).mmDiameter.indData{ee,1});
            indzArea = cat(2,indzArea,Results_BehavData.(animalID).(behavField).zArea.indData{ee,1});
            indzDiameter = cat(2,indzDiameter,Results_BehavData.(animalID).(behavField).zDiameter.indData{ee,1});
        end
        data.(behavField).indArea = cat(2,data.(behavField).indArea,indArea);
        data.(behavField).indDiameter = cat(2,data.(behavField).indDiameter,indDiameter);
        data.(behavField).indzArea = cat(2,data.(behavField).indzArea,indzArea);
        data.(behavField).indzDiameter = cat(2,data.(behavField).indzDiameter,indzDiameter);
    end
end
% mean/std
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    % area
    data.(behavField).meanArea = mean(data.(behavField).indMeanArea,1,'omitnan');
    data.(behavField).stdArea = std(data.(behavField).indMeanArea,0,1,'omitnan');
    realIndex = ~isnan(data.(behavField).indArea);
    data.(behavField).indArea = data.(behavField).indArea(realIndex);
    % diameter
    data.(behavField).meanDiameter = mean(data.(behavField).indMeanDiameter,1,'omitnan');
    data.(behavField).stdDiameter = std(data.(behavField).indMeanDiameter,0,1,'omitnan');
    realIndex = ~isnan(data.(behavField).indDiameter);
    data.(behavField).indDiameter = data.(behavField).indDiameter(realIndex);
    % z area
    data.(behavField).meanzArea = mean(data.(behavField).indMeanzArea,1,'omitnan');
    data.(behavField).stdzArea = std(data.(behavField).indMeanzArea,0,1,'omitnan');
    realIndex = ~isnan(data.(behavField).indzArea);
    data.(behavField).indzArea = data.(behavField).indzArea(realIndex);
    % z diameter
    data.(behavField).meanzDiameter = mean(data.(behavField).indMeanzDiameter,1,'omitnan');
    data.(behavField).stdzDiameter = std(data.(behavField).indMeanzDiameter,0,1,'omitnan');
    realIndex = ~isnan(data.(behavField).indzDiameter);
    data.(behavField).indzDiameter = data.(behavField).indzDiameter(realIndex);
end
%% figures
summaryFigure = figure;
sgtitle('Pupil area/diameter during different arousal-states')
%% mm pupil area scatter
ax1 = subplot(2,4,1);
s1 = scatter(ones(1,length(data.Rest.indMeanArea))*1,data.Rest.indMeanArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.meanArea,data.Rest.stdArea,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(data.Whisk.indMeanArea))*2,data.Whisk.indMeanArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Whisk.meanArea,data.Whisk.stdArea,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(data.Stim.indMeanArea))*3,data.Stim.indMeanArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.Stim.meanArea,data.Stim.stdArea,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(ones(1,length(data.NREM.indMeanArea))*4,data.NREM.indMeanArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.NREM.meanArea,data.NREM.stdArea,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
s5 = scatter(ones(1,length(data.REM.indMeanArea))*5,data.REM.indMeanArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.REM.meanArea,data.REM.stdArea,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title('mm Area')
ylabel('Area (mm^2)')
legend([s1,s2,s3,s4,s5],'Rest','Whisk','Stim','NREM','REM','Location','NorthEast')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% mm pupil area histogram
ax5 = subplot(2,4,5);
[xCurve1,yCurve1] = SmoothHistogramBinsFit(data.Rest.indArea,11,'kernel');
[xCurve2,yCurve2] = SmoothHistogramBinsFit(data.Whisk.indArea,11,'kernel');
[xCurve3,yCurve3] = SmoothHistogramBinsFit(data.Stim.indArea,11,'kernel');
[xCurve4,yCurve4] = SmoothHistogramBinsFit(data.NREM.indArea,31,'kernel');
[xCurve5,yCurve5] = SmoothHistogramBinsFit(data.REM.indArea,11,'kernel');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
title('mm area')
ylabel('Area (mm^2)')
ylabel('Probability')
axis square
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% mm pupil diameter scatter
ax2 = subplot(2,4,2);
scatter(ones(1,length(data.Rest.indMeanDiameter))*1,data.Rest.indMeanDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.meanDiameter,data.Rest.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Whisk.indMeanDiameter))*2,data.Whisk.indMeanDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Whisk.meanDiameter,data.Whisk.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Stim.indMeanDiameter))*3,data.Stim.indMeanDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.Stim.meanDiameter,data.Stim.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.NREM.indMeanDiameter))*4,data.NREM.indMeanDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.NREM.meanDiameter,data.NREM.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.REM.indMeanDiameter))*5,data.REM.indMeanDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.REM.meanDiameter,data.REM.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title('mm Diameter')
ylabel('Diameter (mm)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% mm pupil diameter histogram
ax6 = subplot(2,4,6);
[xCurve1,yCurve1] = SmoothHistogramBinsFit(data.Rest.indDiameter,11,'kernel');
[xCurve2,yCurve2] = SmoothHistogramBinsFit(data.Whisk.indDiameter,11,'kernel');
[xCurve3,yCurve3] = SmoothHistogramBinsFit(data.Stim.indDiameter,11,'kernel');
[xCurve4,yCurve4] = SmoothHistogramBinsFit(data.NREM.indDiameter,31,'kernel');
[xCurve5,yCurve5] = SmoothHistogramBinsFit(data.REM.indDiameter,11,'kernel');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
title('mm diameter')
ylabel('Diameter (mm)')
ylabel('Probability')
axis square
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% mm pupil area scatter
ax3 = subplot(2,4,3);
scatter(ones(1,length(data.Rest.indMeanzArea))*1,data.Rest.indMeanzArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.meanzArea,data.Rest.stdzArea,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Whisk.indMeanzArea))*2,data.Whisk.indMeanzArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Whisk.meanzArea,data.Whisk.stdzArea,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Stim.indMeanzArea))*3,data.Stim.indMeanzArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.Stim.meanzArea,data.Stim.stdzArea,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.NREM.indMeanzArea))*4,data.NREM.indMeanzArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.NREM.meanzArea,data.NREM.stdzArea,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.REM.indMeanzArea))*5,data.REM.indMeanzArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.REM.meanzArea,data.REM.stdzArea,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title('Z-scored area')
ylabel('Z units')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% mm pupil area histogram
ax7 = subplot(2,4,7);
[xCurve1,yCurve1] = SmoothHistogramBinsFit(data.Rest.indzArea,11,'kernel');
[xCurve2,yCurve2] = SmoothHistogramBinsFit(data.Whisk.indzArea,11,'kernel');
[xCurve3,yCurve3] = SmoothHistogramBinsFit(data.Stim.indzArea,11,'kernel');
[xCurve4,yCurve4] = SmoothHistogramBinsFit(data.NREM.indzArea,31,'kernel');
[xCurve5,yCurve5] = SmoothHistogramBinsFit(data.REM.indzArea,11,'kernel');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
title('Z-scored area')
ylabel('Z units')
ylabel('Probability')
axis square
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% mm pupil diameter scatter
ax4 = subplot(2,4,4);
scatter(ones(1,length(data.Rest.indMeanzDiameter))*1,data.Rest.indMeanzDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.meanzDiameter,data.Rest.stdzDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Whisk.indMeanzDiameter))*2,data.Whisk.indMeanzDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Whisk.meanzDiameter,data.Whisk.stdzDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Stim.indMeanzDiameter))*3,data.Stim.indMeanzDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.Stim.meanzDiameter,data.Stim.stdzDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.NREM.indMeanzDiameter))*4,data.NREM.indMeanzDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.NREM.meanzDiameter,data.NREM.stdzDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.REM.indMeanzDiameter))*5,data.REM.indMeanzDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.REM.meanzDiameter,data.REM.stdzDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title('Z-scored diameter')
ylabel('Z units')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% mm pupil diameter histogram
ax8 = subplot(2,4,8);
[xCurve1,yCurve1] = SmoothHistogramBinsFit(data.Rest.indzDiameter,11,'kernel');
[xCurve2,yCurve2] = SmoothHistogramBinsFit(data.Whisk.indzDiameter,11,'kernel');
[xCurve3,yCurve3] = SmoothHistogramBinsFit(data.Stim.indzDiameter,11,'kernel');
[xCurve4,yCurve4] = SmoothHistogramBinsFit(data.NREM.indzDiameter,31,'kernel');
[xCurve5,yCurve5] = SmoothHistogramBinsFit(data.REM.indzDiameter,11,'kernel');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
title('Z-scored diameter')
ylabel('Z units')
ylabel('Probability')
axis square
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Arousal State Pupil Data' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Pupil_BehavData']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Pupil_BehavData'])
end

end
