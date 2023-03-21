function [] = ArousalStateArea_Pupil(rootFolder,saveFigs,delim)
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
resultsStruct = 'Results_BehavArea.mat';
load(resultsStruct);
animalIDs = fieldnames(Results_BehavArea);
behavFields = {'Rest','Whisk','Stim','NREM','REM'};
% mean HbT comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        data.(animalID).(behavField).meanArea = Results_BehavArea.(animalID).(behavField).meanPupilArea;
        data.(animalID).(behavField).indArea = Results_BehavArea.(animalID).(behavField).indPupilArea;
    end
end
% pre-allocate
for cc = 1:length(behavFields)
    behavField = behavFields{1,cc};
    procData.(behavField).meanArea = [];
    procData.(behavField).indArea = [];
end
% concatenate
for cc = 1:length(animalIDs)
    animalID = animalIDs{cc,1};
    for dd = 1:length(behavFields)
        behavField = behavFields{1,dd};
        procData.(behavField).meanArea = cat(1,procData.(behavField).meanArea,nanmean(data.(animalID).(behavField).meanArea));
        indArea = [];
        for ee = 1:length(data.(animalID).(behavField).indArea)
            indArea = cat(2,indArea,data.(animalID).(behavField).indArea{ee,1});
        end
        procData.(behavField).indArea = cat(2,procData.(behavField).indArea,indArea);
    end
end
% mean/std
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    procData.(behavField).mean = nanmean(procData.(behavField).meanArea,1);
    procData.(behavField).stdev = nanstd(procData.(behavField).meanArea,0,1);
    realIndex = ~isnan(procData.(behavField).indArea);
    procData.(behavField).ind = procData.(behavField).indArea == realIndex;
end
% figures
summaryFigure = figure;
ax1 = subplot(1,2,1);
s1 = scatter(ones(1,length(procData.Rest.meanArea))*1,procData.Rest.meanArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,procData.Rest.mean,procData.Rest.stdev,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(procData.Whisk.meanArea))*2,procData.Whisk.meanArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,procData.Whisk.mean,procData.Whisk.stdev,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(procData.Stim.meanArea))*3,procData.Stim.meanArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,procData.Stim.mean,procData.Stim.stdev,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(ones(1,length(procData.NREM.meanArea))*4,procData.NNREM.meanArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,procData.NREM.mean,procData.NREM.stdev,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
s5 = scatter(ones(1,length(procData.REM.meanArea))*4,procData.REM.meanArea,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(4,procData.REM.mean,procData.REM.stdev,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title({'Pupil area during arousal-states'})
ylabel('Area (pixels)')
legend([s1,s2,s3,s4,s5],'Rest','Whisk','Stim','NREM','REM','Location','NorthEast')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ax2 = subplot(1,2,2);
[xCurve1,yCurve1] = SmoothHistogramBinsFit(procData.Rest.indArea,11,'normal');
[xCurve2,yCurve2] = SmoothHistogramBinsFit(procData.Whisk.indArea,11,'normal');
[xCurve3,yCurve3] = SmoothHistogramBinsFit(procData.Stim.indArea,11,'normal');
[xCurve4,yCurve4] = SmoothHistogramBinsFit(procData.NREM.indArea,11,'normal');
[xCurve5,yCurve5] = SmoothHistogramBinsFit(procData.REM.indArea,11,'normal');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
title({'Pupil area distribution during arousal-states'})
ylabel('Area (pixels)')
ylabel('Probability')
axis square
set(gca,'box','off')
% ylim([0,0.6])
ax2.TickLength = [0.03,0.03];
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Arousal State Hemodynamics - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Arousal_HbT']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Arousal_HbT'])
end

end
