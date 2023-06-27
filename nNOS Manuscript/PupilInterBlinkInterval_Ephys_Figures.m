function [] = PupilInterBlinkInterval_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PupilBlinkInterval_Ephys';
load(resultsStruct);
cd(rootFolder)
groups = {'Naive','SSP_SAP','Blank_SAP'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    data.(group).data = [];
    animalIDs = fieldnames(Results_PupilBlinkInterval_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        % concatenate data across animals
        data.(group).data = cat(1,data.(group).data,mean(Results_PupilBlinkInterval_Ephys.(group).(animalID).interBlinkInterval));
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    data.(group).meanData = mean(data.(group).data,1);
    data.(group).stdData = std(data.(group).data,0,1);
end
% figure
summaryFigure = figure;
sgtitle('Inter-blink interval [Ephys]')
xInds = ones(1,length(data.Naive.data));
s1 = scatter(xInds*1,data.Naive.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off','jitterAmount',0.25);
hold on;
e1 = errorbar(1,data.Naive.meanData,data.Naive.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(data.Blank_SAP.data));
s2 = scatter(xInds*2,data.Blank_SAP.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,data.Blank_SAP.meanData,data.Blank_SAP.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(data.SSP_SAP.data));
s3 = scatter(xInds*3,data.SSP_SAP.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(3,data.SSP_SAP.meanData,data.SSP_SAP.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('IBI (s)')
xlim([0,4])
legend([s1,s2,s3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures' delim 'Behavior' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Pupil_IBI_Ephys']);
end