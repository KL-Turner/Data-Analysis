function [] = ModelAccuracy_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_ModelAccuracy_Ephys';
load(resultsStruct);
cd(rootFolder)
% sleep model accuracy using RF and out of bag error
data.physio.holdXlabels = []; data.physio.holdYlabels = []; data.physio.loss = [];
groups = {'Naive','Blank_SAP','SSP_SAP'};
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_ModelAccuracy_Ephys.(group));
    data.(group).holdXlabels = []; data.(group).holdYlabels = []; data.(group).loss = [];
    % extract data from summary structures
    for dd = 1:length(animalIDs)
        animalID = animalIDs{dd,1};
        % physio MDL
        data.(group).holdXlabels = cat(1,data.(group).holdXlabels,Results_ModelAccuracy_Ephys.(group).(animalID).predictedTestingLabels);
        data.(group).holdYlabels = cat(1,data.(group).holdYlabels,Results_ModelAccuracy_Ephys.(group).(animalID).trueTestingLabels);
        data.(group).loss = cat(1,data.(group).loss,Results_ModelAccuracy_Ephys.(group).(animalID).outOfBagError);
    end
    data.(group).meanLoss = mean(data.(group).loss,1);
    data.(group).stdLoss = std(data.(group).loss,0,1);
end
% figure
summaryFigure = figure;
sgtitle('Sleep model accuracy [Ephys]')
for aa = 1:length(groups)
    group = groups{1,aa};
    % sleep model confusion matrix
    subplot(1,4,aa)
    cm = confusionchart(data.(group).holdYlabels,data.(group).holdXlabels);
    cm.ColumnSummary = 'column-normalized';
    cm.RowSummary = 'row-normalized';
    confVals = cm.NormalizedValues;
    totalScores = sum(confVals(:));
    modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
    cm.Title = {strrep(group,'_',' '),['total accuracy: ' num2str(modelAccuracy) ' (%)']};
end
% sleep model 10-fold loss
subplot(1,4,4);
s1 = scatter(ones(1,length(data.Naive.loss))*1,data.Naive.loss,75,'MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('sapphire'),'jitter','on','jitterAmount',0);
hold on
e1 = errorbar(1,data.Naive.meanLoss,data.Naive.stdLoss,'d','MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('black'));
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(data.Blank_SAP.loss))*2,data.Blank_SAP.loss,75,'MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0);
hold on
e2 = errorbar(2,data.Blank_SAP.meanLoss,data.Blank_SAP.stdLoss,'d','MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('black'));
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(data.SSP_SAP.loss))*3,data.SSP_SAP.loss,75,'MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0);
hold on
e3 = errorbar(3,data.SSP_SAP.meanLoss,data.SSP_SAP.stdLoss,'d','MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('black'));
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('Out of bag error')
legend([s1,s2,s3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
set(gca,'box','off')
% save figure
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures' delim 'Model Accuracy' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'ModelAccuracy_Ephys']);
end