function [] = FigS2_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

%% Ephys sleep model accuracy
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_ModelAccuracy_Ephys';
load(resultsStruct);
cd(rootFolder)
% sleep model accuracy using RF and out of bag error
ephysData.physio.holdXlabels = []; ephysData.physio.holdYlabels = []; ephysData.physio.loss = [];
groups = {'Naive','Blank_SAP','SSP_SAP'};
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_ModelAccuracy_Ephys.(group));
    ephysData.(group).holdXlabels = []; ephysData.(group).holdYlabels = []; ephysData.(group).loss = [];
    % extract data from summary structures
    for dd = 1:length(animalIDs)
        animalID = animalIDs{dd,1};
        % physio MDL
        ephysData.(group).holdXlabels = cat(1,ephysData.(group).holdXlabels,Results_ModelAccuracy_Ephys.(group).(animalID).predictedTestingLabels);
        ephysData.(group).holdYlabels = cat(1,ephysData.(group).holdYlabels,Results_ModelAccuracy_Ephys.(group).(animalID).trueTestingLabels);
        ephysData.(group).loss = cat(1,ephysData.(group).loss,Results_ModelAccuracy_Ephys.(group).(animalID).outOfBagError);
    end
    ephysData.(group).meanLoss = mean(ephysData.(group).loss,1);
    ephysData.(group).stdLoss = std(ephysData.(group).loss,0,1);
end

%% GCaMP sleep model accuracy
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_ModelAccuracy_GCaMP';
load(resultsStruct);
cd(rootFolder)
% sleep model accuracy using RF and out of bag error
gcampData.physio.holdXlabels = []; gcampData.physio.holdYlabels = []; gcampData.physio.loss = [];
groups = {'Blank_SAP','SSP_SAP'};
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_ModelAccuracy_GCaMP.(group));
    gcampData.(group).holdXlabels = []; gcampData.(group).holdYlabels = []; gcampData.(group).loss = [];
    % extract data from summary structures
    for dd = 1:length(animalIDs)
        animalID = animalIDs{dd,1};
        % physio MDL
        gcampData.(group).holdXlabels = cat(1,gcampData.(group).holdXlabels,Results_ModelAccuracy_GCaMP.(group).(animalID).predictedTestingLabels);
        gcampData.(group).holdYlabels = cat(1,gcampData.(group).holdYlabels,Results_ModelAccuracy_GCaMP.(group).(animalID).trueTestingLabels);
        gcampData.(group).loss = cat(1,gcampData.(group).loss,Results_ModelAccuracy_GCaMP.(group).(animalID).outOfBagError);
    end
    gcampData.(group).meanLoss = mean(gcampData.(group).loss,1);
    gcampData.(group).stdLoss = std(gcampData.(group).loss,0,1);
end
%% EGFP blue-green
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_GFP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'EGFP'};
hemispheres = {'LH','RH'};
egfpData.blue = [];
egfpData.green = [];
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_GFP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            egfpData.blue = cat(2,egfpData.blue,Results_GFP.(group).(animalID).(hemisphere).blue);
            egfpData.green = cat(2,egfpData.green,Results_GFP.(group).(animalID).(hemisphere).green);
        end
    end
end

%% EGFP cross correlation
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_CrossCorr_EGFP';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'EGFP'};
hemispheres = {'LH','RH'};
behaviors = {'All'};
variables = {'lags','xcVals'};
samplingRate = 30;
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_CrossCorr_EGFP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for ee = 1:length(behaviors)
                behavior = behaviors{1,ee};
                egfpXCorrData.(group).(hemisphere).dummyCheck = 1;
                if isfield(egfpXCorrData.(group).(hemisphere),(behavior)) == false
                    for ff = 1:length(variables)
                        variable = variables{1,ff};
                        egfpXCorrData.(group).(hemisphere).(behavior).(variable) = [];
                    end
                end
                if isempty(Results_CrossCorr_EGFP.(group).(animalID).(hemisphere).(behavior).xcVals) == false
                    egfpXCorrData.(group).(hemisphere).(behavior).lags = cat(1,egfpXCorrData.(group).(hemisphere).(behavior).lags,Results_CrossCorr_EGFP.(group).(animalID).(hemisphere).(behavior).lags/samplingRate);
                    egfpXCorrData.(group).(hemisphere).(behavior).xcVals = cat(1,egfpXCorrData.(group).(hemisphere).(behavior).xcVals,Results_CrossCorr_EGFP.(group).(animalID).(hemisphere).(behavior).xcVals);
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for dd = 1:length(behaviors)
            behavior = behaviors{1,dd};
            egfpXCorrData.(group).(hemisphere).(behavior).mean_xcVals = mean(egfpXCorrData.(group).(hemisphere).(behavior).xcVals,1);
            egfpXCorrData.(group).(hemisphere).(behavior).stdErr_xcVals = std(egfpXCorrData.(group).(hemisphere).(behavior).xcVals,0,1)./sqrt(size(egfpXCorrData.(group).(hemisphere).(behavior).xcVals,1));
            egfpXCorrData.(group).(hemisphere).(behavior).mean_lags = mean(egfpXCorrData.(group).(hemisphere).(behavior).lags,1);
        end
    end
end

%% figure
FigS2 = figure('Name','Figure S2','units','normalized','outerposition',[0 0 1 1]);

% ephys naive confusion chart
subplot(3,4,1)
cm = confusionchart(ephysData.Naive.holdYlabels,ephysData.Naive.holdXlabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = ['Naive total accuracy: ' num2str(modelAccuracy) ' (%)'];

% ephys blank-sap confusion chart
subplot(3,4,2)
cm = confusionchart(ephysData.Blank_SAP.holdYlabels,ephysData.Blank_SAP.holdXlabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = ['Blank-SAP total accuracy: ' num2str(modelAccuracy) ' (%)'];

% ephys ssp-sap confusion chart
subplot(3,4,3)
cm = confusionchart(ephysData.SSP_SAP.holdYlabels,ephysData.SSP_SAP.holdXlabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = ['SSP-SAP total accuracy: ' num2str(modelAccuracy) ' (%)'];

% sleep model OOB-error
subplot(3,4,4);
s1 = scatter(ones(1,length(ephysData.Naive.loss))*1,ephysData.Naive.loss,75,'MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('sapphire'),'jitter','on','jitterAmount',0);
hold on
e1 = errorbar(1,ephysData.Naive.meanLoss,ephysData.Naive.stdLoss,'d','MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('black'));
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(ephysData.Blank_SAP.loss))*2,ephysData.Blank_SAP.loss,75,'MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0);
hold on
e2 = errorbar(2,ephysData.Blank_SAP.meanLoss,ephysData.Blank_SAP.stdLoss,'d','MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('black'));
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(ephysData.SSP_SAP.loss))*3,ephysData.SSP_SAP.loss,75,'MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0);
hold on
e3 = errorbar(3,ephysData.SSP_SAP.meanLoss,ephysData.SSP_SAP.stdLoss,'d','MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('black'));
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

% GCaMP blank-sap confusion chart
subplot(3,3,4)
cm = confusionchart(gcampData.Blank_SAP.holdYlabels,gcampData.Blank_SAP.holdXlabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = ['Blank-SAP total accuracy: ' num2str(modelAccuracy) ' (%)'];

% GCaMP ssp-sap confusion chart
subplot(3,3,5)
cm = confusionchart(ephysData.SSP_SAP.holdYlabels,ephysData.SSP_SAP.holdXlabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = ['SSP-SAP total accuracy: ' num2str(modelAccuracy) ' (%)'];

% sleep model OOB-error
subplot(3,3,6);
s1 = scatter(ones(1,length(ephysData.Blank_SAP.loss))*2,ephysData.Blank_SAP.loss,75,'MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0);
hold on
e1 = errorbar(2,ephysData.Blank_SAP.meanLoss,ephysData.Blank_SAP.stdLoss,'d','MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('black'));
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(ephysData.SSP_SAP.loss))*3,ephysData.SSP_SAP.loss,75,'MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0);
hold on
e2 = errorbar(3,ephysData.SSP_SAP.meanLoss,ephysData.SSP_SAP.stdLoss,'d','MarkerEdgeColor',colors('black'),'MarkerFaceColor',colors('black'));
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('Out of bag error')
legend([s1,s2],'Blank-SAP','SSP-SAP')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
set(gca,'box','off')

subplot(3,2,5)
histogram2(egfpData.blue,egfpData.green,'DisplayStyle','tile','ShowEmptyBins','on')
axis square
set(gca,'box','off')

subplot(3,2,6)
plot(egfpXCorrData.EGFP.RH.All.mean_lags,egfpXCorrData.EGFP.RH.All.mean_xcVals,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(egfpXCorrData.EGFP.RH.All.mean_lags,egfpXCorrData.EGFP.RH.All.mean_xcVals + egfpXCorrData.EGFP.RH.All.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(egfpXCorrData.EGFP.RH.All.mean_lags,egfpXCorrData.EGFP.RH.All.mean_xcVals - egfpXCorrData.EGFP.RH.All.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
title('blue-green XCorr')
ylabel('Corr Coef')
xlabel('Lags (s)')
xlim([-5,5])
set(gca,'box','off')
axis square

%% save figure
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figure Panels' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(FigS2,[dirpath 'FigS2']);
    set(Fig1,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'FigS2'])
end
