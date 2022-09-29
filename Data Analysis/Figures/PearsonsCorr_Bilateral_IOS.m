function [] = PearsonsCorr_Bilateral_IOS(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
% colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
colorRest = [(0/256),(166/256),(81/256)];
colorWhisk = [(31/256),(120/256),(179/256)];
% colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
colorAlert = [(255/256),(191/256),(0/256)];
colorAsleep = [(0/256),(128/256),(255/256)];
colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
resultsStruct = 'Results_PearsonCorr';
load(resultsStruct);
expGroups = {'Naive','SSP-SAP','Blank-SAP'};
setName = 'IOS Set A';
animalIDs.all = {};
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs.all = horzcat(animalIDs.all,{folderList.name});
    animalIDs.(strrep(expGroups{1,aa},'-','_')) = {folderList.name};
end
treatments = {'Naive','SSP_SAP','Blank_SAP'};
behavFields = {'Rest','Whisk','NREM','REM','Awake','Sleep','All'};
dataTypes = {'CBV_HbT','gammaBandPower'};
%% Pearson's correlations during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.Naive = []; data.SSP_SAP = []; data.Blank_SAP = [];
for aa = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,aa};
    % recognize treatment based on animal group
    if ismember(animalID,animalIDs.Naive) == true
        treatment = 'Naive';
    elseif ismember(animalID,animalIDs.SSP_SAP) == true
        treatment = 'SSP_SAP';
    elseif ismember(animalID,animalIDs.Blank_SAP) == true
        treatment = 'Blank_SAP';
    end
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        % create the behavior folder for the first iteration of the loop
        if isfield(data.(treatment),behavField) == false
            data.(treatment).(behavField) = [];
        end
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(Results_PearsonCorr.(animalID).(behavField).(dataType).meanR) == false
                % create the data type folder for the first iteration of the loop
                if isfield(data.(treatment).(behavField),dataType) == false
                    data.(treatment).(behavField).(dataType).meanRs = [];
                    data.(treatment).(behavField).(dataType).animalID = {};
                    data.(treatment).(behavField).(dataType).treatment = {};
                    data.(treatment).(behavField).(dataType).behavior = {};
                end
                % concatenate mean R and the animalID/behavior
                data.(treatment).(behavField).(dataType).meanRs = cat(1,data.(treatment).(behavField).(dataType).meanRs,Results_PearsonCorr.(animalID).(behavField).(dataType).meanR);
                data.(treatment).(behavField).(dataType).animalID = cat(1,data.(treatment).(behavField).(dataType).animalID,animalID);
                data.(treatment).(behavField).(dataType).behavior = cat(1,data.(treatment).(behavField).(dataType).behavior,behavField);
                data.(treatment).(behavField).(dataType).treatment = cat(1,data.(treatment).(behavField).(dataType).treatment,treatment);
            end
        end
    end
end
%% take mean/STD of R
for qq = 1:length(treatments)
    treatment = treatments{1,qq};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        for ff = 1:length(dataTypes)
            dataType = dataTypes{1,ff};
            data.(treatment).(behavField).(dataType).meanR = mean(data.(treatment).(behavField).(dataType).meanRs,1);
            data.(treatment).(behavField).(dataType).stdR = std(data.(treatment).(behavField).(dataType).meanRs,0,1);
        end
    end
end
%% statistics - generalized linear mixed effects model
for qq = 1:length(dataTypes)
    dataType = dataTypes{1,qq};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        % statistics - generalized linear mixed effects model
        Stats.(dataType).(behavField).tableSize = cat(1,data.Blank_SAP.(behavField).(dataType).meanRs,data.SSP_SAP.(behavField).(dataType).meanRs);
        Stats.(dataType).(behavField).Table = table('Size',[size(Stats.(dataType).(behavField).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Treatment','Correlation'});
        Stats.(dataType).(behavField).Table.Mouse = cat(1,data.Blank_SAP.(behavField).(dataType).animalID,data.SSP_SAP.(behavField).(dataType).animalID);
        Stats.(dataType).(behavField).Table.Treatment = cat(1,data.Blank_SAP.(behavField).(dataType).treatment,data.SSP_SAP.(behavField).(dataType).treatment);
        Stats.(dataType).(behavField).Table.Correlation = cat(1,data.Blank_SAP.(behavField).(dataType).meanRs,data.SSP_SAP.(behavField).(dataType).meanRs);
        Stats.(dataType).(behavField).FitFormula = 'Correlation ~ 1 + Treatment + (1|Mouse)';
        Stats.(dataType).(behavField).Stats = fitglme(Stats.(dataType).(behavField).Table,Stats.(dataType).(behavField).FitFormula);
    end
end
%% average correlation coefficient figure
summaryFigure1 = figure;
CC_xInds1A = ones(1,length(data.Naive.Rest.CBV_HbT.meanRs));
CC_xInds1B = ones(1,length(data.SSP_SAP.Rest.CBV_HbT.meanRs));
CC_xInds1C = ones(1,length(data.Blank_SAP.Rest.CBV_HbT.meanRs));
CC_xInds2A = ones(1,length(data.Naive.Awake.CBV_HbT.animalID));
CC_xInds2B = ones(1,length(data.SSP_SAP.Awake.CBV_HbT.animalID));
CC_xInds2C = ones(1,length(data.Blank_SAP.Awake.CBV_HbT.animalID));
CC_xInds3A = ones(1,length(data.Naive.Sleep.CBV_HbT.animalID));
CC_xInds3B = ones(1,length(data.SSP_SAP.Sleep.CBV_HbT.animalID));
CC_xInds3C = ones(1,length(data.Blank_SAP.Sleep.CBV_HbT.animalID));
CC_xInds = ones(1,length(data.Naive.Awake.gammaBandPower.animalID));
%% Pearson's correlations between bilateral HbT during different Rest
b1 = scatter(1,data.Naive.Rest.CBV_HbT.meanR,'d','MarkerFaceColor',colors('sapphire'));
hold on
scatter(CC_xInds1A*1,data.Naive.Rest.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
b2 = scatter(2,data.SSP_SAP.Rest.CBV_HbT.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds1B*2,data.SSP_SAP.Rest.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
b3 = scatter(3,data.Blank_SAP.Rest.CBV_HbT.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds1C*3,data.Blank_SAP.Rest.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
%% Pearson's correlations between bilateral HbT during different Whisk
scatter(5,data.Naive.Whisk.CBV_HbT.meanR,'d','MarkerFaceColor',colors('sapphire'));
scatter(CC_xInds1A*5,data.Naive.Whisk.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on', 'jitterAmount',0.25);
scatter(6,data.SSP_SAP.Whisk.CBV_HbT.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds1B*6,data.SSP_SAP.Whisk.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on', 'jitterAmount',0.25);
scatter(7,data.Blank_SAP.Whisk.CBV_HbT.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds1C*7,data.Blank_SAP.Whisk.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on', 'jitterAmount',0.25);
%% Pearson's correlations between bilateral HbT during different Rest
scatter(9,data.Naive.NREM.CBV_HbT.meanR,'d','MarkerFaceColor',colors('sapphire'));
scatter(CC_xInds1A*9,data.Naive.NREM.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
scatter(10,data.SSP_SAP.NREM.CBV_HbT.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds1B*10,data.SSP_SAP.NREM.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
scatter(11,data.Blank_SAP.NREM.CBV_HbT.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds1C*11,data.Blank_SAP.NREM.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
%% Pearson's correlations between bilateral HbT during different Rest
scatter(13,data.Naive.REM.CBV_HbT.meanR,'d','MarkerFaceColor',colors('sapphire'));
scatter(CC_xInds1A*13,data.Naive.REM.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
scatter(14,data.SSP_SAP.REM.CBV_HbT.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds1B*14,data.SSP_SAP.REM.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
scatter(15,data.Blank_SAP.REM.CBV_HbT.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds1C*15,data.Blank_SAP.REM.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
%% Pearson's correlations between bilateral HbT during different Rest
scatter(17,data.Naive.Awake.CBV_HbT.meanR,'d','MarkerFaceColor',colors('sapphire'));
scatter(CC_xInds2A*17,data.Naive.Awake.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
scatter(18,data.SSP_SAP.Awake.CBV_HbT.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds2B*18,data.SSP_SAP.Awake.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
scatter(19,data.Blank_SAP.Awake.CBV_HbT.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds2C*19,data.Blank_SAP.Awake.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
%% Pearson's correlations between bilateral HbT during different Rest
scatter(21,data.Naive.Sleep.CBV_HbT.meanR,'d','MarkerFaceColor',colors('sapphire'));
scatter(CC_xInds3A*21,data.Naive.Sleep.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
scatter(22,data.SSP_SAP.Sleep.CBV_HbT.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds3B*22,data.SSP_SAP.Sleep.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
scatter(23,data.Blank_SAP.Sleep.CBV_HbT.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds3C*23,data.Blank_SAP.Sleep.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
%% Pearson's correlations between bilateral HbT during different Rest
scatter(25,data.Naive.All.CBV_HbT.meanR,'d','MarkerFaceColor',colors('sapphire'));
scatter(CC_xInds1A*25,data.Naive.All.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
scatter(26,data.SSP_SAP.All.CBV_HbT.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds1B*26,data.SSP_SAP.All.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
scatter(27,data.Blank_SAP.All.CBV_HbT.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds1C*27,data.Blank_SAP.All.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
%% stat lines
plot([2,3],[1,1],'k');
text(2.5,1,'**','FontSize',16)
plot([6,7],[1,1],'k');
text(6.5,1,'*','FontSize',16)
plot([10,11],[1,1],'k');
text(10.5,1,'*','FontSize',16)
plot([14,15],[1,1],'k');
text(14.5,1,'ns','FontSize',16)
plot([18,19],[1,1],'k');
text(18.5,1,'ns','FontSize',16)
plot([22,23],[1,1],'k');
text(22.5,1,'**','FontSize',16)
plot([26,27],[1,1],'k');
text(26.5,1,'*','FontSize',16)
%% figure characteristics
title({'Cortical Pearson''s corr. coef','\DeltaHbT \muM (%)'})
ylabel('Corr. coefficient')
set(gca,'xtick',[2,6,10,14,18,22,26])
set(gca,'xticklabel',{'Rest','Whisk','NREM','REM','Alert','Asleep','All'})
xtickangle(45)
axis square
xlim([0,28])
% ylim([0,1])
set(gca,'box','off')
legend([b1,b2,b3],'Naive','SSP-SAP','Blank-SAP','Location','NorthWest')
axis square
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Pearsons Correlation - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'PearsonsCorrelation_HbT']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'PearsonsCorrelation_HbT'])
end
%% statistical diary
diaryFile = [dirpath 'PearsonsCorrelation_Statistics_HbT.txt'];
if exist(diaryFile,'file') == 2
    delete(diaryFile)
end
diary(diaryFile)
diary on
disp('======================================================================================================================')
disp('GLME statistics for HbT Pearsons Correlations for Rest data')
disp('======================================================================================================================')
disp(Stats.CBV_HbT.Rest.Stats)
disp('======================================================================================================================')
disp('GLME statistics for HbT Pearsons Correlations for Whisk data')
disp('======================================================================================================================')
disp(Stats.CBV_HbT.Whisk.Stats)
disp('======================================================================================================================')
disp('GLME statistics for HbT Pearsons Correlations for NREM data')
disp('======================================================================================================================')
disp(Stats.CBV_HbT.NREM.Stats)
disp('======================================================================================================================')
disp('GLME statistics for HbT Pearsons Correlations for REM data')
disp('======================================================================================================================')
disp(Stats.CBV_HbT.REM.Stats)
disp('======================================================================================================================')
disp('GLME statistics for HbT Pearsons Correlations for Alert data')
disp('======================================================================================================================')
disp(Stats.CBV_HbT.Awake.Stats)
disp('======================================================================================================================')
disp('GLME statistics for HbT Pearsons Correlations for Asleep data')
disp('======================================================================================================================')
disp(Stats.CBV_HbT.Sleep.Stats)
disp('======================================================================================================================')
disp('GLME statistics for HbT Pearsons Correlations for All data')
disp('======================================================================================================================')
disp(Stats.CBV_HbT.All.Stats)
%% average correlation coefficient figure gammma-band
summaryFigure2 = figure;
%% Pearson's correlations between bilateral gamma-band during different Rest
b1 = scatter(1,data.Naive.Rest.gammaBandPower.meanR,'d','MarkerFaceColor',colors('sapphire'));
hold on
scatter(CC_xInds1A*1,data.Naive.Rest.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
b2 = scatter(2,data.SSP_SAP.Rest.gammaBandPower.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds1B*2,data.SSP_SAP.Rest.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
b3 = scatter(3,data.Blank_SAP.Rest.gammaBandPower.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds1C*3,data.Blank_SAP.Rest.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
%% Pearson's correlations between bilateral gamma-band during different Whisk
scatter(5,data.Naive.Whisk.gammaBandPower.meanR,'d','MarkerFaceColor',colors('sapphire'));
scatter(CC_xInds1A*5,data.Naive.Whisk.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on', 'jitterAmount',0.25);
scatter(6,data.SSP_SAP.Whisk.gammaBandPower.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds1B*6,data.SSP_SAP.Whisk.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on', 'jitterAmount',0.25);
scatter(7,data.Blank_SAP.Whisk.gammaBandPower.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds1C*7,data.Blank_SAP.Whisk.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on', 'jitterAmount',0.25);
%% Pearson's correlations between bilateral gamma-band during different Rest
scatter(9,data.Naive.NREM.gammaBandPower.meanR,'d','MarkerFaceColor',colors('sapphire'));
scatter(CC_xInds1A*9,data.Naive.NREM.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
scatter(10,data.SSP_SAP.NREM.gammaBandPower.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds1B*10,data.SSP_SAP.NREM.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
scatter(11,data.Blank_SAP.NREM.gammaBandPower.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds1C*11,data.Blank_SAP.NREM.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
%% Pearson's correlations between bilateral gamma-band during different Rest
scatter(13,data.Naive.REM.gammaBandPower.meanR,'d','MarkerFaceColor',colors('sapphire'));
scatter(CC_xInds1A*13,data.Naive.REM.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
scatter(14,data.SSP_SAP.REM.gammaBandPower.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds1B*14,data.SSP_SAP.REM.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
scatter(15,data.Blank_SAP.REM.gammaBandPower.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds1C*15,data.Blank_SAP.REM.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
%% Pearson's correlations between bilateral gamma-band during different Rest
scatter(17,data.Naive.Awake.gammaBandPower.meanR,'d','MarkerFaceColor',colors('sapphire'));
scatter(CC_xInds*17,data.Naive.Awake.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
scatter(18,data.SSP_SAP.Awake.gammaBandPower.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds2B*18,data.SSP_SAP.Awake.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
scatter(19,data.Blank_SAP.Awake.gammaBandPower.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds2B*19,data.Blank_SAP.Awake.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
%% Pearson's correlations between bilateral gamma-band during different Rest
scatter(21,data.Naive.Sleep.gammaBandPower.meanR,'d','MarkerFaceColor',colors('sapphire'));
scatter(CC_xInds*21,data.Naive.Sleep.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
scatter(22,data.SSP_SAP.Sleep.gammaBandPower.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds3C*22,data.SSP_SAP.Sleep.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
scatter(23,data.Blank_SAP.Sleep.gammaBandPower.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds3C*23,data.Blank_SAP.Sleep.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
%% Pearson's correlations between bilateral gamma-band during different Rest
scatter(25,data.Naive.All.gammaBandPower.meanR,'d','MarkerFaceColor',colors('sapphire'));
scatter(CC_xInds1A*25,data.Naive.All.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
scatter(26,data.SSP_SAP.All.gammaBandPower.meanR,'d','MarkerFaceColor',colors('electric purple'));
scatter(CC_xInds1B*26,data.SSP_SAP.All.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
scatter(27,data.Blank_SAP.All.gammaBandPower.meanR,'d','MarkerFaceColor',colors('north texas green'));
scatter(CC_xInds1C*27,data.Blank_SAP.All.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
%% stat lines
plot([2,3],[1,1],'k');
text(2.5,1,'ns','FontSize',16)
plot([6,7],[1,1],'k');
text(6.5,1,'ns','FontSize',16)
plot([10,11],[1,1],'k');
text(10.5,1,'*','FontSize',16)
plot([14,15],[1,1],'k');
text(14.5,1,'ns','FontSize',16)
plot([18,19],[1,1],'k');
text(18.5,1,'ns','FontSize',16)
plot([22,23],[1,1],'k');
text(22.5,1,'*','FontSize',16)
plot([26,27],[1,1],'k');
text(26.5,1,'ns','FontSize',16)
%% figure characteristics
title({'Cortical Pearson''s corr. coef','Gamma-band [30-100] Hz'})
ylabel('Corr. coefficient')
set(gca,'xtick',[2,6,10,14,18,22,26])
set(gca,'xticklabel',{'Rest','Whisk','NREM','REM','Alert','Asleep','All'})
xtickangle(45)
axis square
xlim([0,28])
% ylim([0,1])
set(gca,'box','off')
legend([b1,b2,b3],'Naive','SSP-SAP','Blank-SAP','Location','NorthWest')
axis square
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Pearsons Correlation - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'PearsonsCorrelation_Gamma']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'PearsonsCorrelation_Gamma'])
end
%% statistical diary
diaryFile = [dirpath 'PearsonsCorrelation_Statistics_gammaBandPower.txt'];
if exist(diaryFile,'file') == 2
    delete(diaryFile)
end
diary(diaryFile)
diary on
disp('======================================================================================================================')
disp('GLME statistics for gamma-band Pearsons Correlations for Rest data')
disp('======================================================================================================================')
disp(Stats.gammaBandPower.Rest.Stats)
disp('======================================================================================================================')
disp('GLME statistics for gamma-band Pearsons Correlations for Whisk data')
disp('======================================================================================================================')
disp(Stats.gammaBandPower.Whisk.Stats)
disp('======================================================================================================================')
disp('GLME statistics for gamma-band Pearsons Correlations for NREM data')
disp('======================================================================================================================')
disp(Stats.gammaBandPower.NREM.Stats)
disp('======================================================================================================================')
disp('GLME statistics for gamma-band Pearsons Correlations for REM data')
disp('======================================================================================================================')
disp(Stats.gammaBandPower.REM.Stats)
disp('======================================================================================================================')
disp('GLME statistics for gamma-band Pearsons Correlations for Alert data')
disp('======================================================================================================================')
disp(Stats.gammaBandPower.Awake.Stats)
disp('======================================================================================================================')
disp('GLME statistics for gamma-band Pearsons Correlations for Asleep data')
disp('======================================================================================================================')
disp(Stats.gammaBandPower.Sleep.Stats)
disp('======================================================================================================================')
disp('GLME statistics for gamma-band Pearsons Correlations for All data')
disp('======================================================================================================================')
disp(Stats.gammaBandPower.All.Stats)

end
