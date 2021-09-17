function [] = WhiskingBehavior_Bilateral_IOS(rootFolder,saveFigs,delim)
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
resultsStruct = 'Results_WhiskBehav';
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
    % don't concatenate empty arrays where there was no data for this behavior
    % create the data type folder for the first iteration of the loop
    if isfield(data.(treatment),'whiskDurations') == false
        data.(treatment).whiskDurations = [];
        data.(treatment).whiskDurationSec = [];
        data.(treatment).whiskDurationPerc = [];
        data.(treatment).treatment = {};
    end
    % concatenate mean R and the animalID/behavior
    data.(treatment).whiskDurations = cat(1,data.(treatment).whiskDurations,Results_WhiskBehav.(animalID).whiskDurations);
    data.(treatment).whiskDurationSec = cat(1,data.(treatment).whiskDurationSec,Results_WhiskBehav.(animalID).whiskDurationSec/60);
    data.(treatment).whiskDurationPerc = cat(1,data.(treatment).whiskDurationPerc,Results_WhiskBehav.(animalID).whiskDurationPerc);
    data.(treatment).treatment = cat(1,data.(treatment).treatment,treatment);
end
%% take mean/STD of R
for qq = 1:length(treatments)
    treatment = treatments{1,qq};
    data.(treatment).meanWhiskDurationSec = mean(data.(treatment).whiskDurationSec,1);
    data.(treatment).stdWhiskDurationSec = std(data.(treatment).whiskDurationSec,0,1);
    data.(treatment).meanWhiskDurationPerc = mean(data.(treatment).whiskDurationPerc,1);
    data.(treatment).stdWhiskDurationPerc = std(data.(treatment).whiskDurationPerc,0,1);
end
%% whisking time and percent stats
Stats.CompA.tableSize = cat(1,data.Blank_SAP.whiskDurationSec,data.SSP_SAP.whiskDurationSec,data.Naive.whiskDurationSec);
Stats.CompA.Table = table('Size',[size(Stats.CompA.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Treatment','Time'});
Stats.CompA.Table.Treatment = cat(1,data.Blank_SAP.treatment,data.SSP_SAP.treatment,data.Naive.treatment);
Stats.CompA.Table.Time = cat(1,data.Blank_SAP.whiskDurationSec,data.SSP_SAP.whiskDurationSec,data.Naive.whiskDurationSec);
Stats.CompA.FitFormula = 'Time ~ 1 + Treatment';
Stats.CompA.Stats = fitglme(Stats.CompA.Table,Stats.CompA.FitFormula);
% percent
Stats.CompB.tableSize = cat(1,data.Blank_SAP.whiskDurationPerc,data.SSP_SAP.whiskDurationPerc,data.Naive.whiskDurationPerc);
Stats.CompB.Table = table('Size',[size(Stats.CompB.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Treatment','Percent'});
Stats.CompB.Table.Treatment = cat(1,data.Blank_SAP.treatment,data.SSP_SAP.treatment,data.Naive.treatment);
Stats.CompB.Table.Percent = cat(1,data.Blank_SAP.whiskDurationPerc,data.SSP_SAP.whiskDurationPerc,data.Naive.whiskDurationPerc);
Stats.CompB.FitFormula = 'Percent ~ 1 + Treatment';
Stats.CompB.Stats = fitglme(Stats.CompB.Table,Stats.CompB.FitFormula);
%% whisking histogram
summaryFigure = figure;
edges = 0:0.25:10;
ax1 = subplot(2,3,1);
histogram(data.Naive.whiskDurations,edges,'Normalization','Probability','EdgeColor','k','FaceColor',colors('sapphire'))
xlabel('Whisk duration (s)') 
ylabel('Probabilty')
title('Naive whisking duration')
set(gca,'box','off')
axis square
ax2 = subplot(2,3,2);
histogram(data.SSP_SAP.whiskDurations,edges,'Normalization','Probability','EdgeColor','k','FaceColor',colors('electric purple'))
xlabel('Whisk duration (s)') 
ylabel('Probabilty')
title('SSP-SAP whisking duration')
set(gca,'box','off')
axis square
ax3 = subplot(2,3,3);
histogram(data.Blank_SAP.whiskDurations,edges,'Normalization','Probability','EdgeColor','k','FaceColor',colors('north texas green'))
xlabel('Whisk duration (s)') 
ylabel('Probabilty')
title('Blank-SAP whisking duration')
set(gca,'box','off')
axis square
% figure characteristics
linkaxes([ax1,ax2,ax3],'xy')
%% whisking duration minutes/percent
subplot(2,3,4)
s1 = scatter(ones(1,length(data.Naive.whiskDurationSec))*1,data.Naive.whiskDurationSec,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Naive.meanWhiskDurationSec,data.Naive.stdWhiskDurationSec,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(data.SSP_SAP.whiskDurationSec))*2,data.SSP_SAP.whiskDurationSec,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,data.SSP_SAP.meanWhiskDurationSec,data.SSP_SAP.stdWhiskDurationSec,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(data.Blank_SAP.whiskDurationSec))*3,data.Blank_SAP.whiskDurationSec,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,data.Blank_SAP.meanWhiskDurationSec,data.Blank_SAP.stdWhiskDurationSec,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel('Whisking time (min)')
title('Total time spent whisking')
legend([s1,s2,s3],'Naive','SSP-SAP','Blank-SAP')
set(gca,'box','off')
axis square
xlim([0,4])
subplot(2,3,5)
s1 = scatter(ones(1,length(data.Naive.whiskDurationPerc))*1,data.Naive.whiskDurationPerc,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Naive.meanWhiskDurationPerc,data.Naive.stdWhiskDurationPerc,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(data.SSP_SAP.whiskDurationPerc))*2,data.SSP_SAP.whiskDurationPerc,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,data.SSP_SAP.meanWhiskDurationPerc,data.SSP_SAP.stdWhiskDurationPerc,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(data.Blank_SAP.whiskDurationPerc))*3,data.Blank_SAP.whiskDurationPerc,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,data.Blank_SAP.meanWhiskDurationPerc,data.Blank_SAP.stdWhiskDurationPerc,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel('Whisking percentage (%)')
title('Percentage of time spent whisking')
legend([s1,s2,s3],'Naive','SSP-SAP','Blank-SAP')
set(gca,'box','off')
axis square
xlim([0,4])
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Whisking Behavior - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'WhiskingBehavior_Bilateral_IOS']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'WhiskingBehavior_Bilateral_IOS'])
    %% statistical diary
    diaryFile = [dirpath 'WhiskingBehavior_Statistics_Bilateral_IOS.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    disp('======================================================================================================================')
    disp('GLME statistics for whisking duration (sec)')
    disp('======================================================================================================================')
    disp(Stats.CompA.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for whisking duration (percent)')
    disp('======================================================================================================================')
    disp(Stats.CompB.Stats)
end

end
