function [] = WhiskingBehavior_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_WhiskBehav_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Naive','Blank_SAP','SSP_SAP'};
dataTypes = {'whiskDurationSec','whiskDurationPerc'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_WhiskBehav_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        data.(group).dummCheck = 1;
        data.(group).group = {};
        data.(group).animalID = {};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            if isfield(data.(group),dataType) == false
                data.(group).(dataType) = [];
            end
        end
        % concatenate data across animals
        data.(group).whiskDurationSec = cat(1,data.(group).whiskDurationSec,Results_WhiskBehav_Ephys.(group).(animalID).whiskDurationSec/60);
        data.(group).whiskDurationPerc = cat(1,data.(group).whiskDurationPerc,Results_WhiskBehav_Ephys.(group).(animalID).whiskDurationPerc);
        data.(group).group = cat(1,data.(group).group,group);
        data.(group).animalID = cat(1,data.(group).animalID,animalID);
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.(group).(['mean_' dataType]) = mean(data.(group).(dataType),1);
        data.(group).(['std_' dataType]) = std(data.(group).(dataType),0,1);
    end
end
% figure
summaryFigure = figure;
sgtitle('Time spent whisking [Ephys]')
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    subplot(1,2,aa);
    s1 = scatter(ones(1,length(data.Naive.(dataType)))*1,data.Naive.(dataType),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off','jitterAmount',0.25);
    hold on;
    e1 = errorbar(1,data.Naive.(['mean_' dataType]),data.Naive.(['std_' dataType]),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
    e1.Color = 'black';
    e1.MarkerSize = 10;
    e1.CapSize = 10;
    s2 = scatter(ones(1,length(data.Blank_SAP.(dataType)))*2,data.Blank_SAP.(dataType),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
    e2 = errorbar(2,data.Blank_SAP.(['mean_' dataType]),data.Blank_SAP.(['std_' dataType]),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
    e2.Color = 'black';
    e2.MarkerSize = 10;
    e2.CapSize = 10;
    s3 = scatter(ones(1,length(data.SSP_SAP.(dataType)))*3,data.SSP_SAP.(dataType),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
    e3 = errorbar(3,data.SSP_SAP.(['mean_' dataType]),data.SSP_SAP.(['std_' dataType]),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
    e3.Color = 'black';
    e3.MarkerSize = 10;
    e3.CapSize = 10;
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    ylabel(dataType)
    title(dataType)
    legend([s1,s2,s3],'Naive','Blank-SAP','SSP-SAP')
    set(gca,'box','off')
    axis square
    xlim([0,4])
end
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures' delim 'Behavior' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'WhiskingBehavior_Ephys']);
end
% statistics - unpaired ttest
[whiskPercStats1.h,whiskPercStats1.p,whiskPercStats1.ci,whiskPercStats1.stats] = ttest2(data.Naive.whiskDurationPerc,data.Blank_SAP.whiskDurationPerc);
[whiskPercStats2.h,whiskPercStats2.p,whiskPercStats2.ci,whiskPercStats2.stats] = ttest2(data.Blank_SAP.whiskDurationPerc,data.SSP_SAP.whiskDurationPerc);
% statistical diary
if saveFigs == true
    % statistical diary
    diaryFile = [dirpath 'WhiskingBehavior_Ephys_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    disp('======================================================================================================================')
    disp('ttest2 statistics:')
    disp('======================================================================================================================')
    disp(['Naive vs. Blank (whiskPerc) p < ' num2str(whiskPercStats1.p)]); disp(' ')
    disp(['Blank vs. SAP (whiskPerc) p < ' num2str(whiskPercStats2.p)]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
end