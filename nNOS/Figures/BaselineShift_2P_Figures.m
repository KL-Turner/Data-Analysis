function [] = BaselineShift_2P_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Baseline_2P';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
variables = {'diameter','baseline'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Baseline_2P.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        vIDs = fieldnames(Results_Baseline_2P.(group).(animalID));
        for cc = 1:length(vIDs)
            vID = vIDs{cc,1};
            data.(group).dummCheck = 1;
            for dd = 1:length(variables)
                if isfield(data.(group),(variables{1,dd})) == false
                    data.(group).(variables{1,dd}) = [];
                end
            end
            data.(group).diameter = cat(1,data.(group).diameter,((Results_Baseline_2P.(group).(animalID).(vID).diameter - Results_Baseline_2P.(group).(animalID).(vID).baseline)/Results_Baseline_2P.(group).(animalID).(vID).baseline)*100);
            data.(group).baseline = cat(1,data.(group).baseline,Results_Baseline_2P.(group).(animalID).(vID).baseline);
        end
    end
end
% mean/std
for ee = 1:length(groups)
    group = groups{1,ee};
    data.(group).meanDiameter = mean(data.(group).diameter,1);
    data.(group).stdDiameter = std(data.(group).diameter,0,1);
    data.(group).meanBaseline = mean(data.(group).baseline,1);
end
% figure
summaryFigure = figure;
sgtitle('Isoflurane % Increase')
s1 = scatter(1,data.Blank_SAP.meanDiameter,'d','MarkerFaceColor','k');
hold on
scatter(ones(1,length(data.Blank_SAP.diameter))*1,data.Blank_SAP.diameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on', 'jitterAmount',0.25);
s2 = scatter(2,data.SSP_SAP.meanDiameter,'d','MarkerFaceColor','k');
scatter(ones(1,length(data.SSP_SAP.diameter))*2,data.SSP_SAP.diameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on', 'jitterAmount',0.25);
ylabel('\DeltaD/D (%)')
set(gca,'xtick',[1,1.2,2,2.2])
set(gca,'xticklabel',{['Avg Baseline: ' num2str(round(data.Blank_SAP.meanBaseline,1)) ' \muM'],['n = ' num2str(length(data.Blank_SAP.diameter)) ' arterioles'],['Avg Baseline: ' num2str(round(data.SSP_SAP.meanBaseline,1)) ' \muM'],['n = ' num2str(length(data.SSP_SAP.diameter)) ' arterioles']})
xtickangle(45)
axis square
xlim([0.5,2.5])
set(gca,'box','off')
legend([s1,s2],'Blank-SAP','SSP-SAP')
axis square
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures' delim 'Arteriole Diameter' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'BaselineShift_2P']);
end