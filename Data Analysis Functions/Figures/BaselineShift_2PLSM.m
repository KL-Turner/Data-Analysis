function [] = BaselineShift_2PLSM(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% set-up
resultsStruct = 'Results_VesselBaselineShift';
load(resultsStruct);
expGroups = {'SSP-SAP','Blank-SAP'};
setName = '2PLSM Set B';
animalIDs.all = {};
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs.all = horzcat(animalIDs.all,{folderList.name});
    animalIDs.(strrep(expGroups{1,aa},'-','_')) = {folderList.name};
end
treatments = {'SSP_SAP','Blank_SAP'};
data = [];
variables = {'diameter','baseline'};
%% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs.all)
    % recognize treatment based on animal group
    if ismember(animalIDs.all{1,aa},animalIDs.SSP_SAP) == true
        treatment = 'SSP_SAP';
    elseif ismember(animalIDs.all{1,aa},animalIDs.Blank_SAP) == true
        treatment = 'Blank_SAP';
    end
    vIDs = fieldnames(Results_VesselBaselineShift.(animalIDs.all{1,aa}));
    for bb = 1:size(vIDs)
        data.(treatment).dummCheck = 1;
        for dd = 1:length(variables)
            if isfield(data.(treatment),(variables{1,dd})) == false
                data.(treatment).(variables{1,dd}) = [];
            end
        end
        data.(treatment).diameter = cat(1,data.(treatment).diameter,((Results_VesselBaselineShift.(animalIDs.all{1,aa}).(vIDs{bb,1}).diameter - Results_VesselBaselineShift.(animalIDs.all{1,aa}).(vIDs{bb,1}).baseline)/Results_VesselBaselineShift.(animalIDs.all{1,aa}).(vIDs{bb,1}).baseline)*100);
        data.(treatment).baseline = cat(1,data.(treatment).baseline,Results_VesselBaselineShift.(animalIDs.all{1,aa}).(vIDs{bb,1}).baseline);
    end
end
%% take the averages of each field through the proper dimension
for ee = 1:length(treatments)
    treatment = treatments{1,ee};
    data.(treatment).meanDiameter = mean(data.(treatment).diameter,1);
    data.(treatment).stdErrDiameter = std(data.(treatment).diameter,0,1)./sqrt(size(data.(treatment).diameter,1));
    data.(treatment).meanBaseline = mean(data.(treatment).baseline,1);
end
%% average stim-evoked figures
summaryFigure = figure;
sgtitle('Isoflurane % Increase')
s1 = scatter(1,data.Blank_SAP.meanDiameter,'d','MarkerFaceColor','k');
hold on
scatter(ones(1,length(data.Blank_SAP.diameter))*1,data.Blank_SAP.diameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on', 'jitterAmount',0.25);
s2 = scatter(2,data.SSP_SAP.meanDiameter,'d','MarkerFaceColor','k');
scatter(ones(1,length(data.SSP_SAP.diameter))*2,data.SSP_SAP.diameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','on', 'jitterAmount',0.25);
ylabel('\DeltaD/D (%)')
set(gca,'xtick',[1,1.2,2,2.2])
set(gca,'xticklabel',{['Avg Baseline: ' num2str(round(data.Blank_SAP.meanBaseline,1)) ' \muM'],['n = ' num2str(length(data.Blank_SAP.diameter)) ' arterioles'],['Avg Baseline: ' num2str(round(data.SSP_SAP.meanBaseline,1)) ' \muM'],['n = ' num2str(length(data.SSP_SAP.diameter)) ' arterioles']})
xtickangle(45)
axis square
xlim([0.5,2.5])
set(gca,'box','off')
legend([s1,s2],'Blank-SAP','SSP-SAP','Location','NorthWest')
axis square
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Baseline Shift - 2PLSM' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'BaselineShift_2PLSM']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'BaselineShift_2PLSM'])
end

end
