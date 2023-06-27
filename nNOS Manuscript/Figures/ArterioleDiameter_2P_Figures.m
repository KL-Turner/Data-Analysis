function [] = ArterioleDiameter_2P_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Diameter_2P';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
behaviors = {'Rest','Whisk'};
variables = {'data'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Diameter_2P.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        vIDs = fieldnames(Results_Diameter_2P.(group).(animalID));
        for cc = 1:length(vIDs)
            vID = vIDs{cc,1};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                data.(group).(behavior).dummCheck = 1;
                for ee = 1:length(variables)
                    if isfield(data.(group).(behavior),(variables{1,ee})) == false
                        data.(group).(behavior).(variables{1,ee}) = [];
                    end
                end
                if isfield(Results_Diameter_2P.(group).(animalID).(vID),behavior) == true
                    data.(group).(behavior).data = cat(1,data.(group).(behavior).data,mean(Results_Diameter_2P.(group).(animalID).(vID).(behavior).mean));
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        data.(group).(behavior).meanData = mean(data.(group).(behavior).data,1);
        data.(group).(behavior).stdData = std(data.(group).(behavior).data,0,1);
    end
end
% figure
summaryFigure = figure;
sgtitle('Arteriole diameter [2P]')
for aa = 1:length(behaviors)
    behavior = behaviors{1,aa};
    subplot(1,2,aa);
    xInds = ones(1,length(data.Blank_SAP.(behavior).data));
    s1 = scatter(xInds*1,data.Blank_SAP.(behavior).data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
    hold on
    e1 = errorbar(1,data.Blank_SAP.(behavior).meanData,data.Blank_SAP.(behavior).stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
    e1.Color = 'black';
    e1.MarkerSize = 10;
    e1.CapSize = 10;
    xInds = ones(1,length(data.SSP_SAP.(behavior).data));
    s2 = scatter(xInds*2,data.SSP_SAP.(behavior).data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
    e2 = errorbar(2,data.SSP_SAP.(behavior).meanData,data.SSP_SAP.(behavior).stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
    e2.Color = 'black';
    e2.MarkerSize = 10;
    e2.CapSize = 10;
    title(behavior)
    ylabel('\DeltaD/D (%)')
    xlim([0,3])
    if aa == 1
        legend([s1,s2],'Blank-SAP','SSP-SAP')
    end
    set(gca,'box','off')
    set(gca,'xtick',[])
    axis square
    % save figure(s)
    if saveFigs == true
        dirpath = [rootFolder delim 'Summary Figures' delim 'Arteriole Diameter' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigure,[dirpath 'ArterioleDiameter_2P']);
    end
end