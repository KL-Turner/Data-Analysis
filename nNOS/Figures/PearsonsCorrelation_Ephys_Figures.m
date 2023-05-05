function [] = PearsonsCorrelation_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PearsonCorr_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','SSP_SAP','Blank_SAP'};
dataTypes = {'HbT','gammaBandPower','deltaBandPower'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PearsonCorr_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            data.(group).(dataType).dummCheck = 1;
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                if isfield(data.(group).(dataType),behavior) == false
                    data.(group).(dataType).(behavior).data = [];
                    data.(group).(dataType).(behavior).group = {};
                    data.(group).(dataType).(behavior).animalID = {};
                end
                if isempty(Results_PearsonCorr_Ephys.(group).(animalID).(dataType).(behavior).R) == false
                    data.(group).(dataType).(behavior).data = cat(1,data.(group).(dataType).(behavior).data,mean(Results_PearsonCorr_Ephys.(group).(animalID).(dataType).(behavior).R));
                    data.(group).(dataType).(behavior).group = cat(1,data.(group).(dataType).(behavior).group,group);
                    data.(group).(dataType).(behavior).animalID = cat(1,data.(group).(dataType).(behavior).animalID,animalID);
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            data.(group).(dataType).(behavior).meanData = mean(data.(group).(dataType).(behavior).data,1);
            data.(group).(dataType).(behavior).stdData = std(data.(group).(dataType).(behavior).data,0,1);
        end
    end
end
% figure
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    summaryFigure = figure;
    sgtitle([dataType ' Pearsons Corr [Ephys]'])
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        subplot(2,3,bb);
        xInds = ones(1,length(data.Naive.(dataType).(behavior).data));
        s1 = scatter(xInds*1,data.Naive.(dataType).(behavior).data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off','jitterAmount',0.25);
        hold on;
        e1 = errorbar(1,data.Naive.(dataType).(behavior).meanData,data.Naive.(dataType).(behavior).stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e1.Color = 'black';
        e1.MarkerSize = 10;
        e1.CapSize = 10;
        xInds = ones(1,length(data.Blank_SAP.(dataType).(behavior).data));
        s2 = scatter(xInds*2,data.Blank_SAP.(dataType).(behavior).data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
        e2 = errorbar(2,data.Blank_SAP.(dataType).(behavior).meanData,data.Blank_SAP.(dataType).(behavior).stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e2.Color = 'black';
        e2.MarkerSize = 10;
        e2.CapSize = 10;
        xInds = ones(1,length(data.SSP_SAP.(dataType).(behavior).data));
        s3 = scatter(xInds*3,data.SSP_SAP.(dataType).(behavior).data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
        e3 = errorbar(3,data.SSP_SAP.(dataType).(behavior).meanData,data.SSP_SAP.(dataType).(behavior).stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e3.Color = 'black';
        e3.MarkerSize = 10;
        e3.CapSize = 10;
        title(behavior)
        ylabel('Corr Coef')
        xlim([0,4])
        if bb == 1
            legend([s1,s2,s3],'Naive','Blank-SAP','SSP-SAP')
        end
        set(gca,'box','off')
        set(gca,'xtick',[])
        axis square
    end
    % save figure(s)
    if saveFigs == true
        dirpath = [rootFolder delim 'Summary Figures' delim 'Pearsons Correlations' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigure,[dirpath 'PearsonsCorr_Ephys_ ' dataType]);
    end
end