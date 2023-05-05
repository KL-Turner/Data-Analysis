function [] = PupilEvoked_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PupilEvoked_Ephys';
load(resultsStruct);
cd(rootFolder)
groups = {'Naive','SSP_SAP','Blank_SAP'};
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
behaviors = {'control','Whisk','Stim','NREM','REM'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PupilEvoked_Ephys.(group));
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
                % concatenate data across animals
                data.(group).(dataType).(behavior).data = cat(1,data.(group).(dataType).(behavior).data,nanmean(Results_PupilArea_Ephys.(group).(animalID).(dataType).(behavior).eventMeans));
                data.(group).(dataType).(behavior).group = cat(1,data.(group).(dataType).(behavior).group,group);
                data.(group).(dataType).(behavior).animalID = cat(1,data.(group).(dataType).(behavior).animalID,animalID);
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
            data.(group).(dataType).(behavior).meanData = nanmean(data.(group).(dataType).(behavior).data,1);
            data.(group).(dataType).(behavior).stdData = nanstd(data.(group).(dataType).(behavior).data,0,1);
        end
    end
end
% figure
labels = {'mm^2','mm','z-units','z-units'};
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    summaryFigure = figure;
    sgtitle(['Pupil ' dataType])
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        subplot(2,3,bb)
        s1 = scatter(ones(1,length(data.Naive.(dataType).(behavior).data))*1,data.Naive.(dataType).(behavior).data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off','jitterAmount',0);
        hold on
        e1 = errorbar(1,data.Naive.(dataType).(behavior).meanData,data.Naive.(dataType).(behavior).stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e1.Color = 'black';
        e1.MarkerSize = 10;
        e1.CapSize = 10;
        s2 = scatter(ones(1,length(data.Blank_SAP.(dataType).(behavior).data))*2,data.Blank_SAP.(dataType).(behavior).data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0);
        hold on
        e2 = errorbar(2,data.Blank_SAP.(dataType).(behavior).meanData,data.Blank_SAP.(dataType).(behavior).stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e2.Color = 'black';
        e2.MarkerSize = 10;
        e2.CapSize = 10;
        s3 = scatter(ones(1,length(data.SSP_SAP.(dataType).(behavior).data))*3,data.SSP_SAP.(dataType).(behavior).data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0);
        hold on
        e3 = errorbar(3,data.SSP_SAP.(dataType).(behavior).meanData,data.SSP_SAP.(dataType).(behavior).stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e3.Color = 'black';
        e3.MarkerSize = 10;
        e3.CapSize = 10;
        ylabel(labels{1,aa})
        title(behavior)
        legend([s1,s2,s3],'Naive','Blank-SAP','SSP-SAP','Location','NorthEast')
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        axis square
        xlim([0,4])
        set(gca,'box','off')
    end
    % save figure(s)
    if saveFigs == true
        dirpath = [rootFolder delim 'Summary Figures' delim 'Pupil Size Ephys' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigure,[dirpath 'PupilSize_' dataType '_Ephys']);
    end
end