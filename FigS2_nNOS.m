function [] = FigS2_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PearsonCorr_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'SSP_SAP','Blank_SAP'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PearsonCorr_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            gcampCorrCoefData.(group).(dataType).dummCheck = 1;
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                if isfield(gcampCorrCoefData.(group).(dataType),behavior) == false
                    gcampCorrCoefData.(group).(dataType).(behavior).data = [];
                    gcampCorrCoefData.(group).(dataType).(behavior).group = {};
                    gcampCorrCoefData.(group).(dataType).(behavior).animalID = {};
                end
                if isempty(Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).(behavior).R) == false
                    gcampCorrCoefData.(group).(dataType).(behavior).data = cat(1,gcampCorrCoefData.(group).(dataType).(behavior).data,mean(Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).(behavior).R));
                    gcampCorrCoefData.(group).(dataType).(behavior).group = cat(1,gcampCorrCoefData.(group).(dataType).(behavior).group,group);
                    gcampCorrCoefData.(group).(dataType).(behavior).animalID = cat(1,gcampCorrCoefData.(group).(dataType).(behavior).animalID,animalID);
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
            gcampCorrCoefData.(group).(dataType).(behavior).meanData = mean(gcampCorrCoefData.(group).(dataType).(behavior).data,1);
            gcampCorrCoefData.(group).(dataType).(behavior).stdData = std(gcampCorrCoefData.(group).(dataType).(behavior).data,0,1);
        end
    end
end
% figure
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    summaryFigure = figure;
    sgtitle([dataType ' Pearsons Corr [GCaMP SI]'])
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        subplot(2,3,bb);
        xInds = ones(1,length(gcampCorrCoefData.Blank_SAP.(dataType).(behavior).data));
        s1 = scatter(xInds*1,gcampCorrCoefData.Blank_SAP.(dataType).(behavior).data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
        hold on;
        e1 = errorbar(1,gcampCorrCoefData.Blank_SAP.(dataType).(behavior).meanData,gcampCorrCoefData.Blank_SAP.(dataType).(behavior).stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e1.Color = 'black';
        e1.MarkerSize = 10;
        e1.CapSize = 10;
        xInds = ones(1,length(gcampCorrCoefData.SSP_SAP.(dataType).(behavior).data));
        s2 = scatter(xInds*2,gcampCorrCoefData.SSP_SAP.(dataType).(behavior).data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
        e2 = errorbar(2,gcampCorrCoefData.SSP_SAP.(dataType).(behavior).meanData,gcampCorrCoefData.SSP_SAP.(dataType).(behavior).stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e2.Color = 'black';
        e2.MarkerSize = 10;
        e2.CapSize = 10;
        title(behavior)
        ylabel('Corr Coef')
        xlim([0,3])
        if bb == 1
            legend([s1,s2],'Blank-SAP','SSP-SAP')
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
        savefig(summaryFigure,[dirpath 'PearsonsCorr_GCaMP_SI_ ' dataType]);
        set(summaryFigure,'PaperPositionMode','auto');
        print('-vector','-dpdf','-fillpage',[dirpath 'PearsonsCorr_GCaMP_SI_ ' dataType])
    end
end