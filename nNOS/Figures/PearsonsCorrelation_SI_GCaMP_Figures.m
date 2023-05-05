function [] = PearsonsCorrelation_SI_GCaMP_Figures(rootFolder,saveFigs,delim)
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
        for dd = 1:length(dataTypes)
            dataType = dataTypes{1,dd};
            data.(group).(dataType).dummCheck = 1;
            for ee = 1:length(behaviors)
                behavior = behaviors{1,ee};
                if isfield(data.(group).(dataType),behavior) == false
                    data.(group).(dataType).(behavior).data = [];
                    data.(group).(dataType).(behavior).group = {};
                    data.(group).(dataType).(behavior).animalID = {};
                end
                data.(group).(dataType).(behavior).data = cat(1,data.(group).(dataType).(behavior).data,mean(Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).(behavior).R));
                data.(group).(dataType).(behavior).group = cat(1,data.(group).(dataType).(behavior).group,group);
                data.(group).(dataType).(behavior).animalID = cat(1,data.(group).(dataType).(behavior).animalID,animalID);
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for dd = 1:length(dataTypes)
        dataType = dataTypes{1,dd};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            data.(group).(dataType).(behavior).meanData = nanmean(data.(group).(dataType).(behavior).data,1);
            data.(group).(dataType).(behavior).stdData = nanstd(data.(group).(dataType).(behavior).data,0,1);
        end
    end
end
% figure
for qq = 1:length(dataTypes)
    dataType = dataTypes{1,qq};
    summaryFigure = figure;
    sgtitle([dataType ' Pearsons Corr [SI GCaMP]'])
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        subplot(2,3,bb);
        xInds = ones(1,length(data.Blank_SAP.(dataType).(behavior).data));
        s1 = scatter(xInds*2,data.Blank_SAP.(dataType).(behavior).data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
        hold on
        e1 = errorbar(2,data.Blank_SAP.(dataType).(behavior).meanData,data.Blank_SAP.(dataType).(behavior).stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e1.Color = 'black';
        e1.MarkerSize = 10;
        e1.CapSize = 10;
        xInds = ones(1,length(data.SSP_SAP.(dataType).(behavior).data));
        s2 = scatter(xInds*3,data.SSP_SAP.(dataType).(behavior).data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
        hold on
        e2 = errorbar(3,data.SSP_SAP.(dataType).(behavior).meanData,data.SSP_SAP.(dataType).(behavior).stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e2.Color = 'black';
        e2.MarkerSize = 10;
        e2.CapSize = 10;
        title(behavior)
        ylabel('Corr Coef')
        xlim([0,4])
        if bb == 1
            legend([s1,s2],'Blank-SAP','SSP-SAP')
        end
        set(gca,'box','off')
        set(gca,'xtick',[])
        axis square
    end
    linkaxes
    % save figure(s)
    if saveFigs == true
        dirpath = [rootFolder delim 'Summary Figures' delim 'Pearsons Corr SI GCaMP' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigure,[dirpath 'PearsonsCorr_GCaMP_SI_ ' dataType]);
    end
end
% % statistics - generalized linear mixed effects model
% for bb = 1:length(hemispheres)
%     hemisphere = hemispheres{1,bb};
%     for aa = 1:length(dataTypes)
%         dataType = dataTypes{1,aa};
%         for dd = 1:length(behaviors)
%             behavior = behaviors{1,dd};
%             for cc = 1:length(variables)
%                 variable = variables{1,cc};
%                 % statistics - unpaired ttest
%                 [stats.(dataType).(behavior).data.h,...
%                     stats.(dataType).(behavior).data.p,...
%                     stats.(dataType).(behavior).data.ci,...
%                     stats.(dataType).(behavior).data.stats] ...
%                     = ttest2(data.Blank_SAP.(dataType).(behavior).data,data.SSP_SAP.(dataType).(behavior).data);
%             end
%         end
%     end
% end
% % statistical diary
% if saveFigs == true
%     % statistical diary
%     diaryFile = [dirpath 'WhiskingBehavior_GCaMP_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     disp('======================================================================================================================')
%     disp('ttest2 statistics: Blank vs. SAP')
%     disp('======================================================================================================================')
%     % statistics - generalized linear mixed effects model
%     for bb = 1:length(hemispheres)
%         hemisphere = hemispheres{1,bb};
%         for aa = 1:length(dataTypes)
%             dataType = dataTypes{1,aa};
%             for dd = 1:length(behaviors)
%                 behavior = behaviors{1,dd};
%                 for cc = 1:length(variables)
%                     variable = variables{1,cc};
%                     disp([hemisphere ' ' dataType ' ' behavior ' ' variable ' p < ' num2str(stats.(dataType).(behavior).data.p)]); disp(' ')
%                     disp('----------------------------------------------------------------------------------------------------------------------')
%                 end
%             end
%         end
%     end
%     diary off
% end