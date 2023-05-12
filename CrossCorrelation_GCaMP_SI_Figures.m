function [] = CrossCorrelation_GCaMP_SI_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_CrossCorr_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
dataTypes = {'HbT','HbO','HbR'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
samplingRate = 10;
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_CrossCorr_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    data.(group).(hemisphere).(dataType).dummyCheck = 1;
                    if isfield(data.(group).(hemisphere).(dataType),(behavior)) == false
                        data.(group).(hemisphere).(dataType).(behavior).lags = [];
                        data.(group).(hemisphere).(dataType).(behavior).xcVals = [];
                    end
                    if isempty(Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals) == false
                        data.(group).(hemisphere).(dataType).(behavior).lags = cat(1,data.(group).(hemisphere).(dataType).(behavior).lags,Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).lags/samplingRate);
                        data.(group).(hemisphere).(dataType).(behavior).xcVals = cat(1,data.(group).(hemisphere).(dataType).(behavior).xcVals,Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals);
                    end
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                data.(group).(hemisphere).(dataType).(behavior).mean_xcVals = mean(data.(group).(hemisphere).(dataType).(behavior).xcVals,1);
                data.(group).(hemisphere).(dataType).(behavior).stdErr_xcVals = std(data.(group).(hemisphere).(dataType).(behavior).xcVals,0,1)./sqrt(size(data.(group).(hemisphere).(dataType).(behavior).xcVals,1));
                data.(group).(hemisphere).(dataType).(behavior).mean_lags = mean(data.(group).(hemisphere).(dataType).(behavior).lags,1);
            end
        end
    end
end
% figure
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        summaryFigure = figure;
        sgtitle([hemisphere ' ' dataType ' cross correlation with GCaMP [GCaMP SI]'])
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            subplot(2,3,cc);
            p1 = plot(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_lags,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_xcVals,'color',colors('north texas green'),'LineWidth',2);
            hold on;
            plot(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_lags,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_xcVals + data.Blank_SAP.(hemisphere).(dataType).(behavior).stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
            plot(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_lags,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_xcVals - data.Blank_SAP.(hemisphere).(dataType).(behavior).stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
            p2 = plot(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_lags,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_xcVals,'color',colors('electric purple'),'LineWidth',2);
            plot(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_lags,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_xcVals + data.SSP_SAP.(hemisphere).(dataType).(behavior).stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
            plot(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_lags,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_xcVals - data.SSP_SAP.(hemisphere).(dataType).(behavior).stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
            ylabel('Corr Coef')
            xlabel('Lags (s)')
            xlim([-5,5])
            title(behavior)
            if cc == 1
                legend([p1,p2],'Blank-SAP','SSP-SAP')
            end
            set(gca,'box','off')
            axis square
        end
        % save figure(s)
        if saveFigs == true
            dirpath = [rootFolder delim 'Summary Figures' delim 'Cross Correlation' delim];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(summaryFigure,[dirpath 'CrossCorr_GCaMP_SI_' hemisphere '_' dataType]);
        end
    end
end