function [] = CrossCorr_GCaMP_Figures(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_CrossCorr_GCaMP';
load(resultsStruct);
cd(rootFolder)
groups = {'SSP_SAP','Blank_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Rest','NREM','REM'};
dataTypes = {'HbTvGCaMPxcVals','GCaMP_lags'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_CrossCorr_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                data.(group).(hemisphere).(behavior).dummCheck = 1;
                for ee = 1:length(dataTypes)
                    dataType = dataTypes{1,ee};
                    if isfield(data.(group).(hemisphere).(behavior),dataType) == false
                        data.(group).(hemisphere).(behavior).(dataType) = [];
                    end
                end
                data.(group).(hemisphere).(behavior).HbTvGCaMPxcVals = cat(1,data.(group).(hemisphere).(behavior).HbTvGCaMPxcVals,Results_CrossCorr_GCaMP.(group).(animalID).(behavior).(hemisphere).HbTvGCaMPxcVals);
                data.(group).(hemisphere).(behavior).GCaMP_lags = cat(1,data.(group).(hemisphere).(behavior).GCaMP_lags,Results_CrossCorr_GCaMP.(group).(animalID).(behavior).(hemisphere).GCaMP_lags);
            end
        end
    end
end
% mean cross correlation per group for each data type
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                data.(group).(hemisphere).(behavior).(dataType) = data.(group).(hemisphere).(behavior).(dataType);
                data.(group).(hemisphere).(behavior).(['mean_' dataType]) = mean(data.(group).(hemisphere).(behavior).(dataType),1);
                data.(group).(hemisphere).(behavior).(['stdErr_' dataType]) = std(data.(group).(hemisphere).(behavior).(dataType),1)./sqrt(size(data.(group).(hemisphere).(behavior).(dataType),1));
            end
        end
    end
end
% figure
lag = 5;
freq = 10;
hemNames = {'LH','RH'};
for aa = 1:length(behaviors)
    figure;
    sgtitle('GCaMP-[HbT] XCorr')
    behavior = behaviors{1,aa};
    for bb = 1:2
        subplot(1,2,bb);
        hemName = hemNames{1,bb};
        p1 = plot(data.Blank_SAP.(hemName).(behavior).mean_GCaMP_lags,data.Blank_SAP.(hemName).(behavior).mean_HbTvGCaMPxcVals,'color',colors('north texas green'),'LineWidth',2);
        hold on
        p2 = plot(data.SSP_SAP.(hemName).(behavior).mean_GCaMP_lags,data.SSP_SAP.(hemName).(behavior).mean_HbTvGCaMPxcVals,'color',colors('electric purple'),'LineWidth',2);
        title([hemName ' ' behavior])
        xticks([-lag*freq,-lag*freq/2,0,lag*freq/2,lag*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lag*freq,lag*freq])
        xlabel('Lags (s)')
        ylabel('Corr')
        if bb == 1
            legend([p1,p2],'Blank-SAP','SSP-SAP')
        end
        axis square
        set(gca,'box','off')
    end
end
