function [] = CrossCorr_Ephys_Figures(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_CrossCorr_Ephys';
load(resultsStruct);
cd(rootFolder)
groups = {'Naive','SSP_SAP','Blank_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Rest','NREM','REM'};
dataTypes = {'HbTvLFPxcVals','LFP_lags','F','HbTvMUAxcVals','MUA_lags'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_CrossCorr_Ephys.(group));
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
                data.(group).(hemisphere).(behavior).HbTvLFPxcVals = cat(3,data.(group).(hemisphere).(behavior).HbTvLFPxcVals,Results_CrossCorr_Ephys.(group).(animalID).(behavior).(hemisphere).HbTvLFPxcVals);
                data.(group).(hemisphere).(behavior).LFP_lags = cat(3,data.(group).(hemisphere).(behavior).LFP_lags,Results_CrossCorr_Ephys.(group).(animalID).(behavior).(hemisphere).LFP_lags);
                data.(group).(hemisphere).(behavior).F = cat(3,data.(group).(hemisphere).(behavior).F,Results_CrossCorr_Ephys.(group).(animalID).(behavior).(hemisphere).F);
                data.(group).(hemisphere).(behavior).HbTvMUAxcVals = cat(1,data.(group).(hemisphere).(behavior).HbTvMUAxcVals,Results_CrossCorr_Ephys.(group).(animalID).(behavior).(hemisphere).HbTvMUAxcVals);
                data.(group).(hemisphere).(behavior).MUA_lags = cat(1,data.(group).(hemisphere).(behavior).MUA_lags,Results_CrossCorr_Ephys.(group).(animalID).(behavior).(hemisphere).LFP_lags);
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
                if any(strcmp(dataType,{'HbTvLFPxcVals','LFP_lags','F'}))
                    data.(group).(hemisphere).(behavior).(dataType) = data.(group).(hemisphere).(behavior).(dataType);
                    data.(group).(hemisphere).(behavior).(['mean_' dataType]) = mean(data.(group).(hemisphere).(behavior).(dataType),3);
                else
                    data.(group).(hemisphere).(behavior).(dataType) = data.(group).(hemisphere).(behavior).(dataType);
                    data.(group).(hemisphere).(behavior).(['mean_' dataType]) = mean(data.(group).(hemisphere).(behavior).(dataType),1);
                    data.(group).(hemisphere).(behavior).(['stdErr_' dataType]) = std(data.(group).(hemisphere).(behavior).(dataType),1)./sqrt(size(data.(group).(hemisphere).(behavior).(dataType),1));
                end
            end
        end
    end
end
% figure
freq = 30;
lag = 5;
groupNames = {'Naive','Naive','Blank_SAP','Blank_SAP','SSP_SAP','SSP_SAP'};
hemNames = {'LH','RH','LH','RH','LH','RH'};
for aa = 1:length(behaviors)
    figure;
    sgtitle('LFP-[HbT] XCorr')
    behavior = behaviors{1,aa};
    for bb = 1:6
        subplot(3,2,bb);
        groupName = groupNames{1,bb};
        hemName = hemNames{1,bb};
        imagesc(data.(groupName).(hemName).(behavior).mean_LFP_lags,data.(groupName).(hemName).(behavior).mean_F,data.(groupName).(hemName).(behavior).mean_HbTvLFPxcVals)
        title([groupName ' ' hemName ' ' behavior])
        xticks([-lag*freq,-lag*freq/2,0,lag*freq/2,lag*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lag*freq,lag*freq])
        xlabel('Lags (s)')
        ylabel('Freq (Hz)')
        ylim([1,100])
        c1 = colorbar;
        ylabel(c1,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
        axis xy
        axis square
        set(gca,'box','off')
    end
end
% figure
hemNames = {'LH','RH'};
for aa = 1:length(behaviors)
    figure;
    sgtitle('MUA-[HbT] XCorr')
    behavior = behaviors{1,aa};
    for bb = 1:2
        subplot(1,2,bb);
        hemName = hemNames{1,bb};
        p1 = plot(data.Naive.(hemName).(behavior).mean_MUA_lags,data.Naive.(hemName).(behavior).mean_HbTvMUAxcVals,'color',colors('sapphire'),'LineWidth',2);
        hold on
        p2 = plot(data.Blank_SAP.(hemName).(behavior).mean_MUA_lags,data.Blank_SAP.(hemName).(behavior).mean_HbTvMUAxcVals,'color',colors('north texas green'),'LineWidth',2);
        p3 = plot(data.SSP_SAP.(hemName).(behavior).mean_MUA_lags,data.SSP_SAP.(hemName).(behavior).mean_HbTvMUAxcVals,'color',colors('electric purple'),'LineWidth',2);
        title([hemName ' ' behavior])
        xticks([-lag*freq,-lag*freq/2,0,lag*freq/2,lag*freq])
        xticklabels({'-5','-2.5','0','2.5','5'})
        xlim([-lag*freq,lag*freq])
        xlabel('Lags (s)')
        ylabel('Corr')
        if bb == 1
            legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
        end
        axis square
        set(gca,'box','off')
    end
end
