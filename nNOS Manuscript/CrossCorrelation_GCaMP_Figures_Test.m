function [] = CrossCorrelation_GCaMP_Figures_Test(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_CrossCorr_EGFP';
load(resultsStruct);
cd(rootFolder)
groups = {'EGFP'};
hemispheres = {'LH','RH'};
behaviors = {'All'};
samplingRate = 10;
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_CrossCorr_EGFP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for ee = 1:length(behaviors)
                behavior = behaviors{1,ee};
                data.(group).(hemisphere).dummyCheck = 1;
                if isfield(data.(group).(hemisphere),(behavior)) == false
                    data.(group).(hemisphere).(behavior).lags = [];
                    data.(group).(hemisphere).(behavior).xcVals = [];
                end
                % concatenate data across animals
                data.(group).(hemisphere).(behavior).lags = cat(1,data.(group).(hemisphere).(behavior).lags,Results_CrossCorr_EGFP.(group).(animalID).(hemisphere).(behavior).lags/samplingRate);
                data.(group).(hemisphere).(behavior).xcVals = cat(1,data.(group).(hemisphere).(behavior).xcVals,Results_CrossCorr_EGFP.(group).(animalID).(hemisphere).(behavior).xcVals);
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for dd = 1:length(behaviors)
            behavior = behaviors{1,dd};
            data.(group).(hemisphere).(behavior).mean_xcVals = mean(data.(group).(hemisphere).(behavior).xcVals,1);
            data.(group).(hemisphere).(behavior).stdErr_xcVals = std(data.(group).(hemisphere).(behavior).xcVals,0,1)./sqrt(size(data.(group).(hemisphere).(behavior).xcVals,1));
            data.(group).(hemisphere).(behavior).mean_lags = mean(data.(group).(hemisphere).(behavior).lags,1);
        end
    end
end
figure; plot(data.EGFP.RH.All.mean_lags,data.EGFP.RH.All.mean_xcVals)
hold on
plot(data.EGFP.RH.All.mean_lags,data.EGFP.RH.All.mean_xcVals + data.EGFP.RH.All.stdErr_xcVals)
plot(data.EGFP.RH.All.mean_lags,data.EGFP.RH.All.mean_xcVals - data.EGFP.RH.All.stdErr_xcVals)