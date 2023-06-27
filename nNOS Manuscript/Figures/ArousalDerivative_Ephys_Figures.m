function [] = ArousalDerivative_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Derivative_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
dataTypes = {'HbT'};
behaviors = {'Alert','Asleep','All'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Derivative_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                data.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(data.(group).(hemisphere).(dataType),behavior) == false
                        data.(group).(hemisphere).(dataType).(behavior).data = [];
                    end
                    for ff = 1:length(Results_Derivative_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).data)
                        data.(group).(hemisphere).(dataType).(behavior).data = cat(2,data.(group).(hemisphere).(dataType).(behavior).data,Results_Derivative_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).data{ff,1});
                    end
                end
            end
        end
    end
end
% figure
edges = -2:0.05:2;
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    summaryFigure = figure;
    sgtitle(['IOS ' hemisphere ' ' dataType ' [Ephys]'])
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        subplot(1,3,bb);
        %         h1 = histogram(data.Naive.(hemisphere).HbT.(behavior).data,edges,'Normalization','probability','FaceColor',colors('sapphire'));
        h3 = histogram(data.SSP_SAP.(hemisphere).HbT.(behavior).data,edges,'Normalization','probability','FaceColor',colors('electric purple'));
        hold on
        h2 = histogram(data.Blank_SAP.(hemisphere).HbT.(behavior).data,edges,'Normalization','probability','FaceColor',colors('north texas green'));
        title(behavior)
        xlabel('diff(\Delta[HbT])')
        ylabel('Probability')
        if bb == 1
            legend([h2,h3],'Blank-SAP','SSP-SAP')
        end
        set(gca,'box','off')
        axis square
    end
    % save figure(s)
    if saveFigs == true
        dirpath = [rootFolder delim 'Summary Figures' delim 'Intrinsic Signals' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigure,[dirpath 'ArousalDerivative_Ephys_ ' hemisphere]);
    end
end