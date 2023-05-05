function [] = NeuralHemoCoherence_GCaMP_SI_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_NeuralHemoCoher_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
dataTypes = {'HbT','HbO','HbR'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'C','f'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_NeuralHemoCoher_GCaMP.(group));
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
                        data.(group).(hemisphere).(dataType).(behavior).C = [];
                        data.(group).(hemisphere).(dataType).(behavior).f = [];
                        data.(group).(hemisphere).(dataType).(behavior).group = {};
                        data.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    if isempty(Results_NeuralHemoCoher_GCaMP.(group).(animalID).(dataType).(hemisphere).(behavior).C) == false
                        data.(group).(hemisphere).(dataType).(behavior).C = cat(1,data.(group).(hemisphere).(dataType).(behavior).C,Results_NeuralHemoCoher_GCaMP.(group).(animalID).(dataType).(hemisphere).(behavior).C.^2');
                        data.(group).(hemisphere).(dataType).(behavior).f = cat(1,data.(group).(hemisphere).(dataType).(behavior).f,Results_NeuralHemoCoher_GCaMP.(group).(animalID).(dataType).(hemisphere).(behavior).f);
                        data.(group).(hemisphere).(dataType).(behavior).group = cat(1,data.(group).(hemisphere).(dataType).(behavior).group,group);
                        data.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,data.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    data.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(data.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    data.(group).(hemisphere).(dataType).(behavior).(['stdErr_' variable]) = std(data.(group).(hemisphere).(dataType).(behavior).(variable),0,1)./sqrt(size(data.(group).(hemisphere).(dataType).(behavior).(variable),1));
                end
            end
        end
    end
end
% coherence figures
xlimits = ({[1/10,0.5],[1/30,0.5],[1/60,0.5],[0.003,0.5],[0.003,0.5],[0.003,0.5]});
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        summaryFigure = figure;
        sgtitle([hemisphere ' ' dataType ' neural-hemo coherence [GCaMP SI]'])
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            subplot(2,3,cc);
            p1 = semilogx(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_f,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_C,'color',colors('north texas green'),'LineWidth',2);
            hold on;
            semilogx(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_f,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_C + data.Blank_SAP.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
            semilogx(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_f,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_C - data.Blank_SAP.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
            p2 = semilogx(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_f,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_C,'color',colors('electric purple'),'LineWidth',2);
            semilogx(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_f,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_C + data.SSP_SAP.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
            semilogx(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_f,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_C - data.SSP_SAP.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
            ylabel('Coherence^2')
            xlabel('Freq (Hz)')
            title(behavior)
            xlim(xlimits{1,cc})
            ylim([0,1])
            if cc == 1
                legend([p1,p2],'Blank-SAP','SSP-SAP')
            end
            set(gca,'box','off')
            axis square
        end
        % save figure(s)
        if saveFigs == true
            dirpath = [rootFolder delim 'Summary Figures' delim 'Neural-Hemo Coherence' delim];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(summaryFigure,[dirpath 'NeuralHemoCoherence_GCaMP_SI_' hemisphere '_' dataType]);
        end
    end
end