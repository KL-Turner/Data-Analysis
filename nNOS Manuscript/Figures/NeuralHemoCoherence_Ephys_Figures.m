function [] = NeuralHemoCoherence_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_NeuralHemoCoher_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
dataTypes = {'gammaBandPower','deltaBandPower'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'C','f'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_NeuralHemoCoher_Ephys.(group));
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
                    if isempty(Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).C) == false
                        data.(group).(hemisphere).(dataType).(behavior).C = cat(1,data.(group).(hemisphere).(dataType).(behavior).C,Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).C');
                        data.(group).(hemisphere).(dataType).(behavior).f = cat(1,data.(group).(hemisphere).(dataType).(behavior).f,Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).f);
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
        sgtitle([hemisphere ' ' dataType ' neural-hemo coherence [Ephys]'])
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            subplot(2,3,cc);
            p1 = plot(data.Naive.(hemisphere).(dataType).(behavior).mean_f,data.Naive.(hemisphere).(dataType).(behavior).mean_C,'color',colors('sapphire'),'LineWidth',2);
            hold on;
            plot(data.Naive.(hemisphere).(dataType).(behavior).mean_f,data.Naive.(hemisphere).(dataType).(behavior).mean_C + data.Naive.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('sapphire'),'LineWidth',0.25);
            plot(data.Naive.(hemisphere).(dataType).(behavior).mean_f,data.Naive.(hemisphere).(dataType).(behavior).mean_C - data.Naive.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('sapphire'),'LineWidth',0.25);
            p2 = plot(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_f,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_C,'color',colors('north texas green'),'LineWidth',2);
            plot(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_f,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_C + data.Blank_SAP.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
            plot(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_f,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_C - data.Blank_SAP.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
            p3 = plot(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_f,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_C,'color',colors('electric purple'),'LineWidth',2);
            plot(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_f,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_C + data.SSP_SAP.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
            plot(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_f,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_C - data.SSP_SAP.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
            ylabel('Coherence')
            xlabel('Freq (Hz)')
            title(behavior)
            xlim(xlimits{1,cc})
            ylim([0,1])
            if cc == 1
                legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
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
            savefig(summaryFigure,[dirpath 'NeuralHemoCoherence_Ephys_' hemisphere '_' dataType]);
            set(summaryFigure,'PaperPositionMode','auto');
            print('-vector','-dpdf','-fillpage',[dirpath 'NeuralHemoCoherence_Ephys_' hemisphere '_' dataType])
        end
    end
end