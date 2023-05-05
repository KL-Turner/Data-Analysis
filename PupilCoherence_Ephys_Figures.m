function [] = PupilCoherence_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PupilCoher_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Naive','Blank_SAP','SSP_SAP'};
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
compDataTypes = {'LH_HbT','RH_HbT','LH_gammaBandPower','RH_gammaBandPower'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PupilCoher_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(compDataTypes)
                compDataType = compDataTypes{1,dd};
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    data.(group).(dataType).(compDataType).dummyCheck = 1;
                    if isfield(data.(group).(dataType).(compDataType),(behavior)) == false
                        data.(group).(dataType).(compDataType).(behavior).C = [];
                        data.(group).(dataType).(compDataType).(behavior).f = [];
                    end
                    if isempty(Results_PupilCoher_Ephys.(group).(animalID).(dataType).(compDataType).(behavior).C) == false
                        data.(group).(dataType).(compDataType).(behavior).C = cat(1,data.(group).(dataType).(compDataType).(behavior).C,Results_PupilCoher_Ephys.(group).(animalID).(dataType).(compDataType).(behavior).C.^2');
                        data.(group).(dataType).(compDataType).(behavior).f = cat(1,data.(group).(dataType).(compDataType).(behavior).f,Results_PupilCoher_Ephys.(group).(animalID).(dataType).(compDataType).(behavior).f);
                    end
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
        for cc = 1:length(compDataTypes)
            compDataType = compDataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                data.(group).(dataType).(compDataType).(behavior).mean_C = mean(data.(group).(dataType).(compDataType).(behavior).C,1);
                data.(group).(dataType).(compDataType).(behavior).stdErr_C = std(data.(group).(dataType).(compDataType).(behavior).C,0,1)./sqrt(size(data.(group).(dataType).(compDataType).(behavior).C,1));
                data.(group).(dataType).(compDataType).(behavior).mean_f = mean(data.(group).(dataType).(compDataType).(behavior).f,1);
            end
        end
    end
end
% figure
xlimits = ({[1/10,0.5],[1/30,0.5],[1/60,0.5],[0.003,0.5],[0.003,0.5],[0.003,0.5]});
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    for bb = 1:length(compDataTypes)
        compDataType = compDataTypes{1,bb};
        summaryFigure = figure;
        sgtitle([dataType ' coherence with ' strrep(compDataType,'_',' ') ' [Ephys]'])
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            subplot(2,3,cc);
            p1 = semilogx(data.Naive.(dataType).(compDataType).(behavior).mean_f,data.Blank_SAP.(dataType).(compDataType).(behavior).mean_C,'color',colors('sapphire'),'LineWidth',2);
            hold on;
            semilogx(data.Naive.(dataType).(compDataType).(behavior).mean_f,data.Blank_SAP.(dataType).(compDataType).(behavior).mean_C + data.Naive.(dataType).(compDataType).(behavior).stdErr_C,'color',colors('sapphire'),'LineWidth',0.25);
            semilogx(data.Naive.(dataType).(compDataType).(behavior).mean_f,data.Blank_SAP.(dataType).(compDataType).(behavior).mean_C - data.Naive.(dataType).(compDataType).(behavior).stdErr_C,'color',colors('sapphire'),'LineWidth',0.25);
            p2 = semilogx(data.Blank_SAP.(dataType).(compDataType).(behavior).mean_f,data.Blank_SAP.(dataType).(compDataType).(behavior).mean_C,'color',colors('north texas green'),'LineWidth',2);
            semilogx(data.Blank_SAP.(dataType).(compDataType).(behavior).mean_f,data.Blank_SAP.(dataType).(compDataType).(behavior).mean_C + data.Blank_SAP.(dataType).(compDataType).(behavior).stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
            semilogx(data.Blank_SAP.(dataType).(compDataType).(behavior).mean_f,data.Blank_SAP.(dataType).(compDataType).(behavior).mean_C - data.Blank_SAP.(dataType).(compDataType).(behavior).stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
            p3 = semilogx(data.SSP_SAP.(dataType).(compDataType).(behavior).mean_f,data.SSP_SAP.(dataType).(compDataType).(behavior).mean_C,'color',colors('electric purple'),'LineWidth',2);
            semilogx(data.SSP_SAP.(dataType).(compDataType).(behavior).mean_f,data.SSP_SAP.(dataType).(compDataType).(behavior).mean_C + data.SSP_SAP.(dataType).(compDataType).(behavior).stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
            semilogx(data.SSP_SAP.(dataType).(compDataType).(behavior).mean_f,data.SSP_SAP.(dataType).(compDataType).(behavior).mean_C - data.SSP_SAP.(dataType).(compDataType).(behavior).stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
            ylabel('Coherence^2')
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
            dirpath = [rootFolder delim 'Summary Figures' delim 'Behavior' delim];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(summaryFigure,[dirpath 'PupilCoherence_Ephys_' dataType '_' compDataType]);
        end
    end
end