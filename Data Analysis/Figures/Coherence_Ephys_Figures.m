function [] = Coherence_Ephys_Figures(rootFolder,~,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Coher_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','SSP_SAP','Blank_SAP'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
dataTypes = {'CBV_HbT','gammaBandPower'};
variables = {'C','f'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Coher_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            data.(group).(behavior).dummCheck = 1;
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if isfield(data.(group).(behavior),(dataType)) == false
                    data.(group).(behavior).(dataType).C = [];
                    data.(group).(behavior).(dataType).f = [];
                end
                data.(group).(behavior).(dataType).C = cat(1,data.(group).(behavior).(dataType).C,Results_Coher_Ephys.(group).(animalID).(behavior).(dataType).C.^2');
                data.(group).(behavior).(dataType).f = cat(1,data.(group).(behavior).(dataType).f,Results_Coher_Ephys.(group).(animalID).(behavior).(dataType).f);
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                data.(group).(behavior).(dataType).(['mean_' variable]) = mean(data.(group).(behavior).(dataType).(variable),1);
                data.(group).(behavior).(dataType).(['stdErr_' variable]) = std(data.(group).(behavior).(dataType).(variable),1)./sqrt(size(data.(group).(behavior).(dataType).(variable),1));
            end
        end
    end
end
% figure
xlimits = ({[1/10,0.5],[1/30,0.5],[1/60,0.5],[0.003,0.5],[0.003,0.5],[0.003,0.5]});
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    figure;
    sgtitle([dataType ' bilateral coherence'])
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        subplot(2,3,bb);
        p1 = semilogx(data.Naive.(behavior).(dataType).mean_f,data.Naive.(behavior).(dataType).mean_C,'color',colors('north texas green'),'LineWidth',2);
        hold on
        semilogx(data.Naive.(behavior).(dataType).mean_f,data.Naive.(behavior).(dataType).mean_C + data.Naive.(behavior).(dataType).stdErr_C,'color',colors('sapphire'),'LineWidth',0.5);
        semilogx(data.Naive.(behavior).(dataType).mean_f,data.Naive.(behavior).(dataType).mean_C - data.Naive.(behavior).(dataType).stdErr_C,'color',colors('sapphire'),'LineWidth',0.5);
        p2 = semilogx(data.Blank_SAP.(behavior).(dataType).mean_f,data.Blank_SAP.(behavior).(dataType).mean_C,'color',colors('north texas green'),'LineWidth',2);
        semilogx(data.Blank_SAP.(behavior).(dataType).mean_f,data.Blank_SAP.(behavior).(dataType).mean_C + data.Blank_SAP.(behavior).(dataType).stdErr_C,'color',colors('north texas green'),'LineWidth',0.5);
        semilogx(data.Blank_SAP.(behavior).(dataType).mean_f,data.Blank_SAP.(behavior).(dataType).mean_C - data.Blank_SAP.(behavior).(dataType).stdErr_C,'color',colors('north texas green'),'LineWidth',0.5);
        p3 = semilogx(data.SSP_SAP.(behavior).(dataType).mean_f,data.SSP_SAP.(behavior).(dataType).mean_C,'color',colors('electric purple'),'LineWidth',2);
        semilogx(data.SSP_SAP.(behavior).(dataType).mean_f,data.SSP_SAP.(behavior).(dataType).mean_C + data.SSP_SAP.(behavior).(dataType).stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
        semilogx(data.SSP_SAP.(behavior).(dataType).mean_f,data.SSP_SAP.(behavior).(dataType).mean_C - data.SSP_SAP.(behavior).(dataType).stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
        ylabel('Coherence^2')
        xlabel('Freq (Hz)')
        title(behavior)
        xlim(xlimits{1,bb})
        ylim([0,1])
        if bb == 1
            legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
        end
        set(gca,'box','off')
        axis square
    end
end
