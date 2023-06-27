function [] = VolitionalWhisk_Pulse_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Evoked_Pulse';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Blank_SAP','SSP_SAP'};
whiskTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
dataTypes = {'HbT','timeVector'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_Pulse.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(whiskTypes)
            whiskType = whiskTypes{1,cc};
            data.(group).(whiskType).dummCheck = 1;
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if isfield(data.(group).(whiskType),dataType) == false
                    data.(group).(whiskType).(dataType) = [];
                end
            end
            data.(group).(whiskType).HbT = cat(1,data.(group).(whiskType).HbT,Results_Evoked_Pulse.(group).(animalID).Whisk.(whiskType).HbT);
            data.(group).(whiskType).timeVector = cat(1,data.(group).(whiskType).timeVector,Results_Evoked_Pulse.(group).(animalID).Whisk.(whiskType).timeVector);
        end
    end
end
% pair stimulation types with hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(whiskTypes)
        whiskType = whiskTypes{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            data.(group).(whiskType).(dataType) = data.(group).(whiskType).(dataType);
            data.(group).(whiskType).(['mean_'  dataType]) = mean(data.(group).(whiskType).(dataType),1);
            data.(group).(whiskType).(['stdErr_' dataType]) = std(data.(group).(whiskType).(dataType),1)./sqrt(size(data.(group).(whiskType).(dataType),1));
        end
    end
end
% figure
summaryFigure = figure;
title('Volitional whisking [Pulse]')
for aa = 1:length(whiskTypes)
    whiskType = whiskTypes{1,aa};
    subplot(1,3,aa)
    p1 = plot(data.Blank_SAP.(whiskType).mean_timeVector,data.Blank_SAP.(whiskType).mean_HbT,'color',colors('north texas green'),'LineWidth',2);
    hold on;
    plot(data.Blank_SAP.(whiskType).mean_timeVector,data.Blank_SAP.(whiskType).mean_HbT + data.Blank_SAP.(whiskType).stdErr_HbT,'color',colors('north texas green'),'LineWidth',0.25)
    plot(data.Blank_SAP.(whiskType).mean_timeVector,data.Blank_SAP.(whiskType).mean_HbT - data.Blank_SAP.(whiskType).stdErr_HbT,'color',colors('north texas green'),'LineWidth',0.25)
    p2 = plot(data.SSP_SAP.(whiskType).mean_timeVector,data.SSP_SAP.(whiskType).mean_HbT,'color',colors('electric purple'),'LineWidth',2);
    plot(data.SSP_SAP.(whiskType).mean_timeVector,data.SSP_SAP.(whiskType).mean_HbT + data.SSP_SAP.(whiskType).stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
    plot(data.SSP_SAP.(whiskType).mean_timeVector,data.SSP_SAP.(whiskType).mean_HbT - data.SSP_SAP.(whiskType).stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
    title(whiskType)
    ylabel('\DeltaHbT (\muM)')
    xlabel('Peri-whisk time (s)')
    if aa == 1
        legend([p1,p2],'Blank-SAP','SSP-SAP')
    end
    set(gca,'box','off')
    xlim([-2,10])
    axis square
end
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures' delim 'Volitional Whisk' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'VolitionalWhisk_Pulse']);
end