function [] = PupilCrossCorrelation_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PupilCrossCorr_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Naive','Blank_SAP','SSP_SAP'};
compDataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
dataTypes = {'LH_HbT','RH_HbT','LH_gammaBandPower','RH_gammaBandPower'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
samplingRate = 30;
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PupilCrossCorr_Ephys.(group));
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
                        data.(group).(dataType).(compDataType).(behavior).lags = [];
                        data.(group).(dataType).(compDataType).(behavior).xcVals = [];
                    end
                    if isempty(Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(compDataType).(behavior).xcVals) == false
                        data.(group).(dataType).(compDataType).(behavior).lags = cat(1,data.(group).(dataType).(compDataType).(behavior).lags,Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(compDataType).(behavior).lags/samplingRate);
                        data.(group).(dataType).(compDataType).(behavior).xcVals = cat(1,data.(group).(dataType).(compDataType).(behavior).xcVals,Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(compDataType).(behavior).xcVals);
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
                data.(group).(dataType).(compDataType).(behavior).mean_xcVals = mean(data.(group).(dataType).(compDataType).(behavior).xcVals,1);
                data.(group).(dataType).(compDataType).(behavior).stdErr_xcVals = std(data.(group).(dataType).(compDataType).(behavior).xcVals,0,1)./sqrt(size(data.(group).(dataType).(compDataType).(behavior).xcVals,1));
                data.(group).(dataType).(compDataType).(behavior).mean_lags = mean(data.(group).(dataType).(compDataType).(behavior).lags,1);
            end
        end
    end
end
% figure
xlimits = ({[-5,5],[-10,10],[-10,10],[-30,30],[-30,30],[-30,30]});
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    for bb = 1:length(compDataTypes)
        compDataType = compDataTypes{1,bb};
        summaryFigure = figure;
        sgtitle([strrep(dataType,'_',' ') ' cross correlation with ' compDataType ' [Ephys]'])
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            subplot(2,3,cc);
            p1 = plot(data.Naive.(dataType).(compDataType).(behavior).mean_lags,data.Naive.(dataType).(compDataType).(behavior).mean_xcVals,'color',colors('sapphire'),'LineWidth',2);
            hold on;
            plot(data.Naive.(dataType).(compDataType).(behavior).mean_lags,data.Naive.(dataType).(compDataType).(behavior).mean_xcVals + data.Naive.(dataType).(compDataType).(behavior).stdErr_xcVals,'color',colors('sapphire'),'LineWidth',0.25);
            plot(data.Naive.(dataType).(compDataType).(behavior).mean_lags,data.Naive.(dataType).(compDataType).(behavior).mean_xcVals - data.Naive.(dataType).(compDataType).(behavior).stdErr_xcVals,'color',colors('sapphire'),'LineWidth',0.25);
            p2 = plot(data.Blank_SAP.(dataType).(compDataType).(behavior).mean_lags,data.Blank_SAP.(dataType).(compDataType).(behavior).mean_xcVals,'color',colors('north texas green'),'LineWidth',2);
            plot(data.Blank_SAP.(dataType).(compDataType).(behavior).mean_lags,data.Blank_SAP.(dataType).(compDataType).(behavior).mean_xcVals + data.Blank_SAP.(dataType).(compDataType).(behavior).stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
            plot(data.Blank_SAP.(dataType).(compDataType).(behavior).mean_lags,data.Blank_SAP.(dataType).(compDataType).(behavior).mean_xcVals - data.Blank_SAP.(dataType).(compDataType).(behavior).stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
            p3 = plot(data.SSP_SAP.(dataType).(compDataType).(behavior).mean_lags,data.SSP_SAP.(dataType).(compDataType).(behavior).mean_xcVals,'color',colors('electric purple'),'LineWidth',2);
            plot(data.SSP_SAP.(dataType).(compDataType).(behavior).mean_lags,data.SSP_SAP.(dataType).(compDataType).(behavior).mean_xcVals + data.SSP_SAP.(dataType).(compDataType).(behavior).stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
            plot(data.SSP_SAP.(dataType).(compDataType).(behavior).mean_lags,data.SSP_SAP.(dataType).(compDataType).(behavior).mean_xcVals - data.SSP_SAP.(dataType).(compDataType).(behavior).stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
            ylabel('Corr Coef')
            xlabel('Lags (s)')
            xlim(xlimits{1,cc})
            title(behavior)
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
            savefig(summaryFigure,[dirpath 'PupilCrossCorr_Ephys_' dataType '_' compDataType]);
        end
    end
end