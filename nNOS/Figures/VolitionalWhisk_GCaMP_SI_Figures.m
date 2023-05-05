function [] = VolitionalWhisk_GCaMP_SI_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Evoked_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Blank_SAP','SSP_SAP'};
whiskTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
hemispheres = {'LH','RH'};
dataTypes = {'HbT','HbO','HbR','GCaMP','timeVector','group','animalID'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(whiskTypes)
                whiskType = whiskTypes{1,dd};
                data.(group).(hemisphere).(whiskType).dummCheck = 1;
                for ee = 1:length(dataTypes)
                    dataType = dataTypes{1,ee};
                    if isfield(data.(group).(hemisphere).(whiskType),dataType) == false
                        if any(strcmp(dataType,{'group','animalID'})) == true
                            data.(group).(hemisphere).(whiskType).(dataType) = {};
                        else
                            data.(group).(hemisphere).(whiskType).(dataType) = [];
                        end
                    end
                end
                data.(group).(hemisphere).(whiskType).HbT = cat(1,data.(group).(hemisphere).(whiskType).HbT,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Whisk.(whiskType).HbT);
                data.(group).(hemisphere).(whiskType).HbO = cat(1,data.(group).(hemisphere).(whiskType).HbO,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Whisk.(whiskType).HbO);
                data.(group).(hemisphere).(whiskType).HbR = cat(1,data.(group).(hemisphere).(whiskType).HbR,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Whisk.(whiskType).HbR);
                data.(group).(hemisphere).(whiskType).GCaMP = cat(1,data.(group).(hemisphere).(whiskType).GCaMP,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Whisk.(whiskType).GCaMP*100);
                data.(group).(hemisphere).(whiskType).timeVector = cat(1,data.(group).(hemisphere).(whiskType).timeVector,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Whisk.(whiskType).timeVector);
                data.(group).(hemisphere).(whiskType).group = cat(1,data.(group).(hemisphere).(whiskType).group,group);
                data.(group).(hemisphere).(whiskType).animalID = cat(1,data.(group).(hemisphere).(whiskType).animalID,animalID);
            end
        end
    end
end
% pair stimulation with hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(whiskTypes)
            whiskType = whiskTypes{1,cc};
            [whiskType] = FindSolenoidComparison(hemisphere,whiskType);
            data.(group).(hemisphere).(whiskType).group = {};
            data.(group).(hemisphere).(whiskType).animalID = {};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if any(strcmp(dataType,{'group','animalID'})) == true
                    data.(group).(hemisphere).(whiskType).(dataType) = data.(group).(hemisphere).(whiskType).(dataType);
                else
                    data.(group).(hemisphere).(whiskType).(dataType) = data.(group).(hemisphere).(whiskType).(dataType);
                    data.(group).(hemisphere).(whiskType).(['mean_' dataType]) = mean(data.(group).(hemisphere).(whiskType).(dataType),1);
                    data.(group).(hemisphere).(whiskType).(['stdErr_' dataType]) = std(data.(group).(hemisphere).(whiskType).(dataType),1)./sqrt(size(data.(group).(hemisphere).(whiskType).(dataType),1));
                    if strcmp(dataType,'timeVector') == false
                        for ee = 1:size(data.(group).(hemisphere).(whiskType).(dataType),1)
                            startIdx = find(data.(group).(hemisphere).(whiskType).timeVector(ee,:) == 2);
                            endIdx =  find(data.(group).(hemisphere).(whiskType).timeVector(ee,:) == 4);
                            timeSnip = data.(group).(hemisphere).(whiskType).timeVector(ee,:);
                            dataSnip = data.(group).(hemisphere).(whiskType).(dataType)(ee,:);
                            data.(group).(hemisphere).(whiskType).(['AUC_' dataType])(ee,1) = trapz(timeSnip(startIdx:endIdx),dataSnip(startIdx:endIdx));
                            [maxData,maxIdx] = max(dataSnip);
                            data.(group).(hemisphere).(whiskType).(['Peak_' dataType])(ee,1) = maxData;
                            data.(group).(hemisphere).(whiskType).(['TTP_' dataType])(ee,1) = timeSnip(maxIdx);
                        end
                    end
                end
            end
        end
    end
end
% figure
dataTypes = {'HbT','HbO','HbR','GCaMP'};
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        summaryFigure = figure;
        sgtitle([hemisphere ' ' dataType ' whisker stimlation [GCaMP SI]']);
        for cc = 1:length(whiskTypes)
            whiskType = whiskTypes{1,cc};
            subplot(1,3,cc);
            p1 = plot(data.Blank_SAP.(hemisphere).(whiskType).mean_timeVector,data.Blank_SAP.(hemisphere).(whiskType).(['mean_' dataType]),'color',colors('north texas green'),'LineWidth',2);
            hold on;
            plot(data.Blank_SAP.(hemisphere).(whiskType).mean_timeVector,data.Blank_SAP.(hemisphere).(whiskType).(['mean_' dataType]) + data.Blank_SAP.(hemisphere).(whiskType).(['stdErr_' dataType]),'color',colors('north texas green'),'LineWidth',0.25)
            plot(data.Blank_SAP.(hemisphere).(whiskType).mean_timeVector,data.Blank_SAP.(hemisphere).(whiskType).(['mean_' dataType]) - data.Blank_SAP.(hemisphere).(whiskType).(['stdErr_' dataType]),'color',colors('north texas green'),'LineWidth',0.25)
            p2 = plot(data.SSP_SAP.(hemisphere).(whiskType).mean_timeVector,data.SSP_SAP.(hemisphere).(whiskType).(['mean_' dataType]),'color',colors('electric purple'),'LineWidth',2);
            plot(data.SSP_SAP.(hemisphere).(whiskType).mean_timeVector,data.SSP_SAP.(hemisphere).(whiskType).(['mean_' dataType]) + data.SSP_SAP.(hemisphere).(whiskType).(['stdErr_' dataType]),'color',colors('electric purple'),'LineWidth',0.25)
            plot(data.SSP_SAP.(hemisphere).(whiskType).mean_timeVector,data.SSP_SAP.(hemisphere).(whiskType).(['mean_' dataType]) - data.SSP_SAP.(hemisphere).(whiskType).(['stdErr_' dataType]),'color',colors('electric purple'),'LineWidth',0.25)
            title([hemisphere ' ' dataType])
            if strcmp(dataType,'GCaMP') == true
                label = '\DeltaF/F (%)';
            else
                label = ['\Delta' dataType ' (\muM)'];
            end
            ylabel(label)
            xlabel('Peri-stimulus time (s)')
            if cc == 1
                legend([p1,p2],'Blank-SAP','SSP-SAP')
            end
            set(gca,'box','off')
            xlim([-2,10])
            axis square
            linkaxes
            % save figure(s)
            if saveFigs == true
                dirpath = [rootFolder delim 'Summary Figures' delim 'Volitional Whisk' delim];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(summaryFigure,[dirpath 'VolitionalWhisk_GCaMP_SI_' hemisphere '_' dataType '_' whiskType]);
            end
        end
    end
end