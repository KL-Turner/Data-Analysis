function [] = StimEvoked_GCaMP_Figures(rootFolder,saveFigs,delim)
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
solenoids = {'LPadSol','RPadSol','AudSol'};
hemispheres = {'LH','RH','fLH','fRH'};
dataTypes = {'HbT','HbO','HbR','GCaMP','timeVector','group','animalID'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(solenoids)
                solenoid = solenoids{1,dd};
                data.(group).(hemisphere).(solenoid).dummCheck = 1;
                for ee = 1:length(dataTypes)
                    dataType = dataTypes{1,ee};
                    if isfield(data.(group).(hemisphere).(solenoid),dataType) == false
                        if any(strcmp(dataType,{'group','animalID'})) == true
                            data.(group).(hemisphere).(solenoid).(dataType) = {};
                        else
                            data.(group).(hemisphere).(solenoid).(dataType) = [];
                        end
                    end
                end
                data.(group).(hemisphere).(solenoid).HbT = cat(1,data.(group).(hemisphere).(solenoid).HbT,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).HbT);
                data.(group).(hemisphere).(solenoid).HbO = cat(1,data.(group).(hemisphere).(solenoid).HbO,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).HbO);
                data.(group).(hemisphere).(solenoid).HbR = cat(1,data.(group).(hemisphere).(solenoid).HbR,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).HbR);
                data.(group).(hemisphere).(solenoid).GCaMP = cat(1,data.(group).(hemisphere).(solenoid).GCaMP,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).GCaMP*100);
                data.(group).(hemisphere).(solenoid).timeVector = cat(1,data.(group).(hemisphere).(solenoid).timeVector,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).timeVector);
                data.(group).(hemisphere).(solenoid).group = cat(1,data.(group).(hemisphere).(solenoid).group,group);
                data.(group).(hemisphere).(solenoid).animalID = cat(1,data.(group).(hemisphere).(solenoid).animalID,animalID);
            end
        end
    end
end
% pair stimulation with hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(solenoids)
            solenoid = solenoids{1,cc};
            [comparison] = FindSolenoidComparison(hemisphere,solenoid);
            data.(group).(hemisphere).(comparison).group = {};
            data.(group).(hemisphere).(comparison).animalID = {};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if any(strcmp(dataType,{'group','animalID'})) == true
                    data.(group).(hemisphere).(comparison).(dataType) = data.(group).(hemisphere).(solenoid).(dataType);
                else
                    data.(group).(hemisphere).(comparison).(dataType) = data.(group).(hemisphere).(solenoid).(dataType);
                    data.(group).(hemisphere).(comparison).(['mean_' dataType]) = mean(data.(group).(hemisphere).(solenoid).(dataType),1);
                    data.(group).(hemisphere).(comparison).(['stdErr_' dataType]) = std(data.(group).(hemisphere).(solenoid).(dataType),1)./sqrt(size(data.(group).(hemisphere).(solenoid).(dataType),1));
                    if strcmp(dataType,'timeVector') == false
                        for ee = 1:size(data.(group).(hemisphere).(comparison).(dataType),1)
                            startIdx = find(data.(group).(hemisphere).(solenoid).timeVector(ee,:) == 2);
                            endIdx =  find(data.(group).(hemisphere).(solenoid).timeVector(ee,:) == 4);
                            timeSnip = data.(group).(hemisphere).(solenoid).timeVector(ee,:);
                            dataSnip = data.(group).(hemisphere).(comparison).(dataType)(ee,:);
                            data.(group).(hemisphere).(comparison).(['AUC_' dataType])(ee,1) = trapz(timeSnip(startIdx:endIdx),dataSnip(startIdx:endIdx));
                            [maxData,maxIdx] = max(dataSnip);
                            data.(group).(hemisphere).(comparison).(['Peak_' dataType])(ee,1) = maxData;
                            data.(group).(hemisphere).(comparison).(['TTP_' dataType])(ee,1) = timeSnip(maxIdx);
                        end
                    end
                end
            end
        end
    end
end
% figure
dataTypes = {'HbT','HbO','HbR','GCaMP'};
for bb = 1:length(dataTypes)
    dataType = dataTypes{1,bb};
    summaryFigure = figure;
    sgtitle(['Contralateral whisker stimlation ' dataType ' [GCaMP]'])
    for aa = 1:length(hemispheres)
        subplot(2,2,aa);
        hemisphere = hemispheres{1,aa};
        hold on
        % Blank-SAP
        p1 = plot(data.Blank_SAP.(hemisphere).contra.mean_timeVector,data.Blank_SAP.(hemisphere).contra.(['mean_' dataType]),'color',colors('north texas green'),'LineWidth',2);
        plot(data.Blank_SAP.(hemisphere).contra.mean_timeVector,data.Blank_SAP.(hemisphere).contra.(['mean_' dataType]) + data.Blank_SAP.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('north texas green'),'LineWidth',0.5)
        plot(data.Blank_SAP.(hemisphere).contra.mean_timeVector,data.Blank_SAP.(hemisphere).contra.(['mean_' dataType]) - data.Blank_SAP.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('north texas green'),'LineWidth',0.5)
        % SSP-SAP
        p2 = plot(data.SSP_SAP.(hemisphere).contra.mean_timeVector,data.SSP_SAP.(hemisphere).contra.(['mean_' dataType]),'color',colors('electric purple'),'LineWidth',2);
        plot(data.SSP_SAP.(hemisphere).contra.mean_timeVector,data.SSP_SAP.(hemisphere).contra.(['mean_' dataType]) + data.SSP_SAP.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('electric purple'),'LineWidth',0.5)
        plot(data.SSP_SAP.(hemisphere).contra.mean_timeVector,data.SSP_SAP.(hemisphere).contra.(['mean_' dataType]) - data.SSP_SAP.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('electric purple'),'LineWidth',0.5)
        title([hemisphere ' ' dataType])
        if strcmp(dataType,'GCaMP') == true
            label = '\DeltaF/F (%)';
        else
            label = ['\Delta' dataType ' (\muM)'];
        end
        ylabel(label)
        xlabel('Peri-stimulus time (s)')
        if aa == 1
            legend([p1,p2],'Blank-SAP','SSP-SAP')
        end
        set(gca,'box','off')
        xlim([-2,10])
        axis square
        linkaxes
        % save figure(s)
        if saveFigs == true
            dirpath = [rootFolder delim 'Summary Figures' delim 'Stim Evoked GCaMP' delim];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(summaryFigure,[dirpath 'StimEvoked_GCaMP_' dataType]);
        end
    end
end
% figure
statsVariables = {'AUC','Peak','TTP'};
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    summaryFigure = figure;
    sgtitle([dataType ' stats [GCaMP]'])
    for bb = 1:length(statsVariables)
        statsVariable = statsVariables{1,bb};
        subplot(1,3,bb);
        xInds = ones(1,length(data.Blank_SAP.RH.contra.([statsVariable '_' dataType])));
        s1 = scatter(xInds*1,data.Blank_SAP.RH.contra.(([statsVariable '_' dataType])),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
        hold on
        e1 = errorbar(1,mean(data.Blank_SAP.RH.contra.(([statsVariable '_' dataType]))),std(data.Blank_SAP.RH.contra.(([statsVariable '_' dataType]))),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e1.Color = 'black';
        e1.MarkerSize = 10;
        e1.CapSize = 10;
        xInds = ones(1,length(data.SSP_SAP.RH.contra.(([statsVariable '_' dataType]))));
        s2 = scatter(xInds*2,data.SSP_SAP.RH.contra.(([statsVariable '_' dataType])),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
        hold on
        e2 = errorbar(2,mean(data.SSP_SAP.RH.contra.(([statsVariable '_' dataType]))),std(data.SSP_SAP.RH.contra.(([statsVariable '_' dataType]))),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e2.Color = 'black';
        e2.MarkerSize = 10;
        e2.CapSize = 10;
        ylabel(strrep(statsVariable,'_',' '))
        title(strrep(statsVariable,'_',' '))
        xlim([0,3])
        if aa == 1
            legend([s1,s2],'Blank-SAP','SSP-SAP')
        end
        set(gca,'box','off')
        set(gca,'xtick',[])
        axis square
    end
    % save figure(s)
    if saveFigs == true
        savefig(summaryFigure,[dirpath 'StimEvoked_GCaMP_HbT_StatsFigure']);
    end
end
% statistics - generalized linear mixed effects model
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    for bb = 1:length(statsVariables)
        statsVariable = statsVariables{1,bb};
        Stats.(dataType).(statsVariable).tableSize = cat(1,data.Blank_SAP.RH.contra.([statsVariable '_' dataType]),data.SSP_SAP.RH.contra.([statsVariable '_' dataType]));
        Stats.(dataType).(statsVariable).Table = table('Size',[size(Stats.(dataType).(statsVariable).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'AnimalID','Group','Data'});
        Stats.(dataType).(statsVariable).Table.AnimalID = cat(1,data.Blank_SAP.RH.contra.animalID,data.SSP_SAP.RH.contra.animalID);
        Stats.(dataType).(statsVariable).Table.Group = cat(1,data.Blank_SAP.RH.contra.group,data.SSP_SAP.RH.contra.group);
        Stats.(dataType).(statsVariable).Table.Data = cat(1,data.Blank_SAP.RH.contra.([statsVariable '_' dataType]),data.SSP_SAP.RH.contra.([statsVariable '_' dataType]));
        Stats.(dataType).(statsVariable).FitFormula = 'Data ~ 1 + Group + (1|AnimalID)';
        Stats.(dataType).(statsVariable).Stats = fitglme(Stats.(dataType).(statsVariable).Table,Stats.(dataType).(statsVariable).FitFormula);
    end
end
% statistical diary
if saveFigs == true
    % statistical diary
    diaryFile = [dirpath 'StimEvoked_GCaMP_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        for bb = 1:length(statsVariables)
            statsVariable = statsVariables{1,bb};
            disp('======================================================================================================================')
            disp(['GLME statistics: ' dataType ' ' statsVariable])
            disp('======================================================================================================================')
            disp(Stats.(dataType).(statsVariable).Stats)
            disp('----------------------------------------------------------------------------------------------------------------------')
        end
    end
    diary off
end