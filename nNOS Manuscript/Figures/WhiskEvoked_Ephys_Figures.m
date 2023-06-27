function [] = WhiskEvoked_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Evoked_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Naive','Blank_SAP','SSP_SAP'};
whiskTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
hemispheres = {'LH','RH'};
dataTypes = {'HbT','cortMUA','cortGam','cortS','T','F','timeVector','group','animalID'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_Ephys.(group));
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
                data.(group).(hemisphere).(whiskType).HbT = cat(1,data.(group).(hemisphere).(whiskType).HbT,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(whiskType).HbT);
                data.(group).(hemisphere).(whiskType).cortMUA = cat(1,data.(group).(hemisphere).(whiskType).cortMUA,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(whiskType).cortMUA);
                data.(group).(hemisphere).(whiskType).cortGam = cat(1,data.(group).(hemisphere).(whiskType).cortGam,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(whiskType).cortGam);
                data.(group).(hemisphere).(whiskType).cortS = cat(3,data.(group).(hemisphere).(whiskType).cortS,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(whiskType).cortLFP);
                data.(group).(hemisphere).(whiskType).T = cat(1,data.(group).(hemisphere).(whiskType).T,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(whiskType).T);
                data.(group).(hemisphere).(whiskType).F = cat(1,data.(group).(hemisphere).(whiskType).F,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(whiskType).F);
                data.(group).(hemisphere).(whiskType).timeVector = cat(1,data.(group).(hemisphere).(whiskType).timeVector,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(whiskType).timeVector);
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
            data.(group).(hemisphere).(whiskType).group = {};
            data.(group).(hemisphere).(whiskType).animalID = {};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if any(strcmp(dataType,{'group','animalID'})) == true
                    data.(group).(hemisphere).(whiskType).(dataType) = data.(group).(hemisphere).(whiskType).(dataType);
                elseif strcmp(dataType,'cortS')
                    data.(group).(hemisphere).(whiskType).(dataType) = data.(group).(hemisphere).(whiskType).(dataType);
                    data.(group).(hemisphere).(whiskType).(['mean_' dataType]) = mean(data.(group).(hemisphere).(whiskType).(dataType),3).*100;
                else
                    data.(group).(hemisphere).(whiskType).(dataType) = data.(group).(hemisphere).(whiskType).(dataType);
                    data.(group).(hemisphere).(whiskType).(['mean_' dataType]) = mean(data.(group).(hemisphere).(whiskType).(dataType),1);
                    data.(group).(hemisphere).(whiskType).(['stdErr_' dataType]) = std(data.(group).(hemisphere).(whiskType).(dataType),1)./sqrt(size(data.(group).(hemisphere).(whiskType).(dataType),1));
                    if strcmp(dataType,'HbT') == true
                        for ee = 1:size(data.(group).(hemisphere).(whiskType).HbT,1)
                            startIdx = find(data.(group).(hemisphere).(whiskType).timeVector(ee,:) == 2);
                            endIdx =  find(data.(group).(hemisphere).(whiskType).timeVector(ee,:) == 4);
                            timeSnip = data.(group).(hemisphere).(whiskType).timeVector(ee,:);
                            hbtSnip = data.(group).(hemisphere).(whiskType).HbT(ee,:);
                            data.(group).(hemisphere).(whiskType).AUC_HbT(ee,1) = trapz(timeSnip(startIdx:endIdx),hbtSnip(startIdx:endIdx));
                            [maxHbT,maxIdx] = max(hbtSnip);
                            data.(group).(hemisphere).(whiskType).Peak_HbT(ee,1) = maxHbT;
                            data.(group).(hemisphere).(whiskType).TTP_HbT(ee,1) = timeSnip(maxIdx);
                        end
                    end
                end
            end
        end
    end
end
% figure
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(whiskTypes)
        whiskType = whiskTypes{1,bb};
        summaryFigure = figure;
        sgtitle([hemisphere ' ' whiskType ' (HbT) [Ephys]'])
        dataType = 'HbT';
        subplot(1,3,bb)
        p1 = plot(data.Naive.(hemisphere).(whiskType).mean_timeVector,data.Naive.(hemisphere).(whiskType).(['mean_' dataType]),'color',colors('sapphire'),'LineWidth',2);
        hold on;
        plot(data.Naive.(hemisphere).(whiskType).mean_timeVector,data.Naive.(hemisphere).(whiskType).(['mean_' dataType]) + data.Naive.(hemisphere).(whiskType).(['stdErr_' dataType]),'color',colors('sapphire'),'LineWidth',0.25)
        plot(data.Naive.(hemisphere).(whiskType).mean_timeVector,data.Naive.(hemisphere).(whiskType).(['mean_' dataType]) - data.Naive.(hemisphere).(whiskType).(['stdErr_' dataType]),'color',colors('sapphire'),'LineWidth',0.25)
        p2 = plot(data.Blank_SAP.(hemisphere).(whiskType).mean_timeVector,data.Blank_SAP.(hemisphere).(whiskType).(['mean_' dataType]),'color',colors('north texas green'),'LineWidth',2);
        plot(data.Blank_SAP.(hemisphere).(whiskType).mean_timeVector,data.Blank_SAP.(hemisphere).(whiskType).(['mean_' dataType]) + data.Blank_SAP.(hemisphere).(whiskType).(['stdErr_' dataType]),'color',colors('north texas green'),'LineWidth',0.25)
        plot(data.Blank_SAP.(hemisphere).(whiskType).mean_timeVector,data.Blank_SAP.(hemisphere).(whiskType).(['mean_' dataType]) - data.Blank_SAP.(hemisphere).(whiskType).(['stdErr_' dataType]),'color',colors('north texas green'),'LineWidth',0.25)
        p3 = plot(data.SSP_SAP.(hemisphere).(whiskType).mean_timeVector,data.SSP_SAP.(hemisphere).(whiskType).(['mean_' dataType]),'color',colors('electric purple'),'LineWidth',2);
        plot(data.SSP_SAP.(hemisphere).(whiskType).mean_timeVector,data.SSP_SAP.(hemisphere).(whiskType).(['mean_' dataType]) + data.SSP_SAP.(hemisphere).(whiskType).(['stdErr_' dataType]),'color',colors('electric purple'),'LineWidth',0.25)
        plot(data.SSP_SAP.(hemisphere).(whiskType).mean_timeVector,data.SSP_SAP.(hemisphere).(whiskType).(['mean_' dataType]) - data.SSP_SAP.(hemisphere).(whiskType).(['stdErr_' dataType]),'color',colors('electric purple'),'LineWidth',0.25)
        title([hemisphere ' ' dataType])
        label = '\DeltaHbT (\muM)';
        ylabel(label)
        xlabel('Peri-whisk time (s)')
        if aa == 1
            legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
        end
        set(gca,'box','off')
        xlim([-2,10])
        axis square
    end
    % save figure(s)
    if saveFigs == true
        dirpath = [rootFolder delim 'Summary Figures' delim 'Volitional Whisking' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigure,[dirpath 'VolitionalWhisks_Ephys_' hemisphere]);
    end
end
% figure
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(whiskTypes)
            whiskType = whiskTypes{1,cc};
            summaryFigure = figure;
            sgtitle([group ' ' hemisphere ' ' whiskType ' whisker stimlation (LFP) [Ephys]'])
            subplot(1,3,cc);
            imagesc(data.(group).(hemisphere).(whiskType).mean_T,data.(group).(hemisphere).(whiskType).mean_F,data.(group).(hemisphere).(whiskType).mean_cortS)
            title(whiskType)
            ylabel('Freq (Hz)')
            xlabel('Peri-stimulus time (s)')
            c1 = colorbar;
            ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
            caxis([-100,100])
            set(gca,'Ticklength',[0,0])
            axis xy
            set(gca,'box','off')
            axis square
        end
        % save figure(s)
        if saveFigs == true
            savefig(summaryFigure,[dirpath 'VolitionalWhisks_Ephys_LFP_' group '_' hemisphere]);
        end
    end
end
% statistics - generalized linear mixed effects model
for aa = 1:length(statsVariables)
    statsVariable = statsVariables{1,aa};
    % statistics - generalized linear mixed effects model
    Stats.(statsVariable).tableSize = cat(1,data.Blank_SAP.RH.contra.(statsVariable),data.SSP_SAP.RH.contra.(statsVariable));
    Stats.(statsVariable).Table = table('Size',[size(Stats.(statsVariable).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'AnimalID','Group','Coherence'});
    Stats.(statsVariable).Table.AnimalID = cat(1,data.Blank_SAP.RH.contra.animalID,data.SSP_SAP.RH.contra.animalID);
    Stats.(statsVariable).Table.Group = cat(1,data.Blank_SAP.RH.contra.group,data.SSP_SAP.RH.contra.group);
    Stats.(statsVariable).Table.Data = cat(1,data.Blank_SAP.RH.contra.(statsVariable),data.SSP_SAP.RH.contra.(statsVariable));
    Stats.(statsVariable).FitFormula = 'Data ~ 1 + Group + (1|AnimalID)';
    Stats.(statsVariable).Stats = fitglme(Stats.(statsVariable).Table,Stats.(statsVariable).FitFormula);
end
% statistical diary
if saveFigs == true
    % statistical diary
    diaryFile = [dirpath 'StimEvoked_Ephys_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    for aa = 1:length(statsVariables)
        statsVariable = statsVariables{1,aa};
        disp('======================================================================================================================')
        disp(['GLME statistics: ' statsVariable])
        disp('======================================================================================================================')
        disp(Stats.(statsVariable).Stats)
        disp('----------------------------------------------------------------------------------------------------------------------')
    end
    diary off
end