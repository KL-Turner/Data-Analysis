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
dataTypes = {'gammaBandPower'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'C','f'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_NeuralHemoCoher_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for qq = 1:length(hemispheres)
            hemisphere = hemispheres{1,qq};
            for cc = 1:length(dataTypes)
                dataType = dataTypes{1,cc};
                data.(group).(hemisphere).(dataType).dummCheck = 1;
                for dd = 1:length(behaviors)
                    behavior = behaviors{1,dd};
                    if isfield(data.(group).(hemisphere).(dataType),behavior) == false
                        data.(group).(hemisphere).(dataType).(behavior).C = [];
                        data.(group).(hemisphere).(dataType).(behavior).f = [];
                        data.(group).(hemisphere).(dataType).(behavior).group = {};
                        data.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    if isempty(Results_NeuralHemoCoher_Ephys.(group).(animalID).(dataType).(hemisphere).(behavior).C) == false
                        data.(group).(hemisphere).(dataType).(behavior).C = cat(1,data.(group).(hemisphere).(dataType).(behavior).C,Results_NeuralHemoCoher_Ephys.(group).(animalID).(dataType).(hemisphere).(behavior).C.^2');
                        data.(group).(hemisphere).(dataType).(behavior).f = cat(1,data.(group).(hemisphere).(dataType).(behavior).f,Results_NeuralHemoCoher_Ephys.(group).(animalID).(dataType).(hemisphere).(behavior).f);
                        data.(group).(hemisphere).(dataType).(behavior).group = cat(1,data.(group).(hemisphere).(dataType).(behavior).group,group);
                        data.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,data.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
                    end
                end
            end
        end
    end
end
% mean/stdanimalID
for aa = 1:length(groups)
    group = groups{1,aa};
    for qq = 1:length(hemispheres)
        hemisphere = hemispheres{1,qq};
        for bb = 1:length(behaviors)
            behavior = behaviors{1,bb};
            for cc = 1:length(dataTypes)
                dataType = dataTypes{1,cc};
                for dd = 1:length(variables)
                    variable = variables{1,dd};
                    data.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(data.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    data.(group).(hemisphere).(dataType).(behavior).(['stdErr_' variable]) = std(data.(group).(hemisphere).(dataType).(behavior).(variable),1)./sqrt(size(data.(group).(hemisphere).(dataType).(behavior).(variable),1));
                end
            end
        end
    end
end
% coherence figures
xlimits = ({[1/10,0.5],[1/30,0.5],[1/60,0.5],[0.003,0.5],[0.003,0.5],[0.003,0.5]});
for qq = 1:length(hemispheres)
    hemisphere = hemispheres{1,qq};
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        summaryFigure = figure;
        sgtitle([hemisphere ' ' dataType ' neural-hemo coherence [Ephys]'])
        for bb = 1:length(behaviors)
            behavior = behaviors{1,bb};
            subplot(2,3,bb);
            p1 = semilogx(data.Naive.(hemisphere).(dataType).(behavior).mean_f,data.Naive.(hemisphere).(dataType).(behavior).mean_C,'color',colors('sapphire'),'LineWidth',2);
            hold on
            semilogx(data.Naive.(hemisphere).(dataType).(behavior).mean_f,data.Naive.(hemisphere).(dataType).(behavior).mean_C + data.Naive.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('north texas green'),'LineWidth',0.5);
            semilogx(data.Naive.(hemisphere).(dataType).(behavior).mean_f,data.Naive.(hemisphere).(dataType).(behavior).mean_C - data.Naive.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('north texas green'),'LineWidth',0.5);
            p2 = semilogx(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_f,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_C,'color',colors('north texas green'),'LineWidth',2);
            semilogx(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_f,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_C + data.Blank_SAP.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('north texas green'),'LineWidth',0.5);
            semilogx(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_f,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_C - data.Blank_SAP.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('north texas green'),'LineWidth',0.5);
            p3 = semilogx(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_f,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_C,'color',colors('electric purple'),'LineWidth',2);
            semilogx(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_f,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_C + data.SSP_SAP.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
            semilogx(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_f,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_C - data.SSP_SAP.(hemisphere).(dataType).(behavior).stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
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
        % save figure(s)
        if saveFigs == true
            dirpath = [rootFolder delim 'Summary Figures' delim 'Neural-Hemo Coherence Ephys' delim];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(summaryFigure,[dirpath 'NeuralHemoCoherence_Ephys_' hemisphere '_' dataType]);
        end
    end
end
% % find Hz peaks in coherence
% arousalStates = {'Alert','Asleep','All'};
% for aa = 1:length(groups)
%     group = groups{1,aa};
%     for bb = 1:length(dataTypes)
%         dataType = dataTypes{1,bb};
%         for cc = 1:length(arousalStates)
%             arousalState = arousalStates{1,cc};
%             for dd = 1:size(data.(group).(hemisphere).(dataType).(arousalState).C,1)
%                 F = round(data.(group).(hemisphere).(dataType).(arousalState).f(dd,:),3);
%                 C = data.(group).(hemisphere).(dataType).(arousalState).C(dd,:);
%                 index001 = find(F == 0.01);
%                 index01 = find(F == 0.1);
%                 index05 = find(F == 0.5);
%                 data.(group).(hemisphere).(dataType).(arousalState).C001(dd,1) = mean(C(1:index001(1)));
%                 data.(group).(hemisphere).(dataType).(arousalState).C01(dd,1) = mean(C(index001(1) + 1:index01(1)));
%                 data.(group).(hemisphere).(dataType).(arousalState).C05(dd,1) = mean(C(index01(1) + 1:index05(1)));
%             end
%         end
%     end
% end
% % peak coherence figures
% freqs = {'C001','C01','C05'};
% for aa = 1:length(dataTypes)
%     dataType = dataTypes{1,aa};
%     figure;
%     sgtitle([dataType ' peak coherence [Ephys]'])
%     for bb = 1:length(arousalStates)
%         arousalState = arousalStates{1,bb};
%         subplot(1,3,bb); zz = 1;
%         for cc = 1:length(freqs)
%             freq = freqs{1,cc};
%             xInds = ones(1,length(data.Blank_SAP.(hemisphere).(dataType).(arousalState).(freq)));
%             s1 = scatter(xInds*zz,data.Blank_SAP.(hemisphere).(dataType).(arousalState).(freq),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
%             hold on
%             e1 = errorbar(zz,mean(data.Blank_SAP.(hemisphere).(dataType).(arousalState).(freq)),std(data.Blank_SAP.(hemisphere).(dataType).(arousalState).(freq)),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
%             e1.Color = 'black';
%             e1.MarkerSize = 10;
%             e1.CapSize = 10;
%             zz = zz + 1;
%             xInds = ones(1,length(data.SSP_SAP.(hemisphere).(dataType).(arousalState).(freq)));
%             s2 = scatter(xInds*zz,data.SSP_SAP.(hemisphere).(dataType).(arousalState).(freq),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
%             hold on
%             e2 = errorbar(zz,mean(data.SSP_SAP.(hemisphere).(dataType).(arousalState).(freq)),std(data.SSP_SAP.(hemisphere).(dataType).(arousalState).(freq)),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
%             e2.Color = 'black';
%             e2.MarkerSize = 10;
%             e2.CapSize = 10;
%             zz = zz + 1;
%             ylabel('Coherence^2')
%             xticks([.5,3.5,5.5])
%             xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
%             title(arousalState)
%             xlim([0,7])
%             ylim([0,1])
%         end
%         if bb == 1
%             legend([s1,s2],'Blank-SAP','SSP-SAP')
%         end
%         set(gca,'box','off')
%         axis square
%     end
%     % save figure(s)
%     if saveFigs == true
%         savefig(summaryFigure,[dirpath 'BilateralCoherence_Ephys_' dataType '_FreqPeaks']);
%     end
% end
% % statistics - generalized linear mixed effects model
% for aa = 1:length(freqs)
%     freq = freqs{1,aa};
%     for bb = 1:length(dataTypes)
%         dataType = dataTypes{1,bb};
%         for cc = 1:length(arousalStates)
%             arousalState = arousalStates{1,cc};
%             % statistics - generalized linear mixed effects model
%             Stats.(dataType).(arousalState).(freq).tableSize = cat(1,data.Blank_SAP.(hemisphere).(dataType).(arousalState).(freq),data.SSP_SAP.(hemisphere).(dataType).(arousalState).(freq));
%             Stats.(dataType).(arousalState).(freq).Table = table('Size',[size(Stats.(dataType).(arousalState).(freq).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'AnimalID','Group','Coherence'});
%             Stats.(dataType).(arousalState).(freq).Table.AnimalID = cat(1,data.Blank_SAP.(hemisphere).(dataType).(arousalState).animalID,data.SSP_SAP.(hemisphere).(dataType).(arousalState).animalID);
%             Stats.(dataType).(arousalState).(freq).Table.Group = cat(1,data.Blank_SAP.(hemisphere).(dataType).(arousalState).group,data.SSP_SAP.(hemisphere).(dataType).(arousalState).group);
%             Stats.(dataType).(arousalState).(freq).Table.Coherence = cat(1,data.Blank_SAP.(hemisphere).(dataType).(arousalState).(freq),data.SSP_SAP.(hemisphere).(dataType).(arousalState).(freq));
%             Stats.(dataType).(arousalState).(freq).FitFormula = 'Coherence ~ 1 + Group + (1|AnimalID)';
%             Stats.(dataType).(arousalState).(freq).Stats = fitglme(Stats.(dataType).(arousalState).(freq).Table,Stats.(dataType).(arousalState).(freq).FitFormula);
%         end
%     end
% end
% % statistical diary
% if saveFigs == true
%     % statistical diary
%     diaryFile = [dirpath 'BilateralCoherence_Ephys_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     for aa = 1:length(dataTypes)
%         dataType = dataTypes{1,aa};
%         for bb = 1:length(arousalStates)
%             arousalState = arousalStates{1,bb};
%             disp('======================================================================================================================')
%             disp(['GLME statistics: ' dataType ' ' arousalState ' Coherence^2'])
%             disp('======================================================================================================================')
%             disp('0 -> 0.01 Hz')
%             disp(Stats.(dataType).(arousalState).C001.Stats)
%             disp('----------------------------------------------------------------------------------------------------------------------')
%             disp('0.01 -> 0.1 Hz')
%             disp(Stats.(dataType).(arousalState).C01.Stats)
%             disp('----------------------------------------------------------------------------------------------------------------------')
%             disp('0.1 -> 0.5 Hz')
%             disp(Stats.(dataType).(arousalState).C05.Stats)
%         end
%     end
%     diary off
% end