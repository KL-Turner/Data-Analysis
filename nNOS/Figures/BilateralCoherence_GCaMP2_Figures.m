function [] = BilateralCoherence_GCaMP2_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_BilatCoher_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
variables = {'C','f'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_BilatCoher_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            data.(group).(dataType).dummCheck = 1;
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                if isfield(data.(group).(dataType),behavior) == false
                    data.(group).(dataType).(behavior).C = [];
                    data.(group).(dataType).(behavior).f = [];
                    data.(group).(dataType).(behavior).group = {};
                    data.(group).(dataType).(behavior).animalID = {};
                end
                if isempty(Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).C) == false
                    data.(group).(dataType).(behavior).C = cat(1,data.(group).(dataType).(behavior).C,Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).fC.^2');
                    data.(group).(dataType).(behavior).f = cat(1,data.(group).(dataType).(behavior).f,Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).f);
                    data.(group).(dataType).(behavior).group = cat(1,data.(group).(dataType).(behavior).group,group);
                    data.(group).(dataType).(behavior).animalID = cat(1,data.(group).(dataType).(behavior).animalID,animalID);
                end
            end
        end
    end
end
% mean/stdanimalID
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                data.(group).(dataType).(behavior).(['mean_' variable]) = mean(data.(group).(dataType).(behavior).(variable),1);
                data.(group).(dataType).(behavior).(['stdErr_' variable]) = std(data.(group).(dataType).(behavior).(variable),1)./sqrt(size(data.(group).(dataType).(behavior).(variable),1));
            end
        end
    end
end
% coherence figures
xlimits = ({[1/10,0.5],[1/30,0.5],[1/60,0.5],[0.003,0.5],[0.003,0.5],[0.003,0.5]});
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    summaryFigure = figure;
    sgtitle([dataType ' bilateral coherence [GCaMP]'])
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        subplot(2,3,bb);
        p1 = semilogx(data.Blank_SAP.(dataType).(behavior).mean_f,data.Blank_SAP.(dataType).(behavior).mean_C,'color',colors('north texas green'),'LineWidth',2);
        hold on
        semilogx(data.Blank_SAP.(dataType).(behavior).mean_f,data.Blank_SAP.(dataType).(behavior).mean_C + data.Blank_SAP.(dataType).(behavior).stdErr_C,'color',colors('north texas green'),'LineWidth',0.5);
        semilogx(data.Blank_SAP.(dataType).(behavior).mean_f,data.Blank_SAP.(dataType).(behavior).mean_C - data.Blank_SAP.(dataType).(behavior).stdErr_C,'color',colors('north texas green'),'LineWidth',0.5);
        p2 = semilogx(data.SSP_SAP.(dataType).(behavior).mean_f,data.SSP_SAP.(dataType).(behavior).mean_C,'color',colors('electric purple'),'LineWidth',2);
        semilogx(data.SSP_SAP.(dataType).(behavior).mean_f,data.SSP_SAP.(dataType).(behavior).mean_C + data.SSP_SAP.(dataType).(behavior).stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
        semilogx(data.SSP_SAP.(dataType).(behavior).mean_f,data.SSP_SAP.(dataType).(behavior).mean_C - data.SSP_SAP.(dataType).(behavior).stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
        ylabel('Coherence^2')
        xlabel('Freq (Hz)')
        title(behavior)
        xlim(xlimits{1,bb})
        ylim([0,1])
        if bb == 1
            legend([p1,p2],'Blank-SAP','SSP-SAP')
        end
        set(gca,'box','off')
        axis square
        % save figure(s)
        if saveFigs == true
            dirpath = [rootFolder delim 'Summary Figures' delim 'Bilateral Coherence GCaMP' delim];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(summaryFigure,[dirpath 'BilateralCoherence_GCaMP_' dataType]);
        end
    end
end
% find Hz peaks in coherence
arousalStates = {'Alert','Asleep','All'};
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(arousalStates)
            arousalState = arousalStates{1,cc};
            for dd = 1:size(data.(group).(dataType).(arousalState).C,1)
                F = round(data.(group).(dataType).(arousalState).f(dd,:),3);
                C = data.(group).(dataType).(arousalState).C(dd,:);
                index001 = find(F == 0.01);
                index01 = find(F == 0.1);
                index05 = find(F == 0.5);
                data.(group).(dataType).(arousalState).C001(dd,1) = mean(C(1:index001(1)));
                data.(group).(dataType).(arousalState).C01(dd,1) = mean(C(index001(1) + 1:index01(1)));
                data.(group).(dataType).(arousalState).C05(dd,1) = mean(C(index01(1) + 1:index05(1)));
            end
        end
    end
end
% peak coherence figures
freqs = {'C001','C01','C05'};
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    figure;
    sgtitle([dataType ' peak coherence [GCaMP]'])
    for bb = 1:length(arousalStates)
        arousalState = arousalStates{1,bb};
        subplot(1,3,bb); zz = 1;
        for cc = 1:length(freqs)
            freq = freqs{1,cc};
            xInds = ones(1,length(data.Blank_SAP.(dataType).(arousalState).(freq)));
            s1 = scatter(xInds*zz,data.Blank_SAP.(dataType).(arousalState).(freq),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
            hold on
            e1 = errorbar(zz,mean(data.Blank_SAP.(dataType).(arousalState).(freq)),std(data.Blank_SAP.(dataType).(arousalState).(freq)),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
            e1.Color = 'black';
            e1.MarkerSize = 10;
            e1.CapSize = 10;
            zz = zz + 1;
            xInds = ones(1,length(data.SSP_SAP.(dataType).(arousalState).(freq)));
            s2 = scatter(xInds*zz,data.SSP_SAP.(dataType).(arousalState).(freq),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
            hold on
            e2 = errorbar(zz,mean(data.SSP_SAP.(dataType).(arousalState).(freq)),std(data.SSP_SAP.(dataType).(arousalState).(freq)),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
            e2.Color = 'black';
            e2.MarkerSize = 10;
            e2.CapSize = 10;
            zz = zz + 1;
            ylabel('Coherence^2')
            xticks([.5,3.5,5.5])
            xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
            title(arousalState)
            xlim([0,7])
            ylim([0,1])
        end
        if bb == 1
            legend([s1,s2],'Blank-SAP','SSP-SAP')
        end
        set(gca,'box','off')
        axis square
    end
    % save figure(s)
    if saveFigs == true
        savefig(summaryFigure,[dirpath 'BilateralCoherence_GCaMP_' dataType '_FreqPeaks']);
    end
end
% statistics - generalized linear mixed effects model
for aa = 1:length(freqs)
    freq = freqs{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(arousalStates)
            arousalState = arousalStates{1,cc};
            % statistics - generalized linear mixed effects model
            Stats.(dataType).(arousalState).(freq).tableSize = cat(1,data.Blank_SAP.(dataType).(arousalState).(freq),data.SSP_SAP.(dataType).(arousalState).(freq));
            Stats.(dataType).(arousalState).(freq).Table = table('Size',[size(Stats.(dataType).(arousalState).(freq).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'AnimalID','Group','Coherence'});
            Stats.(dataType).(arousalState).(freq).Table.AnimalID = cat(1,data.Blank_SAP.(dataType).(arousalState).animalID,data.SSP_SAP.(dataType).(arousalState).animalID);
            Stats.(dataType).(arousalState).(freq).Table.Group = cat(1,data.Blank_SAP.(dataType).(arousalState).group,data.SSP_SAP.(dataType).(arousalState).group);
            Stats.(dataType).(arousalState).(freq).Table.Coherence = cat(1,data.Blank_SAP.(dataType).(arousalState).(freq),data.SSP_SAP.(dataType).(arousalState).(freq));
            Stats.(dataType).(arousalState).(freq).FitFormula = 'Coherence ~ 1 + Group + (1|AnimalID)';
            Stats.(dataType).(arousalState).(freq).Stats = fitglme(Stats.(dataType).(arousalState).(freq).Table,Stats.(dataType).(arousalState).(freq).FitFormula);
        end
    end
end
% statistical diary
if saveFigs == true
    % statistical diary
    diaryFile = [dirpath 'BilateralCoherence_GCaMP_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        for bb = 1:length(arousalStates)
            arousalState = arousalStates{1,bb};
            disp('======================================================================================================================')
            disp(['GLME statistics: ' dataType ' ' arousalState ' Coherence^2'])
            disp('======================================================================================================================')
            disp('0 -> 0.01 Hz')
            disp(Stats.(dataType).(arousalState).C001.Stats)
            disp('----------------------------------------------------------------------------------------------------------------------')
            disp('0.01 -> 0.1 Hz')
            disp(Stats.(dataType).(arousalState).C01.Stats)
            disp('----------------------------------------------------------------------------------------------------------------------')
            disp('0.1 -> 0.5 Hz')
            disp(Stats.(dataType).(arousalState).C05.Stats)
        end
    end
    diary off
end