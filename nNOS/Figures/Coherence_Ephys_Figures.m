function [] = Coherence_Ephys_Figures(rootFolder,~,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_BilatCoher_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','SSP_SAP','Blank_SAP'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
dataTypes = {'HbT','gammaBandPower','deltaBandPower'};
variables = {'C','f'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_BilatCoher_Ephys.(group));
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
                data.(group).(behavior).(dataType).C = cat(1,data.(group).(behavior).(dataType).C,Results_BilatCoher_Ephys.(group).(animalID).(behavior).(dataType).C.^2');
                data.(group).(behavior).(dataType).f = cat(1,data.(group).(behavior).(dataType).f,Results_BilatCoher_Ephys.(group).(animalID).(behavior).(dataType).f);
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
    sgtitle([dataType ' bilateral coherence [Ephys]'])
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
% find Hz peaks in coherence
arousalStates = {'Alert','Asleep','All'};
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(arousalStates)
            arousalState = arousalStates{1,cc};
            for dd = 1:size(data.(group).(arousalState).(dataType).C,2)
                F = round(data.(group).(dataType).(arousalState).f(dd,:),3);
                C = data.(group).(dataType).(arousalState).C(:,dd);
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
% figure
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    figure;
    sgtitle([dataType ' peak coherence [Ephys]'])
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
% xInds = ones(1,length(animalIDs.Blank_SAP));
% s1 = scatter(xInds*1,data.Blank_SAP.CBV_HbT.Alert.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e1 = errorbar(1,mean(data.Blank_SAP.CBV_HbT.Alert.C001),std(data.Blank_SAP.CBV_HbT.Alert.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'black';
% e1.MarkerSize = 10;
% e1.CapSize = 10;
% s2 = scatter(xInds*2,data.SSP_SAP.CBV_HbT.Alert.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e2 = errorbar(2,mean(data.SSP_SAP.CBV_HbT.Alert.C001),std(data.Blank_SAP.CBV_HbT.Alert.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e2.Color = 'black';
% e2.MarkerSize = 10;
% e2.CapSize = 10;
% scatter(xInds*3,data.Blank_SAP.CBV_HbT.Alert.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e3 = errorbar(3,mean(data.Blank_SAP.CBV_HbT.Alert.C01),std(data.Blank_SAP.CBV_HbT.Alert.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e3.Color = 'black';
% e3.MarkerSize = 10;
% e3.CapSize = 10;
% scatter(xInds*4,data.SSP_SAP.CBV_HbT.Alert.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e4 = errorbar(4,mean(data.SSP_SAP.CBV_HbT.Alert.C01),std(data.Blank_SAP.CBV_HbT.Alert.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'black';
% e4.MarkerSize = 10;
% e4.CapSize = 10;
% scatter(xInds*5,data.Blank_SAP.CBV_HbT.Alert.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e5 = errorbar(5,mean(data.Blank_SAP.CBV_HbT.Alert.C05),std(data.Blank_SAP.CBV_HbT.Alert.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e5.Color = 'black';
% e5.MarkerSize = 10;
% e5.CapSize = 10;
% scatter(xInds*6,data.SSP_SAP.CBV_HbT.Alert.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e6 = errorbar(6,mean(data.SSP_SAP.CBV_HbT.Alert.C05),std(data.Blank_SAP.CBV_HbT.Alert.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e6.Color = 'black';
% e6.MarkerSize = 10;
% e6.CapSize = 10;
% statistics - generalized linear mixed effects model
freqs = {'C001','C01','C05'};
for aa = 1:length(freqs)
    freq = freqs{1,aa};
    for bb = 1:length(behavFields2)
        behavField = behavFields2{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % statistics - generalized linear mixed effects model
            Stats.(dataType).(behavField).(freq).tableSize = cat(1,data.Blank_SAP.(behavField).(dataType).(freq),data.SSP_SAP.(behavField).(dataType).(freq));
            Stats.(dataType).(behavField).(freq).Table = table('Size',[size(Stats.(dataType).(behavField).(freq).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Treatment','Coherence'});
            Stats.(dataType).(behavField).(freq).Table.Mouse = cat(1,data.Blank_SAP.(behavField).(dataType).animalID,data.SSP_SAP.(behavField).(dataType).animalID);
            Stats.(dataType).(behavField).(freq).Table.Treatment = cat(1,data.Blank_SAP.(behavField).(dataType).treatment,data.SSP_SAP.(behavField).(dataType).treatment);
            Stats.(dataType).(behavField).(freq).Table.Coherence = cat(1,data.Blank_SAP.(behavField).(dataType).(freq),data.SSP_SAP.(behavField).(dataType).(freq));
            Stats.(dataType).(behavField).(freq).FitFormula = 'Coherence ~ 1 + Treatment + (1|Mouse)';
            Stats.(dataType).(behavField).(freq).Stats = fitglme(Stats.(dataType).(behavField).(freq).Table,Stats.(dataType).(behavField).(freq).FitFormula);
        end
    end
end

%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Bilateral Coherence - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure5,[dirpath 'AverageBilateralCoherence_Statistics']);
    set(summaryFigure5,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'AverageBilateralCoherence_Statistics'])
    %% statistical diary
    diaryFile = [dirpath 'AverageBilateralCoherence_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % Awake stats
    disp('======================================================================================================================')
    disp('GLME statistics for gamma Coherence^2 for Awake data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.gammaBandPower.Awake.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.gammaBandPower.Awake.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.gammaBandPower.Awake.C05.Stats)
    % Sleep stats
    disp('======================================================================================================================')
    disp('GLME statistics for gamma Coherence^2 for Sleep data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.gammaBandPower.Sleep.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.gammaBandPower.Sleep.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.gammaBandPower.Sleep.C05.Stats)
    % All stats
    disp('======================================================================================================================')
    disp('GLME statistics for gamma Coherence^2 for All data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.gammaBandPower.All.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.gammaBandPower.All.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.gammaBandPower.All.C05.Stats)
    % Awake stats
    disp('======================================================================================================================')
    disp('GLME statistics for CBV_HbT Coherence^2 for Awake data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.CBV_HbT.Awake.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.CBV_HbT.Awake.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.CBV_HbT.Awake.C05.Stats)
    % Sleep stats
    disp('======================================================================================================================')
    disp('GLME statistics for CBV_HbT Coherence^2 for Sleep data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.CBV_HbT.Sleep.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.CBV_HbT.Sleep.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.CBV_HbT.Sleep.C05.Stats)
    % All stats
    disp('======================================================================================================================')
    disp('GLME statistics for CBV_HbT Coherence^2 for All data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.CBV_HbT.All.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.CBV_HbT.All.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.CBV_HbT.All.C05.Stats)
    diary off
end

end
