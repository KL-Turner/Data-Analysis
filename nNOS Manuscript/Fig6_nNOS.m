function [] = Fig6_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PowerSpec_LFP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Alert','Asleep','All'};
variables = {'S','f','deltaS'};
dimensions = [2,1,1];
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PowerSpec_LFP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                data.(group).(hemisphere).(behavior).dummCheck = 1;
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    dimension = dimensions(ee);
                    if isfield(data.(group).(hemisphere).(behavior),(variable)) == false
                        data.(group).(hemisphere).(behavior).(variable) = [];
                        data.(group).(hemisphere).(behavior).group = {};
                        data.(group).(hemisphere).(behavior).animalID = {};
                        data.(group).(hemisphere).(behavior).hemisphere = {};
                        data.(group).(hemisphere).(behavior).behavior = {};
                    end
                    % pull data if field isn't empty
                    if isempty(Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).S) == false
                        if strcmp(variable,'deltaS') == true
                            index = find(round(Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).f,2) == 4);
                            deltaIndex = index(end);
                            data.(group).(hemisphere).(behavior).(variable) = cat(dimension,data.(group).(hemisphere).(behavior).(variable),mean(Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).S(1:deltaIndex)));
                            % for stats
                            data.(group).(hemisphere).(behavior).group = cat(1,data.(group).(hemisphere).(behavior).group,group);
                            data.(group).(hemisphere).(behavior).animalID = cat(1,data.(group).(hemisphere).(behavior).animalID,animalID);
                            data.(group).(hemisphere).(behavior).hemisphere = cat(1,data.(group).(hemisphere).(behavior).hemisphere,hemisphere);
                            data.(group).(hemisphere).(behavior).behavior = cat(1,data.(group).(hemisphere).(behavior).behavior,behavior);
                        else
                            data.(group).(hemisphere).(behavior).(variable) = cat(dimension,data.(group).(hemisphere).(behavior).(variable),Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).(variable));
                        end
                    end
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                dimension = dimensions(dd);
                data.(group).(hemisphere).(behavior).(['mean_' variable]) = mean(data.(group).(hemisphere).(behavior).(variable),dimension);
                data.(group).(hemisphere).(behavior).(['stdErr_' variable]) = std(data.(group).(hemisphere).(behavior).(variable),0,dimension)./sqrt(size(data.(group).(hemisphere).(behavior).(variable),dimension));
            end
        end
    end
end
% figure
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    summaryFigure = figure;
    sgtitle([hemisphere ' LFP power spectra [1-100 Hz]'])
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        subplot(1,3,bb);
        L1 = loglog(data.Naive.(hemisphere).(behavior).mean_f,data.Naive.(hemisphere).(behavior).mean_S,'color',colors('sapphire'),'LineWidth',2);
        hold on;
        loglog(data.Naive.(hemisphere).(behavior).mean_f,data.Naive.(hemisphere).(behavior).mean_S + data.Naive.(hemisphere).(behavior).stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
        loglog(data.Naive.(hemisphere).(behavior).mean_f,data.Naive.(hemisphere).(behavior).mean_S - data.Naive.(hemisphere).(behavior).stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
        L2 = loglog(data.Blank_SAP.(hemisphere).(behavior).mean_f,data.Blank_SAP.(hemisphere).(behavior).mean_S,'color',colors('north texas green'),'LineWidth',2);
        loglog(data.Blank_SAP.(hemisphere).(behavior).mean_f,data.Blank_SAP.(hemisphere).(behavior).mean_S + data.Blank_SAP.(hemisphere).(behavior).stdErr_S,'color',colors('north texas green'),'LineWidth',0.25);
        loglog(data.Blank_SAP.(hemisphere).(behavior).mean_f,data.Blank_SAP.(hemisphere).(behavior).mean_S - data.Blank_SAP.(hemisphere).(behavior).stdErr_S,'color',colors('north texas green'),'LineWidth',0.25);
        L3 = loglog(data.SSP_SAP.(hemisphere).(behavior).mean_f,data.SSP_SAP.(hemisphere).(behavior).mean_S,'color',colors('electric purple'),'LineWidth',2);
        loglog(data.SSP_SAP.(hemisphere).(behavior).mean_f,data.SSP_SAP.(hemisphere).(behavior).mean_S + data.SSP_SAP.(hemisphere).(behavior).stdErr_S,'color',colors('electric purple'),'LineWidth',0.25);
        loglog(data.SSP_SAP.(hemisphere).(behavior).mean_f,data.SSP_SAP.(hemisphere).(behavior).mean_S - data.SSP_SAP.(hemisphere).(behavior).stdErr_S,'color',colors('electric purple'),'LineWidth',0.25);
        title(behavior)
        ylabel('Power (a.u.)')
        xlabel('Freq (Hz)')
        xlim([1,100])
        if bb == 1
            legend([L1,L2,L3],'Naive','Blank-SAP','SSP-SAP')
        end
        set(gca,'box','off')
        axis square
    end
    linkaxes
    % save figure(s)
    if saveFigs == true
        dirpath = [rootFolder delim 'Summary Figures' delim 'Power Spectrum' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigure,[dirpath 'PowerSpectrumLFP_Ephys_' hemisphere]);
        set(summaryFigure,'PaperPositionMode','auto');
        print('-vector','-dpdf','-fillpage',[dirpath 'PowerSpectrumLFP_Ephys_' hemisphere])
    end
end
% figure
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    summaryFigure = figure;
    sgtitle([hemisphere ' LFP delta power'])
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        subplot(1,3,bb);
        s1 = scatter(ones(1,length(data.Naive.(hemisphere).(behavior).deltaS))*1,data.Naive.(hemisphere).(behavior).deltaS,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off','jitterAmount',0.25);
        hold on;
        e1 = errorbar(1,mean(data.Naive.(hemisphere).(behavior).deltaS),std(data.Naive.(hemisphere).(behavior).deltaS,0,1),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e1.Color = 'black';
        e1.MarkerSize = 10;
        e1.CapSize = 10;
        s2 = scatter(ones(1,length(data.Blank_SAP.(hemisphere).(behavior).deltaS))*2,data.Blank_SAP.(hemisphere).(behavior).deltaS,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
        e2 = errorbar(2,mean(data.Blank_SAP.(hemisphere).(behavior).deltaS),std(data.Blank_SAP.(hemisphere).(behavior).deltaS,0,1),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e2.Color = 'black';
        e2.MarkerSize = 10;
        e2.CapSize = 10;
        s3 = scatter(ones(1,length(data.SSP_SAP.(hemisphere).(behavior).deltaS))*3,data.SSP_SAP.(hemisphere).(behavior).deltaS,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
        e3 = errorbar(3,mean(data.SSP_SAP.(hemisphere).(behavior).deltaS),std(data.SSP_SAP.(hemisphere).(behavior).deltaS,0,1),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
        e3.Color = 'black';
        e3.MarkerSize = 10;
        e3.CapSize = 10;
        title(behavior)
        ylabel('Power (a.u.)')
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        xlim([0,4])
        if bb == 1
            legend([s1,s2,s3],'Naive','Blank-SAP','SSP-SAP')
        end
        axis square
        % set(gca,'YScale','log')
    end
    % save figure(s)
    if saveFigs == true
        savefig(summaryFigure,[dirpath 'PowerSpectrumLFP_Ephys_' hemisphere '_deltaPower']);
        set(summaryFigure,'PaperPositionMode','auto');
        print('-vector','-dpdf','-fillpage',[dirpath 'PowerSpectrumLFP_Ephys_' hemisphere '_deltaPower'])
    end
end
% statistics - generalized linear mixed effects model
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        stats.(hemisphere).(behavior).tableSize = cat(1,data.Blank_SAP.(hemisphere).(behavior).deltaS,data.SSP_SAP.(hemisphere).(behavior).deltaS,data.Naive.(hemisphere).(behavior).deltaS);
        stats.(hemisphere).(behavior).Table = table('Size',[size(stats.(hemisphere).(behavior).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'group','animalID','behavior','deltaS'});
        stats.(hemisphere).(behavior).Table.group = cat(1,data.Blank_SAP.(hemisphere).(behavior).group,data.SSP_SAP.(hemisphere).(behavior).group,data.Naive.(hemisphere).(behavior).group);
        stats.(hemisphere).(behavior).Table.animalID = cat(1,data.Blank_SAP.(hemisphere).(behavior).animalID,data.SSP_SAP.(hemisphere).(behavior).animalID,data.Naive.(hemisphere).(behavior).animalID);
        stats.(hemisphere).(behavior).Table.behavior = cat(1,data.Blank_SAP.(hemisphere).(behavior).behavior,data.SSP_SAP.(hemisphere).(behavior).behavior,data.Naive.(hemisphere).(behavior).behavior);
        stats.(hemisphere).(behavior).Table.deltaS = cat(1,data.Blank_SAP.(hemisphere).(behavior).deltaS,data.SSP_SAP.(hemisphere).(behavior).deltaS,data.Naive.(hemisphere).(behavior).deltaS);
        stats.(hemisphere).(behavior).FitFormula = 'deltaS ~ 1 + group + behavior + (1|animalID)';
        stats.(hemisphere).(behavior).Stats = fitglme(stats.(hemisphere).(behavior).Table,stats.(hemisphere).(behavior).FitFormula);
    end
end
% statistical diary
if saveFigs == true
    % statistical diary
    diaryFile = [dirpath 'PowerSpectrumLFP_Ephys_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    for aa = 1:length(hemispheres)
        hemisphere = hemispheres{1,aa};
        for bb = 1:length(behaviors)
            behavior = behaviors{1,bb};
            disp('======================================================================================================================')
            disp(['GLME statistics: ' hemisphere ' ' behavior])
            disp('======================================================================================================================')
            disp(stats.(hemisphere).(behavior).Stats)
            disp('----------------------------------------------------------------------------------------------------------------------')
        end
    end
    diary off
end