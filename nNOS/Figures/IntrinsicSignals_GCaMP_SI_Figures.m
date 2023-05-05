function [] = IntrinsicSignals_GCaMP_SI_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_IntSig_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Rest','Whisk','Stim','NREM','REM'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
variables = {'avg','p2p','vari'};
fs = 30;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_IntSig_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                data.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(data.(group).(hemisphere).(dataType),behavior) == false
                        data.(group).(hemisphere).(dataType).(behavior).avg = [];
                        data.(group).(hemisphere).(dataType).(behavior).p2p = [];
                        data.(group).(hemisphere).(dataType).(behavior).vari = [];
                        data.(group).(hemisphere).(dataType).(behavior).group = {};
                        data.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    animalVar = [];
                    animalP2P = [];
                    for ff = 1:length(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(behavior).indHbT)
                        if strcmp(dataType,'Rest') == true
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(behavior).indHbT{ff,1}(2*fs:end);
                        elseif strcmp(dataType,'Stim') == true
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(behavior).indHbT{ff,1} - mean(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(behavior).HbT);
                        else
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(behavior).indHbT{ff,1};
                        end
                        animalVar(ff,1) = var(dataArray);
                        animalP2P(ff,1) = max(dataArray) - min(dataArray);
                    end
                    data.(group).(hemisphere).(dataType).(behavior).avg = cat(1,data.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(behavior).HbT));
                    data.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,data.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                    data.(group).(hemisphere).(dataType).(behavior).vari = cat(1,data.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                    data.(group).(hemisphere).(dataType).(behavior).group = cat(1,data.(group).(hemisphere).(dataType).(behavior).group,group);
                    data.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,data.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    data.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(data.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    data.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(data.(group).(hemisphere).(dataType).(behavior).(variable),1);
                end
            end
        end
    end
end
% figure
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(variables)
            variable = variables{1,cc};
            summaryFigure = figure;
            sgtitle(['IOS ' hemisphere ' ' dataType ' ' variable ' [GCaMP SI]'])
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                subplot(2,3,dd);
                xInds = ones(1,length(data.Blank_SAP.(hemisphere).(dataType).(behavior).(variable)));
                s1 = scatter(xInds*1,data.Blank_SAP.(hemisphere).(dataType).(behavior).(variable),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
                e1 = errorbar(1,data.Blank_SAP.(hemisphere).(dataType).(behavior).(['mean_' variable]),data.Blank_SAP.(hemisphere).(dataType).(behavior).(['std_' variable]),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
                e1.Color = 'black';
                e1.MarkerSize = 10;
                e1.CapSize = 10;
                xInds = ones(1,length(data.SSP_SAP.(hemisphere).(dataType).(behavior).(variable)));
                s2 = scatter(xInds*2,data.SSP_SAP.(hemisphere).(dataType).(behavior).(variable),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
                e2 = errorbar(2,data.SSP_SAP.(hemisphere).(dataType).(behavior).(['mean_' variable]),data.SSP_SAP.(hemisphere).(dataType).(behavior).(['std_' variable]),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
                e2.Color = 'black';
                e2.MarkerSize = 10;
                e2.CapSize = 10;
                title(behavior)
                if strcmp(dataType,'GCaMP') == true
                    ylabel('\DeltaF/F (%)')
                else
                    ylabel('\Delta[Hb] (\muM)')
                end
                xlim([0,3])
                if dd == 1
                    legend([s1,s2],'Blank-SAP','SSP-SAP')
                end
                set(gca,'box','off')
                set(gca,'xtick',[])
                axis square
                % save figure(s)
                if saveFigs == true
                    dirpath = [rootFolder delim 'Summary Figures' delim 'Intrinsic Signals' delim];
                    if ~exist(dirpath,'dir')
                        mkdir(dirpath);
                    end
                    savefig(summaryFigure,[dirpath 'IntrinsicSignals_GCaMP_SI_ ' hemisphere '_' dataType '_' variable]);
                end
            end
        end
    end
end
% statistics - generalized linear mixed effects model
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                % statistics - unpaired ttest
                [stats.(hemisphere).(dataType).(behavior).(variable).h,...
                    stats.(hemisphere).(dataType).(behavior).(variable).p,...
                    stats.(hemisphere).(dataType).(behavior).(variable).ci,...
                    stats.(hemisphere).(dataType).(behavior).(variable).stats] ...
                    = ttest2(data.Blank_SAP.(hemisphere).(dataType).(behavior).(variable),data.SSP_SAP.(hemisphere).(dataType).(behavior).(variable));
            end
        end
    end
end
% statistical diary
if saveFigs == true
    % statistical diary
    diaryFile = [dirpath 'IntrinsicSignals_GCaMP_SI_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    disp('======================================================================================================================')
    disp('ttest2 statistics: Blank vs. SAP')
    disp('======================================================================================================================')
    % statistics - generalized linear mixed effects model
    for aa = 1:length(hemispheres)
        hemisphere = hemispheres{1,aa};
        for bb = 1:length(dataTypes)
            dataType = dataTypes{1,bb};
            for cc = 1:length(behaviors)
                behavior = behaviors{1,cc};
                for dd = 1:length(variables)
                    variable = variables{1,dd};
                    disp([hemisphere ' ' dataType ' ' behavior ' ' variable ' p < ' num2str(stats.(hemisphere).(dataType).(behavior).(variable).p)]); disp(' ')
                    disp('----------------------------------------------------------------------------------------------------------------------')
                end
            end
        end
    end
    diary off
end