function [] = IntrinsicSignals_GCaMP_Figures(rootFolder,saveFigs,delim)
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
groups = {'SSP_SAP','Blank_SAP'};
hemispheres = {'LH','RH','fLH','fRH'};
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
                    for zz = 1:length(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(behavior).indHbT)
                        if strcmp(dataType,'Rest') == true
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(behavior).indHbT{zz,1}(2*fs:end);
                        else
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(behavior).indHbT{zz,1};
                        end
                        animalVar(zz,1) = var(dataArray);
                        animalP2P(zz,1) = max(dataArray) - min(dataArray);
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
        for dd = 1:length(dataTypes)
            dataType = dataTypes{1,dd};
            for cc = 1:length(behaviors)
                behavior = behaviors{1,cc};
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
for qq = 1:length(dataTypes)
    dataType = dataTypes{1,qq};
    for zz = 1:length(hemispheres)
        hemisphere = hemispheres{1,zz};
        for aa = 1:length(variables)
            variable = variables{1,aa};
            summaryFigure = figure;
            sgtitle(['IOS ' hemisphere ' ' dataType ' ' variable ' [GCaMP]'])
            for bb = 1:length(behaviors)
                behavior = behaviors{1,bb};
                subplot(2,3,bb);
                xInds = ones(1,length(data.Blank_SAP.(hemisphere).(dataType).(behavior).(variable)));
                s1 = scatter(xInds*1,data.Blank_SAP.(hemisphere).(dataType).(behavior).(variable),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
                hold on
                e1 = errorbar(1,data.Blank_SAP.(hemisphere).(dataType).(behavior).(['mean_' variable]),data.Blank_SAP.(hemisphere).(dataType).(behavior).(['std_' variable]),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
                e1.Color = 'black';
                e1.MarkerSize = 10;
                e1.CapSize = 10;
                xInds = ones(1,length(data.SSP_SAP.(hemisphere).(dataType).(behavior).(variable)));
                s2 = scatter(xInds*2,data.SSP_SAP.(hemisphere).(dataType).(behavior).(variable),75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
                hold on
                e2 = errorbar(2,data.SSP_SAP.(hemisphere).(dataType).(behavior).(['mean_' variable]),data.SSP_SAP.(hemisphere).(dataType).(behavior).(['std_' variable]),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
                e2.Color = 'black';
                e2.MarkerSize = 10;
                e2.CapSize = 10;
                title(behavior)
                if strcmp(dataType,'GCaMP') == true
                    ylabel('\DeltaF/F (%)')
                else
                ylabel('\Delta[HbT] (\muM)')
                end
                xlim([0,4])
                if bb == 1
                    legend([s1,s2],'Blank-SAP','SSP-SAP')
                end
                set(gca,'box','off')
                set(gca,'xtick',[])
                axis square
                % save figure(s)
                if saveFigs == true
                    dirpath = [rootFolder delim 'Summary Figures' delim 'Intrinsic Signals GCaMP' delim];
                    if ~exist(dirpath,'dir')
                        mkdir(dirpath);
                    end
                    savefig(summaryFigure,[dirpath 'IOS_ ' hemisphere '_' dataType '_' variable '_GCaMP']);
                end
            end
        end
    end
end
% statistics - generalized linear mixed effects model
for bb = 1:length(hemispheres)
    hemisphere = hemispheres{1,bb};
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        for dd = 1:length(behaviors)
            behavior = behaviors{1,dd};
            for cc = 1:length(variables)
                variable = variables{1,cc};
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
    diaryFile = [dirpath 'WhiskingBehavior_GCaMP_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    disp('======================================================================================================================')
    disp('ttest2 statistics: Blank vs. SAP')
    disp('======================================================================================================================')
    % statistics - generalized linear mixed effects model
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for aa = 1:length(dataTypes)
            dataType = dataTypes{1,aa};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                for cc = 1:length(variables)
                    variable = variables{1,cc};
                    disp([hemisphere ' ' dataType ' ' behavior ' ' variable ' p < ' num2str(stats.(hemisphere).(dataType).(behavior).(variable).p)]); disp(' ')
                    disp('----------------------------------------------------------------------------------------------------------------------')
                end
            end
        end
    end
    diary off
end