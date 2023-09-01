function [] = FigS3_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

%% arousal-state hemodnyamics [Ephys]
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_IntSig_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Rest','Whisk','Stim','NREM','REM','Iso'};
dataTypes = {'HbT'};
variables = {'avg','p2p','vari'};
fs = 30;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_IntSig_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                ephysData.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(ephysData.(group).(hemisphere).(dataType),behavior) == false
                        ephysData.(group).(hemisphere).(dataType).(behavior).avg = [];
                        ephysData.(group).(hemisphere).(dataType).(behavior).p2p = [];
                        ephysData.(group).(hemisphere).(dataType).(behavior).vari = [];
                        ephysData.(group).(hemisphere).(dataType).(behavior).group = {};
                        ephysData.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    animalVar = [];
                    animalP2P = [];
                    for ff = 1:length(Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT)
                        if strcmp(behavior,'Rest') == true
                            dataArray = Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT{ff,1}(2*fs:end);
                        elseif strcmp(behavior,'Stim') == true
                            dataArray = Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT{ff,1} - mean(cell2mat(Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT),1);
                        else
                            dataArray = Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT{ff,1};
                        end
                        animalVar(ff,1) = var(dataArray);
                        animalP2P(ff,1) = max(dataArray) - min(dataArray);
                    end
                    ephysData.(group).(hemisphere).(dataType).(behavior).avg = cat(1,ephysData.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).HbT));
                    ephysData.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,ephysData.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                    ephysData.(group).(hemisphere).(dataType).(behavior).vari = cat(1,ephysData.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                    ephysData.(group).(hemisphere).(dataType).(behavior).group = cat(1,ephysData.(group).(hemisphere).(dataType).(behavior).group,group);
                    ephysData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,ephysData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
                    ephysData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(ephysData.(group).(hemisphere).(dataType).(behavior).(variable),1,'omitnan');
                    ephysData.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(ephysData.(group).(hemisphere).(dataType).(behavior).(variable),1,'omitnan');
                end
            end
        end
    end
end
% statistics - ttest
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                % statistics - unpaired ttest
                [EphysStats.(hemisphere).(dataType).(behavior).(variable).h,EphysStats.(hemisphere).(dataType).(behavior).(variable).p] = ttest2(ephysData.Blank_SAP.(hemisphere).(dataType).(behavior).(variable),ephysData.SSP_SAP.(hemisphere).(dataType).(behavior).(variable));
            end
        end
    end
end

%% arousal-state hemodnyamics [GCaMP]
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
fs = 10;
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
                gcampData.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(gcampData.(group).(hemisphere).(dataType),behavior) == false
                        gcampData.(group).(hemisphere).(dataType).(behavior).avg = [];
                        gcampData.(group).(hemisphere).(dataType).(behavior).p2p = [];
                        gcampData.(group).(hemisphere).(dataType).(behavior).vari = [];
                        gcampData.(group).(hemisphere).(dataType).(behavior).group = {};
                        gcampData.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    animalVar = [];
                    animalP2P = [];
                    for ff = 1:length(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData)
                        if strcmp(behavior,'Rest') == true
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData{ff,1}(2*fs:end);
                        elseif strcmp(behavior,'Stim') == true
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData{ff,1} - mean(cell2mat(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData),1);
                        else
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData{ff,1};
                        end
                        animalVar(ff,1) = var(dataArray);
                        animalP2P(ff,1) = max(dataArray) - min(dataArray);
                    end
                    try
                        gcampData.(group).(hemisphere).(dataType).(behavior).avg = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).mean));
                    catch
                        gcampData.(group).(hemisphere).(dataType).(behavior).avg = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).avg,mean(cell2mat(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).mean)));
                    end
                    gcampData.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                    gcampData.(group).(hemisphere).(dataType).(behavior).vari = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                    gcampData.(group).(hemisphere).(dataType).(behavior).group = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).group,group);
                    gcampData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
                    gcampData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(gcampData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    gcampData.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(gcampData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                end
            end
        end
    end
end
% statistics - ttest
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                % statistics - unpaired ttest
                [GCaMPStats.(hemisphere).(dataType).(behavior).(variable).h,GCaMPStats.(hemisphere).(dataType).(behavior).(variable).p] = ttest2(gcampData.Blank_SAP.(hemisphere).(dataType).(behavior).(variable),gcampData.SSP_SAP.(hemisphere).(dataType).(behavior).(variable));
            end
        end
    end
end

%% two photo diameter
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Diameter_2P';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
behaviors = {'Rest','Whisk'};
variables = {'data'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Diameter_2P.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        vIDs = fieldnames(Results_Diameter_2P.(group).(animalID));
        for cc = 1:length(vIDs)
            vID = vIDs{cc,1};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                diameterData.(group).(behavior).dummCheck = 1;
                for ee = 1:length(variables)
                    if isfield(diameterData.(group).(behavior),(variables{1,ee})) == false
                        diameterData.(group).(behavior).(variables{1,ee}) = [];
                    end
                end
                if isfield(Results_Diameter_2P.(group).(animalID).(vID),behavior) == true
                    diameterData.(group).(behavior).data = cat(1,diameterData.(group).(behavior).data,mean(Results_Diameter_2P.(group).(animalID).(vID).(behavior).mean));
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        diameterData.(group).(behavior).meanData = mean(diameterData.(group).(behavior).data,1);
        diameterData.(group).(behavior).stdData = std(diameterData.(group).(behavior).data,0,1);
    end
end
% statistics - ttest
for cc = 1:length(behaviors)
    behavior = behaviors{1,cc};
    [DiameterStats.(behavior).h,DiameterStats.(behavior).p] = ttest2(diameterData.Blank_SAP.(behavior).data,diameterData.SSP_SAP.(behavior).data);
end

%% two photon isoflurane shift
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Baseline_2P';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
variables = {'diameter','baseline'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Baseline_2P.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        vIDs = fieldnames(Results_Baseline_2P.(group).(animalID));
        for cc = 1:length(vIDs)
            vID = vIDs{cc,1};
            diameterShiftData.(group).dummCheck = 1;
            for dd = 1:length(variables)
                if isfield(diameterShiftData.(group),(variables{1,dd})) == false
                    diameterShiftData.(group).(variables{1,dd}) = [];
                end
            end
            diameterShiftData.(group).diameter = cat(1,diameterShiftData.(group).diameter,((Results_Baseline_2P.(group).(animalID).(vID).diameter - Results_Baseline_2P.(group).(animalID).(vID).baseline)/Results_Baseline_2P.(group).(animalID).(vID).baseline)*100);
            diameterShiftData.(group).baseline = cat(1,diameterShiftData.(group).baseline,Results_Baseline_2P.(group).(animalID).(vID).baseline);
        end
    end
end
% mean/std
for ee = 1:length(groups)
    group = groups{1,ee};
    diameterShiftData.(group).meanDiameter = mean(diameterShiftData.(group).diameter,1);
    diameterShiftData.(group).stdDiameter = std(diameterShiftData.(group).diameter,0,1);
    diameterShiftData.(group).meanBaseline = mean(diameterShiftData.(group).baseline,1);
end
[TwoPIsoStats.h,TwoPIsoStats.p] = ttest2(diameterShiftData.Blank_SAP.diameter,diameterShiftData.SSP_SAP.diameter);

%% figure
FigS3 = figure('Name','Figure S3','units','normalized','outerposition',[0 0 1 1]);

subplot(2,3,1);
hold on
xInds = ones(1,length(ephysData.Blank_SAP.RH.HbT.Rest.avg));
s2 = scatter(xInds*1,ephysData.Blank_SAP.RH.HbT.Rest.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(1,ephysData.Blank_SAP.RH.HbT.Rest.mean_avg,ephysData.Blank_SAP.RH.HbT.Rest.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysData.SSP_SAP.RH.HbT.Rest.avg));
s3 = scatter(xInds*2,ephysData.SSP_SAP.RH.HbT.Rest.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,ephysData.SSP_SAP.RH.HbT.Rest.mean_avg,ephysData.SSP_SAP.RH.HbT.Rest.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;

xInds = ones(1,length(ephysData.Blank_SAP.RH.HbT.NREM.avg));
scatter(xInds*4,ephysData.Blank_SAP.RH.HbT.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(4,ephysData.Blank_SAP.RH.HbT.NREM.mean_avg,ephysData.Blank_SAP.RH.HbT.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysData.SSP_SAP.RH.HbT.NREM.avg));
scatter(xInds*5,ephysData.SSP_SAP.RH.HbT.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,ephysData.SSP_SAP.RH.HbT.NREM.mean_avg,ephysData.SSP_SAP.RH.HbT.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;

xInds = ones(1,length(ephysData.Blank_SAP.RH.HbT.REM.avg));
scatter(xInds*7,ephysData.Blank_SAP.RH.HbT.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(7,ephysData.Blank_SAP.RH.HbT.REM.mean_avg,ephysData.Blank_SAP.RH.HbT.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysData.SSP_SAP.RH.HbT.REM.avg));
scatter(xInds*8,ephysData.SSP_SAP.RH.HbT.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(8,ephysData.SSP_SAP.RH.HbT.REM.mean_avg,ephysData.SSP_SAP.RH.HbT.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(ephysData.Blank_SAP.RH.HbT.Iso.avg));
scatter(xInds*10,ephysData.Blank_SAP.RH.HbT.Iso.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(10,ephysData.Blank_SAP.RH.HbT.Iso.mean_avg,ephysData.Blank_SAP.RH.HbT.Iso.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysData.SSP_SAP.RH.HbT.Iso.avg));
scatter(xInds*11,ephysData.SSP_SAP.RH.HbT.Iso.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(11,ephysData.SSP_SAP.RH.HbT.Iso.mean_avg,ephysData.SSP_SAP.RH.HbT.Iso.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('\Delta[HbT] (\muM)')
xlim([0,12])
legend([s2,s3],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(2,3,2)
xInds = ones(1,length(diameterData.Blank_SAP.Rest.data));
s1 = scatter(xInds*1,diameterData.Blank_SAP.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,diameterData.Blank_SAP.Rest.meanData,diameterData.Blank_SAP.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(diameterData.SSP_SAP.Rest.data));
s2 = scatter(xInds*2,diameterData.SSP_SAP.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,diameterData.SSP_SAP.Rest.meanData,diameterData.SSP_SAP.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
title(behavior)
ylabel('\DeltaD/D (%)')
xlim([0,3])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(2,3,3)
xInds = ones(1,length(diameterShiftData.Blank_SAP.diameter));
scatter(xInds*1,diameterShiftData.Blank_SAP.diameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,diameterShiftData.Blank_SAP.meanDiameter,diameterShiftData.Blank_SAP.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(diameterShiftData.SSP_SAP.diameter));
scatter(xInds*2,diameterShiftData.SSP_SAP.diameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,diameterShiftData.SSP_SAP.meanDiameter,diameterShiftData.SSP_SAP.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('\DeltaD/D (%)')
set(gca,'xtick',[1,1.2,2,2.2])
set(gca,'xticklabel',{['Avg Baseline: ' num2str(round(diameterShiftData.Blank_SAP.meanBaseline,1)) ' \muM'],['n = ' num2str(length(diameterShiftData.Blank_SAP.diameter)) ' arterioles'],['Avg Baseline: ' num2str(round(diameterShiftData.SSP_SAP.meanBaseline,1)) ' \muM'],['n = ' num2str(length(diameterShiftData.SSP_SAP.diameter)) ' arterioles']})
xtickangle(45)
axis square
xlim([0.5,2.5])
set(gca,'box','off')
legend([s1,s2],'Blank-SAP','SSP-SAP')
axis square

subplot(2,3,4);
hold on
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbT.Rest.avg));
scatter(xInds*1,gcampData.Blank_SAP.RH.HbT.Rest.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(1,gcampData.Blank_SAP.RH.HbT.Rest.mean_avg,gcampData.Blank_SAP.RH.HbT.Rest.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbT.Rest.avg));
scatter(xInds*2,gcampData.SSP_SAP.RH.HbT.Rest.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,gcampData.SSP_SAP.RH.HbT.Rest.mean_avg,gcampData.SSP_SAP.RH.HbT.Rest.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbT.NREM.avg));
scatter(xInds*4,gcampData.Blank_SAP.RH.HbT.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(4,gcampData.Blank_SAP.RH.HbT.NREM.mean_avg,gcampData.Blank_SAP.RH.HbT.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbT.NREM.avg));
scatter(xInds*5,gcampData.SSP_SAP.RH.HbT.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,gcampData.SSP_SAP.RH.HbT.NREM.mean_avg,gcampData.SSP_SAP.RH.HbT.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbT.REM.avg));
scatter(xInds*7,gcampData.Blank_SAP.RH.HbT.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(7,gcampData.Blank_SAP.RH.HbT.REM.mean_avg,gcampData.Blank_SAP.RH.HbT.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbT.REM.avg));
scatter(xInds*8,gcampData.SSP_SAP.RH.HbT.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(8,gcampData.SSP_SAP.RH.HbT.REM.mean_avg,gcampData.SSP_SAP.RH.HbT.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('\Delta[HbT] (\muM)')
xlim([0,9])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(2,3,5);
hold on
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbO.Rest.avg));
scatter(xInds*1,gcampData.Blank_SAP.RH.HbO.Rest.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(1,gcampData.Blank_SAP.RH.HbO.Rest.mean_avg,gcampData.Blank_SAP.RH.HbO.Rest.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbO.Rest.avg));
scatter(xInds*2,gcampData.SSP_SAP.RH.HbO.Rest.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,gcampData.SSP_SAP.RH.HbO.Rest.mean_avg,gcampData.SSP_SAP.RH.HbO.Rest.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbO.NREM.avg));
scatter(xInds*4,gcampData.Blank_SAP.RH.HbO.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(4,gcampData.Blank_SAP.RH.HbO.NREM.mean_avg,gcampData.Blank_SAP.RH.HbO.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbO.NREM.avg));
scatter(xInds*5,gcampData.SSP_SAP.RH.HbO.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,gcampData.SSP_SAP.RH.HbO.NREM.mean_avg,gcampData.SSP_SAP.RH.HbO.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbO.REM.avg));
scatter(xInds*7,gcampData.Blank_SAP.RH.HbO.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(7,gcampData.Blank_SAP.RH.HbO.REM.mean_avg,gcampData.Blank_SAP.RH.HbO.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbO.REM.avg));
scatter(xInds*8,gcampData.SSP_SAP.RH.HbO.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(8,gcampData.SSP_SAP.RH.HbO.REM.mean_avg,gcampData.SSP_SAP.RH.HbO.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('\Delta[HbO] (\muM)')
xlim([0,9])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(2,3,6);
hold on
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbR.Rest.avg));
scatter(xInds*1,gcampData.Blank_SAP.RH.HbR.Rest.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(1,gcampData.Blank_SAP.RH.HbR.Rest.mean_avg,gcampData.Blank_SAP.RH.HbR.Rest.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbR.Rest.avg));
scatter(xInds*2,gcampData.SSP_SAP.RH.HbR.Rest.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,gcampData.SSP_SAP.RH.HbR.Rest.mean_avg,gcampData.SSP_SAP.RH.HbR.Rest.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbR.NREM.avg));
scatter(xInds*4,gcampData.Blank_SAP.RH.HbR.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(4,gcampData.Blank_SAP.RH.HbR.NREM.mean_avg,gcampData.Blank_SAP.RH.HbR.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbR.NREM.avg));
scatter(xInds*5,gcampData.SSP_SAP.RH.HbR.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,gcampData.SSP_SAP.RH.HbR.NREM.mean_avg,gcampData.SSP_SAP.RH.HbR.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbR.REM.avg));
scatter(xInds*7,gcampData.Blank_SAP.RH.HbR.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(7,gcampData.Blank_SAP.RH.HbR.REM.mean_avg,gcampData.Blank_SAP.RH.HbR.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbR.REM.avg));
scatter(xInds*8,gcampData.SSP_SAP.RH.HbR.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(8,gcampData.SSP_SAP.RH.HbR.REM.mean_avg,gcampData.SSP_SAP.RH.HbR.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('\Delta[HbR] (\muM)')
xlim([0,9])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

%% save figure
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figure Panels' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(FigS3,[dirpath 'FigS3']);
    set(FigS3,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'FigS3'])
    % statistical diary
    diaryFile = [dirpath 'FigS3_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on

    % [HbT] during each arousal-state
    disp('======================================================================================================================')
    disp('[HbT] during each arousal-state (ephys), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(ephysData.Blank_SAP.RH.HbT.Rest.mean_avg) ' +/- ' num2str(ephysData.Blank_SAP.RH.HbT.Rest.std_avg)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(ephysData.SSP_SAP.RH.HbT.Rest.mean_avg) ' +/- ' num2str(ephysData.SSP_SAP.RH.HbT.Rest.std_avg)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(EphysStats.RH.HbT.Rest.avg.p)]); disp(' ')
    disp(['Blank-SAP NREM: ' num2str(ephysData.Blank_SAP.RH.HbT.NREM.mean_avg) ' +/- ' num2str(ephysData.Blank_SAP.RH.HbT.NREM.std_avg)]); disp(' ')
    disp(['SSP-SAP NREM: ' num2str(ephysData.SSP_SAP.RH.HbT.NREM.mean_avg) ' +/- ' num2str(ephysData.SSP_SAP.RH.HbT.NREM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP NREM ttest p = ' num2str(EphysStats.RH.HbT.NREM.avg.p)]); disp(' ')
    disp(['Blank-SAP REM: ' num2str(ephysData.Blank_SAP.RH.HbT.REM.mean_avg) ' +/- ' num2str(ephysData.Blank_SAP.RH.HbT.REM.std_avg)]); disp(' ')
    disp(['SSP-SAP REM: ' num2str(ephysData.SSP_SAP.RH.HbT.REM.mean_avg) ' +/- ' num2str(ephysData.SSP_SAP.RH.HbT.REM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP REM ttest p = ' num2str(EphysStats.RH.HbT.REM.avg.p)]); disp(' ')
    disp(['Blank-SAP Iso: ' num2str(ephysData.Blank_SAP.RH.HbT.Iso.mean_avg) ' +/- ' num2str(ephysData.Blank_SAP.RH.HbT.Iso.std_avg)]); disp(' ')
    disp(['SSP-SAP Iso: ' num2str(ephysData.SSP_SAP.RH.HbT.Iso.mean_avg) ' +/- ' num2str(ephysData.SSP_SAP.RH.HbT.Iso.std_avg)]); disp(' ')
    disp(['Blank vs. SAP Iso ttest p = ' num2str(EphysStats.RH.HbT.Iso.avg.p)]); disp(' ')

    % resting-state 2P diameter
    disp('======================================================================================================================')
    disp('Resting-state 2P diameter, n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(diameterData.Blank_SAP.Rest.meanData) ' +/- ' num2str(diameterData.Blank_SAP.Rest.stdData)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(diameterData.SSP_SAP.Rest.meanData) ' +/- ' num2str(diameterData.SSP_SAP.Rest.stdData)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(DiameterStats.Rest.p)]); disp(' ')

    % isoflurane shift 2P diameter
    disp('======================================================================================================================')
    disp('isoflurane diameter shift, n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(diameterShiftData.Blank_SAP.meanDiameter) ' +/- ' num2str(diameterShiftData.Blank_SAP.stdDiameter) ' baseline diameter: ' num2str(diameterShiftData.Blank_SAP.meanBaseline) ' \muM']); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(diameterShiftData.SSP_SAP.meanDiameter) ' +/- ' num2str(diameterShiftData.SSP_SAP.stdDiameter) ' baseline diameter: ' num2str(diameterShiftData.SSP_SAP.meanBaseline) ' \muM']); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(TwoPIsoStats.p)]); disp(' ')

    % [HbT] during each arousal-state
    disp('======================================================================================================================')
    disp('[HbT] during each arousal-state (GCaMP), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(gcampData.Blank_SAP.RH.HbT.Rest.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbT.Rest.std_avg)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(gcampData.SSP_SAP.RH.HbT.Rest.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbT.Rest.std_avg)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(GCaMPStats.RH.HbT.Rest.avg.p)]); disp(' ')
    disp(['Blank-SAP NREM: ' num2str(gcampData.Blank_SAP.RH.HbT.NREM.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbT.NREM.std_avg)]); disp(' ')
    disp(['SSP-SAP NREM: ' num2str(gcampData.SSP_SAP.RH.HbT.NREM.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbT.NREM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP NREM ttest p = ' num2str(GCaMPStats.RH.HbT.NREM.avg.p)]); disp(' ')
    disp(['Blank-SAP REM: ' num2str(gcampData.Blank_SAP.RH.HbT.REM.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbT.REM.std_avg)]); disp(' ')
    disp(['SSP-SAP REM: ' num2str(gcampData.SSP_SAP.RH.HbT.REM.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbT.REM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP REM ttest p = ' num2str(GCaMPStats.RH.HbT.REM.avg.p)]); disp(' ')

    % [HbO] during each arousal-state
    disp('======================================================================================================================')
    disp('[HbO] during each arousal-state (GCaMP), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(gcampData.Blank_SAP.RH.HbO.Rest.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbO.Rest.std_avg)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(gcampData.SSP_SAP.RH.HbO.Rest.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbO.Rest.std_avg)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(GCaMPStats.RH.HbO.Rest.avg.p)]); disp(' ')
    disp(['Blank-SAP NREM: ' num2str(gcampData.Blank_SAP.RH.HbO.NREM.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbO.NREM.std_avg)]); disp(' ')
    disp(['SSP-SAP NREM: ' num2str(gcampData.SSP_SAP.RH.HbO.NREM.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbO.NREM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP NREM ttest p = ' num2str(GCaMPStats.RH.HbO.NREM.avg.p)]); disp(' ')
    disp(['Blank-SAP REM: ' num2str(gcampData.Blank_SAP.RH.HbO.REM.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbO.REM.std_avg)]); disp(' ')
    disp(['SSP-SAP REM: ' num2str(gcampData.SSP_SAP.RH.HbO.REM.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbO.REM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP REM ttest p = ' num2str(GCaMPStats.RH.HbO.REM.avg.p)]); disp(' ')

    % [HbR] during each arousal-state
    disp('======================================================================================================================')
    disp('[HbR] during each arousal-state (GCaMP), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(gcampData.Blank_SAP.RH.HbR.Rest.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbR.Rest.std_avg)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(gcampData.SSP_SAP.RH.HbR.Rest.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbR.Rest.std_avg)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(GCaMPStats.RH.HbR.Rest.avg.p)]); disp(' ')
    disp(['Blank-SAP NREM: ' num2str(gcampData.Blank_SAP.RH.HbR.NREM.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbR.NREM.std_avg)]); disp(' ')
    disp(['SSP-SAP NREM: ' num2str(gcampData.SSP_SAP.RH.HbR.NREM.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbR.NREM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP NREM ttest p = ' num2str(GCaMPStats.RH.HbR.NREM.avg.p)]); disp(' ')
    disp(['Blank-SAP REM: ' num2str(gcampData.Blank_SAP.RH.HbR.REM.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbR.REM.std_avg)]); disp(' ')
    disp(['SSP-SAP REM: ' num2str(gcampData.SSP_SAP.RH.HbR.REM.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbR.REM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP REM ttest p = ' num2str(GCaMPStats.RH.HbR.REM.avg.p)]); disp(' ')

    diary off
end
