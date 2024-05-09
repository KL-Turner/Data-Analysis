function [] = Fig4_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

%% IOS variance signals
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_IntSig_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
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
                iosSigData.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(iosSigData.(group).(hemisphere).(dataType),behavior) == false
                        iosSigData.(group).(hemisphere).(dataType).(behavior).avg = [];
                        iosSigData.(group).(hemisphere).(dataType).(behavior).p2p = [];
                        iosSigData.(group).(hemisphere).(dataType).(behavior).vari = [];
                        iosSigData.(group).(hemisphere).(dataType).(behavior).group = {};
                        iosSigData.(group).(hemisphere).(dataType).(behavior).animalID = {};
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
                    iosSigData.(group).(hemisphere).(dataType).(behavior).avg = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).HbT));
                    iosSigData.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                    iosSigData.(group).(hemisphere).(dataType).(behavior).vari = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                    iosSigData.(group).(hemisphere).(dataType).(behavior).group = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).group,group);
                    iosSigData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
                    iosSigData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(iosSigData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    iosSigData.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(iosSigData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                end
            end
        end
    end
end

%% IOS pulse variance
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_IntSig_Pulse';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
variables = {'avg','p2p','vari'};
fs = 30;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_IntSig_Pulse.(group));
    pulseSigData.(group).dummCheck = 1;
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        if isfield(pulseSigData.(group),'avg') == false
            pulseSigData.(group).avg = [];
            pulseSigData.(group).p2p = [];
            pulseSigData.(group).vari = [];
            pulseSigData.(group).group = {};
            pulseSigData.(group).animalID = {};
        end
        animalVar = [];
        animalP2P = [];
        if isfield(Results_IntSig_Pulse.(group).(animalID),'Rest') == true
            for ff = 1:length(Results_IntSig_Pulse.(group).(animalID).Rest.indHbT)
                dataArray = Results_IntSig_Pulse.(group).(animalID).Rest.indHbT{ff,1}(2*fs:end);
                animalVar(ff,1) = var(dataArray);
                animalP2P(ff,1) = max(dataArray) - min(dataArray);
            end
            pulseSigData.(group).avg = cat(1,pulseSigData.(group).avg,mean(Results_IntSig_Pulse.(group).(animalID).Rest.HbT,'omitnan'));
            pulseSigData.(group).p2p = cat(1,pulseSigData.(group).p2p,mean(animalP2P,'omitnan'));
            pulseSigData.(group).vari = cat(1,pulseSigData.(group).vari,mean(animalVar,'omitnan'));
            pulseSigData.(group).group = cat(1,pulseSigData.(group).group,group);
            pulseSigData.(group).animalID = cat(1,pulseSigData.(group).animalID,animalID);
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for ee = 1:length(variables)
        variable = variables{1,ee};
        pulseSigData.(group).(['mean_' variable]) = mean(pulseSigData.(group).(variable),1);
        pulseSigData.(group).(['std_' variable]) = std(pulseSigData.(group).(variable),1);
    end
end

%% GCaMP variance signals
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
                gcampSigdata.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(gcampSigdata.(group).(hemisphere).(dataType),behavior) == false
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg = [];
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).p2p = [];
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).vari = [];
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).group = {};
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).animalID = {};
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
                        if strcmp(dataType,'GCaMP') == true
                            dataArray = (dataArray - 1)*100;
                        end
                        animalVar(ff,1) = var(dataArray);
                        animalP2P(ff,1) = max(dataArray) - min(dataArray);
                    end
                    try
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).mean));
                    catch
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg,mean(cell2mat(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).mean)));
                    end
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).vari = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).group = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).group,group);
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(gcampSigdata.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(gcampSigdata.(group).(hemisphere).(dataType).(behavior).(variable),1);
                end
            end
        end
    end
end

%% 2P variance
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Diameter_2P';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
variables = {'avg','p2p','vari'};
fs = 5;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Diameter_2P.(group));
    diameterData.(group).dummCheck = 1;
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        vIDs = fieldnames(Results_Diameter_2P.(group).(animalID));
        for cc = 1:length(vIDs)
            vID = vIDs{cc,1};
            if isfield(diameterData.(group),'avg') == false
                diameterData.(group).avg = [];
                diameterData.(group).p2p = [];
                diameterData.(group).vari = [];
                diameterData.(group).group = {};
                diameterData.(group).animalID = {};
            end
            animalVar = [];
            animalP2P = [];
            if isfield(Results_Diameter_2P.(group).(animalID).(vID),'Rest') == true
                for ff = 1:length(Results_Diameter_2P.(group).(animalID).(vID).Rest.indEvents)
                    dataArray = Results_Diameter_2P.(group).(animalID).(vID).Rest.indEvents{ff,1}(2*fs:end);
                    animalVar(ff,1) = var(dataArray);
                    animalP2P(ff,1) = max(dataArray) - min(dataArray);
                end
                diameterData.(group).avg = cat(1,diameterData.(group).avg,mean(Results_Diameter_2P.(group).(animalID).(vID).Rest.mean,'omitnan'));
                diameterData.(group).p2p = cat(1,diameterData.(group).p2p,mean(animalP2P,'omitnan'));
                diameterData.(group).vari = cat(1,diameterData.(group).vari,mean(animalVar,'omitnan'));
                diameterData.(group).group = cat(1,diameterData.(group).group,group);
                diameterData.(group).animalID = cat(1,diameterData.(group).animalID,animalID);
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for ee = 1:length(variables)
        variable = variables{1,ee};
        diameterData.(group).(['mean_' variable]) = mean(diameterData.(group).(variable),1);
        diameterData.(group).(['std_' variable]) = std(diameterData.(group).(variable),1);
    end
end

%% HbT variance statistics
blankHbTVarData = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.vari,gcampSigdata.Blank_SAP.RH.HbT.Rest.vari,pulseSigData.Blank_SAP.vari);
sspHbTVarData = cat(1,iosSigData.SSP_SAP.RH.HbT.Rest.vari,gcampSigdata.SSP_SAP.RH.HbT.Rest.vari,pulseSigData.SSP_SAP.vari);
restHbTVarStats.tableSize = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.vari,iosSigData.SSP_SAP.RH.HbT.Rest.vari,gcampSigdata.Blank_SAP.RH.HbT.Rest.vari,gcampSigdata.SSP_SAP.RH.HbT.Rest.vari,pulseSigData.Blank_SAP.vari,pulseSigData.SSP_SAP.vari);
restHbTVarStats.Table = table('Size',[size(restHbTVarStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','Variance'});
restHbTVarStats.Table.Mouse = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.animalID,iosSigData.SSP_SAP.RH.HbT.Rest.animalID,gcampSigdata.Blank_SAP.RH.HbT.Rest.animalID,gcampSigdata.SSP_SAP.RH.HbT.Rest.animalID,pulseSigData.Blank_SAP.animalID,pulseSigData.SSP_SAP.animalID);
restHbTVarStats.Table.Group = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.group,iosSigData.SSP_SAP.RH.HbT.Rest.group,gcampSigdata.Blank_SAP.RH.HbT.Rest.group,gcampSigdata.SSP_SAP.RH.HbT.Rest.group,pulseSigData.Blank_SAP.group,pulseSigData.SSP_SAP.group);
restHbTVarStats.Table.Variance = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.vari,iosSigData.SSP_SAP.RH.HbT.Rest.vari,gcampSigdata.Blank_SAP.RH.HbT.Rest.vari,gcampSigdata.SSP_SAP.RH.HbT.Rest.vari,pulseSigData.Blank_SAP.vari,pulseSigData.SSP_SAP.vari);
restHbTVarStats.FitFormula = 'Variance ~ 1 + Group + (1|Mouse)';
restHbTVarStats.Stats = fitglme(restHbTVarStats.Table,restHbTVarStats.FitFormula);

%% Diameter variance statistics
blankDiameterVarData = cat(1,diameterData.Blank_SAP.vari);
sspDiameterVarData = cat(1,diameterData.SSP_SAP.vari);
restDiameterVarStats.tableSize = cat(1,diameterData.Blank_SAP.vari,diameterData.SSP_SAP.vari);
restDiameterVarStats.Table = table('Size',[size(restDiameterVarStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','Variance'});
restDiameterVarStats.Table.Mouse = cat(1,diameterData.Blank_SAP.animalID,diameterData.SSP_SAP.animalID);
restDiameterVarStats.Table.Group = cat(1,diameterData.Blank_SAP.group,diameterData.SSP_SAP.group);
% restDiameterVarStats.Table.Variance = cat(1,diameterData.Blank_SAP.vari,diameterData.SSP_SAP.vari);
restDiameterVarStats.Table.Variance = cat(1,diameterData.Blank_SAP.vari,sspDiameterVarData);
restDiameterVarStats.FitFormula = 'Variance ~ 1 + Group + (1|Mouse)';
restDiameterVarStats.Stats = fitglme(restDiameterVarStats.Table,restDiameterVarStats.FitFormula);

%% GCaMP variance statistics
blankGCaMPVarData = cat(1,gcampSigdata.Blank_SAP.RH.GCaMP.Rest.vari);
sspGCaMPVarData = cat(1,gcampSigdata.SSP_SAP.RH.GCaMP.Rest.vari);
restGCaMPVarStats.tableSize = cat(1,gcampSigdata.Blank_SAP.RH.GCaMP.Rest.vari,gcampSigdata.SSP_SAP.RH.GCaMP.Rest.vari);
restGCaMPVarStats.Table = table('Size',[size(restGCaMPVarStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','Variance'});
restGCaMPVarStats.Table.Mouse = cat(1,gcampSigdata.Blank_SAP.RH.GCaMP.Rest.animalID,gcampSigdata.SSP_SAP.RH.GCaMP.Rest.animalID);
restGCaMPVarStats.Table.Group = cat(1,gcampSigdata.Blank_SAP.RH.GCaMP.Rest.group,gcampSigdata.SSP_SAP.RH.GCaMP.Rest.group);
restGCaMPVarStats.Table.Variance = cat(1,gcampSigdata.Blank_SAP.RH.GCaMP.Rest.vari,gcampSigdata.SSP_SAP.RH.GCaMP.Rest.vari);
restGCaMPVarStats.FitFormula = 'Variance ~ 1 + Group + (1|Mouse)';
restGCaMPVarStats.Stats = fitglme(restGCaMPVarStats.Table,restGCaMPVarStats.FitFormula);

%% LFP power
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

%% statistics - generalized linear mixed effects model
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        lfpStats.(hemisphere).(behavior).tableSize = cat(1,data.Blank_SAP.(hemisphere).(behavior).deltaS,data.SSP_SAP.(hemisphere).(behavior).deltaS,data.Naive.(hemisphere).(behavior).deltaS);
        lfpStats.(hemisphere).(behavior).Table = table('Size',[size(lfpStats.(hemisphere).(behavior).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'group','animalID','behavior','deltaS'});
        lfpStats.(hemisphere).(behavior).Table.group = cat(1,data.Blank_SAP.(hemisphere).(behavior).group,data.SSP_SAP.(hemisphere).(behavior).group,data.Naive.(hemisphere).(behavior).group);
        lfpStats.(hemisphere).(behavior).Table.animalID = cat(1,data.Blank_SAP.(hemisphere).(behavior).animalID,data.SSP_SAP.(hemisphere).(behavior).animalID,data.Naive.(hemisphere).(behavior).animalID);
        lfpStats.(hemisphere).(behavior).Table.behavior = cat(1,data.Blank_SAP.(hemisphere).(behavior).behavior,data.SSP_SAP.(hemisphere).(behavior).behavior,data.Naive.(hemisphere).(behavior).behavior);
        lfpStats.(hemisphere).(behavior).Table.deltaS = cat(1,data.Blank_SAP.(hemisphere).(behavior).deltaS,data.SSP_SAP.(hemisphere).(behavior).deltaS,data.Naive.(hemisphere).(behavior).deltaS);
        lfpStats.(hemisphere).(behavior).FitFormula = 'deltaS ~ 1 + group + behavior + (1|animalID)';
        lfpStats.(hemisphere).(behavior).Stats = fitglme(lfpStats.(hemisphere).(behavior).Table,lfpStats.(hemisphere).(behavior).FitFormula);
    end
end

%% figure
Fig4 = figure('Name','Figure 4','units','normalized','outerposition',[0 0 1 1]);

% HbT rest variance
subplot(1,5,1)
xInds = ones(1,length(blankHbTVarData));
scatter(xInds*1,blankHbTVarData,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(blankHbTVarData,'omitnan'),std(blankHbTVarData,0,1,'omitnan'),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(sspHbTVarData));
scatter(xInds*2,sspHbTVarData,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,mean(sspHbTVarData,'omitnan'),std(sspHbTVarData,0,1,'omitnan'),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('\Delta[HbT]^2 (\muM)')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
axis tight
xlim([0,3]);
ylim([0,90])

% Diameter rest variance
subplot(1,5,2)
xInds = ones(1,length(blankDiameterVarData));
scatter(xInds*1,blankDiameterVarData,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(blankDiameterVarData,'omitnan'),std(blankDiameterVarData,0,1,'omitnan'),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(sspDiameterVarData));
scatter(xInds*2,sspDiameterVarData,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,mean(sspDiameterVarData,'omitnan'),std(sspDiameterVarData,0,1,'omitnan'),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('\Delta[Diameter]^2 (\muM)')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
axis tight
xlim([0,3]);
ylim([-5,60])

% Ephys alert LFP
ax1 = subplot(1,5,3);
loglog(data.Blank_SAP.RH.Alert.mean_f,data.Blank_SAP.RH.Alert.mean_S,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.RH.Alert.mean_f,data.Blank_SAP.RH.Alert.mean_S + data.Blank_SAP.RH.Alert.stdErr_S,'color',colors('north texas green'),'LineWidth',0.25);
loglog(data.Blank_SAP.RH.Alert.mean_f,data.Blank_SAP.RH.Alert.mean_S - data.Blank_SAP.RH.Alert.stdErr_S,'color',colors('north texas green'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.Alert.mean_f,data.SSP_SAP.RH.Alert.mean_S,'color',colors('electric purple'),'LineWidth',2);
loglog(data.SSP_SAP.RH.Alert.mean_f,data.SSP_SAP.RH.Alert.mean_S + data.SSP_SAP.RH.Alert.stdErr_S,'color',colors('electric purple'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.Alert.mean_f,data.SSP_SAP.RH.Alert.mean_S - data.SSP_SAP.RH.Alert.stdErr_S,'color',colors('electric purple'),'LineWidth',0.25);
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
axis square
axis tight
xlim([1,100]);

% Ephys asleep LFP
ax2 = subplot(1,5,4);
loglog(data.Blank_SAP.RH.Asleep.mean_f,data.Blank_SAP.RH.Asleep.mean_S,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.RH.Asleep.mean_f,data.Blank_SAP.RH.Asleep.mean_S + data.Blank_SAP.RH.Asleep.stdErr_S,'color',colors('north texas green'),'LineWidth',0.25);
loglog(data.Blank_SAP.RH.Asleep.mean_f,data.Blank_SAP.RH.Asleep.mean_S - data.Blank_SAP.RH.Asleep.stdErr_S,'color',colors('north texas green'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.Asleep.mean_f,data.SSP_SAP.RH.Asleep.mean_S,'color',colors('electric purple'),'LineWidth',2);
loglog(data.SSP_SAP.RH.Asleep.mean_f,data.SSP_SAP.RH.Asleep.mean_S + data.SSP_SAP.RH.Asleep.stdErr_S,'color',colors('electric purple'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.Asleep.mean_f,data.SSP_SAP.RH.Asleep.mean_S - data.SSP_SAP.RH.Asleep.stdErr_S,'color',colors('electric purple'),'LineWidth',0.25);
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
axis square
axis tight
xlim([1,100]);

% Ephys all LFP
ax3 = subplot(1,5,5);
loglog(data.Blank_SAP.RH.All.mean_f,data.Blank_SAP.RH.All.mean_S,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.RH.All.mean_f,data.Blank_SAP.RH.All.mean_S + data.Blank_SAP.RH.All.stdErr_S,'color',colors('north texas green'),'LineWidth',0.25);
loglog(data.Blank_SAP.RH.All.mean_f,data.Blank_SAP.RH.All.mean_S - data.Blank_SAP.RH.All.stdErr_S,'color',colors('north texas green'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.All.mean_f,data.SSP_SAP.RH.All.mean_S,'color',colors('electric purple'),'LineWidth',2);
loglog(data.SSP_SAP.RH.All.mean_f,data.SSP_SAP.RH.All.mean_S + data.SSP_SAP.RH.All.stdErr_S,'color',colors('electric purple'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.All.mean_f,data.SSP_SAP.RH.All.mean_S - data.SSP_SAP.RH.All.stdErr_S,'color',colors('electric purple'),'LineWidth',0.25);
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
axis square
axis tight
xlim([1,100]);
linkaxes([ax1,ax2,ax3],'xy')

%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figure Panels' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig4,[dirpath 'Fig4']);
    set(Fig4,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'Fig3'])
    diaryFile = [dirpath 'Fig3_Readout.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    % statistical diary
    diary(diaryFile)
    diary on

    % IOS [HbT] resting variance
    disp('IOS resting variance, n = 24 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(blankHbTVarData)) ' +/- ' num2str(std(blankHbTVarData,0,1)./sqrt(size(blankHbTVarData,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(sspHbTVarData)) ' +/- ' num2str(std(sspHbTVarData,0,1)./sqrt(size(sspHbTVarData,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for resting HbT variance')
    disp('======================================================================================================================')
    disp(restHbTVarStats.Stats)

    % IOS [HbT] resting variance
    disp('IOS resting variance, n = 24 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(blankHbTVarData)) ' +/- ' num2str(std(blankHbTVarData,0,1)./sqrt(size(blankHbTVarData,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(sspHbTVarData)) ' +/- ' num2str(std(sspHbTVarData,0,1)./sqrt(size(sspHbTVarData,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for resting HbT variance')
    disp('======================================================================================================================')
    disp(restHbTVarStats.Stats)

    % LFP delta power
    disp('LFP delta power, n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(data.Blank_SAP.RH.All.deltaS)) ' +/- ' num2str(std(data.Blank_SAP.RH.All.deltaS,0,1)./sqrt(size(data.Blank_SAP.RH.All.deltaS,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(data.SSP_SAP.RH.All.deltaS)) ' +/- ' num2str(std(data.SSP_SAP.RH.All.deltaS,0,1)./sqrt(size(data.SSP_SAP.RH.All.deltaS,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for LFP delta (1:4 Hz)')
    disp('======================================================================================================================')
    disp(lfpStats.RH.All.Stats)

    diary off
end