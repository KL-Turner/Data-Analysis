function [] = BilateralCoherence_Figures(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

% set-up and process data
resultsStruct = 'Results_BilatCoher';
load(resultsStruct);
% experimental groups and sets for analysis
expGroups = {'SSP_SAP','Blank_SAP'};
sets = {'IOS_Ephys','IOS_GCaMP7s';'IOS_Ephys','IOS_GCaMP7s'};
behavFields = {'Rest','NREM','REM','Alert','Asleep','All'};
% extract animal IDs for each group
for aa = 1:length(expGroups)
    setNames = sets(aa,:);
    for bb = 1:length(setNames)
        folderList = dir([expGroups{1,aa} delim setNames{1,bb}]);
        folderList = folderList(~startsWith({folderList.name}, '.'));
        data.(setNames{1,bb}).(expGroups{1,aa}).animalIDs = {folderList.name};
    end
end
% concatenate data from each data type
for aa = 1:length(expGroups)
    expGroup = expGroups{1,aa};
    setNames = sets(aa,:);
    for bb = 1:length(setNames)
        setName = setNames{1,bb};
        if strcmp(setName,'IOS_GCaMP7s') == true
            dataTypes = {'CBV_HbT','GCaMP7s','Deoxy'};
        elseif strcmp(setName,'IOS_Ephys') == true
            dataTypes = {'CBV_HbT','gammaBandPower'};
        end
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behavFields)
                behavField = behavFields{1,dd};
                % pre-allocate necessary variable fields
                data.(setName).(expGroup).(dataType).dummCheck = 1;
                if isfield(data.(setName).(expGroup).(dataType),behavField) == false
                    data.(setName).(expGroup).(dataType).(behavField).C = [];
                    data.(setName).(expGroup).(dataType).(behavField).f = [];
                end
                % go through each animal and concatenate data
                for ee = 1:length(data.(setName).(expGroup).animalIDs)
                    animalID = data.(setName).(expGroup).animalIDs{1,ee};
                    if isempty(Results_BilatCoher.(animalID).(behavField).(dataType).C) == false
                        data.(setName).(expGroup).(dataType).(behavField).C = cat(2,data.(setName).(expGroup).(dataType).(behavField).C,Results_BilatCoher.(animalID).(behavField).(dataType).C.^2);
                        data.(setName).(expGroup).(dataType).(behavField).f = cat(1,data.(setName).(expGroup).(dataType).(behavField).f,Results_BilatCoher.(animalID).(behavField).(dataType).f);
                    end
                end
            end
        end
    end
end
% average data from each data type
for aa = 1:length(expGroups)
    expGroup = expGroups{1,aa};
    setNames = sets(aa,:);
    for bb = 1:length(setNames)
        setName = setNames{1,bb};
        if strcmp(setName,'IOS_GCaMP7s') == true
            dataTypes = {'CBV_HbT','GCaMP7s','Deoxy'};
        elseif strcmp(setName,'IOS_Ephys') == true
            dataTypes = {'CBV_HbT','gammaBandPower'};
        end
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behavFields)
                behavField = behavFields{1,dd};
                data.(setName).(expGroup).(dataType).(behavField).meanC = mean(data.(setName).(expGroup).(dataType).(behavField).C,2);
                data.(setName).(expGroup).(dataType).(behavField).stdErrC = std(data.(setName).(expGroup).(dataType).(behavField).C,0,2)./sqrt(size(data.(setName).(expGroup).(dataType).(behavField).C,2));
                data.(setName).(expGroup).(dataType).(behavField).meanf = mean(data.(setName).(expGroup).(dataType).(behavField).f,1);
            end
        end
    end
end
%% average HbT coherence
summaryFigure = figure;
sgtitle('Ephys Bilateral Coherence \DeltaHbT')
%% coherence^2 between bilateral HbT during rest
subplot(2,3,1);
p1 = semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanC + data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanC - data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
p2 = semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanC + data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanC - data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Rest] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/10,0.5])
ylim([0.01,1])
set(gca,'box','off')
legend([p1,p2],'Blank-SAP','SSP-SAP')
%% coherence^2 between bilateral HbT during NREM
subplot(2,3,2);
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanC + data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanC - data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanC + data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanC - data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[NREM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/30,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during REM
subplot(2,3,3);
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanC + data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanC - data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanC + data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanC - data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[REM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/60,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Awake
subplot(2,3,4);
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanC + data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanC - data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanC + data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanC - data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(2,3,5);
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanC + data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanC - data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanC + data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanC - data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(2,3,6);
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanC + data.IOS_Ephys.Blank_SAP.CBV_HbT.All.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanC - data.IOS_Ephys.Blank_SAP.CBV_HbT.All.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanC + data.IOS_Ephys.SSP_SAP.CBV_HbT.All.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanC - data.IOS_Ephys.SSP_SAP.CBV_HbT.All.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0.01,1])
set(gca,'box','off')


%%
%% average HbT coherence
summaryFigure = figure;
sgtitle('Ephys Bilateral Coherence \DeltaHbT')
%% coherence^2 between bilateral HbT during rest
subplot(2,3,1);
p1 = semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.Rest.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.Rest.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.Rest.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.Rest.meanC + data.IOS_Ephys.Blank_SAP.gammaBandPower.Rest.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.Rest.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.Rest.meanC - data.IOS_Ephys.Blank_SAP.gammaBandPower.Rest.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
p2 = semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.Rest.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.Rest.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.Rest.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.Rest.meanC + data.IOS_Ephys.SSP_SAP.gammaBandPower.Rest.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.Rest.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.Rest.meanC - data.IOS_Ephys.SSP_SAP.gammaBandPower.Rest.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Rest] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/10,0.5])
ylim([0.01,1])
set(gca,'box','off')
legend([p1,p2],'Blank-SAP','SSP-SAP')
%% coherence^2 between bilateral HbT during NREM
subplot(2,3,2);
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.NREM.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.NREM.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.NREM.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.NREM.meanC + data.IOS_Ephys.Blank_SAP.gammaBandPower.NREM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.NREM.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.NREM.meanC - data.IOS_Ephys.Blank_SAP.gammaBandPower.NREM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.NREM.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.NREM.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.NREM.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.NREM.meanC + data.IOS_Ephys.SSP_SAP.gammaBandPower.NREM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.NREM.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.NREM.meanC - data.IOS_Ephys.SSP_SAP.gammaBandPower.NREM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[NREM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/30,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during REM
subplot(2,3,3);
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.REM.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.REM.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.REM.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.REM.meanC + data.IOS_Ephys.Blank_SAP.gammaBandPower.REM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.REM.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.REM.meanC - data.IOS_Ephys.Blank_SAP.gammaBandPower.REM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.REM.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.REM.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.REM.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.REM.meanC + data.IOS_Ephys.SSP_SAP.gammaBandPower.REM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.REM.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.REM.meanC - data.IOS_Ephys.SSP_SAP.gammaBandPower.REM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[REM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/60,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Awake
subplot(2,3,4);
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.Alert.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.Alert.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.Alert.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.Alert.meanC + data.IOS_Ephys.Blank_SAP.gammaBandPower.Alert.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.Alert.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.Alert.meanC - data.IOS_Ephys.Blank_SAP.gammaBandPower.Alert.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.Alert.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.Alert.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.Alert.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.Alert.meanC + data.IOS_Ephys.SSP_SAP.gammaBandPower.Alert.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.Alert.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.Alert.meanC - data.IOS_Ephys.SSP_SAP.gammaBandPower.Alert.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(2,3,5);
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.Asleep.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.Asleep.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.Asleep.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.Asleep.meanC + data.IOS_Ephys.Blank_SAP.gammaBandPower.Asleep.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.Asleep.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.Asleep.meanC - data.IOS_Ephys.Blank_SAP.gammaBandPower.Asleep.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.Asleep.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.Asleep.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.Asleep.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.Asleep.meanC + data.IOS_Ephys.SSP_SAP.gammaBandPower.Asleep.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.Asleep.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.Asleep.meanC - data.IOS_Ephys.SSP_SAP.gammaBandPower.Asleep.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(2,3,6);
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.All.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.All.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.All.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.All.meanC + data.IOS_Ephys.Blank_SAP.gammaBandPower.All.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.Blank_SAP.gammaBandPower.All.meanf,data.IOS_Ephys.Blank_SAP.gammaBandPower.All.meanC - data.IOS_Ephys.Blank_SAP.gammaBandPower.All.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.All.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.All.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.All.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.All.meanC + data.IOS_Ephys.SSP_SAP.gammaBandPower.All.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_Ephys.SSP_SAP.gammaBandPower.All.meanf,data.IOS_Ephys.SSP_SAP.gammaBandPower.All.meanC - data.IOS_Ephys.SSP_SAP.gammaBandPower.All.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0.01,1])
set(gca,'box','off')

%%
%% average HbT coherence
summaryFigure = figure;
sgtitle('Ephys Bilateral Coherence \DeltaHbT')
%% coherence^2 between bilateral HbT during rest
subplot(2,3,1);
p1 = semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Rest.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Rest.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Rest.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Rest.meanC + data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Rest.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Rest.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Rest.meanC - data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Rest.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
p2 = semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Rest.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Rest.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Rest.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Rest.meanC + data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Rest.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Rest.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Rest.meanC - data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Rest.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Rest] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/10,0.5])
ylim([0.01,1])
set(gca,'box','off')
legend([p1,p2],'Blank-SAP','SSP-SAP')
%% coherence^2 between bilateral HbT during NREM
subplot(2,3,2);
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.NREM.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.NREM.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.NREM.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.NREM.meanC + data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.NREM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.NREM.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.NREM.meanC - data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.NREM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.NREM.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.NREM.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.NREM.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.NREM.meanC + data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.NREM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.NREM.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.NREM.meanC - data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.NREM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[NREM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/30,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during REM
subplot(2,3,3);
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.REM.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.REM.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.REM.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.REM.meanC + data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.REM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.REM.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.REM.meanC - data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.REM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.REM.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.REM.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.REM.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.REM.meanC + data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.REM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.REM.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.REM.meanC - data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.REM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[REM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/60,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Awake
subplot(2,3,4);
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Alert.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Alert.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Alert.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Alert.meanC + data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Alert.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Alert.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Alert.meanC - data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Alert.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Alert.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Alert.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Alert.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Alert.meanC + data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Alert.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Alert.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Alert.meanC - data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Alert.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(2,3,5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Asleep.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Asleep.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Asleep.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Asleep.meanC + data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Asleep.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Asleep.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Asleep.meanC - data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.Asleep.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Asleep.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Asleep.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Asleep.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Asleep.meanC + data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Asleep.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Asleep.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Asleep.meanC - data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.Asleep.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(2,3,6);
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.All.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.All.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.All.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.All.meanC + data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.All.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.All.meanf,data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.All.meanC - data.IOS_GCaMP7s.Blank_SAP.CBV_HbT.All.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.All.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.All.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.All.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.All.meanC + data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.All.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.All.meanf,data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.All.meanC - data.IOS_GCaMP7s.SSP_SAP.CBV_HbT.All.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0.01,1])
set(gca,'box','off')

%%
%% average HbT coherence
summaryFigure = figure;
sgtitle('Ephys Bilateral Coherence \DeltaHbT')
%% coherence^2 between bilateral HbT during rest
subplot(2,3,1);
p1 = semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Rest.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Rest.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Rest.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Rest.meanC + data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Rest.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Rest.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Rest.meanC - data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Rest.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
p2 = semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Rest.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Rest.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Rest.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Rest.meanC + data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Rest.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Rest.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Rest.meanC - data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Rest.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Rest] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/10,0.5])
ylim([0.01,1])
set(gca,'box','off')
legend([p1,p2],'Blank-SAP','SSP-SAP')
%% coherence^2 between bilateral HbT during NREM
subplot(2,3,2);
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.NREM.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.NREM.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.NREM.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.NREM.meanC + data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.NREM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.NREM.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.NREM.meanC - data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.NREM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.NREM.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.NREM.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.NREM.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.NREM.meanC + data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.NREM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.NREM.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.NREM.meanC - data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.NREM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[NREM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/30,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during REM
subplot(2,3,3);
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.REM.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.REM.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.REM.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.REM.meanC + data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.REM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.REM.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.REM.meanC - data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.REM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.REM.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.REM.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.REM.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.REM.meanC + data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.REM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.REM.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.REM.meanC - data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.REM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[REM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/60,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Awake
subplot(2,3,4);
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Alert.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Alert.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Alert.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Alert.meanC + data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Alert.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Alert.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Alert.meanC - data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Alert.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Alert.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Alert.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Alert.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Alert.meanC + data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Alert.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Alert.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Alert.meanC - data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Alert.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(2,3,5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Asleep.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Asleep.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Asleep.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Asleep.meanC + data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Asleep.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Asleep.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Asleep.meanC - data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.Asleep.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Asleep.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Asleep.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Asleep.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Asleep.meanC + data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Asleep.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Asleep.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Asleep.meanC - data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.Asleep.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0.01,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(2,3,6);
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.All.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.All.meanC,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.All.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.All.meanC + data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.All.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.All.meanf,data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.All.meanC - data.IOS_GCaMP7s.Blank_SAP.GCaMP7s.All.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.All.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.All.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.All.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.All.meanC + data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.All.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.All.meanf,data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.All.meanC - data.IOS_GCaMP7s.SSP_SAP.GCaMP7s.All.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0.01,1])
set(gca,'box','off')
% 
% 
% 
% 
% 
% %% find Hz peaks in coherence
% treatments2 = {'SSP_SAP','Blank_SAP'};
% behavFields2 = {'Awake','Sleep','All'};
% for qq = 1:length(treatments2)
%     treatment = treatments2{1,qq};
%     for ee = 1:length(behavFields2)
%         behavField = behavFields2{1,ee};
%         for ff = 1:length(dataTypes)
%             dataType = dataTypes{1,ff};
%             for gg = 1:size(data.(treatment).(behavField).(dataType).C,2)
%                 F = round(data.(treatment).(behavField).(dataType).f(gg,:),3);
%                 C = data.(treatment).(behavField).(dataType).C(:,gg);
%                 index001 = find(F == 0.01);
%                 index01 = find(F == 0.1);
%                 index05 = find(F == 0.5);
%                 data.(treatment).(behavField).(dataType).C001(gg,1) = mean(C(1:index001(1)));
%                 data.(treatment).(behavField).(dataType).C01(gg,1) = mean(C(index001(1) + 1:index01(1)));
%                 data.(treatment).(behavField).(dataType).C05(gg,1) = mean(C(index01(1) + 1:index05(1)));
%             end
%         end
%     end
% end
% %% statistics - generalized linear mixed effects model
% freqs = {'C001','C01','C05'};
% for aa = 1:length(freqs)
%     freq = freqs{1,aa};
%     for bb = 1:length(behavFields2)
%         behavField = behavFields2{1,bb};
%         for cc = 1:length(dataTypes)
%             dataType = dataTypes{1,cc};
%             % statistics - generalized linear mixed effects model
%             Stats.(dataType).(behavField).(freq).tableSize = cat(1,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).(freq),data.IOS_Ephys.SSP_SAP.(behavField).(dataType).(freq));
%             Stats.(dataType).(behavField).(freq).Table = table('Size',[size(Stats.(dataType).(behavField).(freq).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Treatment','Coherence'});
%             Stats.(dataType).(behavField).(freq).Table.Mouse = cat(1,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).animalID,data.IOS_Ephys.SSP_SAP.(behavField).(dataType).animalID);
%             Stats.(dataType).(behavField).(freq).Table.Treatment = cat(1,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).treatment,data.IOS_Ephys.SSP_SAP.(behavField).(dataType).treatment);
%             Stats.(dataType).(behavField).(freq).Table.Coherence = cat(1,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).(freq),data.IOS_Ephys.SSP_SAP.(behavField).(dataType).(freq));
%             Stats.(dataType).(behavField).(freq).FitFormula = 'Coherence ~ 1 + Treatment + (1|Mouse)';
%             Stats.(dataType).(behavField).(freq).Stats = fitglme(Stats.(dataType).(behavField).(freq).Table,Stats.(dataType).(behavField).(freq).FitFormula);
%         end
%     end
% end
% %%
% % for bb = 1:length(behavFields)
% %     behavField = behavFields{1,bb};
% %     for cc = 1:length(dataTypes)
% %         dataType = dataTypes{1,cc};
% %         % statistics - generalized linear mixed effects model
% %         Stats.(dataType).(behavField).tableSize = cat(1,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).RH.C001,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).RH.C01,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).RH.C05,...
% %             data.IOS_Ephys.SSP_SAP.(behavField).(dataType).RH.C001,data.IOS_Ephys.SSP_SAP.(behavField).(dataType).RH.C01,data.IOS_Ephys.SSP_SAP.(behavField).(dataType).RH.C05);
% %         Stats.(dataType).(behavField).Table = table('Size',[size(Stats.(dataType).(behavField).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Treatment','Frequency','Coherence'});
% %         Stats.(dataType).(behavField).Table.Mouse = cat(1,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).animalID,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).animalID,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).animalID,...
% %             data.IOS_Ephys.SSP_SAP.(behavField).(dataType).animalID,data.IOS_Ephys.SSP_SAP.(behavField).(dataType).animalID,data.IOS_Ephys.SSP_SAP.(behavField).(dataType).animalID);
% %         Stats.(dataType).(behavField).Table.Treatment = cat(1,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).treatment,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).treatment,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).treatment,...,
% %             data.IOS_Ephys.SSP_SAP.(behavField).(dataType).treatment,data.IOS_Ephys.SSP_SAP.(behavField).(dataType).treatment,data.IOS_Ephys.SSP_SAP.(behavField).(dataType).treatment);
% %         Stats.(dataType).(behavField).Table.Frequency = cat(1,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).freqC001,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).freqC01,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).freqC05,...
% %             data.IOS_Ephys.SSP_SAP.(behavField).(dataType).freqC001,data.IOS_Ephys.SSP_SAP.(behavField).(dataType).freqC01,data.IOS_Ephys.SSP_SAP.(behavField).(dataType).freqC05);
% %         Stats.(dataType).(behavField).Table.Coherence = cat(1,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).RH.C001,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).RH.C01,data.IOS_Ephys.Blank_SAP.(behavField).(dataType).RH.C05,...
% %             data.IOS_Ephys.SSP_SAP.(behavField).(dataType).RH.C001,data.IOS_Ephys.SSP_SAP.(behavField).(dataType).RH.C01,data.IOS_Ephys.SSP_SAP.(behavField).(dataType).RH.C05);
% %         Stats.(dataType).(behavField).FitFormula = 'Coherence ~ 1 + Treatment + Frequency + Frequency*Treatment + (1|Mouse)';
% %         Stats.(dataType).(behavField).Stats = fitglme(Stats.(dataType).(behavField).Table,Stats.(dataType).(behavField).FitFormula);
% %     end
% % end
% %% average HbT coherence
% summaryFigure1 = figure;
% sgtitle('Bilateral Coherence \DeltaHbT')
% %% coherence^2 between bilateral HbT during rest
% subplot(2,3,1);
% p1 = semilogx(data.Naive.CBV_HbT.Rest.meanf,data.Naive.CBV_HbT.Rest.meanC,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Naive.CBV_HbT.Rest.meanf,data.Naive.CBV_HbT.Rest.meanC + data.Naive.CBV_HbT.Rest.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.Naive.CBV_HbT.Rest.meanf,data.Naive.CBV_HbT.Rest.meanC - data.Naive.CBV_HbT.Rest.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% p1 = semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanC,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanC + data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanC - data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% p2 = semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanC,'color',colors('electric purple'),'LineWidth',2);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanC + data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanC - data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Rest] Bilateral coherence^2','\Delta[HbT] (\muM)'})
% xlim([1/10,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% legend([p1,p1,p2],'Naive','Blank-SAP','SSP-SAP')
% %% coherence^2 between bilateral HbT during NREM
% subplot(2,3,2);
% semilogx(data.Naive.CBV_HbT.NREM.meanf,data.Naive.CBV_HbT.NREM.meanC,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Naive.CBV_HbT.NREM.meanf,data.Naive.CBV_HbT.NREM.meanC + data.Naive.CBV_HbT.NREM.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.Naive.CBV_HbT.NREM.meanf,data.Naive.CBV_HbT.NREM.meanC - data.Naive.CBV_HbT.NREM.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanC,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanC + data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanC - data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanC,'color',colors('electric purple'),'LineWidth',2);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanC + data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanC - data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[NREM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
% xlim([1/30,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral HbT during REM
% subplot(2,3,3);
% semilogx(data.Naive.CBV_HbT.REM.meanf,data.Naive.CBV_HbT.REM.meanC,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Naive.CBV_HbT.REM.meanf,data.Naive.CBV_HbT.REM.meanC + data.Naive.CBV_HbT.REM.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.Naive.CBV_HbT.REM.meanf,data.Naive.CBV_HbT.REM.meanC - data.Naive.CBV_HbT.REM.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanC,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanC + data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanC - data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanC,'color',colors('electric purple'),'LineWidth',2);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanC + data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanC - data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[REM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
% xlim([1/60,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral HbT during Awake
% subplot(2,3,4);
% semilogx(data.Naive.CBV_HbT.Alert.meanf,data.Naive.CBV_HbT.Alert.meanC,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Naive.CBV_HbT.Alert.meanf,data.Naive.CBV_HbT.Alert.meanC + data.Naive.CBV_HbT.Alert.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.Naive.CBV_HbT.Alert.meanf,data.Naive.CBV_HbT.Alert.meanC - data.Naive.CBV_HbT.Alert.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanC,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanC + data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanC - data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanC,'color',colors('electric purple'),'LineWidth',2);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanC + data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanC - data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Alert] Bilateral coherence^2','\Delta[HbT] (\muM)'})
% xlim([0.003,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral HbT during Sleep
% subplot(2,3,5);
% semilogx(data.Naive.CBV_HbT.Asleep.meanf,data.Naive.CBV_HbT.Asleep.meanC,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Naive.CBV_HbT.Asleep.meanf,data.Naive.CBV_HbT.Asleep.meanC + data.Naive.CBV_HbT.Asleep.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.Naive.CBV_HbT.Asleep.meanf,data.Naive.CBV_HbT.Asleep.meanC - data.Naive.CBV_HbT.Asleep.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanC,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanC + data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanC - data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanC,'color',colors('electric purple'),'LineWidth',2);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanC + data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanC - data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Asleep] Bilateral coherence^2','\Delta[HbT] (\muM)'})
% xlim([0.003,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral HbT during All data
% subplot(2,3,6);
% semilogx(data.Naive.CBV_HbT.All.meanf,data.Naive.CBV_HbT.All.meanC,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Naive.CBV_HbT.All.meanf,data.Naive.CBV_HbT.All.meanC + data.Naive.CBV_HbT.All.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.Naive.CBV_HbT.All.meanf,data.Naive.CBV_HbT.All.meanC - data.Naive.CBV_HbT.All.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanC,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanC + data.IOS_Ephys.Blank_SAP.CBV_HbT.All.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanC - data.IOS_Ephys.Blank_SAP.CBV_HbT.All.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanC,'color',colors('electric purple'),'LineWidth',2);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanC + data.IOS_Ephys.SSP_SAP.CBV_HbT.All.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanC - data.IOS_Ephys.SSP_SAP.CBV_HbT.All.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[All] Bilateral coherence^2','\Delta[HbT] (\muM)'})
% xlim([0.003,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% save figure(s)
% if saveFigs == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Bilateral Coherence - Bilateral IOS' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure1,[dirpath 'AverageBilateralCoherence_HbT']);
%     set(summaryFigure1,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'AverageBilateralCoherence_HbT'])
% end
% %% individual HbT coherence
% summaryFigure2 = figure;
% sgtitle('Bilateral Coherence \DeltaHbT - individual animals')
% %% coherence^2 between bilateral HbT during rest
% subplot(2,3,1);
% % Naive
% for aa = 1:size(data.Naive.CBV_HbT.Rest.C,2)
%     semilogx(data.Naive.CBV_HbT.Rest.meanf,data.Naive.CBV_HbT.Rest.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.C,2)
%     semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Rest.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.C,2)
%     semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Rest.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Rest] Bilateral coherence^2','\Delta[HbT] (\muM)'})
% xlim([1/10,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral HbT during NREM
% subplot(2,3,2);
% % Naive
% for aa = 1:size(data.Naive.CBV_HbT.NREM.C,2)
%     semilogx(data.Naive.CBV_HbT.NREM.meanf,data.Naive.CBV_HbT.NREM.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.C,2)
%     semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.NREM.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.C,2)
%     semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.NREM.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[NREM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
% xlim([1/30,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral HbT during REM
% subplot(2,3,3);
% % Naive
% for aa = 1:size(data.Naive.CBV_HbT.REM.C,2)
%     semilogx(data.Naive.CBV_HbT.REM.meanf,data.Naive.CBV_HbT.REM.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.C,2)
%     semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.REM.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.C,2)
%     semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.REM.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[REM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
% xlim([0.01,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral HbT during Awake
% subplot(2,3,4);
% % Naive
% for aa = 1:size(data.Naive.CBV_HbT.Alert.C,2)
%     semilogx(data.Naive.CBV_HbT.Alert.meanf,data.Naive.CBV_HbT.Alert.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C,2)
%     semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.C,2)
%     semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Alert] Bilateral coherence^2','\Delta[HbT] (\muM)'})
% xlim([0.003,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral HbT during Sleep
% subplot(2,3,5);
% % Naive
% for aa = 1:size(data.Naive.CBV_HbT.Asleep.C,2)
%     semilogx(data.Naive.CBV_HbT.Asleep.meanf,data.Naive.CBV_HbT.Asleep.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C,2)
%     semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.C,2)
%     semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Asleep] Bilateral coherence^2','\Delta[HbT] (\muM)'})
% xlim([0.003,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral HbT during All data
% subplot(2,3,6);
% % Naive
% for aa = 1:size(data.Naive.CBV_HbT.All.C,2)
%     semilogx(data.Naive.CBV_HbT.All.meanf,data.Naive.CBV_HbT.All.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C,2)
%     semilogx(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.IOS_Ephys.SSP_SAP.CBV_HbT.All.C,2)
%     semilogx(data.IOS_Ephys.SSP_SAP.CBV_HbT.All.meanf,data.IOS_Ephys.SSP_SAP.CBV_HbT.All.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[All] Bilateral coherence^2','\Delta[HbT] (\muM)'})
% xlim([0.003,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% save figure(s)
% if strcmp(saveFigs,'y') == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Bilateral Coherence - Bilateral IOS' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure2,[dirpath 'IndividualCoherence_HbT']);
%     set(summaryFigure2,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'IndividualCoherence_HbT'])
% end
% %% average gamma-band coherence
% summaryFigure3 = figure;
% sgtitle('Bilateral Coherence Gamma-band [30-100 Hz] (envelope)')
% %% coherence^2 between bilateral gamma-band during Rest
% subplot(2,3,1);
% p1 = semilogx(data.Naive.Rest.gammaBandPower.meanf,data.Naive.Rest.gammaBandPower.meanC,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Naive.Rest.gammaBandPower.meanf,data.Naive.Rest.gammaBandPower.meanC + data.Naive.Rest.gammaBandPower.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.Naive.Rest.gammaBandPower.meanf,data.Naive.Rest.gammaBandPower.meanC - data.Naive.Rest.gammaBandPower.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% p1 = semilogx(data.IOS_Ephys.Blank_SAP.Rest.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.Rest.gammaBandPower.meanC,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.IOS_Ephys.Blank_SAP.Rest.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.Rest.gammaBandPower.meanC + data.IOS_Ephys.Blank_SAP.Rest.gammaBandPower.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.Rest.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.Rest.gammaBandPower.meanC - data.IOS_Ephys.Blank_SAP.Rest.gammaBandPower.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% p2 = semilogx(data.IOS_Ephys.SSP_SAP.Rest.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.Rest.gammaBandPower.meanC,'color',colors('electric purple'),'LineWidth',2);
% semilogx(data.IOS_Ephys.SSP_SAP.Rest.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.Rest.gammaBandPower.meanC + data.IOS_Ephys.SSP_SAP.Rest.gammaBandPower.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.Rest.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.Rest.gammaBandPower.meanC - data.IOS_Ephys.SSP_SAP.Rest.gammaBandPower.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Rest] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
% xlim([1/10,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% legend([p1,p1,p2],'Naive','Blank-SAP','SSP-SAP')
% %% coherence^2 between bilateral gamma-band during NREM
% subplot(2,3,2);
% semilogx(data.Naive.NREM.gammaBandPower.meanf,data.Naive.NREM.gammaBandPower.meanC,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Naive.NREM.gammaBandPower.meanf,data.Naive.NREM.gammaBandPower.meanC + data.Naive.NREM.gammaBandPower.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.Naive.NREM.gammaBandPower.meanf,data.Naive.NREM.gammaBandPower.meanC - data.Naive.NREM.gammaBandPower.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.NREM.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.NREM.gammaBandPower.meanC,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.IOS_Ephys.Blank_SAP.NREM.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.NREM.gammaBandPower.meanC + data.IOS_Ephys.Blank_SAP.NREM.gammaBandPower.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.NREM.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.NREM.gammaBandPower.meanC - data.IOS_Ephys.Blank_SAP.NREM.gammaBandPower.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.NREM.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.NREM.gammaBandPower.meanC,'color',colors('electric purple'),'LineWidth',2);
% semilogx(data.IOS_Ephys.SSP_SAP.NREM.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.NREM.gammaBandPower.meanC + data.IOS_Ephys.SSP_SAP.NREM.gammaBandPower.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.NREM.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.NREM.gammaBandPower.meanC - data.IOS_Ephys.SSP_SAP.NREM.gammaBandPower.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[NREM] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
% xlim([1/30,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral gamma-band during REM
% subplot(2,3,3);
% semilogx(data.Naive.REM.gammaBandPower.meanf,data.Naive.REM.gammaBandPower.meanC,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Naive.REM.gammaBandPower.meanf,data.Naive.REM.gammaBandPower.meanC + data.Naive.REM.gammaBandPower.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.Naive.REM.gammaBandPower.meanf,data.Naive.REM.gammaBandPower.meanC - data.Naive.REM.gammaBandPower.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.REM.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.REM.gammaBandPower.meanC,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.IOS_Ephys.Blank_SAP.REM.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.REM.gammaBandPower.meanC + data.IOS_Ephys.Blank_SAP.REM.gammaBandPower.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.REM.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.REM.gammaBandPower.meanC - data.IOS_Ephys.Blank_SAP.REM.gammaBandPower.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.REM.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.REM.gammaBandPower.meanC,'color',colors('electric purple'),'LineWidth',2);
% semilogx(data.IOS_Ephys.SSP_SAP.REM.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.REM.gammaBandPower.meanC + data.IOS_Ephys.SSP_SAP.REM.gammaBandPower.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.REM.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.REM.gammaBandPower.meanC - data.IOS_Ephys.SSP_SAP.REM.gammaBandPower.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[REM] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
% xlim([1/60,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral gamma-band during Awake
% subplot(2,3,4);
% semilogx(data.Naive.Awake.gammaBandPower.meanf,data.Naive.Awake.gammaBandPower.meanC,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Naive.Awake.gammaBandPower.meanf,data.Naive.Awake.gammaBandPower.meanC + data.Naive.Awake.gammaBandPower.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.Naive.Awake.gammaBandPower.meanf,data.Naive.Awake.gammaBandPower.meanC - data.Naive.Awake.gammaBandPower.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.meanC,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.meanC + data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.meanC - data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.meanC,'color',colors('electric purple'),'LineWidth',2);
% semilogx(data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.meanC + data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.meanC - data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Alert] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
% xlim([0.003,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral gamma-band during Sleep
% subplot(2,3,5);
% semilogx(data.Naive.Sleep.gammaBandPower.meanf,data.Naive.Sleep.gammaBandPower.meanC,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Naive.Sleep.gammaBandPower.meanf,data.Naive.Sleep.gammaBandPower.meanC + data.Naive.Sleep.gammaBandPower.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.Naive.Sleep.gammaBandPower.meanf,data.Naive.Sleep.gammaBandPower.meanC - data.Naive.Sleep.gammaBandPower.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.meanC,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.meanC + data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.meanC - data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.meanC,'color',colors('electric purple'),'LineWidth',2);
% semilogx(data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.meanC + data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.meanC - data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Asleep] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
% xlim([0.003,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral gamma-band during All data
% subplot(2,3,6);
% semilogx(data.Naive.All.gammaBandPower.meanf,data.Naive.All.gammaBandPower.meanC,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Naive.All.gammaBandPower.meanf,data.Naive.All.gammaBandPower.meanC + data.Naive.All.gammaBandPower.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.Naive.All.gammaBandPower.meanf,data.Naive.All.gammaBandPower.meanC - data.Naive.All.gammaBandPower.stdErrC,'color',colors('sapphire'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.All.gammaBandPower.meanC,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.All.gammaBandPower.meanC + data.IOS_Ephys.Blank_SAP.All.gammaBandPower.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.All.gammaBandPower.meanC - data.IOS_Ephys.Blank_SAP.All.gammaBandPower.stdErrC,'color',colors('north texas green'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.All.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.All.gammaBandPower.meanC,'color',colors('electric purple'),'LineWidth',2);
% semilogx(data.IOS_Ephys.SSP_SAP.All.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.All.gammaBandPower.meanC + data.IOS_Ephys.SSP_SAP.All.gammaBandPower.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% semilogx(data.IOS_Ephys.SSP_SAP.All.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.All.gammaBandPower.meanC - data.IOS_Ephys.SSP_SAP.All.gammaBandPower.stdErrC,'color',colors('electric purple'),'LineWidth',0.5);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[All] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
% xlim([0.003,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% save figure(s)
% if saveFigs == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Bilateral Coherence - Bilateral IOS' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure3,[dirpath 'AverageBilateralCoherence_Gamma']);
%     set(summaryFigure3,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'AverageBilateralCoherence_Gamma'])
% end
% %% individual gamma-band coherence
% summaryFigure4 = figure;
% sgtitle('Bilateral Coherence Gamma-band [30-100 Hz] - individual animals')
% %% coherence^2 between bilateral gamma-band during rest
% subplot(2,3,1);
% % Naive
% for aa = 1:size(data.Naive.Rest.gammaBandPower.C,2)
%     semilogx(data.Naive.Rest.gammaBandPower.meanf,data.Naive.Rest.gammaBandPower.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.IOS_Ephys.Blank_SAP.Rest.gammaBandPower.C,2)
%     semilogx(data.IOS_Ephys.Blank_SAP.Rest.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.Rest.gammaBandPower.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.IOS_Ephys.SSP_SAP.Rest.gammaBandPower.C,2)
%     semilogx(data.IOS_Ephys.SSP_SAP.Rest.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.Rest.gammaBandPower.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Rest] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
% xlim([1/10,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral gamma-band during NREM
% subplot(2,3,2);
% % Naive
% for aa = 1:size(data.Naive.NREM.gammaBandPower.C,2)
%     semilogx(data.Naive.NREM.gammaBandPower.meanf,data.Naive.NREM.gammaBandPower.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.IOS_Ephys.Blank_SAP.NREM.gammaBandPower.C,2)
%     semilogx(data.IOS_Ephys.Blank_SAP.NREM.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.NREM.gammaBandPower.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.IOS_Ephys.SSP_SAP.NREM.gammaBandPower.C,2)
%     semilogx(data.IOS_Ephys.SSP_SAP.NREM.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.NREM.gammaBandPower.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[NREM] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
% xlim([1/30,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral gamma-band during REM
% subplot(2,3,3);
% % Naive
% for aa = 1:size(data.Naive.REM.gammaBandPower.C,2)
%     semilogx(data.Naive.REM.gammaBandPower.meanf,data.Naive.REM.gammaBandPower.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.IOS_Ephys.Blank_SAP.REM.gammaBandPower.C,2)
%     semilogx(data.IOS_Ephys.Blank_SAP.REM.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.REM.gammaBandPower.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.IOS_Ephys.SSP_SAP.REM.gammaBandPower.C,2)
%     semilogx(data.IOS_Ephys.SSP_SAP.REM.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.REM.gammaBandPower.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[REM] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
% xlim([1/60,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral gamma-band during Awake
% subplot(2,3,4);
% % Naive
% for aa = 1:size(data.Naive.Awake.gammaBandPower.C,2)
%     semilogx(data.Naive.Awake.gammaBandPower.meanf,data.Naive.Awake.gammaBandPower.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C,2)
%     semilogx(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.C,2)
%     semilogx(data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Alert] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
% xlim([0.003,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral gamma-band during Sleep
% subplot(2,3,5);
% % Naive
% for aa = 1:size(data.Naive.Sleep.gammaBandPower.C,2)
%     semilogx(data.Naive.Sleep.gammaBandPower.meanf,data.Naive.Sleep.gammaBandPower.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C,2)
%     semilogx(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.C,2)
%     semilogx(data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Asleep] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
% xlim([0.003,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral gamma-band during All data
% subplot(2,3,6);
% % Naive
% for aa = 1:size(data.Naive.All.gammaBandPower.C,2)
%     semilogx(data.Naive.All.gammaBandPower.meanf,data.Naive.All.gammaBandPower.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C,2)
%     semilogx(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.meanf,data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.IOS_Ephys.SSP_SAP.All.gammaBandPower.C,2)
%     semilogx(data.IOS_Ephys.SSP_SAP.All.gammaBandPower.meanf,data.IOS_Ephys.SSP_SAP.All.gammaBandPower.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[All] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
% xlim([0.003,0.5])
% ylim([0.01,1])
% set(gca,'box','off')
% %% save figure(s)
% if strcmp(saveFigs,'y') == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Bilateral Coherence - Bilateral IOS' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure4,[dirpath 'IndividualBilateralCoherence_Gamma']);
%     set(summaryFigure4,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'IndividualBilateralCoherence_Gamma'])
% end
% %% HbT and gamma-band coherence stats
% summaryFigure5 = figure;
% sgtitle('Bilateral Coherence Statics')
% %% Alert HbT Stats
% ax1 = subplot(2,3,1);
% xInds = ones(1,length(animalIDs.Blank_SAP));
% s1 = scatter(xInds*1,data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e1 = errorbar(1,mean(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C001),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'black';
% e1.MarkerSize = 10;
% e1.CapSize = 10;
% s2 = scatter(xInds*2,data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e2 = errorbar(2,mean(data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.C001),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e2.Color = 'black';
% e2.MarkerSize = 10;
% e2.CapSize = 10;
% scatter(xInds*3,data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e3 = errorbar(3,mean(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C01),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e3.Color = 'black';
% e3.MarkerSize = 10;
% e3.CapSize = 10;
% scatter(xInds*4,data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e4 = errorbar(4,mean(data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.C01),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'black';
% e4.MarkerSize = 10;
% e4.CapSize = 10;
% scatter(xInds*5,data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e5 = errorbar(5,mean(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C05),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e5.Color = 'black';
% e5.MarkerSize = 10;
% e5.CapSize = 10;
% scatter(xInds*6,data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e6 = errorbar(6,mean(data.IOS_Ephys.SSP_SAP.CBV_HbT.Alert.C05),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.Alert.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e6.Color = 'black';
% e6.MarkerSize = 10;
% e6.CapSize = 10;
% % stat lines
% plot([1,2],[1,1],'k');
% text(1.5,1,'ns','FontSize',16)
% plot([3,4],[1,1],'k');
% text(3.5,1,'ns','FontSize',16)
% plot([5,6],[1,1],'k');
% text(5.5,1,'*','FontSize',16)
% title({'Alert HbT',''})
% ylabel('Average coherence^2')
% legend([s1,s2],'Blank-SAP','SSP-SAP')
% set(gca,'xtick',[1.5,3.5,5.5])
% xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
% axis square
% axis tight
% xlim([0,7])
% ylim([0.01,1])
% set(gca,'box','off')
% ax1.TickLength = [0.03,0.03];
% %% Asleep HbT stats
% ax2 = subplot(2,3,2);
% xInds2 = ones(1,length(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C001));
% scatter(xInds2*1,data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e1 = errorbar(1,mean(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C001),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'black';
% e1.MarkerSize = 10;
% e1.CapSize = 10;
% scatter(xInds2*2,data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e2 = errorbar(2,mean(data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.C001),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e2.Color = 'black';
% e2.MarkerSize = 10;
% e2.CapSize = 10;
% scatter(xInds2*3,data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e3 = errorbar(3,mean(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C01),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e3.Color = 'black';
% e3.MarkerSize = 10;
% e3.CapSize = 10;
% scatter(xInds2*4,data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e4 = errorbar(4,mean(data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.C01),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'black';
% e4.MarkerSize = 10;
% e4.CapSize = 10;
% scatter(xInds2*5,data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e5 = errorbar(5,mean(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C05),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e5.Color = 'black';
% e5.MarkerSize = 10;
% e5.CapSize = 10;
% scatter(xInds2*6,data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e6 = errorbar(6,mean(data.IOS_Ephys.SSP_SAP.CBV_HbT.Asleep.C05),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.Asleep.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e6.Color = 'black';
% e6.MarkerSize = 10;
% e6.CapSize = 10;
% % stat lines
% plot([1,2],[1,1],'k');
% text(1.5,1,'ns','FontSize',16)
% plot([3,4],[1,1],'k');
% text(3.5,1,'*','FontSize',16)
% plot([5,6],[1,1],'k');
% text(5.5,1,'*','FontSize',16)
% title({'Asleep HbT',''})
% ylabel('Average coherence^2')
% set(gca,'xtick',[1.5,3.5,5.5])
% xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
% axis square
% axis tight
% xlim([0,7])
% ylim([0.01,1])
% set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
% %% All HbT stats
% ax3 = subplot(2,3,3);
% xInds = ones(1,length(animalIDs.Blank_SAP));
% scatter(xInds*1,data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e1 = errorbar(1,mean(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C001),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'black';
% e1.MarkerSize = 10;
% e1.CapSize = 10;
% scatter(xInds*2,data.IOS_Ephys.SSP_SAP.CBV_HbT.All.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e2 = errorbar(2,mean(data.IOS_Ephys.SSP_SAP.CBV_HbT.All.C001),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e2.Color = 'black';
% e2.MarkerSize = 10;
% e2.CapSize = 10;
% scatter(xInds*3,data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e3 = errorbar(3,mean(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C01),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e3.Color = 'black';
% e3.MarkerSize = 10;
% e3.CapSize = 10;
% scatter(xInds*4,data.IOS_Ephys.SSP_SAP.CBV_HbT.All.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e4 = errorbar(4,mean(data.IOS_Ephys.SSP_SAP.CBV_HbT.All.C01),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'black';
% e4.MarkerSize = 10;
% e4.CapSize = 10;
% scatter(xInds*5,data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e5 = errorbar(5,mean(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C05),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e5.Color = 'black';
% e5.MarkerSize = 10;
% e5.CapSize = 10;
% scatter(xInds*6,data.IOS_Ephys.SSP_SAP.CBV_HbT.All.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e6 = errorbar(6,mean(data.IOS_Ephys.SSP_SAP.CBV_HbT.All.C05),std(data.IOS_Ephys.Blank_SAP.CBV_HbT.All.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e6.Color = 'black';
% e6.MarkerSize = 10;
% e6.CapSize = 10;
% % stat lines
% plot([1,2],[1,1],'k');
% text(1.5,1,'ns','FontSize',16)
% plot([3,4],[1,1],'k');
% text(3.5,1,'ns','FontSize',16)
% plot([5,6],[1,1],'k');
% text(5.5,1,'***','FontSize',16)
% title({'All HbT',''})
% ylabel('Average coherence^2')
% set(gca,'xtick',[1.5,3.5,5.5])
% xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
% axis square
% axis tight
% xlim([0,7])
% ylim([0.01,1])
% set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
% %% Alert Gamma-band power Stats
% ax4 = subplot(2,3,4);
% scatter(xInds*1,data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e1 = errorbar(1,mean(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C001),std(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'black';
% e1.MarkerSize = 10;
% e1.CapSize = 10;
% scatter(xInds*2,data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e2 = errorbar(2,mean(data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.C001),std(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e2.Color = 'black';
% e2.MarkerSize = 10;
% e2.CapSize = 10;
% scatter(xInds*3,data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e3 = errorbar(3,mean(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C01),std(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e3.Color = 'black';
% e3.MarkerSize = 10;
% e3.CapSize = 10;
% scatter(xInds*4,data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e4 = errorbar(4,mean(data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.C01),std(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'black';
% e4.MarkerSize = 10;
% e4.CapSize = 10;
% scatter(xInds*5,data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e5 = errorbar(5,mean(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C05),std(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e5.Color = 'black';
% e5.MarkerSize = 10;
% e5.CapSize = 10;
% scatter(xInds*6,data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e6 = errorbar(6,mean(data.IOS_Ephys.SSP_SAP.Awake.gammaBandPower.C05),std(data.IOS_Ephys.Blank_SAP.Awake.gammaBandPower.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e6.Color = 'black';
% e6.MarkerSize = 10;
% e6.CapSize = 10;
% % stat lines
% plot([1,2],[1,1],'k');
% text(1.5,1,'ns','FontSize',16)
% plot([3,4],[1,1],'k');
% text(3.5,1,'ns','FontSize',16)
% plot([5,6],[1,1],'k');
% text(5.5,1,'ns','FontSize',16)
% title({'Alert Gamma-band power',''})
% ylabel('Average coherence^2')
% set(gca,'xtick',[1.5,3.5,5.5])
% xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
% axis square
% axis tight
% xlim([0,7])
% ylim([0.01,1])
% set(gca,'box','off')
% ax4.TickLength = [0.03,0.03];
% %% Asleep Gamma-band power stats
% ax5 = subplot(2,3,5);
% scatter(xInds2*1,data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e1 = errorbar(1,mean(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C001),std(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'black';
% e1.MarkerSize = 10;
% e1.CapSize = 10;
% scatter(xInds2*2,data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e2 = errorbar(2,mean(data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.C001),std(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e2.Color = 'black';
% e2.MarkerSize = 10;
% e2.CapSize = 10;
% scatter(xInds2*3,data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e3 = errorbar(3,mean(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C01),std(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e3.Color = 'black';
% e3.MarkerSize = 10;
% e3.CapSize = 10;
% scatter(xInds2*4,data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e4 = errorbar(4,mean(data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.C01),std(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'black';
% e4.MarkerSize = 10;
% e4.CapSize = 10;
% scatter(xInds2*5,data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e5 = errorbar(5,mean(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C05),std(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e5.Color = 'black';
% e5.MarkerSize = 10;
% e5.CapSize = 10;
% scatter(xInds2*6,data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e6 = errorbar(6,mean(data.IOS_Ephys.SSP_SAP.Sleep.gammaBandPower.C05),std(data.IOS_Ephys.Blank_SAP.Sleep.gammaBandPower.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e6.Color = 'black';
% e6.MarkerSize = 10;
% e6.CapSize = 10;
% % stat lines
% plot([1,2],[1,1],'k');
% text(1.5,1,'ns','FontSize',16)
% plot([3,4],[1,1],'k');
% text(3.5,1,'*','FontSize',16)
% plot([5,6],[1,1],'k');
% text(5.5,1,'ns','FontSize',16)
% title({'Asleep Gamma-band power',''})
% ylabel('Average coherence^2')
% set(gca,'xtick',[1.5,3.5,5.5])
% xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
% axis square
% axis tight
% xlim([0,7])
% ylim([0.01,1])
% set(gca,'box','off')
% ax5.TickLength = [0.03,0.03];
% %% All Gamma-band power stats
% ax6 = subplot(2,3,6);
% scatter(xInds*1,data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e1 = errorbar(1,mean(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C001),std(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'black';
% e1.MarkerSize = 10;
% e1.CapSize = 10;
% scatter(xInds*2,data.IOS_Ephys.SSP_SAP.All.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e2 = errorbar(2,mean(data.IOS_Ephys.SSP_SAP.All.gammaBandPower.C001),std(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e2.Color = 'black';
% e2.MarkerSize = 10;
% e2.CapSize = 10;
% scatter(xInds*3,data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e3 = errorbar(3,mean(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C01),std(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e3.Color = 'black';
% e3.MarkerSize = 10;
% e3.CapSize = 10;
% scatter(xInds*4,data.IOS_Ephys.SSP_SAP.All.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e4 = errorbar(4,mean(data.IOS_Ephys.SSP_SAP.All.gammaBandPower.C01),std(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'black';
% e4.MarkerSize = 10;
% e4.CapSize = 10;
% scatter(xInds*5,data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
% hold on
% e5 = errorbar(5,mean(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C05),std(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e5.Color = 'black';
% e5.MarkerSize = 10;
% e5.CapSize = 10;
% scatter(xInds*6,data.IOS_Ephys.SSP_SAP.All.gammaBandPower.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
% hold on
% e6 = errorbar(6,mean(data.IOS_Ephys.SSP_SAP.All.gammaBandPower.C05),std(data.IOS_Ephys.Blank_SAP.All.gammaBandPower.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e6.Color = 'black';
% e6.MarkerSize = 10;
% e6.CapSize = 10;
% % stat lines
% plot([1,2],[1,1],'k');
% text(1.5,1,'ns','FontSize',16)
% plot([3,4],[1,1],'k');
% text(3.5,1,'ns','FontSize',16)
% plot([5,6],[1,1],'k');
% text(5.5,1,'ns','FontSize',16)
% title({'All Gamma-band power',''})
% ylabel('Average coherence^2')
% set(gca,'xtick',[1.5,3.5,5.5])
% xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
% axis square
% axis tight
% xlim([0,7])
% ylim([0.01,1])
% set(gca,'box','off')
% ax6.TickLength = [0.03,0.03];
% %% save figure(s)
% if saveFigs == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Bilateral Coherence - Bilateral IOS' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure5,[dirpath 'AverageBilateralCoherence_Statistics']);
%     set(summaryFigure5,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'AverageBilateralCoherence_Statistics'])
%     %% statistical diary
%     diaryFile = [dirpath 'AverageBilateralCoherence_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     % Awake stats
%     disp('======================================================================================================================')
%     disp('GLME statistics for gamma Coherence^2 for Awake data')
%     disp('======================================================================================================================')
%     disp('0 -> 0.01 Hz')
%     disp(Stats.gammaBandPower.Awake.C001.Stats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp('0.01 -> 0.1 Hz')
%     disp(Stats.gammaBandPower.Awake.C01.Stats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp('0.1 -> 0.5 Hz')
%     disp(Stats.gammaBandPower.Awake.C05.Stats)
%     % Sleep stats
%     disp('======================================================================================================================')
%     disp('GLME statistics for gamma Coherence^2 for Sleep data')
%     disp('======================================================================================================================')
%     disp('0 -> 0.01 Hz')
%     disp(Stats.gammaBandPower.Sleep.C001.Stats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp('0.01 -> 0.1 Hz')
%     disp(Stats.gammaBandPower.Sleep.C01.Stats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp('0.1 -> 0.5 Hz')
%     disp(Stats.gammaBandPower.Sleep.C05.Stats)
%     % All stats
%     disp('======================================================================================================================')
%     disp('GLME statistics for gamma Coherence^2 for All data')
%     disp('======================================================================================================================')
%     disp('0 -> 0.01 Hz')
%     disp(Stats.gammaBandPower.All.C001.Stats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp('0.01 -> 0.1 Hz')
%     disp(Stats.gammaBandPower.All.C01.Stats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp('0.1 -> 0.5 Hz')
%     disp(Stats.gammaBandPower.All.C05.Stats)
%     % Awake stats
%     disp('======================================================================================================================')
%     disp('GLME statistics for CBV_HbT Coherence^2 for Awake data')
%     disp('======================================================================================================================')
%     disp('0 -> 0.01 Hz')
%     disp(Stats.CBV_HbT.Awake.C001.Stats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp('0.01 -> 0.1 Hz')
%     disp(Stats.CBV_HbT.Awake.C01.Stats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp('0.1 -> 0.5 Hz')
%     disp(Stats.CBV_HbT.Awake.C05.Stats)
%     % Sleep stats
%     disp('======================================================================================================================')
%     disp('GLME statistics for CBV_HbT Coherence^2 for Sleep data')
%     disp('======================================================================================================================')
%     disp('0 -> 0.01 Hz')
%     disp(Stats.CBV_HbT.Sleep.C001.Stats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp('0.01 -> 0.1 Hz')
%     disp(Stats.CBV_HbT.Sleep.C01.Stats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp('0.1 -> 0.5 Hz')
%     disp(Stats.CBV_HbT.Sleep.C05.Stats)
%     % All stats
%     disp('======================================================================================================================')
%     disp('GLME statistics for CBV_HbT Coherence^2 for All data')
%     disp('======================================================================================================================')
%     disp('0 -> 0.01 Hz')
%     disp(Stats.CBV_HbT.All.C001.Stats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp('0.01 -> 0.1 Hz')
%     disp(Stats.CBV_HbT.All.C01.Stats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp('0.1 -> 0.5 Hz')
%     disp(Stats.CBV_HbT.All.C05.Stats)
%     diary off
% end
% 
% end
