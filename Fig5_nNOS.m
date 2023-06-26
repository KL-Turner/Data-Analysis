function [] = Fig5_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
%% Ephys bilateral coherence
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_BilatCoher_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
dataTypes = {'HbT','gammaBandPower','deltaBandPower'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'C','f'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_BilatCoher_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        if any(strcmp(animalID,{'T142','T172'})) == false
            for cc = 1:length(dataTypes)
                dataType = dataTypes{1,cc};
                ephysBilatCoherData.(group).(dataType).dummCheck = 1;
                for dd = 1:length(behaviors)
                    behavior = behaviors{1,dd};
                    if isfield(ephysBilatCoherData.(group).(dataType),behavior) == false
                        ephysBilatCoherData.(group).(dataType).(behavior).C = [];
                        ephysBilatCoherData.(group).(dataType).(behavior).f = [];
                        ephysBilatCoherData.(group).(dataType).(behavior).group = {};
                        ephysBilatCoherData.(group).(dataType).(behavior).animalID = {};
                    end
                    if isempty(Results_BilatCoher_Ephys.(group).(animalID).(dataType).(behavior).C) == false
                        ephysBilatCoherData.(group).(dataType).(behavior).C = cat(1,ephysBilatCoherData.(group).(dataType).(behavior).C,Results_BilatCoher_Ephys.(group).(animalID).(dataType).(behavior).C');
                        ephysBilatCoherData.(group).(dataType).(behavior).f = cat(1,ephysBilatCoherData.(group).(dataType).(behavior).f,Results_BilatCoher_Ephys.(group).(animalID).(dataType).(behavior).f);
                        ephysBilatCoherData.(group).(dataType).(behavior).group = cat(1,ephysBilatCoherData.(group).(dataType).(behavior).group,group);
                        ephysBilatCoherData.(group).(dataType).(behavior).animalID = cat(1,ephysBilatCoherData.(group).(dataType).(behavior).animalID,animalID);
                    end
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
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                ephysBilatCoherData.(group).(dataType).(behavior).(['mean_' variable]) = mean(ephysBilatCoherData.(group).(dataType).(behavior).(variable),1);
                ephysBilatCoherData.(group).(dataType).(behavior).(['stdErr_' variable]) = std(ephysBilatCoherData.(group).(dataType).(behavior).(variable),0,1)./sqrt(size(ephysBilatCoherData.(group).(dataType).(behavior).(variable),1));
            end
        end
    end
end
%% GCaMP bilateral coherence
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_BilatCoher_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'C','f'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_BilatCoher_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            gcampBilatCoherData.(group).(dataType).dummCheck = 1;
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                if isfield(gcampBilatCoherData.(group).(dataType),behavior) == false
                    gcampBilatCoherData.(group).(dataType).(behavior).C = [];
                    gcampBilatCoherData.(group).(dataType).(behavior).f = [];
                    gcampBilatCoherData.(group).(dataType).(behavior).group = {};
                    gcampBilatCoherData.(group).(dataType).(behavior).animalID = {};
                end
                if isempty(Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).C) == false
                    gcampBilatCoherData.(group).(dataType).(behavior).C = cat(1,gcampBilatCoherData.(group).(dataType).(behavior).C,Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).C.^2');
                    gcampBilatCoherData.(group).(dataType).(behavior).f = cat(1,gcampBilatCoherData.(group).(dataType).(behavior).f,Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).f);
                    gcampBilatCoherData.(group).(dataType).(behavior).group = cat(1,gcampBilatCoherData.(group).(dataType).(behavior).group,group);
                    gcampBilatCoherData.(group).(dataType).(behavior).animalID = cat(1,gcampBilatCoherData.(group).(dataType).(behavior).animalID,animalID);
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
                gcampBilatCoherData.(group).(dataType).(behavior).(['mean_' variable]) = mean(gcampBilatCoherData.(group).(dataType).(behavior).(variable),1);
                gcampBilatCoherData.(group).(dataType).(behavior).(['stdErr_' variable]) = std(gcampBilatCoherData.(group).(dataType).(behavior).(variable),0,1)./sqrt(size(gcampBilatCoherData.(group).(dataType).(behavior).(variable),1));
            end
        end
    end
end
%% figure
figure;
% Ephys bilateral HbT - Rest
subplot(2,3,1);
p1 = semilogx(ephysBilatCoherData.Blank_SAP.HbT.Rest.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Rest.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Rest.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Rest.mean_C + ephysBilatCoherData.Blank_SAP.HbT.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Rest.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Rest.mean_C - ephysBilatCoherData.Blank_SAP.HbT.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
p2 = semilogx(ephysBilatCoherData.SSP_SAP.HbT.Rest.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Rest.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Rest.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Rest.mean_C + ephysBilatCoherData.SSP_SAP.HbT.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Rest.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Rest.mean_C - ephysBilatCoherData.SSP_SAP.HbT.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbT Rest')
xlabel('Freq (Hz)')
ylabel('Coherence')
legend([p1,p2],'Blank-SAP','SSP-SAP')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% Ephys bilateral HbT - Alert
subplot(2,3,2);
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Alert.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Alert.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Alert.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Alert.mean_C + ephysBilatCoherData.Blank_SAP.HbT.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Alert.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Alert.mean_C - ephysBilatCoherData.Blank_SAP.HbT.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Alert.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Alert.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Alert.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Alert.mean_C + ephysBilatCoherData.SSP_SAP.HbT.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Alert.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Alert.mean_C - ephysBilatCoherData.SSP_SAP.HbT.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbT Alert')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% Ephys bilateral HbT - Asleep
subplot(2,3,3);
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Asleep.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Asleep.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Asleep.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Asleep.mean_C + ephysBilatCoherData.Blank_SAP.HbT.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Asleep.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Asleep.mean_C - ephysBilatCoherData.Blank_SAP.HbT.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Asleep.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Asleep.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Asleep.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Asleep.mean_C + ephysBilatCoherData.SSP_SAP.HbT.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Asleep.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Asleep.mean_C - ephysBilatCoherData.SSP_SAP.HbT.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbT Asleep')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% Ephys bilateral gammaBandPower - Rest
subplot(2,3,4);
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.mean_C + ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.mean_C - ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.mean_C + ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.mean_C - ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat gammaBandPower Rest')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% Ephys bilateral gammaBandPower - Alert
subplot(2,3,5);
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.mean_C + ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.mean_C - ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.mean_C + ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.mean_C - ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat gammaBandPower Alert')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% Ephys bilateral gammaBandPower - Asleep
subplot(2,3,6);
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.mean_C + ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.mean_C - ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.mean_C + ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.mean_C - ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat gammaBandPower Asleep')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
%% figure
figure
% GCaMP bilateral HbT - Rest
subplot(4,3,1);
p1 = semilogx(gcampBilatCoherData.Blank_SAP.HbT.Rest.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Rest.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Rest.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Rest.mean_C + gcampBilatCoherData.Blank_SAP.HbT.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Rest.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Rest.mean_C - gcampBilatCoherData.Blank_SAP.HbT.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
p2 = semilogx(gcampBilatCoherData.SSP_SAP.HbT.Rest.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Rest.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Rest.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Rest.mean_C + gcampBilatCoherData.SSP_SAP.HbT.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Rest.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Rest.mean_C - gcampBilatCoherData.SSP_SAP.HbT.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbT Rest')
xlabel('Freq (Hz)')
ylabel('Coherence')
legend([p1,p2],'Blank-SAP','SSP-SAP')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% GCaMP bilateral HbT - Alert
subplot(4,3,2);
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Alert.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Alert.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Alert.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Alert.mean_C + gcampBilatCoherData.Blank_SAP.HbT.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Alert.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Alert.mean_C - gcampBilatCoherData.Blank_SAP.HbT.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Alert.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Alert.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Alert.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Alert.mean_C + gcampBilatCoherData.SSP_SAP.HbT.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Alert.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Alert.mean_C - gcampBilatCoherData.SSP_SAP.HbT.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbT Alert')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% GCaMP bilateral HbT - Asleep
subplot(4,3,3);
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Asleep.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Asleep.mean_C + gcampBilatCoherData.Blank_SAP.HbT.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Asleep.mean_C - gcampBilatCoherData.Blank_SAP.HbT.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Asleep.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Asleep.mean_C + gcampBilatCoherData.SSP_SAP.HbT.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Asleep.mean_C - gcampBilatCoherData.SSP_SAP.HbT.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbT Asleep')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% GCaMP bilateral GCaMP - Rest
subplot(4,3,4);
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Rest.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Rest.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Rest.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Rest.mean_C + gcampBilatCoherData.Blank_SAP.GCaMP.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Rest.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Rest.mean_C - gcampBilatCoherData.Blank_SAP.GCaMP.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Rest.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Rest.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Rest.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Rest.mean_C + gcampBilatCoherData.SSP_SAP.GCaMP.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Rest.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Rest.mean_C - gcampBilatCoherData.SSP_SAP.GCaMP.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat GCaMP Rest')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% GCaMP bilateral GCaMP - Alert
subplot(4,3,5);
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Alert.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Alert.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Alert.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Alert.mean_C + gcampBilatCoherData.Blank_SAP.GCaMP.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Alert.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Alert.mean_C - gcampBilatCoherData.Blank_SAP.GCaMP.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Alert.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Alert.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Alert.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Alert.mean_C + gcampBilatCoherData.SSP_SAP.GCaMP.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Alert.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Alert.mean_C - gcampBilatCoherData.SSP_SAP.GCaMP.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat GCaMP Alert')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% GCaMP bilateral GCaMP - Asleep
subplot(4,3,6);
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.mean_C + gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.mean_C - gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.mean_C + gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.mean_C - gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat GCaMP Asleep')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% HbR bilateral HbO - Rest
subplot(4,3,7);
p1 = semilogx(gcampBilatCoherData.Blank_SAP.HbO.Rest.mean_f,gcampBilatCoherData.Blank_SAP.HbO.Rest.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.HbO.Rest.mean_f,gcampBilatCoherData.Blank_SAP.HbO.Rest.mean_C + gcampBilatCoherData.Blank_SAP.HbO.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.HbO.Rest.mean_f,gcampBilatCoherData.Blank_SAP.HbO.Rest.mean_C - gcampBilatCoherData.Blank_SAP.HbO.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
p2 = semilogx(gcampBilatCoherData.SSP_SAP.HbO.Rest.mean_f,gcampBilatCoherData.SSP_SAP.HbO.Rest.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.HbO.Rest.mean_f,gcampBilatCoherData.SSP_SAP.HbO.Rest.mean_C + gcampBilatCoherData.SSP_SAP.HbO.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbO.Rest.mean_f,gcampBilatCoherData.SSP_SAP.HbO.Rest.mean_C - gcampBilatCoherData.SSP_SAP.HbO.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbO Rest')
xlabel('Freq (Hz)')
ylabel('Coherence')
legend([p1,p2],'Blank-SAP','SSP-SAP')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% HbR bilateral HbO - Alert
subplot(4,3,8);
semilogx(gcampBilatCoherData.Blank_SAP.HbO.Alert.mean_f,gcampBilatCoherData.Blank_SAP.HbO.Alert.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.HbO.Alert.mean_f,gcampBilatCoherData.Blank_SAP.HbO.Alert.mean_C + gcampBilatCoherData.Blank_SAP.HbO.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.HbO.Alert.mean_f,gcampBilatCoherData.Blank_SAP.HbO.Alert.mean_C - gcampBilatCoherData.Blank_SAP.HbO.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbO.Alert.mean_f,gcampBilatCoherData.SSP_SAP.HbO.Alert.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.HbO.Alert.mean_f,gcampBilatCoherData.SSP_SAP.HbO.Alert.mean_C + gcampBilatCoherData.SSP_SAP.HbO.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbO.Alert.mean_f,gcampBilatCoherData.SSP_SAP.HbO.Alert.mean_C - gcampBilatCoherData.SSP_SAP.HbO.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbO Alert')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% HbR bilateral HbO - Asleep
subplot(4,3,9);
semilogx(gcampBilatCoherData.Blank_SAP.HbO.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.HbO.Asleep.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.HbO.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.HbO.Asleep.mean_C + gcampBilatCoherData.Blank_SAP.HbO.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.HbO.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.HbO.Asleep.mean_C - gcampBilatCoherData.Blank_SAP.HbO.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbO.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.HbO.Asleep.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.HbO.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.HbO.Asleep.mean_C + gcampBilatCoherData.SSP_SAP.HbO.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbO.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.HbO.Asleep.mean_C - gcampBilatCoherData.SSP_SAP.HbO.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbO Asleep')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% HbR bilateral HbR - Rest
subplot(4,3,10);
semilogx(gcampBilatCoherData.Blank_SAP.HbR.Rest.mean_f,gcampBilatCoherData.Blank_SAP.HbR.Rest.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.HbR.Rest.mean_f,gcampBilatCoherData.Blank_SAP.HbR.Rest.mean_C + gcampBilatCoherData.Blank_SAP.HbR.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.HbR.Rest.mean_f,gcampBilatCoherData.Blank_SAP.HbR.Rest.mean_C - gcampBilatCoherData.Blank_SAP.HbR.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbR.Rest.mean_f,gcampBilatCoherData.SSP_SAP.HbR.Rest.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.HbR.Rest.mean_f,gcampBilatCoherData.SSP_SAP.HbR.Rest.mean_C + gcampBilatCoherData.SSP_SAP.HbR.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbR.Rest.mean_f,gcampBilatCoherData.SSP_SAP.HbR.Rest.mean_C - gcampBilatCoherData.SSP_SAP.HbR.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbR Rest')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% HbR bilateral HbR - Alert
subplot(4,3,11);
semilogx(gcampBilatCoherData.Blank_SAP.HbR.Alert.mean_f,gcampBilatCoherData.Blank_SAP.HbR.Alert.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.HbR.Alert.mean_f,gcampBilatCoherData.Blank_SAP.HbR.Alert.mean_C + gcampBilatCoherData.Blank_SAP.HbR.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.HbR.Alert.mean_f,gcampBilatCoherData.Blank_SAP.HbR.Alert.mean_C - gcampBilatCoherData.Blank_SAP.HbR.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbR.Alert.mean_f,gcampBilatCoherData.SSP_SAP.HbR.Alert.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.HbR.Alert.mean_f,gcampBilatCoherData.SSP_SAP.HbR.Alert.mean_C + gcampBilatCoherData.SSP_SAP.HbR.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbR.Alert.mean_f,gcampBilatCoherData.SSP_SAP.HbR.Alert.mean_C - gcampBilatCoherData.SSP_SAP.HbR.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbR Alert')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% HbR bilateral HbR - Asleep
subplot(4,3,12);
semilogx(gcampBilatCoherData.Blank_SAP.HbR.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.HbR.Asleep.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.HbR.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.HbR.Asleep.mean_C + gcampBilatCoherData.Blank_SAP.HbR.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.HbR.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.HbR.Asleep.mean_C - gcampBilatCoherData.Blank_SAP.HbR.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbR.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.HbR.Asleep.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.HbR.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.HbR.Asleep.mean_C + gcampBilatCoherData.SSP_SAP.HbR.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbR.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.HbR.Asleep.mean_C - gcampBilatCoherData.SSP_SAP.HbR.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbR Asleep')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
