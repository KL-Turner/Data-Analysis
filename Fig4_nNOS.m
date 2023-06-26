function [] = Fig4_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
%% Ephys cross correlation
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_CrossCorr_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Naive','Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
dataTypes = {'gammaBandPower','muaPower'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
samplingRate = 30;
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_CrossCorr_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        if any(strcmp(animalID,{'T142','T172'})) == false
            for cc = 1:length(hemispheres)
                hemisphere = hemispheres{1,cc};
                for dd = 1:length(dataTypes)
                    dataType = dataTypes{1,dd};
                    for ee = 1:length(behaviors)
                        behavior = behaviors{1,ee};
                        ephysXCorrData.(group).(hemisphere).(dataType).dummyCheck = 1;
                        if isfield(ephysXCorrData.(group).(hemisphere).(dataType),(behavior)) == false
                            ephysXCorrData.(group).(hemisphere).(dataType).(behavior).lags = [];
                            ephysXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals = [];
                        end
                        if isempty(Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals) == false
                            ephysXCorrData.(group).(hemisphere).(dataType).(behavior).lags = cat(1,ephysXCorrData.(group).(hemisphere).(dataType).(behavior).lags,Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).lags/samplingRate);
                            ephysXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals = cat(1,ephysXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals,Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals);
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
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                ephysXCorrData.(group).(hemisphere).(dataType).(behavior).mean_xcVals = mean(ephysXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals,1);
                ephysXCorrData.(group).(hemisphere).(dataType).(behavior).stdErr_xcVals = std(ephysXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals,0,1)./sqrt(size(ephysXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals,1));
                ephysXCorrData.(group).(hemisphere).(dataType).(behavior).mean_lags = mean(ephysXCorrData.(group).(hemisphere).(dataType).(behavior).lags,1);
            end
        end
    end
end
%% Ephys neural-hemo coherence
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_NeuralHemoCoher_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
dataTypes = {'gammaBandPower','deltaBandPower'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'C','f'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_NeuralHemoCoher_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        if any(strcmp(animalID,{'T142','T172'})) == false
            for cc = 1:length(hemispheres)
                hemisphere = hemispheres{1,cc};
                for dd = 1:length(dataTypes)
                    dataType = dataTypes{1,dd};
                    ephysNHCoherData.(group).(hemisphere).(dataType).dummCheck = 1;
                    for ee = 1:length(behaviors)
                        behavior = behaviors{1,ee};
                        if isfield(ephysNHCoherData.(group).(hemisphere).(dataType),behavior) == false
                            ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).C = [];
                            ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).f = [];
                            ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).group = {};
                            ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).animalID = {};
                        end
                        if isempty(Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).C) == false
                            ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).C = cat(1,ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).C,Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).C');
                            ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).f = cat(1,ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).f,Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).f);
                            ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).group = cat(1,ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).group,group);
                            ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).(['stdErr_' variable]) = std(ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable),0,1)./sqrt(size(ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable),1));
                end
            end
        end
    end
end
%% GCaMP cross correlation
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_CrossCorr_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
dataTypes = {'HbT','HbO','HbR'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
samplingRate = 10;
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_CrossCorr_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    gcampXCorrData.(group).(hemisphere).(dataType).dummyCheck = 1;
                    if isfield(gcampXCorrData.(group).(hemisphere).(dataType),(behavior)) == false
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).lags = [];
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals = [];
                    end
                    if isempty(Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals) == false
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).lags = cat(1,gcampXCorrData.(group).(hemisphere).(dataType).(behavior).lags,Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).lags/samplingRate);
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals = cat(1,gcampXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals,Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals);
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
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                gcampXCorrData.(group).(hemisphere).(dataType).(behavior).mean_xcVals = mean(gcampXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals,1);
                gcampXCorrData.(group).(hemisphere).(dataType).(behavior).stdErr_xcVals = std(gcampXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals,0,1)./sqrt(size(gcampXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals,1));
                gcampXCorrData.(group).(hemisphere).(dataType).(behavior).mean_lags = mean(gcampXCorrData.(group).(hemisphere).(dataType).(behavior).lags,1);
            end
        end
    end
end
%% GCaMP neural-hemo coherence
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_NeuralHemoCoher_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
dataTypes = {'HbT','HbO','HbR'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'C','f'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_NeuralHemoCoher_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                gcampNHCoherData.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(gcampNHCoherData.(group).(hemisphere).(dataType),behavior) == false
                        gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).C = [];
                        gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).f = [];
                        gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).group = {};
                        gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    if isempty(Results_NeuralHemoCoher_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).C) == false
                        gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).C = cat(1,gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).C,Results_NeuralHemoCoher_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).C');
                        gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).f = cat(1,gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).f,Results_NeuralHemoCoher_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).f);
                        gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).group = cat(1,gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).group,group);
                        gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).(['stdErr_' variable]) = std(gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable),0,1)./sqrt(size(gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable),1));
                end
            end
        end
    end
end
%% figure panel 4
figure;
% Ephys cross correlation - Rest
subplot(4,3,1);
p1 = plot(ephysXCorrData.Blank_SAP.RH.gammaBandPower.Rest.mean_lags,ephysXCorrData.Blank_SAP.RH.gammaBandPower.Rest.mean_xcVals,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(ephysXCorrData.Blank_SAP.RH.gammaBandPower.Rest.mean_lags,ephysXCorrData.Blank_SAP.RH.gammaBandPower.Rest.mean_xcVals + ephysXCorrData.Blank_SAP.RH.gammaBandPower.Rest.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(ephysXCorrData.Blank_SAP.RH.gammaBandPower.Rest.mean_lags,ephysXCorrData.Blank_SAP.RH.gammaBandPower.Rest.mean_xcVals - ephysXCorrData.Blank_SAP.RH.gammaBandPower.Rest.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
p2 = plot(ephysXCorrData.SSP_SAP.RH.gammaBandPower.Rest.mean_lags,ephysXCorrData.SSP_SAP.RH.gammaBandPower.Rest.mean_xcVals,'color',colors('electric purple'),'LineWidth',2);
plot(ephysXCorrData.SSP_SAP.RH.gammaBandPower.Rest.mean_lags,ephysXCorrData.SSP_SAP.RH.gammaBandPower.Rest.mean_xcVals + ephysXCorrData.SSP_SAP.RH.gammaBandPower.Rest.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
plot(ephysXCorrData.SSP_SAP.RH.gammaBandPower.Rest.mean_lags,ephysXCorrData.SSP_SAP.RH.gammaBandPower.Rest.mean_xcVals - ephysXCorrData.SSP_SAP.RH.gammaBandPower.Rest.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
title('Gamma-HbT Rest')
ylabel('Corr Coef')
xlabel('Lags (s)')
xlim([-5,5])
legend([p1,p2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
axis square
% Ephys cross correlation - Alert
subplot(4,3,2)
plot(ephysXCorrData.Blank_SAP.RH.gammaBandPower.Alert.mean_lags,ephysXCorrData.Blank_SAP.RH.gammaBandPower.Alert.mean_xcVals,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(ephysXCorrData.Blank_SAP.RH.gammaBandPower.Alert.mean_lags,ephysXCorrData.Blank_SAP.RH.gammaBandPower.Alert.mean_xcVals + ephysXCorrData.Blank_SAP.RH.gammaBandPower.Alert.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(ephysXCorrData.Blank_SAP.RH.gammaBandPower.Alert.mean_lags,ephysXCorrData.Blank_SAP.RH.gammaBandPower.Alert.mean_xcVals - ephysXCorrData.Blank_SAP.RH.gammaBandPower.Alert.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(ephysXCorrData.SSP_SAP.RH.gammaBandPower.Alert.mean_lags,ephysXCorrData.SSP_SAP.RH.gammaBandPower.Alert.mean_xcVals,'color',colors('electric purple'),'LineWidth',2);
plot(ephysXCorrData.SSP_SAP.RH.gammaBandPower.Alert.mean_lags,ephysXCorrData.SSP_SAP.RH.gammaBandPower.Alert.mean_xcVals + ephysXCorrData.SSP_SAP.RH.gammaBandPower.Alert.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
plot(ephysXCorrData.SSP_SAP.RH.gammaBandPower.Alert.mean_lags,ephysXCorrData.SSP_SAP.RH.gammaBandPower.Alert.mean_xcVals - ephysXCorrData.SSP_SAP.RH.gammaBandPower.Alert.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
title('Gamma-HbT Alert')
ylabel('Corr Coef')
xlabel('Lags (s)')
xlim([-5,5])
legend([p1,p2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
axis square
% Ephys cross correlation - Asleep
subplot(4,3,3)
plot(ephysXCorrData.Blank_SAP.RH.gammaBandPower.Asleep.mean_lags,ephysXCorrData.Blank_SAP.RH.gammaBandPower.Asleep.mean_xcVals,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(ephysXCorrData.Blank_SAP.RH.gammaBandPower.Asleep.mean_lags,ephysXCorrData.Blank_SAP.RH.gammaBandPower.Asleep.mean_xcVals + ephysXCorrData.Blank_SAP.RH.gammaBandPower.Asleep.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(ephysXCorrData.Blank_SAP.RH.gammaBandPower.Asleep.mean_lags,ephysXCorrData.Blank_SAP.RH.gammaBandPower.Asleep.mean_xcVals - ephysXCorrData.Blank_SAP.RH.gammaBandPower.Asleep.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(ephysXCorrData.SSP_SAP.RH.gammaBandPower.Asleep.mean_lags,ephysXCorrData.SSP_SAP.RH.gammaBandPower.Asleep.mean_xcVals,'color',colors('electric purple'),'LineWidth',2);
plot(ephysXCorrData.SSP_SAP.RH.gammaBandPower.Asleep.mean_lags,ephysXCorrData.SSP_SAP.RH.gammaBandPower.Asleep.mean_xcVals + ephysXCorrData.SSP_SAP.RH.gammaBandPower.Asleep.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
plot(ephysXCorrData.SSP_SAP.RH.gammaBandPower.Asleep.mean_lags,ephysXCorrData.SSP_SAP.RH.gammaBandPower.Asleep.mean_xcVals - ephysXCorrData.SSP_SAP.RH.gammaBandPower.Asleep.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
title('Gamma-HbT Asleep')
ylabel('Corr Coef')
xlabel('Lags (s)')
xlim([-5,5])
set(gca,'box','off')
axis square
% Ephys neural-hemo coherence - Rest
subplot(4,3,4);
semilogx(ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Rest.mean_f,ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Rest.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Rest.mean_f,ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Rest.mean_C + ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Rest.mean_f,ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Rest.mean_C - ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Rest.mean_f,ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Rest.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Rest.mean_f,ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Rest.mean_C + ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Rest.mean_f,ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Rest.mean_C - ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('Gamma-HbT Rest')
ylabel('Coherence')
xlabel('Freq (Hz)')
xlim([0.1,0.35])
set(gca,'box','off')
axis square
% Ephys neural-hemo coherence - Alert
subplot(4,3,5);
semilogx(ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Alert.mean_f,ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Alert.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Alert.mean_f,ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Alert.mean_C + ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Alert.mean_f,ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Alert.mean_C - ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Alert.mean_f,ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Alert.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Alert.mean_f,ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Alert.mean_C + ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Alert.mean_f,ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Alert.mean_C - ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('Gamma-HbT Alert')
ylabel('Coherence')
xlabel('Freq (Hz)')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% Ephys neural-hemo coherence - Asleep
subplot(4,3,6);
semilogx(ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Asleep.mean_f,ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Asleep.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Asleep.mean_f,ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Asleep.mean_C + ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Asleep.mean_f,ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Asleep.mean_C - ephysNHCoherData.Blank_SAP.RH.gammaBandPower.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Asleep.mean_f,ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Asleep.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Asleep.mean_f,ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Asleep.mean_C + ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Asleep.mean_f,ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Asleep.mean_C - ephysNHCoherData.SSP_SAP.RH.gammaBandPower.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('Gamma-HbT Asleep')
ylabel('Coherence')
xlabel('Freq (Hz)')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% GCaMP cross correlation - Rest
subplot(4,3,7);
plot(gcampXCorrData.Blank_SAP.RH.HbT.Rest.mean_lags,gcampXCorrData.Blank_SAP.RH.HbT.Rest.mean_xcVals,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(gcampXCorrData.Blank_SAP.RH.HbT.Rest.mean_lags,gcampXCorrData.Blank_SAP.RH.HbT.Rest.mean_xcVals + gcampXCorrData.Blank_SAP.RH.HbT.Rest.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(gcampXCorrData.Blank_SAP.RH.HbT.Rest.mean_lags,gcampXCorrData.Blank_SAP.RH.HbT.Rest.mean_xcVals - gcampXCorrData.Blank_SAP.RH.HbT.Rest.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(gcampXCorrData.SSP_SAP.RH.HbT.Rest.mean_lags,gcampXCorrData.SSP_SAP.RH.HbT.Rest.mean_xcVals,'color',colors('electric purple'),'LineWidth',2);
plot(gcampXCorrData.SSP_SAP.RH.HbT.Rest.mean_lags,gcampXCorrData.SSP_SAP.RH.HbT.Rest.mean_xcVals + gcampXCorrData.SSP_SAP.RH.HbT.Rest.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
plot(gcampXCorrData.SSP_SAP.RH.HbT.Rest.mean_lags,gcampXCorrData.SSP_SAP.RH.HbT.Rest.mean_xcVals - gcampXCorrData.SSP_SAP.RH.HbT.Rest.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
title('GCaMP-HbT Rest')
ylabel('Corr Coef')
xlabel('Lags (s)')
xlim([-5,5])
set(gca,'box','off')
axis square
% GCaMP cross correlation - Alert
subplot(4,3,8)
plot(gcampXCorrData.Blank_SAP.RH.HbT.Alert.mean_lags,gcampXCorrData.Blank_SAP.RH.HbT.Alert.mean_xcVals,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(gcampXCorrData.Blank_SAP.RH.HbT.Alert.mean_lags,gcampXCorrData.Blank_SAP.RH.HbT.Alert.mean_xcVals + gcampXCorrData.Blank_SAP.RH.HbT.Alert.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(gcampXCorrData.Blank_SAP.RH.HbT.Alert.mean_lags,gcampXCorrData.Blank_SAP.RH.HbT.Alert.mean_xcVals - gcampXCorrData.Blank_SAP.RH.HbT.Alert.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(gcampXCorrData.SSP_SAP.RH.HbT.Alert.mean_lags,gcampXCorrData.SSP_SAP.RH.HbT.Alert.mean_xcVals,'color',colors('electric purple'),'LineWidth',2);
plot(gcampXCorrData.SSP_SAP.RH.HbT.Alert.mean_lags,gcampXCorrData.SSP_SAP.RH.HbT.Alert.mean_xcVals + gcampXCorrData.SSP_SAP.RH.HbT.Alert.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
plot(gcampXCorrData.SSP_SAP.RH.HbT.Alert.mean_lags,gcampXCorrData.SSP_SAP.RH.HbT.Alert.mean_xcVals - gcampXCorrData.SSP_SAP.RH.HbT.Alert.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
title('GCaMP-HbT Alert')
ylabel('Corr Coef')
xlabel('Lags (s)')
xlim([-5,5])
legend([p1,p2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
axis square
% GCaMP cross correlation - Asleep
subplot(4,3,9)
plot(gcampXCorrData.Blank_SAP.RH.HbT.Asleep.mean_lags,gcampXCorrData.Blank_SAP.RH.HbT.Asleep.mean_xcVals,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(gcampXCorrData.Blank_SAP.RH.HbT.Asleep.mean_lags,gcampXCorrData.Blank_SAP.RH.HbT.Asleep.mean_xcVals + gcampXCorrData.Blank_SAP.RH.HbT.Asleep.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(gcampXCorrData.Blank_SAP.RH.HbT.Asleep.mean_lags,gcampXCorrData.Blank_SAP.RH.HbT.Asleep.mean_xcVals - gcampXCorrData.Blank_SAP.RH.HbT.Asleep.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(gcampXCorrData.SSP_SAP.RH.HbT.Asleep.mean_lags,gcampXCorrData.SSP_SAP.RH.HbT.Asleep.mean_xcVals,'color',colors('electric purple'),'LineWidth',2);
plot(gcampXCorrData.SSP_SAP.RH.HbT.Asleep.mean_lags,gcampXCorrData.SSP_SAP.RH.HbT.Asleep.mean_xcVals + gcampXCorrData.SSP_SAP.RH.HbT.Asleep.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
plot(gcampXCorrData.SSP_SAP.RH.HbT.Asleep.mean_lags,gcampXCorrData.SSP_SAP.RH.HbT.Asleep.mean_xcVals - gcampXCorrData.SSP_SAP.RH.HbT.Asleep.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
title('GCaMP-HbT Asleep')
ylabel('Corr Coef')
xlabel('Lags (s)')
xlim([-5,5])
set(gca,'box','off')
axis square
% GCaMP neural-hemo coherence - Rest
subplot(4,3,10);
semilogx(gcampNHCoherData.Blank_SAP.RH.HbT.Rest.mean_f,gcampNHCoherData.Blank_SAP.RH.HbT.Rest.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampNHCoherData.Blank_SAP.RH.HbT.Rest.mean_f,gcampNHCoherData.Blank_SAP.RH.HbT.Rest.mean_C + gcampNHCoherData.Blank_SAP.RH.HbT.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampNHCoherData.Blank_SAP.RH.HbT.Rest.mean_f,gcampNHCoherData.Blank_SAP.RH.HbT.Rest.mean_C - gcampNHCoherData.Blank_SAP.RH.HbT.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampNHCoherData.SSP_SAP.RH.HbT.Rest.mean_f,gcampNHCoherData.SSP_SAP.RH.HbT.Rest.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampNHCoherData.SSP_SAP.RH.HbT.Rest.mean_f,gcampNHCoherData.SSP_SAP.RH.HbT.Rest.mean_C + gcampNHCoherData.SSP_SAP.RH.HbT.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(gcampNHCoherData.SSP_SAP.RH.HbT.Rest.mean_f,gcampNHCoherData.SSP_SAP.RH.HbT.Rest.mean_C - gcampNHCoherData.SSP_SAP.RH.HbT.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('GCaMP-HbT Rest')
ylabel('Coherence')
xlabel('Freq (Hz)')
xlim([0.1,0.35])
set(gca,'box','off')
axis square
% GCaMP neural-hemo coherence - Alert
subplot(4,3,11);
semilogx(gcampNHCoherData.Blank_SAP.RH.HbT.Alert.mean_f,gcampNHCoherData.Blank_SAP.RH.HbT.Alert.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampNHCoherData.Blank_SAP.RH.HbT.Alert.mean_f,gcampNHCoherData.Blank_SAP.RH.HbT.Alert.mean_C + gcampNHCoherData.Blank_SAP.RH.HbT.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampNHCoherData.Blank_SAP.RH.HbT.Alert.mean_f,gcampNHCoherData.Blank_SAP.RH.HbT.Alert.mean_C - gcampNHCoherData.Blank_SAP.RH.HbT.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampNHCoherData.SSP_SAP.RH.HbT.Alert.mean_f,gcampNHCoherData.SSP_SAP.RH.HbT.Alert.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampNHCoherData.SSP_SAP.RH.HbT.Alert.mean_f,gcampNHCoherData.SSP_SAP.RH.HbT.Alert.mean_C + gcampNHCoherData.SSP_SAP.RH.HbT.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(gcampNHCoherData.SSP_SAP.RH.HbT.Alert.mean_f,gcampNHCoherData.SSP_SAP.RH.HbT.Alert.mean_C - gcampNHCoherData.SSP_SAP.RH.HbT.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('GCaMP-HbT Alert')
ylabel('Coherence')
xlabel('Freq (Hz)')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
% GCaMP neural-hemo coherence - Asleep
subplot(4,3,12);
semilogx(gcampNHCoherData.Blank_SAP.RH.HbT.Asleep.mean_f,gcampNHCoherData.Blank_SAP.RH.HbT.Asleep.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on
semilogx(gcampNHCoherData.Blank_SAP.RH.HbT.Asleep.mean_f,gcampNHCoherData.Blank_SAP.RH.HbT.Asleep.mean_C + gcampNHCoherData.Blank_SAP.RH.HbT.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampNHCoherData.Blank_SAP.RH.HbT.Asleep.mean_f,gcampNHCoherData.Blank_SAP.RH.HbT.Asleep.mean_C - gcampNHCoherData.Blank_SAP.RH.HbT.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampNHCoherData.SSP_SAP.RH.HbT.Asleep.mean_f,gcampNHCoherData.SSP_SAP.RH.HbT.Asleep.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(gcampNHCoherData.SSP_SAP.RH.HbT.Asleep.mean_f,gcampNHCoherData.SSP_SAP.RH.HbT.Asleep.mean_C + gcampNHCoherData.SSP_SAP.RH.HbT.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(gcampNHCoherData.SSP_SAP.RH.HbT.Asleep.mean_f,gcampNHCoherData.SSP_SAP.RH.HbT.Asleep.mean_C - gcampNHCoherData.SSP_SAP.RH.HbT.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('GCaMP-HbT Asleep')
ylabel('Coherence')
xlabel('Freq (Hz)')
xlim([0.02,0.35])
set(gca,'box','off')
axis square
