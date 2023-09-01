function [] = FigS5_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

%% Ephys cross correlation
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_HbTCrossCorr_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Blank_SAP','SSP_SAP'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'lags','xcVals','peak','ttp','group','animalID'};
samplingRate = 30;
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_HbTCrossCorr_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for ee = 1:length(behaviors)
            behavior = behaviors{1,ee};
            hbtXCorrData.(group).dummyCheck = 1;
            if isfield(hbtXCorrData.(group),(behavior)) == false
                for ff = 1:length(variables)
                    variable = variables{1,ff};
                    if any(strcmp(variable,{'group','animalID'})) == true
                        hbtXCorrData.(group).(behavior).(variable) = {};
                    else
                        hbtXCorrData.(group).(behavior).(variable) = [];
                    end
                end
            end
            if isempty(Results_HbTCrossCorr_Ephys.(group).(animalID).(behavior).xcVals) == false
                hbtXCorrData.(group).(behavior).lags = cat(1,hbtXCorrData.(group).(behavior).lags,Results_HbTCrossCorr_Ephys.(group).(animalID).(behavior).lags/samplingRate);
                hbtXCorrData.(group).(behavior).xcVals = cat(1,hbtXCorrData.(group).(behavior).xcVals,Results_HbTCrossCorr_Ephys.(group).(animalID).(behavior).xcVals);
                halfIdx = ceil(length(Results_HbTCrossCorr_Ephys.(group).(animalID).(behavior).xcVals)/2);
                [peak,ttp] = max(Results_HbTCrossCorr_Ephys.(group).(animalID).(behavior).xcVals(halfIdx:end));
                hbtXCorrData.(group).(behavior).peak = cat(1,hbtXCorrData.(group).(behavior).peak,peak);
                hbtXCorrData.(group).(behavior).ttp = cat(1,hbtXCorrData.(group).(behavior).ttp,ttp/samplingRate);
                hbtXCorrData.(group).(behavior).group = cat(1,hbtXCorrData.(group).(behavior).group,group);
                hbtXCorrData.(group).(behavior).animalID = cat(1,hbtXCorrData.(group).(behavior).animalID,animalID);
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for dd = 1:length(behaviors)
        behavior = behaviors{1,dd};
        hbtXCorrData.(group).(behavior).mean_xcVals = mean(hbtXCorrData.(group).(behavior).xcVals,1);
        hbtXCorrData.(group).(behavior).stdErr_xcVals = std(hbtXCorrData.(group).(behavior).xcVals,0,1)./sqrt(size(hbtXCorrData.(group).(behavior).xcVals,1));
        hbtXCorrData.(group).(behavior).mean_lags = mean(hbtXCorrData.(group).(behavior).lags,1);
    end
end
% GLME comparing peak correlation
for bb = 1:length(behaviors)
    behavior = behaviors{1,bb};
    hbtXCorrStats.(behavior).tableSize = cat(1,hbtXCorrData.Blank_SAP.(behavior).peak,hbtXCorrData.SSP_SAP.(behavior).peak);
    hbtXCorrStats.(behavior).Table = table('Size',[size(hbtXCorrStats.(behavior).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'AnimalID','Treatment','Correlation'});
    hbtXCorrStats.(behavior).Table.AnimalID = cat(1,hbtXCorrData.Blank_SAP.(behavior).animalID,hbtXCorrData.SSP_SAP.(behavior).animalID);
    hbtXCorrStats.(behavior).Table.Treatment = cat(1,hbtXCorrData.Blank_SAP.(behavior).group,hbtXCorrData.SSP_SAP.(behavior).group);
    hbtXCorrStats.(behavior).Table.Correlation = cat(1,hbtXCorrData.Blank_SAP.(behavior).peak,hbtXCorrData.SSP_SAP.(behavior).peak);
    hbtXCorrStats.(behavior).FitFormula = 'Correlation ~ 1 + Treatment + (1|AnimalID)';
    hbtXCorrStats.(behavior).Stats = fitglme(hbtXCorrStats.(behavior).Table,hbtXCorrStats.(behavior).FitFormula);
end

%% frontal cortex coherence
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_BilatCoher_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'C','f','binC','binf','group','animalID'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_BilatCoher_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            coherData.(group).(dataType).dummCheck = 1;
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                if isfield(coherData.(group).(dataType),behavior) == false
                    for ff = 1:length(variables)
                        variable = variables{1,ff};
                        if any(strcmp(variable,{'animalID','group','binf'})) == true
                            coherData.(group).(dataType).(behavior).(variable) = {};
                        else
                            coherData.(group).(dataType).(behavior).(variable) = [];
                        end
                    end
                end
                if isempty(Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).C) == false
                    coherData.(group).(dataType).(behavior).C = cat(1,coherData.(group).(dataType).(behavior).C,Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).fC');
                    coherData.(group).(dataType).(behavior).f = cat(1,coherData.(group).(dataType).(behavior).f,Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).ff);
                    freqBand = round(Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).f,2);
                    frequencyList = unique(freqBand);
                    for qq = 1:length(frequencyList)/2
                        freqIdx = find(freqBand == frequencyList(1,qq));
                        coherData.(group).(dataType).(behavior).binC = cat(1,coherData.(group).(dataType).(behavior).binC,mean(Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).C(freqIdx)));
                        coherData.(group).(dataType).(behavior).binf = cat(1,coherData.(group).(dataType).(behavior).binf,num2str(mean(freqBand(freqIdx))));
                        coherData.(group).(dataType).(behavior).group = cat(1,coherData.(group).(dataType).(behavior).group,group);
                        coherData.(group).(dataType).(behavior).animalID = cat(1,coherData.(group).(dataType).(behavior).animalID,animalID);
                    end
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
                if any(strcmp(variable,{'C','f'})) == true
                    coherData.(group).(dataType).(behavior).(['mean_' variable]) = mean(coherData.(group).(dataType).(behavior).(variable),1);
                    coherData.(group).(dataType).(behavior).(['stdErr_' variable]) = std(coherData.(group).(dataType).(behavior).(variable),0,1)./sqrt(size(coherData.(group).(dataType).(behavior).(variable),1));
                end
            end
        end
    end
end
% GLME comparing peak correlation
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        FrontalCoherStats.(dataType).(behavior).tableSize = cat(1,coherData.Blank_SAP.(dataType).(behavior).binC,coherData.SSP_SAP.(dataType).(behavior).binC);
        FrontalCoherStats.(dataType).(behavior).Table = table('Size',[size(FrontalCoherStats.(dataType).(behavior).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'AnimalID','Treatment','Frequency','Coherence'});
        FrontalCoherStats.(dataType).(behavior).Table.AnimalID = cat(1,coherData.Blank_SAP.(dataType).(behavior).animalID,coherData.SSP_SAP.(dataType).(behavior).animalID);
        FrontalCoherStats.(dataType).(behavior).Table.Treatment = cat(1,coherData.Blank_SAP.(dataType).(behavior).group,coherData.SSP_SAP.(dataType).(behavior).group);
        FrontalCoherStats.(dataType).(behavior).Table.Frequency = cat(1,coherData.Blank_SAP.(dataType).(behavior).binf,coherData.SSP_SAP.(dataType).(behavior).binf);
        FrontalCoherStats.(dataType).(behavior).Table.Coherence = cat(1,coherData.Blank_SAP.(dataType).(behavior).binC,coherData.SSP_SAP.(dataType).(behavior).binC);
        FrontalCoherStats.(dataType).(behavior).FitFormula = 'Coherence ~ 1 + Treatment + (1|Frequency) + (1|AnimalID)';
        FrontalCoherStats.(dataType).(behavior).Stats = fitglme(FrontalCoherStats.(dataType).(behavior).Table,FrontalCoherStats.(dataType).(behavior).FitFormula);
    end
end

%% Ephys Pearsons correlations
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PearsonCorr_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','SSP_SAP','Blank_SAP'};
dataTypes = {'HbT','gammaBandPower','deltaBandPower'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PearsonCorr_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            ephysCorrData.(group).(dataType).dummCheck = 1;
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                if isfield(ephysCorrData.(group).(dataType),behavior) == false
                    ephysCorrData.(group).(dataType).(behavior).data = [];
                    ephysCorrData.(group).(dataType).(behavior).group = {};
                    ephysCorrData.(group).(dataType).(behavior).animalID = {};
                end
                if isempty(Results_PearsonCorr_Ephys.(group).(animalID).(dataType).(behavior).R) == false
                    ephysCorrData.(group).(dataType).(behavior).data = cat(1,ephysCorrData.(group).(dataType).(behavior).data,mean(Results_PearsonCorr_Ephys.(group).(animalID).(dataType).(behavior).R));
                    ephysCorrData.(group).(dataType).(behavior).group = cat(1,ephysCorrData.(group).(dataType).(behavior).group,group);
                    ephysCorrData.(group).(dataType).(behavior).animalID = cat(1,ephysCorrData.(group).(dataType).(behavior).animalID,animalID);
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            ephysCorrData.(group).(dataType).(behavior).meanData = mean(ephysCorrData.(group).(dataType).(behavior).data,1);
            ephysCorrData.(group).(dataType).(behavior).stdData = std(ephysCorrData.(group).(dataType).(behavior).data,0,1);
        end
    end
end
% statistics - ttest
for bb = 1:length(dataTypes)
    dataType = dataTypes{1,bb};
    for cc = 1:length(behaviors)
        behavior = behaviors{1,cc};
        % statistics - unpaired ttest
        [EphysPearsonStats.(dataType).(behavior).h,EphysPearsonStats.(dataType).(behavior).p] = ttest2(ephysCorrData.Blank_SAP.(dataType).(behavior).data,ephysCorrData.SSP_SAP.(dataType).(behavior).data);
    end
end

%% GCaMP Pearsons correlations
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PearsonCorr_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'SSP_SAP','Blank_SAP'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PearsonCorr_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            gcampCorrData.(group).(dataType).dummCheck = 1;
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                if isfield(gcampCorrData.(group).(dataType),behavior) == false
                    gcampCorrData.(group).(dataType).(behavior).data = [];
                    gcampCorrData.(group).(dataType).(behavior).group = {};
                    gcampCorrData.(group).(dataType).(behavior).animalID = {};
                end
                if isempty(Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).(behavior).R) == false
                    gcampCorrData.(group).(dataType).(behavior).data = cat(1,gcampCorrData.(group).(dataType).(behavior).data,mean(Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).(behavior).R));
                    gcampCorrData.(group).(dataType).(behavior).group = cat(1,gcampCorrData.(group).(dataType).(behavior).group,group);
                    gcampCorrData.(group).(dataType).(behavior).animalID = cat(1,gcampCorrData.(group).(dataType).(behavior).animalID,animalID);
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            gcampCorrData.(group).(dataType).(behavior).meanData = mean(gcampCorrData.(group).(dataType).(behavior).data,1);
            gcampCorrData.(group).(dataType).(behavior).stdData = std(gcampCorrData.(group).(dataType).(behavior).data,0,1);
        end
    end
end
% statistics - ttest
for bb = 1:length(dataTypes)
    dataType = dataTypes{1,bb};
    for cc = 1:length(behaviors)
        behavior = behaviors{1,cc};
        % statistics - unpaired ttest
        [GCaMPPearsonStats.(dataType).(behavior).h,GCaMPPearsonStats.(dataType).(behavior).p] = ttest2(gcampCorrData.Blank_SAP.(dataType).(behavior).data,gcampCorrData.SSP_SAP.(dataType).(behavior).data);
    end
end

%% figure
FigS5 = figure('Name','Figure S5','units','normalized','outerposition',[0 0 1 1]);

% Ephys cross correlation - Rest
subplot(3,3,1);
plot(hbtXCorrData.Blank_SAP.Rest.mean_lags,hbtXCorrData.Blank_SAP.Rest.mean_xcVals,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(hbtXCorrData.Blank_SAP.Rest.mean_lags,hbtXCorrData.Blank_SAP.Rest.mean_xcVals + hbtXCorrData.Blank_SAP.Rest.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(hbtXCorrData.Blank_SAP.Rest.mean_lags,hbtXCorrData.Blank_SAP.Rest.mean_xcVals - hbtXCorrData.Blank_SAP.Rest.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(hbtXCorrData.SSP_SAP.Rest.mean_lags,hbtXCorrData.SSP_SAP.Rest.mean_xcVals,'color',colors('electric purple'),'LineWidth',2);
plot(hbtXCorrData.SSP_SAP.Rest.mean_lags,hbtXCorrData.SSP_SAP.Rest.mean_xcVals + hbtXCorrData.SSP_SAP.Rest.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
plot(hbtXCorrData.SSP_SAP.Rest.mean_lags,hbtXCorrData.SSP_SAP.Rest.mean_xcVals - hbtXCorrData.SSP_SAP.Rest.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
title('Gamma-HbT Rest')
ylabel('Corr Coef')
xlabel('Lags (s)')
xlim([-5,5])
set(gca,'box','off')
axis square

% Ephys cross correlation - Alert
subplot(3,3,2)
plot(hbtXCorrData.Blank_SAP.Alert.mean_lags,hbtXCorrData.Blank_SAP.Alert.mean_xcVals,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(hbtXCorrData.Blank_SAP.Alert.mean_lags,hbtXCorrData.Blank_SAP.Alert.mean_xcVals + hbtXCorrData.Blank_SAP.Alert.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(hbtXCorrData.Blank_SAP.Alert.mean_lags,hbtXCorrData.Blank_SAP.Alert.mean_xcVals - hbtXCorrData.Blank_SAP.Alert.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(hbtXCorrData.SSP_SAP.Alert.mean_lags,hbtXCorrData.SSP_SAP.Alert.mean_xcVals,'color',colors('electric purple'),'LineWidth',2);
plot(hbtXCorrData.SSP_SAP.Alert.mean_lags,hbtXCorrData.SSP_SAP.Alert.mean_xcVals + hbtXCorrData.SSP_SAP.Alert.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
plot(hbtXCorrData.SSP_SAP.Alert.mean_lags,hbtXCorrData.SSP_SAP.Alert.mean_xcVals - hbtXCorrData.SSP_SAP.Alert.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
title('Gamma-HbT Alert')
ylabel('Corr Coef')
xlabel('Lags (s)')
xlim([-5,5])
set(gca,'box','off')
axis square

% Ephys cross correlation - Asleep
subplot(3,3,3)
plot(hbtXCorrData.Blank_SAP.Asleep.mean_lags,hbtXCorrData.Blank_SAP.Asleep.mean_xcVals,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(hbtXCorrData.Blank_SAP.Asleep.mean_lags,hbtXCorrData.Blank_SAP.Asleep.mean_xcVals + hbtXCorrData.Blank_SAP.Asleep.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(hbtXCorrData.Blank_SAP.Asleep.mean_lags,hbtXCorrData.Blank_SAP.Asleep.mean_xcVals - hbtXCorrData.Blank_SAP.Asleep.stdErr_xcVals,'color',colors('north texas green'),'LineWidth',0.25);
plot(hbtXCorrData.SSP_SAP.Asleep.mean_lags,hbtXCorrData.SSP_SAP.Asleep.mean_xcVals,'color',colors('electric purple'),'LineWidth',2);
plot(hbtXCorrData.SSP_SAP.Asleep.mean_lags,hbtXCorrData.SSP_SAP.Asleep.mean_xcVals + hbtXCorrData.SSP_SAP.Asleep.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
plot(hbtXCorrData.SSP_SAP.Asleep.mean_lags,hbtXCorrData.SSP_SAP.Asleep.mean_xcVals - hbtXCorrData.SSP_SAP.Asleep.stdErr_xcVals,'color',colors('electric purple'),'LineWidth',0.25);
title('Gamma-HbT Asleep')
ylabel('Corr Coef')
xlabel('Lags (s)')
xlim([-5,5])
set(gca,'box','off')
axis square

subplot(3,3,4)
semilogx(coherData.Blank_SAP.HbT.Rest.mean_f,coherData.Blank_SAP.HbT.Rest.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on;
semilogx(coherData.Blank_SAP.HbT.Rest.mean_f,coherData.Blank_SAP.HbT.Rest.mean_C + coherData.Blank_SAP.HbT.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(coherData.Blank_SAP.HbT.Rest.mean_f,coherData.Blank_SAP.HbT.Rest.mean_C - coherData.Blank_SAP.HbT.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(coherData.SSP_SAP.HbT.Rest.mean_f,coherData.SSP_SAP.HbT.Rest.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(coherData.SSP_SAP.HbT.Rest.mean_f,coherData.SSP_SAP.HbT.Rest.mean_C + coherData.SSP_SAP.HbT.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(coherData.SSP_SAP.HbT.Rest.mean_f,coherData.SSP_SAP.HbT.Rest.mean_C - coherData.SSP_SAP.HbT.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
ylabel('Coherence')
xlabel('Freq (Hz)')
ylim([0.7,1])
xlim([0.1,0.5])
set(gca,'box','off')
axis square

subplot(3,3,5)
semilogx(coherData.Blank_SAP.HbT.Alert.mean_f,coherData.Blank_SAP.HbT.Alert.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on;
semilogx(coherData.Blank_SAP.HbT.Alert.mean_f,coherData.Blank_SAP.HbT.Alert.mean_C + coherData.Blank_SAP.HbT.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(coherData.Blank_SAP.HbT.Alert.mean_f,coherData.Blank_SAP.HbT.Alert.mean_C - coherData.Blank_SAP.HbT.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(coherData.SSP_SAP.HbT.Alert.mean_f,coherData.SSP_SAP.HbT.Alert.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(coherData.SSP_SAP.HbT.Alert.mean_f,coherData.SSP_SAP.HbT.Alert.mean_C + coherData.SSP_SAP.HbT.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(coherData.SSP_SAP.HbT.Alert.mean_f,coherData.SSP_SAP.HbT.Alert.mean_C - coherData.SSP_SAP.HbT.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
ylabel('Coherence')
xlabel('Freq (Hz)')
ylim([0.7,1])
xlim([0.01,0.5])
set(gca,'box','off')
axis square

subplot(3,3,6)
semilogx(coherData.Blank_SAP.HbT.Asleep.mean_f,coherData.Blank_SAP.HbT.Asleep.mean_C,'color',colors('north texas green'),'LineWidth',2);
hold on;
semilogx(coherData.Blank_SAP.HbT.Asleep.mean_f,coherData.Blank_SAP.HbT.Asleep.mean_C + coherData.Blank_SAP.HbT.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(coherData.Blank_SAP.HbT.Asleep.mean_f,coherData.Blank_SAP.HbT.Asleep.mean_C - coherData.Blank_SAP.HbT.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(coherData.SSP_SAP.HbT.Asleep.mean_f,coherData.SSP_SAP.HbT.Asleep.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(coherData.SSP_SAP.HbT.Asleep.mean_f,coherData.SSP_SAP.HbT.Asleep.mean_C + coherData.SSP_SAP.HbT.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(coherData.SSP_SAP.HbT.Asleep.mean_f,coherData.SSP_SAP.HbT.Asleep.mean_C - coherData.SSP_SAP.HbT.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
ylabel('Coherence')
xlabel('Freq (Hz)')
ylim([0.7,1])
xlim([0.01,0.5])
set(gca,'box','off')
axis square

subplot(3,4,9)
xInds = ones(1,length(ephysCorrData.Blank_SAP.HbT.Rest.data));
scatter(xInds*1,ephysCorrData.Blank_SAP.HbT.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(1,ephysCorrData.Blank_SAP.HbT.Rest.meanData,ephysCorrData.Blank_SAP.HbT.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysCorrData.SSP_SAP.HbT.Rest.data));
scatter(xInds*2,ephysCorrData.SSP_SAP.HbT.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,ephysCorrData.SSP_SAP.HbT.Rest.meanData,ephysCorrData.SSP_SAP.HbT.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(ephysCorrData.Blank_SAP.HbT.Alert.data));
scatter(xInds*4,ephysCorrData.Blank_SAP.HbT.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(4,ephysCorrData.Blank_SAP.HbT.Alert.meanData,ephysCorrData.Blank_SAP.HbT.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysCorrData.SSP_SAP.HbT.Alert.data));
scatter(xInds*5,ephysCorrData.SSP_SAP.HbT.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,ephysCorrData.SSP_SAP.HbT.Alert.meanData,ephysCorrData.SSP_SAP.HbT.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(ephysCorrData.Blank_SAP.HbT.Asleep.data));
scatter(xInds*6,ephysCorrData.Blank_SAP.HbT.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(6,ephysCorrData.Blank_SAP.HbT.Asleep.meanData,ephysCorrData.Blank_SAP.HbT.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysCorrData.SSP_SAP.HbT.Asleep.data));
scatter(xInds*7,ephysCorrData.SSP_SAP.HbT.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(7,ephysCorrData.SSP_SAP.HbT.Asleep.meanData,ephysCorrData.SSP_SAP.HbT.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('Corr Coef')
xlim([0,8])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(3,4,10)
xInds = ones(1,length(ephysCorrData.Blank_SAP.gammaBandPower.Rest.data));
scatter(xInds*1,ephysCorrData.Blank_SAP.gammaBandPower.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(1,ephysCorrData.Blank_SAP.gammaBandPower.Rest.meanData,ephysCorrData.Blank_SAP.gammaBandPower.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysCorrData.SSP_SAP.gammaBandPower.Rest.data));
scatter(xInds*2,ephysCorrData.SSP_SAP.gammaBandPower.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,ephysCorrData.SSP_SAP.gammaBandPower.Rest.meanData,ephysCorrData.SSP_SAP.gammaBandPower.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(ephysCorrData.Blank_SAP.gammaBandPower.Alert.data));
scatter(xInds*4,ephysCorrData.Blank_SAP.gammaBandPower.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(4,ephysCorrData.Blank_SAP.gammaBandPower.Alert.meanData,ephysCorrData.Blank_SAP.gammaBandPower.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysCorrData.SSP_SAP.gammaBandPower.Alert.data));
scatter(xInds*5,ephysCorrData.SSP_SAP.gammaBandPower.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,ephysCorrData.SSP_SAP.gammaBandPower.Alert.meanData,ephysCorrData.SSP_SAP.gammaBandPower.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(ephysCorrData.Blank_SAP.gammaBandPower.Asleep.data));
scatter(xInds*6,ephysCorrData.Blank_SAP.gammaBandPower.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(6,ephysCorrData.Blank_SAP.gammaBandPower.Asleep.meanData,ephysCorrData.Blank_SAP.gammaBandPower.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysCorrData.SSP_SAP.gammaBandPower.Asleep.data));
scatter(xInds*7,ephysCorrData.SSP_SAP.gammaBandPower.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(7,ephysCorrData.SSP_SAP.gammaBandPower.Asleep.meanData,ephysCorrData.SSP_SAP.gammaBandPower.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('Corr Coef')
xlim([0,8])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(3,4,11)
xInds = ones(1,length(gcampCorrData.Blank_SAP.HbT.Rest.data));
scatter(xInds*1,gcampCorrData.Blank_SAP.HbT.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(1,gcampCorrData.Blank_SAP.HbT.Rest.meanData,gcampCorrData.Blank_SAP.HbT.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampCorrData.SSP_SAP.HbT.Rest.data));
scatter(xInds*2,gcampCorrData.SSP_SAP.HbT.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,gcampCorrData.SSP_SAP.HbT.Rest.meanData,gcampCorrData.SSP_SAP.HbT.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampCorrData.Blank_SAP.HbT.Alert.data));
scatter(xInds*4,gcampCorrData.Blank_SAP.HbT.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(4,gcampCorrData.Blank_SAP.HbT.Alert.meanData,gcampCorrData.Blank_SAP.HbT.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampCorrData.SSP_SAP.HbT.Alert.data));
scatter(xInds*5,gcampCorrData.SSP_SAP.HbT.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,gcampCorrData.SSP_SAP.HbT.Alert.meanData,gcampCorrData.SSP_SAP.HbT.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampCorrData.Blank_SAP.HbT.Asleep.data));
scatter(xInds*6,gcampCorrData.Blank_SAP.HbT.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(6,gcampCorrData.Blank_SAP.HbT.Asleep.meanData,gcampCorrData.Blank_SAP.HbT.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampCorrData.SSP_SAP.HbT.Asleep.data));
scatter(xInds*7,gcampCorrData.SSP_SAP.HbT.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(7,gcampCorrData.SSP_SAP.HbT.Asleep.meanData,gcampCorrData.SSP_SAP.HbT.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('Corr Coef')
xlim([0,8])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(3,4,12)
xInds = ones(1,length(gcampCorrData.Blank_SAP.GCaMP.Rest.data));
scatter(xInds*1,gcampCorrData.Blank_SAP.GCaMP.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(1,gcampCorrData.Blank_SAP.GCaMP.Rest.meanData,gcampCorrData.Blank_SAP.GCaMP.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampCorrData.SSP_SAP.GCaMP.Rest.data));
scatter(xInds*2,gcampCorrData.SSP_SAP.GCaMP.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,gcampCorrData.SSP_SAP.GCaMP.Rest.meanData,gcampCorrData.SSP_SAP.GCaMP.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampCorrData.Blank_SAP.GCaMP.Alert.data));
scatter(xInds*4,gcampCorrData.Blank_SAP.GCaMP.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(4,gcampCorrData.Blank_SAP.GCaMP.Alert.meanData,gcampCorrData.Blank_SAP.GCaMP.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampCorrData.SSP_SAP.GCaMP.Alert.data));
scatter(xInds*5,gcampCorrData.SSP_SAP.GCaMP.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,gcampCorrData.SSP_SAP.GCaMP.Alert.meanData,gcampCorrData.SSP_SAP.GCaMP.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampCorrData.Blank_SAP.GCaMP.Asleep.data));
scatter(xInds*6,gcampCorrData.Blank_SAP.GCaMP.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(6,gcampCorrData.Blank_SAP.GCaMP.Asleep.meanData,gcampCorrData.Blank_SAP.GCaMP.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampCorrData.SSP_SAP.GCaMP.Asleep.data));
scatter(xInds*7,gcampCorrData.SSP_SAP.GCaMP.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(7,gcampCorrData.SSP_SAP.GCaMP.Asleep.meanData,gcampCorrData.SSP_SAP.GCaMP.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('Corr Coef')
xlim([0,8])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

%% save figure
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figure Panels' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(FigS5,[dirpath 'FigS5']);
    set(FigS5,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'FigS5'])
    % statistical diary
    diaryFile = [dirpath 'FigS5_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on

    % HbT-HbT XCorr (Rest)
    disp('======================================================================================================================')
    disp('HbT-HbT XCorr (Rest), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(hbtXCorrData.Blank_SAP.Rest.peak)) ' +/- ' num2str(std(hbtXCorrData.Blank_SAP.Rest.peak,0,1)./sqrt(size(hbtXCorrData.Blank_SAP.Rest.peak,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(hbtXCorrData.SSP_SAP.Rest.peak)) ' +/- ' num2str(std(hbtXCorrData.SSP_SAP.Rest.peak,0,1)./sqrt(size(hbtXCorrData.SSP_SAP.Rest.peak,1)))]); disp(' ')
    disp(hbtXCorrStats.Rest.Stats)

    % HbT-HbT XCorr (Alert)
    disp('======================================================================================================================')
    disp('HbT-HbT XCorr (Alert), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(hbtXCorrData.Blank_SAP.Alert.peak)) ' +/- ' num2str(std(hbtXCorrData.Blank_SAP.Alert.peak,0,1)./sqrt(size(hbtXCorrData.Blank_SAP.Alert.peak,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(hbtXCorrData.SSP_SAP.Alert.peak)) ' +/- ' num2str(std(hbtXCorrData.SSP_SAP.Alert.peak,0,1)./sqrt(size(hbtXCorrData.SSP_SAP.Alert.peak,1)))]); disp(' ')
    disp(hbtXCorrStats.Alert.Stats)

    % HbT-HbT XCorr (Asleep)
    disp('======================================================================================================================')
    disp('HbT-HbT XCorr (Asleep), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(hbtXCorrData.Blank_SAP.Asleep.peak)) ' +/- ' num2str(std(hbtXCorrData.Blank_SAP.Asleep.peak,0,1)./sqrt(size(hbtXCorrData.Blank_SAP.Asleep.peak,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(hbtXCorrData.SSP_SAP.Asleep.peak)) ' +/- ' num2str(std(hbtXCorrData.SSP_SAP.Asleep.peak,0,1)./sqrt(size(hbtXCorrData.SSP_SAP.Asleep.peak,1)))]); disp(' ')
    disp(hbtXCorrStats.Asleep.Stats)

    % GCaMP-HbT frontal coherence(Rest)
    disp('======================================================================================================================')
    disp('GCaMP-HbT frontal coherence (Rest), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(coherData.Blank_SAP.HbT.Rest.binC)) ' +/- ' num2str(std(coherData.Blank_SAP.HbT.Rest.binC,0,1)./sqrt(size(coherData.Blank_SAP.HbT.Rest.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(coherData.SSP_SAP.HbT.Rest.binC)) ' +/- ' num2str(std(coherData.SSP_SAP.HbT.Rest.binC,0,1)./sqrt(size(coherData.SSP_SAP.HbT.Rest.binC,1)))]); disp(' ')
    disp(FrontalCoherStats.HbT.Rest.Stats)

    % GCaMP-HbT frontal coherence (Alert)
    disp('======================================================================================================================')
    disp('GCaMP-HbT frontal coherence (Alert), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(coherData.Blank_SAP.HbT.Alert.binC)) ' +/- ' num2str(std(coherData.Blank_SAP.HbT.Alert.binC,0,1)./sqrt(size(coherData.Blank_SAP.HbT.Alert.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(coherData.SSP_SAP.HbT.Alert.binC)) ' +/- ' num2str(std(coherData.SSP_SAP.HbT.Alert.binC,0,1)./sqrt(size(coherData.SSP_SAP.HbT.Alert.binC,1)))]); disp(' ')
    disp(FrontalCoherStats.HbT.Alert.Stats)

    % GCaMP-HbT frontal coherence (Asleep)
    disp('======================================================================================================================')
    disp('GCaMP-HbT frontal coherence (Asleep), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(coherData.Blank_SAP.HbT.Asleep.binC)) ' +/- ' num2str(std(coherData.Blank_SAP.HbT.Asleep.binC,0,1)./sqrt(size(coherData.Blank_SAP.HbT.Asleep.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(coherData.SSP_SAP.HbT.Asleep.binC)) ' +/- ' num2str(std(coherData.SSP_SAP.HbT.Asleep.binC,0,1)./sqrt(size(coherData.SSP_SAP.HbT.Asleep.binC,1)))]); disp(' ')
    disp(FrontalCoherStats.HbT.Asleep.Stats)

    % [HbT] Pearson's correlation during each arousal-state
    disp('======================================================================================================================')
    disp('[HbT] Pearson''s correlation during each arousal-state (ephys), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(ephysCorrData.Blank_SAP.HbT.Rest.meanData) ' +/- ' num2str(ephysCorrData.Blank_SAP.HbT.Rest.stdData)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(ephysCorrData.SSP_SAP.HbT.Rest.meanData) ' +/- ' num2str(ephysCorrData.SSP_SAP.HbT.Rest.stdData)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(EphysPearsonStats.HbT.Rest.p)]); disp(' ')
    disp(['Blank-SAP Alert: ' num2str(ephysCorrData.Blank_SAP.HbT.Alert.meanData) ' +/- ' num2str(ephysCorrData.Blank_SAP.HbT.Alert.stdData)]); disp(' ')
    disp(['SSP-SAP Alert: ' num2str(ephysCorrData.SSP_SAP.HbT.Alert.meanData) ' +/- ' num2str(ephysCorrData.SSP_SAP.HbT.Alert.stdData)]); disp(' ')
    disp(['Blank vs. SAP Alert ttest p = ' num2str(EphysPearsonStats.HbT.Alert.p)]); disp(' ')
    disp(['Blank-SAP Asleep: ' num2str(ephysCorrData.Blank_SAP.HbT.Asleep.meanData) ' +/- ' num2str(ephysCorrData.Blank_SAP.HbT.Asleep.stdData)]); disp(' ')
    disp(['SSP-SAP Asleep: ' num2str(ephysCorrData.SSP_SAP.HbT.Asleep.meanData) ' +/- ' num2str(ephysCorrData.SSP_SAP.HbT.Asleep.stdData)]); disp(' ')
    disp(['Blank vs. SAP Asleep ttest p = ' num2str(EphysPearsonStats.HbT.Asleep.p)]); disp(' ')

    % [gammaBandPower] Pearson's correlation during each arousal-state
    disp('======================================================================================================================')
    disp('[gammaBandPower] Pearson''s correlation during each arousal-state (ephys), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(ephysCorrData.Blank_SAP.gammaBandPower.Rest.meanData) ' +/- ' num2str(ephysCorrData.Blank_SAP.gammaBandPower.Rest.stdData)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(ephysCorrData.SSP_SAP.gammaBandPower.Rest.meanData) ' +/- ' num2str(ephysCorrData.SSP_SAP.gammaBandPower.Rest.stdData)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(EphysPearsonStats.gammaBandPower.Rest.p)]); disp(' ')
    disp(['Blank-SAP Alert: ' num2str(ephysCorrData.Blank_SAP.gammaBandPower.Alert.meanData) ' +/- ' num2str(ephysCorrData.Blank_SAP.gammaBandPower.Alert.stdData)]); disp(' ')
    disp(['SSP-SAP Alert: ' num2str(ephysCorrData.SSP_SAP.gammaBandPower.Alert.meanData) ' +/- ' num2str(ephysCorrData.SSP_SAP.gammaBandPower.Alert.stdData)]); disp(' ')
    disp(['Blank vs. SAP Alert ttest p = ' num2str(EphysPearsonStats.gammaBandPower.Alert.p)]); disp(' ')
    disp(['Blank-SAP Asleep: ' num2str(ephysCorrData.Blank_SAP.gammaBandPower.Asleep.meanData) ' +/- ' num2str(ephysCorrData.Blank_SAP.gammaBandPower.Asleep.stdData)]); disp(' ')
    disp(['SSP-SAP Asleep: ' num2str(ephysCorrData.SSP_SAP.gammaBandPower.Asleep.meanData) ' +/- ' num2str(ephysCorrData.SSP_SAP.gammaBandPower.Asleep.stdData)]); disp(' ')
    disp(['Blank vs. SAP Asleep ttest p = ' num2str(EphysPearsonStats.gammaBandPower.Asleep.p)]); disp(' ')

    % [HbT] Pearson's correlation during each arousal-state
    disp('======================================================================================================================')
    disp('[HbT] Pearson''s correlation during each arousal-state (ephys), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(gcampCorrData.Blank_SAP.HbT.Rest.meanData) ' +/- ' num2str(gcampCorrData.Blank_SAP.HbT.Rest.stdData)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(gcampCorrData.SSP_SAP.HbT.Rest.meanData) ' +/- ' num2str(gcampCorrData.SSP_SAP.HbT.Rest.stdData)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(GCaMPPearsonStats.HbT.Rest.p)]); disp(' ')
    disp(['Blank-SAP Alert: ' num2str(gcampCorrData.Blank_SAP.HbT.Alert.meanData) ' +/- ' num2str(gcampCorrData.Blank_SAP.HbT.Alert.stdData)]); disp(' ')
    disp(['SSP-SAP Alert: ' num2str(gcampCorrData.SSP_SAP.HbT.Alert.meanData) ' +/- ' num2str(gcampCorrData.SSP_SAP.HbT.Alert.stdData)]); disp(' ')
    disp(['Blank vs. SAP Alert ttest p = ' num2str(GCaMPPearsonStats.HbT.Alert.p)]); disp(' ')
    disp(['Blank-SAP Asleep: ' num2str(gcampCorrData.Blank_SAP.HbT.Asleep.meanData) ' +/- ' num2str(gcampCorrData.Blank_SAP.HbT.Asleep.stdData)]); disp(' ')
    disp(['SSP-SAP Asleep: ' num2str(gcampCorrData.SSP_SAP.HbT.Asleep.meanData) ' +/- ' num2str(gcampCorrData.SSP_SAP.HbT.Asleep.stdData)]); disp(' ')
    disp(['Blank vs. SAP Asleep ttest p = ' num2str(GCaMPPearsonStats.HbT.Asleep.p)]); disp(' ')

    % [GCaMP] Pearson's correlation during each arousal-state
    disp('======================================================================================================================')
    disp('[GCaMP] Pearson''s correlation during each arousal-state (ephys), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(gcampCorrData.Blank_SAP.GCaMP.Rest.meanData) ' +/- ' num2str(gcampCorrData.Blank_SAP.GCaMP.Rest.stdData)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(gcampCorrData.SSP_SAP.GCaMP.Rest.meanData) ' +/- ' num2str(gcampCorrData.SSP_SAP.GCaMP.Rest.stdData)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(GCaMPPearsonStats.GCaMP.Rest.p)]); disp(' ')
    disp(['Blank-SAP Alert: ' num2str(gcampCorrData.Blank_SAP.GCaMP.Alert.meanData) ' +/- ' num2str(gcampCorrData.Blank_SAP.GCaMP.Alert.stdData)]); disp(' ')
    disp(['SSP-SAP Alert: ' num2str(gcampCorrData.SSP_SAP.GCaMP.Alert.meanData) ' +/- ' num2str(gcampCorrData.SSP_SAP.GCaMP.Alert.stdData)]); disp(' ')
    disp(['Blank vs. SAP Alert ttest p = ' num2str(GCaMPPearsonStats.GCaMP.Alert.p)]); disp(' ')
    disp(['Blank-SAP Asleep: ' num2str(gcampCorrData.Blank_SAP.GCaMP.Asleep.meanData) ' +/- ' num2str(gcampCorrData.Blank_SAP.GCaMP.Asleep.stdData)]); disp(' ')
    disp(['SSP-SAP Asleep: ' num2str(gcampCorrData.SSP_SAP.GCaMP.Asleep.meanData) ' +/- ' num2str(gcampCorrData.SSP_SAP.GCaMP.Asleep.stdData)]); disp(' ')
    disp(['Blank vs. SAP Asleep ttest p = ' num2str(GCaMPPearsonStats.GCaMP.Asleep.p)]); disp(' ')

    diary off
end