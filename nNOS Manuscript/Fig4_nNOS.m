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
variables = {'lags','xcVals','peak','ttp','group','animalID'};
samplingRate = 30;
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_CrossCorr_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    ephysXCorrData.(group).(hemisphere).(dataType).dummyCheck = 1;
                    if isfield(ephysXCorrData.(group).(hemisphere).(dataType),(behavior)) == false
                        for ff = 1:length(variables)
                            variable = variables{1,ff};
                            if any(strcmp(variable,{'group','animalID'})) == true
                                ephysXCorrData.(group).(hemisphere).(dataType).(behavior).(variable) = {};
                            else
                                ephysXCorrData.(group).(hemisphere).(dataType).(behavior).(variable) = [];
                            end
                        end
                    end
                    if isempty(Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals) == false
                        ephysXCorrData.(group).(hemisphere).(dataType).(behavior).lags = cat(1,ephysXCorrData.(group).(hemisphere).(dataType).(behavior).lags,Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).lags/samplingRate);
                        ephysXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals = cat(1,ephysXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals,Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals);
                        halfIdx = ceil(length(Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals)/2);
                        [peak,ttp] = max(Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals(halfIdx:end));
                        ephysXCorrData.(group).(hemisphere).(dataType).(behavior).peak = cat(1,ephysXCorrData.(group).(hemisphere).(dataType).(behavior).peak,peak);
                        ephysXCorrData.(group).(hemisphere).(dataType).(behavior).ttp = cat(1,ephysXCorrData.(group).(hemisphere).(dataType).(behavior).ttp,ttp/samplingRate);
                        ephysXCorrData.(group).(hemisphere).(dataType).(behavior).group = cat(1,ephysXCorrData.(group).(hemisphere).(dataType).(behavior).group,group);
                        ephysXCorrData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,ephysXCorrData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
% GLME comparing peak correlation
for cc = 1:length(hemispheres)
    hemisphere = hemispheres{1,cc};
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        for bb = 1:length(behaviors)
            behavior = behaviors{1,bb};
            ephysXCorrStats.(hemisphere).(dataType).(behavior).tableSize = cat(1,ephysXCorrData.Blank_SAP.(hemisphere).(dataType).(behavior).peak,ephysXCorrData.SSP_SAP.(hemisphere).(dataType).(behavior).peak);
            ephysXCorrStats.(hemisphere).(dataType).(behavior).Table = table('Size',[size(ephysXCorrStats.(hemisphere).(dataType).(behavior).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'AnimalID','Treatment','Correlation'});
            ephysXCorrStats.(hemisphere).(dataType).(behavior).Table.AnimalID = cat(1,ephysXCorrData.Blank_SAP.(hemisphere).(dataType).(behavior).animalID,ephysXCorrData.SSP_SAP.(hemisphere).(dataType).(behavior).animalID);
            ephysXCorrStats.(hemisphere).(dataType).(behavior).Table.Treatment = cat(1,ephysXCorrData.Blank_SAP.(hemisphere).(dataType).(behavior).group,ephysXCorrData.SSP_SAP.(hemisphere).(dataType).(behavior).group);
            ephysXCorrStats.(hemisphere).(dataType).(behavior).Table.Correlation = cat(1,ephysXCorrData.Blank_SAP.(hemisphere).(dataType).(behavior).peak,ephysXCorrData.SSP_SAP.(hemisphere).(dataType).(behavior).peak);
            ephysXCorrStats.(hemisphere).(dataType).(behavior).FitFormula = 'Correlation ~ 1 + Treatment + (1|AnimalID)';
            ephysXCorrStats.(hemisphere).(dataType).(behavior).Stats = fitglme(ephysXCorrStats.(hemisphere).(dataType).(behavior).Table,ephysXCorrStats.(hemisphere).(dataType).(behavior).FitFormula);
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
variables = {'C','f','binC','binf','group','animalID'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_NeuralHemoCoher_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                ephysNHCoherData.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(ephysNHCoherData.(group).(hemisphere).(dataType),behavior) == false
                        for ff = 1:length(variables)
                            variable = variables{1,ff};
                            if any(strcmp(variable,{'animalID','group','binf'})) == true
                                ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable) = {};
                            else
                                ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable) = [];
                            end
                        end
                    end
                    if isempty(Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).C) == false
                        ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).C = cat(1,ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).C,Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).C');
                        ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).f = cat(1,ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).f,Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).f);
                        freqBand = round(Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).f,2);
                        frequencyList = unique(freqBand);
                        for qq = 1:length(frequencyList)/2
                            freqIdx = find(freqBand == frequencyList(1,qq));
                            ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).binC = cat(1,ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).binC,mean(Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).C(freqIdx)));
                            ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).binf = cat(1,ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).binf,num2str(mean(freqBand(freqIdx))));
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
                    if any(strcmp(variable,{'C','f'})) == true
                        ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                        ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).(['stdErr_' variable]) = std(ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable),0,1)./sqrt(size(ephysNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable),1));
                    end
                end
            end
        end
    end
end
% GLME comparing peak correlation
for cc = 1:length(hemispheres)
    hemisphere = hemispheres{1,cc};
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        for bb = 1:length(behaviors)
            behavior = behaviors{1,bb};
            ephysNHCoherStats.(hemisphere).(dataType).(behavior).tableSize = cat(1,ephysNHCoherData.Blank_SAP.(hemisphere).(dataType).(behavior).binC,ephysNHCoherData.SSP_SAP.(hemisphere).(dataType).(behavior).binC);
            ephysNHCoherStats.(hemisphere).(dataType).(behavior).Table = table('Size',[size(ephysNHCoherStats.(hemisphere).(dataType).(behavior).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'AnimalID','Treatment','Frequency','Coherence'});
            ephysNHCoherStats.(hemisphere).(dataType).(behavior).Table.AnimalID = cat(1,ephysNHCoherData.Blank_SAP.(hemisphere).(dataType).(behavior).animalID,ephysNHCoherData.SSP_SAP.(hemisphere).(dataType).(behavior).animalID);
            ephysNHCoherStats.(hemisphere).(dataType).(behavior).Table.Treatment = cat(1,ephysNHCoherData.Blank_SAP.(hemisphere).(dataType).(behavior).group,ephysNHCoherData.SSP_SAP.(hemisphere).(dataType).(behavior).group);
            ephysNHCoherStats.(hemisphere).(dataType).(behavior).Table.Frequency = cat(1,ephysNHCoherData.Blank_SAP.(hemisphere).(dataType).(behavior).binf,ephysNHCoherData.SSP_SAP.(hemisphere).(dataType).(behavior).binf);
            ephysNHCoherStats.(hemisphere).(dataType).(behavior).Table.Coherence = cat(1,ephysNHCoherData.Blank_SAP.(hemisphere).(dataType).(behavior).binC,ephysNHCoherData.SSP_SAP.(hemisphere).(dataType).(behavior).binC);
            ephysNHCoherStats.(hemisphere).(dataType).(behavior).FitFormula = 'Coherence ~ 1 + Treatment + (1|Frequency) + (1|AnimalID)';
            ephysNHCoherStats.(hemisphere).(dataType).(behavior).Stats = fitglme(ephysNHCoherStats.(hemisphere).(dataType).(behavior).Table,ephysNHCoherStats.(hemisphere).(dataType).(behavior).FitFormula);
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
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).peak = [];
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).ttp = [];
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).group = {};
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    if isempty(Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals) == false
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).lags = cat(1,gcampXCorrData.(group).(hemisphere).(dataType).(behavior).lags,Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).lags/samplingRate);
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals = cat(1,gcampXCorrData.(group).(hemisphere).(dataType).(behavior).xcVals,Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals);
                        halfIdx = ceil(length(Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals)/2);
                        [peak,ttp] = max(Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).xcVals(halfIdx:end));
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).peak = cat(1,gcampXCorrData.(group).(hemisphere).(dataType).(behavior).peak,peak);
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).ttp = cat(1,gcampXCorrData.(group).(hemisphere).(dataType).(behavior).ttp,ttp/samplingRate);
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).group = cat(1,gcampXCorrData.(group).(hemisphere).(dataType).(behavior).group,group);
                        gcampXCorrData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,gcampXCorrData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
% GLME comparing peak correlation
for cc = 1:length(hemispheres)
    hemisphere = hemispheres{1,cc};
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        for bb = 1:length(behaviors)
            behavior = behaviors{1,bb};
            gcampXCorrStats.(hemisphere).(dataType).(behavior).tableSize = cat(1,gcampXCorrData.Blank_SAP.(hemisphere).(dataType).(behavior).peak,gcampXCorrData.SSP_SAP.(hemisphere).(dataType).(behavior).peak);
            gcampXCorrStats.(hemisphere).(dataType).(behavior).Table = table('Size',[size(gcampXCorrStats.(hemisphere).(dataType).(behavior).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'AnimalID','Treatment','Correlation'});
            gcampXCorrStats.(hemisphere).(dataType).(behavior).Table.AnimalID = cat(1,gcampXCorrData.Blank_SAP.(hemisphere).(dataType).(behavior).animalID,gcampXCorrData.SSP_SAP.(hemisphere).(dataType).(behavior).animalID);
            gcampXCorrStats.(hemisphere).(dataType).(behavior).Table.Treatment = cat(1,gcampXCorrData.Blank_SAP.(hemisphere).(dataType).(behavior).group,gcampXCorrData.SSP_SAP.(hemisphere).(dataType).(behavior).group);
            gcampXCorrStats.(hemisphere).(dataType).(behavior).Table.Correlation = cat(1,gcampXCorrData.Blank_SAP.(hemisphere).(dataType).(behavior).peak,gcampXCorrData.SSP_SAP.(hemisphere).(dataType).(behavior).peak);
            gcampXCorrStats.(hemisphere).(dataType).(behavior).FitFormula = 'Correlation ~ 1 + Treatment + (1|AnimalID)';
            gcampXCorrStats.(hemisphere).(dataType).(behavior).Stats = fitglme(gcampXCorrStats.(hemisphere).(dataType).(behavior).Table,gcampXCorrStats.(hemisphere).(dataType).(behavior).FitFormula);
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
variables = {'C','f','binC','binf','group','animalID'};
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
                        for ff = 1:length(variables)
                            variable = variables{1,ff};
                            if any(strcmp(variable,{'animalID','group','binf'})) == true
                                gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable) = {};
                            else
                                gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable) = [];
                            end
                        end
                    end
                    if isempty(Results_NeuralHemoCoher_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).C) == false
                        gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).C = cat(1,gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).C,Results_NeuralHemoCoher_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).C');
                        gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).f = cat(1,gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).f,Results_NeuralHemoCoher_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).f);
                        freqBand = round(Results_NeuralHemoCoher_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).f,2);
                        frequencyList = unique(freqBand);
                        for qq = 1:length(frequencyList)/2
                            freqIdx = find(freqBand == frequencyList(1,qq));
                            gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).binC = cat(1,gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).binC,mean(Results_NeuralHemoCoher_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).C(freqIdx)));
                            gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).binf = cat(1,gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).binf,num2str(mean(freqBand(freqIdx))));
                            gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).group = cat(1,gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).group,group);
                            gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
                    if any(strcmp(variable,{'C','f'})) == true
                        gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                        gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).(['stdErr_' variable]) = std(gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable),0,1)./sqrt(size(gcampNHCoherData.(group).(hemisphere).(dataType).(behavior).(variable),1));
                    end
                end
            end
        end
    end
end
% GLME comparing peak correlation
for cc = 1:length(hemispheres)
    hemisphere = hemispheres{1,cc};
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        for bb = 1:length(behaviors)
            behavior = behaviors{1,bb};
            gcampNHCoherStats.(hemisphere).(dataType).(behavior).tableSize = cat(1,gcampNHCoherData.Blank_SAP.(hemisphere).(dataType).(behavior).binC,gcampNHCoherData.SSP_SAP.(hemisphere).(dataType).(behavior).binC);
            gcampNHCoherStats.(hemisphere).(dataType).(behavior).Table = table('Size',[size(gcampNHCoherStats.(hemisphere).(dataType).(behavior).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'AnimalID','Treatment','Frequency','Coherence'});
            gcampNHCoherStats.(hemisphere).(dataType).(behavior).Table.AnimalID = cat(1,gcampNHCoherData.Blank_SAP.(hemisphere).(dataType).(behavior).animalID,gcampNHCoherData.SSP_SAP.(hemisphere).(dataType).(behavior).animalID);
            gcampNHCoherStats.(hemisphere).(dataType).(behavior).Table.Treatment = cat(1,gcampNHCoherData.Blank_SAP.(hemisphere).(dataType).(behavior).group,gcampNHCoherData.SSP_SAP.(hemisphere).(dataType).(behavior).group);
            gcampNHCoherStats.(hemisphere).(dataType).(behavior).Table.Frequency = cat(1,gcampNHCoherData.Blank_SAP.(hemisphere).(dataType).(behavior).binf,gcampNHCoherData.SSP_SAP.(hemisphere).(dataType).(behavior).binf);
            gcampNHCoherStats.(hemisphere).(dataType).(behavior).Table.Coherence = cat(1,gcampNHCoherData.Blank_SAP.(hemisphere).(dataType).(behavior).binC,gcampNHCoherData.SSP_SAP.(hemisphere).(dataType).(behavior).binC);
            gcampNHCoherStats.(hemisphere).(dataType).(behavior).FitFormula = 'Coherence ~ 1 + Treatment + (1|Frequency) + (1|AnimalID)';
            gcampNHCoherStats.(hemisphere).(dataType).(behavior).Stats = fitglme(gcampNHCoherStats.(hemisphere).(dataType).(behavior).Table,gcampNHCoherStats.(hemisphere).(dataType).(behavior).FitFormula);
        end
    end
end
%% figure panel 4
Fig4 = figure('Name','Figure 4','units','normalized','outerposition',[0 0 1 1]);
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
ylim([-0.05,0.45])
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
ylim([-0.05,0.45])
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
ylim([-0.05,0.45])
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
xlim([0.1,0.5])
ylim([0.1,0.8])
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
xlim([0.01,0.5])
ylim([0.1,0.8])
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
xlim([0.01,0.5])
ylim([0.1,0.8])
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
ylim([-0.1,0.9])
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
ylim([-0.1,0.9])
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
ylim([-0.1,0.9])
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
xlim([0.1,0.5])
ylim([0.2,1])
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
xlim([0.01,0.5])
ylim([0.2,1])
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
xlim([0.01,0.5])
ylim([0.2,1])
set(gca,'box','off')
axis square
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figure Panels' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig4,[dirpath 'Fig4']);
    set(Fig4,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'Fig4'])
    diaryFile = [dirpath 'Fig4_Readout.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    % statistical diary
    diary(diaryFile)
    diary on
    disp('======================================================================================================================')
    disp('GLME statistics for gamma-HbT XCorr during Rest [Ephys]')
    disp('======================================================================================================================')
    disp(ephysXCorrStats.RH.gammaBandPower.Rest.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for gamma-HbT XCorr during Alert [Ephys]')
    disp('======================================================================================================================')
    disp(ephysXCorrStats.RH.gammaBandPower.Alert.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for gamma-HbT XCorr during Asleep [Ephys]')
    disp('======================================================================================================================')
    disp(ephysXCorrStats.RH.gammaBandPower.Asleep.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for gamma-HbT Coherence during Rest [Ephys]')
    disp('======================================================================================================================')
    disp(ephysNHCoherStats.RH.gammaBandPower.Rest.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for gamma-HbT Coherence during Alert [Ephys]')
    disp('======================================================================================================================')
    disp(ephysNHCoherStats.RH.gammaBandPower.Alert.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for gamma-HbT Coherence during Asleep [Ephys]')
    disp('======================================================================================================================')
    disp(ephysNHCoherStats.RH.gammaBandPower.Asleep.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for GCaMP-HbT XCorr during Rest [GCaMP]')
    disp('======================================================================================================================')
    disp(gcampXCorrStats.RH.HbT.Rest.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for GCaMP-HbT XCorr during Alert [GCaMP]')
    disp('======================================================================================================================')
    disp(gcampXCorrStats.RH.HbT.Alert.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for GCaMP-HbT XCorr during Asleep [GCaMP]')
    disp('======================================================================================================================')
    disp(gcampXCorrStats.RH.HbT.Asleep.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for GCaMP-HbT Coherence during Rest [GCaMP]')
    disp('======================================================================================================================')
    disp(gcampNHCoherStats.RH.HbT.Rest.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for GCaMP-HbT Coherence during Alert [GCaMP]')
    disp('======================================================================================================================')
    disp(gcampNHCoherStats.RH.HbT.Alert.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for GCaMP-HbT Coherence during Asleep [GCaMP]')
    disp('======================================================================================================================')
    disp(gcampNHCoherStats.RH.HbT.Asleep.Stats)
    diary off
end