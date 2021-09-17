function [] = NeuralHemoCoherence_Saporin(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________


%% set-up and process data
expGroups = {'C57BL6J','SSP-SAP','Blank-SAP'};
setName = 'IOS Set A';
animalIDs.all = {};
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs.all = horzcat(animalIDs.all,{folderList.name});
    animalIDs.(strrep(expGroups{1,aa},'-','_')) = {folderList.name};
end
treatments = {'C57BL6J','SSP_SAP','Blank_SAP'};
behavFields = {'Awake','Sleep','All'};
dataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
%% average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,aa};
    % recognize treatment based on animal group
    if ismember(animalID,animalIDs.C57BL6J) == true
        treatment = 'C57BL6J';
    elseif ismember(animalID,animalIDs.SSP_SAP) == true
        treatment = 'SSP_SAP';
    elseif ismember(animalID,animalIDs.Blank_SAP) == true
        treatment = 'Blank_SAP';
    end
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % pre-allocate necessary variable fields
            data.(treatment).(behavField).dummCheck = 1;
            if isfield(data.(treatment).(behavField),dataType) == false
                data.(treatment).(behavField).(dataType).LH.C = [];
                data.(treatment).(behavField).(dataType).LH.f = [];
                data.(treatment).(behavField).(dataType).LH.confC = [];
                data.(treatment).(behavField).(dataType).RH.C = [];
                data.(treatment).(behavField).(dataType).RH.f = [];
                data.(treatment).(behavField).(dataType).RH.confC = [];
                data.(treatment).(behavField).(dataType).LH.animalID = {};
                data.(treatment).(behavField).(dataType).RH.animalID = {};
                data.(treatment).(behavField).(dataType).LH.treatment = {};
                data.(treatment).(behavField).(dataType).RH.treatment = {};
            end
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C) == false
                % concatenate C/f for existing data - exclude any empty sets
                data.(treatment).(behavField).(dataType).LH.C = cat(2,data.(treatment).(behavField).(dataType).LH.C,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C.^2);
                data.(treatment).(behavField).(dataType).LH.f = cat(1,data.(treatment).(behavField).(dataType).LH.f,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.f);
                data.(treatment).(behavField).(dataType).LH.confC = cat(1,data.(treatment).(behavField).(dataType).LH.confC,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.confC);
                data.(treatment).(behavField).(dataType).RH.C = cat(2,data.(treatment).(behavField).(dataType).RH.C,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.C.^2);
                data.(treatment).(behavField).(dataType).RH.f = cat(1,data.(treatment).(behavField).(dataType).RH.f,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.f);
                data.(treatment).(behavField).(dataType).RH.confC = cat(1,data.(treatment).(behavField).(dataType).RH.confC,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.confC);
                data.(treatment).(behavField).(dataType).LH.animalID = cat(1,data.(treatment).(behavField).(dataType).LH.animalID,animalID);
                data.(treatment).(behavField).(dataType).RH.animalID = cat(1,data.(treatment).(behavField).(dataType).RH.animalID,animalID);
                data.(treatment).(behavField).(dataType).LH.treatment = cat(1,data.(treatment).(behavField).(dataType).LH.treatment,treatment);
                data.(treatment).(behavField).(dataType).RH.treatment = cat(1,data.(treatment).(behavField).(dataType).RH.treatment,treatment);
            end
        end
    end
end
%% take mean/StD of C/f and determine confC line
for qq = 1:length(treatments)
    treatment = treatments{1,qq};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        for ff = 1:length(dataTypes)
            dataType = dataTypes{1,ff};
            data.(treatment).(behavField).(dataType).LH.meanC = mean(data.(treatment).(behavField).(dataType).LH.C,2);
            data.(treatment).(behavField).(dataType).LH.stdC = std(data.(treatment).(behavField).(dataType).LH.C,0,2);
            data.(treatment).(behavField).(dataType).LH.meanf = mean(data.(treatment).(behavField).(dataType).LH.f,1);
            data.(treatment).(behavField).(dataType).LH.maxConfC = geomean(data.(treatment).(behavField).(dataType).LH.confC);
            data.(treatment).(behavField).(dataType).LH.maxConfC_Y = ones(length(data.(treatment).(behavField).(dataType).LH.meanf),1)*data.(treatment).(behavField).(dataType).LH.maxConfC;
            data.(treatment).(behavField).(dataType).RH.meanC = mean(data.(treatment).(behavField).(dataType).RH.C,2);
            data.(treatment).(behavField).(dataType).RH.stdC = std(data.(treatment).(behavField).(dataType).RH.C,0,2);
            data.(treatment).(behavField).(dataType).RH.meanf = mean(data.(treatment).(behavField).(dataType).RH.f,1);
            data.(treatment).(behavField).(dataType).RH.maxConfC = geomean(data.(treatment).(behavField).(dataType).RH.confC);
            data.(treatment).(behavField).(dataType).RH.maxConfC_Y = ones(length(data.(treatment).(behavField).(dataType).RH.meanf),1)*data.(treatment).(behavField).(dataType).RH.maxConfC;
        end
    end
end
%% find Hz peaks in coherence
for qq = 1:length(treatments)
    treatment = treatments{1,qq};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        for ff = 1:length(dataTypes)
            dataType = dataTypes{1,ff};
            for gg = 1:size(data.(treatment).(behavField).(dataType).RH.C,2)
                F = round(data.(treatment).(behavField).(dataType).RH.f(gg,:),3);
                RH_C = data.(treatment).(behavField).(dataType).RH.C(:,gg);
                index001 = find(F == 0.01);
                index01 = find(F == 0.1);
                index05 = find(F == 0.5);
                data.(treatment).(behavField).(dataType).RH.C001(gg,1) = mean(RH_C(1:index001(1)));
                data.(treatment).(behavField).(dataType).RH.C01(gg,1) = mean(RH_C(index001(1) + 1:index01(1)));
                data.(treatment).(behavField).(dataType).RH.C05(gg,1) = mean(RH_C(index01(1) + 1:index05(1)));
            end
        end
    end
end
%% statistics - generalized linear mixed effects model
freqBands = {'C001','C01','C05'};
for aa = 1:length(freqBands)
    freqBand = freqBands{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        % statistics - generalized linear mixed effects model
        Stats.(freqBand).(behavField).tableSize = cat(1,data.Blank_SAP.(behavField).gammaBandPower.RH.(freqBand),data.SSP_SAP.(behavField).gammaBandPower.RH.(freqBand));
        Stats.(freqBand).(behavField).Table = table('Size',[size(Stats.(freqBand).(behavField).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Treatment',freqBand});
        Stats.(freqBand).(behavField).Table.Mouse = cat(1,data.Blank_SAP.(behavField).gammaBandPower.RH.animalID,data.SSP_SAP.(behavField).gammaBandPower.RH.animalID);
        Stats.(freqBand).(behavField).Table.Treatment = cat(1,data.Blank_SAP.(behavField).gammaBandPower.RH.treatment,data.SSP_SAP.(behavField).gammaBandPower.RH.treatment);
        Stats.(freqBand).(behavField).Table.(freqBand) = cat(1,data.Blank_SAP.(behavField).gammaBandPower.RH.(freqBand),data.SSP_SAP.(behavField).gammaBandPower.RH.(freqBand));
        Stats.(freqBand).(behavField).FitFormula = [freqBand ' ~ 1 + Treatment + (1|Mouse)'];
        Stats.(freqBand).(behavField).Stats = fitglme(Stats.(freqBand).(behavField).Table,Stats.(freqBand).(behavField).FitFormula,'Link','probit');
    end
end
%% average HbT coherence
summaryFigure1 = figure;
sgtitle('Gamma-band - \DeltaHbT Coherence')
%% coherence^2 between bilateral HbT during Awake
subplot(3,2,1);
s1 = semilogx(data.C57BL6J.Awake.gammaBandPower.LH.meanf,data.C57BL6J.Awake.gammaBandPower.LH.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.C57BL6J.Awake.gammaBandPower.LH.meanf,data.C57BL6J.Awake.gammaBandPower.LH.meanC + data.C57BL6J.Awake.gammaBandPower.LH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.C57BL6J.Awake.gammaBandPower.LH.meanf,data.C57BL6J.Awake.gammaBandPower.LH.meanC - data.C57BL6J.Awake.gammaBandPower.LH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
s2 = semilogx(data.Blank_SAP.Awake.gammaBandPower.LH.meanf,data.Blank_SAP.Awake.gammaBandPower.LH.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.Awake.gammaBandPower.LH.meanf,data.Blank_SAP.Awake.gammaBandPower.LH.meanC + data.Blank_SAP.Awake.gammaBandPower.LH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Awake.gammaBandPower.LH.meanf,data.Blank_SAP.Awake.gammaBandPower.LH.meanC - data.Blank_SAP.Awake.gammaBandPower.LH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
s3 = semilogx(data.SSP_SAP.Awake.gammaBandPower.LH.meanf,data.SSP_SAP.Awake.gammaBandPower.LH.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.Awake.gammaBandPower.LH.meanf,data.SSP_SAP.Awake.gammaBandPower.LH.meanC + data.SSP_SAP.Awake.gammaBandPower.LH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Awake.gammaBandPower.LH.meanf,data.SSP_SAP.Awake.gammaBandPower.LH.meanC - data.SSP_SAP.Awake.gammaBandPower.LH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] LH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
% ylim([0,1])
legend([s1,s2,s3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Awake
subplot(3,2,2);
semilogx(data.C57BL6J.Awake.gammaBandPower.RH.meanf,data.C57BL6J.Awake.gammaBandPower.RH.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.C57BL6J.Awake.gammaBandPower.RH.meanf,data.C57BL6J.Awake.gammaBandPower.RH.meanC + data.C57BL6J.Awake.gammaBandPower.RH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.C57BL6J.Awake.gammaBandPower.RH.meanf,data.C57BL6J.Awake.gammaBandPower.RH.meanC - data.C57BL6J.Awake.gammaBandPower.RH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Awake.gammaBandPower.RH.meanf,data.Blank_SAP.Awake.gammaBandPower.RH.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.Awake.gammaBandPower.RH.meanf,data.Blank_SAP.Awake.gammaBandPower.RH.meanC + data.Blank_SAP.Awake.gammaBandPower.RH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Awake.gammaBandPower.RH.meanf,data.Blank_SAP.Awake.gammaBandPower.RH.meanC - data.Blank_SAP.Awake.gammaBandPower.RH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Awake.gammaBandPower.RH.meanf,data.SSP_SAP.Awake.gammaBandPower.RH.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.Awake.gammaBandPower.RH.meanf,data.SSP_SAP.Awake.gammaBandPower.RH.meanC + data.SSP_SAP.Awake.gammaBandPower.RH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Awake.gammaBandPower.RH.meanf,data.SSP_SAP.Awake.gammaBandPower.RH.meanC - data.SSP_SAP.Awake.gammaBandPower.RH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] RH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(3,2,3);
semilogx(data.C57BL6J.Sleep.gammaBandPower.LH.meanf,data.C57BL6J.Sleep.gammaBandPower.LH.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.C57BL6J.Sleep.gammaBandPower.LH.meanf,data.C57BL6J.Sleep.gammaBandPower.LH.meanC + data.C57BL6J.Sleep.gammaBandPower.LH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.C57BL6J.Sleep.gammaBandPower.LH.meanf,data.C57BL6J.Sleep.gammaBandPower.LH.meanC - data.C57BL6J.Sleep.gammaBandPower.LH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Sleep.gammaBandPower.LH.meanf,data.Blank_SAP.Sleep.gammaBandPower.LH.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.Sleep.gammaBandPower.LH.meanf,data.Blank_SAP.Sleep.gammaBandPower.LH.meanC + data.Blank_SAP.Sleep.gammaBandPower.LH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Sleep.gammaBandPower.LH.meanf,data.Blank_SAP.Sleep.gammaBandPower.LH.meanC - data.Blank_SAP.Sleep.gammaBandPower.LH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Sleep.gammaBandPower.LH.meanf,data.SSP_SAP.Sleep.gammaBandPower.LH.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.Sleep.gammaBandPower.LH.meanf,data.SSP_SAP.Sleep.gammaBandPower.LH.meanC + data.SSP_SAP.Sleep.gammaBandPower.LH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Sleep.gammaBandPower.LH.meanf,data.SSP_SAP.Sleep.gammaBandPower.LH.meanC - data.SSP_SAP.Sleep.gammaBandPower.LH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] LH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(3,2,4);
semilogx(data.C57BL6J.Sleep.gammaBandPower.RH.meanf,data.C57BL6J.Sleep.gammaBandPower.RH.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.C57BL6J.Sleep.gammaBandPower.RH.meanf,data.C57BL6J.Sleep.gammaBandPower.RH.meanC + data.C57BL6J.Sleep.gammaBandPower.RH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.C57BL6J.Sleep.gammaBandPower.RH.meanf,data.C57BL6J.Sleep.gammaBandPower.RH.meanC - data.C57BL6J.Sleep.gammaBandPower.RH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Sleep.gammaBandPower.RH.meanf,data.Blank_SAP.Sleep.gammaBandPower.RH.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.Sleep.gammaBandPower.RH.meanf,data.Blank_SAP.Sleep.gammaBandPower.RH.meanC + data.Blank_SAP.Sleep.gammaBandPower.RH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Sleep.gammaBandPower.RH.meanf,data.Blank_SAP.Sleep.gammaBandPower.RH.meanC - data.Blank_SAP.Sleep.gammaBandPower.RH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Sleep.gammaBandPower.RH.meanf,data.SSP_SAP.Sleep.gammaBandPower.RH.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.Sleep.gammaBandPower.RH.meanf,data.SSP_SAP.Sleep.gammaBandPower.RH.meanC + data.SSP_SAP.Sleep.gammaBandPower.RH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Sleep.gammaBandPower.RH.meanf,data.SSP_SAP.Sleep.gammaBandPower.RH.meanC - data.SSP_SAP.Sleep.gammaBandPower.RH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] RH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(3,2,5);
semilogx(data.C57BL6J.All.gammaBandPower.LH.meanf,data.C57BL6J.All.gammaBandPower.LH.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.C57BL6J.All.gammaBandPower.LH.meanf,data.C57BL6J.All.gammaBandPower.LH.meanC + data.C57BL6J.All.gammaBandPower.LH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.C57BL6J.All.gammaBandPower.LH.meanf,data.C57BL6J.All.gammaBandPower.LH.meanC - data.C57BL6J.All.gammaBandPower.LH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.All.gammaBandPower.LH.meanf,data.Blank_SAP.All.gammaBandPower.LH.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.All.gammaBandPower.LH.meanf,data.Blank_SAP.All.gammaBandPower.LH.meanC + data.Blank_SAP.All.gammaBandPower.LH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.All.gammaBandPower.LH.meanf,data.Blank_SAP.All.gammaBandPower.LH.meanC - data.Blank_SAP.All.gammaBandPower.LH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.All.gammaBandPower.LH.meanf,data.SSP_SAP.All.gammaBandPower.LH.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.All.gammaBandPower.LH.meanf,data.SSP_SAP.All.gammaBandPower.LH.meanC + data.SSP_SAP.All.gammaBandPower.LH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.All.gammaBandPower.LH.meanf,data.SSP_SAP.All.gammaBandPower.LH.meanC - data.SSP_SAP.All.gammaBandPower.LH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] LH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(3,2,6);
semilogx(data.C57BL6J.All.gammaBandPower.RH.meanf,data.C57BL6J.All.gammaBandPower.RH.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.C57BL6J.All.gammaBandPower.RH.meanf,data.C57BL6J.All.gammaBandPower.RH.meanC + data.C57BL6J.All.gammaBandPower.RH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.C57BL6J.All.gammaBandPower.RH.meanf,data.C57BL6J.All.gammaBandPower.RH.meanC - data.C57BL6J.All.gammaBandPower.RH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.All.gammaBandPower.RH.meanf,data.Blank_SAP.All.gammaBandPower.RH.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.All.gammaBandPower.RH.meanf,data.Blank_SAP.All.gammaBandPower.RH.meanC + data.Blank_SAP.All.gammaBandPower.RH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.All.gammaBandPower.RH.meanf,data.Blank_SAP.All.gammaBandPower.RH.meanC - data.Blank_SAP.All.gammaBandPower.RH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.All.gammaBandPower.RH.meanf,data.SSP_SAP.All.gammaBandPower.RH.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.All.gammaBandPower.RH.meanf,data.SSP_SAP.All.gammaBandPower.RH.meanC + data.SSP_SAP.All.gammaBandPower.RH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.All.gammaBandPower.RH.meanf,data.SSP_SAP.All.gammaBandPower.RH.meanC - data.SSP_SAP.All.gammaBandPower.RH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] RH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'NeuralHemoCoherence_Gamma']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'NeuralHemoCoherence_Gamma'])
    %% statistical diary
    diaryFile = [dirpath 'NeuralHemoCoherence_Gamma_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % Awake stats
    disp('======================================================================================================================')
    disp('GLME statistics for gamma-HbT Coherence^2 for Awake data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.C001.Awake.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.C01.Awake.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.C05.Awake.Stats)
    % Sleep stats
    disp('======================================================================================================================')
    disp('GLME statistics for gamma-HbT Coherence^2 for Sleep data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.C001.Sleep.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.C01.Sleep.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.C05.Sleep.Stats)
    % All stats
    disp('======================================================================================================================')
    disp('GLME statistics for gamma-HbT Coherence^2 for All data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.C001.All.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.C01.All.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.C05.All.Stats)
    diary off
end

end
