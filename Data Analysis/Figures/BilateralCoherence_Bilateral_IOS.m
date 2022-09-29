function [] = BilateralCoherence_Bilateral_IOS(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% set-up and process data
resultsStruct = 'Results_BilatCoher';
load(resultsStruct);
expGroups = {'Naive','SSP-SAP','Blank-SAP'};
setName = 'IOS Set A';
animalIDs.all = {};
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs.all = horzcat(animalIDs.all,{folderList.name});
    animalIDs.(strrep(expGroups{1,aa},'-','_')) = {folderList.name};
end
treatments = {'Naive','SSP_SAP','Blank_SAP'};
behavFields = {'Rest','NREM','REM','Awake','Sleep','All'};
dataTypes = {'CBV_HbT','gammaBandPower'};
%% average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,aa};
    % recognize treatment based on animal group
    if ismember(animalID,animalIDs.Naive) == true
        treatment = 'Naive';
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
                data.(treatment).(behavField).(dataType).C = [];
                data.(treatment).(behavField).(dataType).f = [];
                data.(treatment).(behavField).(dataType).confC = [];
                data.(treatment).(behavField).(dataType).animalID = {};
                data.(treatment).(behavField).(dataType).treatment = {};
                data.(treatment).(behavField).(dataType).freqC001 = {};
                data.(treatment).(behavField).(dataType).freqC01 = {};
                data.(treatment).(behavField).(dataType).freqC05 = {};
            end
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(Results_BilatCoher.(animalID).(behavField).(dataType).C) == false
                % concatenate C/f for existing data - exclude any empty sets
                data.(treatment).(behavField).(dataType).C = cat(2,data.(treatment).(behavField).(dataType).C,Results_BilatCoher.(animalID).(behavField).(dataType).C.^2);
                data.(treatment).(behavField).(dataType).f = cat(1,data.(treatment).(behavField).(dataType).f,Results_BilatCoher.(animalID).(behavField).(dataType).f);
                data.(treatment).(behavField).(dataType).confC = cat(1,data.(treatment).(behavField).(dataType).confC,Results_BilatCoher.(animalID).(behavField).(dataType).confC);
                data.(treatment).(behavField).(dataType).animalID = cat(1,data.(treatment).(behavField).(dataType).animalID,animalID);
                data.(treatment).(behavField).(dataType).treatment = cat(1,data.(treatment).(behavField).(dataType).treatment,treatment);
                data.(treatment).(behavField).(dataType).freqC001 = cat(1,data.(treatment).(behavField).(dataType).freqC001,'C001');
                data.(treatment).(behavField).(dataType).freqC01 = cat(1,data.(treatment).(behavField).(dataType).freqC01,'C01');
                data.(treatment).(behavField).(dataType).freqC05 = cat(1,data.(treatment).(behavField).(dataType).freqC05,'C05');
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
            data.(treatment).(behavField).(dataType).meanC = mean(data.(treatment).(behavField).(dataType).C,2);
            data.(treatment).(behavField).(dataType).stdC = std(data.(treatment).(behavField).(dataType).C,0,2);
            data.(treatment).(behavField).(dataType).meanf = mean(data.(treatment).(behavField).(dataType).f,1);
            data.(treatment).(behavField).(dataType).maxConfC = geomean(data.(treatment).(behavField).(dataType).confC);
            data.(treatment).(behavField).(dataType).maxConfC_Y = ones(length(data.(treatment).(behavField).(dataType).meanf),1)*data.(treatment).(behavField).(dataType).maxConfC;
        end
    end
end
%% find Hz peaks in coherence
treatments2 = {'SSP_SAP','Blank_SAP'};
behavFields2 = {'Awake','Sleep','All'};
for qq = 1:length(treatments2)
    treatment = treatments2{1,qq};
    for ee = 1:length(behavFields2)
        behavField = behavFields2{1,ee};
        for ff = 1:length(dataTypes)
            dataType = dataTypes{1,ff};
            for gg = 1:size(data.(treatment).(behavField).(dataType).C,2)
                F = round(data.(treatment).(behavField).(dataType).f(gg,:),3);
                C = data.(treatment).(behavField).(dataType).C(:,gg);
                index001 = find(F == 0.01);
                index01 = find(F == 0.1);
                index05 = find(F == 0.5);
                data.(treatment).(behavField).(dataType).C001(gg,1) = mean(C(1:index001(1)));
                data.(treatment).(behavField).(dataType).C01(gg,1) = mean(C(index001(1) + 1:index01(1)));
                data.(treatment).(behavField).(dataType).C05(gg,1) = mean(C(index01(1) + 1:index05(1)));
            end
        end
    end
end
%% statistics - generalized linear mixed effects model
freqs = {'C001','C01','C05'};
for aa = 1:length(freqs)
    freq = freqs{1,aa};
    for bb = 1:length(behavFields2)
        behavField = behavFields2{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % statistics - generalized linear mixed effects model
            Stats.(dataType).(behavField).(freq).tableSize = cat(1,data.Blank_SAP.(behavField).(dataType).(freq),data.SSP_SAP.(behavField).(dataType).(freq));
            Stats.(dataType).(behavField).(freq).Table = table('Size',[size(Stats.(dataType).(behavField).(freq).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Treatment','Coherence'});
            Stats.(dataType).(behavField).(freq).Table.Mouse = cat(1,data.Blank_SAP.(behavField).(dataType).animalID,data.SSP_SAP.(behavField).(dataType).animalID);
            Stats.(dataType).(behavField).(freq).Table.Treatment = cat(1,data.Blank_SAP.(behavField).(dataType).treatment,data.SSP_SAP.(behavField).(dataType).treatment);
            Stats.(dataType).(behavField).(freq).Table.Coherence = cat(1,data.Blank_SAP.(behavField).(dataType).(freq),data.SSP_SAP.(behavField).(dataType).(freq));
            Stats.(dataType).(behavField).(freq).FitFormula = 'Coherence ~ 1 + Treatment + (1|Mouse)';
            Stats.(dataType).(behavField).(freq).Stats = fitglme(Stats.(dataType).(behavField).(freq).Table,Stats.(dataType).(behavField).(freq).FitFormula);
        end
    end
end
%%
% for bb = 1:length(behavFields)
%     behavField = behavFields{1,bb};
%     for cc = 1:length(dataTypes)
%         dataType = dataTypes{1,cc};
%         % statistics - generalized linear mixed effects model
%         Stats.(dataType).(behavField).tableSize = cat(1,data.Blank_SAP.(behavField).(dataType).RH.C001,data.Blank_SAP.(behavField).(dataType).RH.C01,data.Blank_SAP.(behavField).(dataType).RH.C05,...
%             data.SSP_SAP.(behavField).(dataType).RH.C001,data.SSP_SAP.(behavField).(dataType).RH.C01,data.SSP_SAP.(behavField).(dataType).RH.C05);
%         Stats.(dataType).(behavField).Table = table('Size',[size(Stats.(dataType).(behavField).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Treatment','Frequency','Coherence'});
%         Stats.(dataType).(behavField).Table.Mouse = cat(1,data.Blank_SAP.(behavField).(dataType).animalID,data.Blank_SAP.(behavField).(dataType).animalID,data.Blank_SAP.(behavField).(dataType).animalID,...
%             data.SSP_SAP.(behavField).(dataType).animalID,data.SSP_SAP.(behavField).(dataType).animalID,data.SSP_SAP.(behavField).(dataType).animalID);
%         Stats.(dataType).(behavField).Table.Treatment = cat(1,data.Blank_SAP.(behavField).(dataType).treatment,data.Blank_SAP.(behavField).(dataType).treatment,data.Blank_SAP.(behavField).(dataType).treatment,...,
%             data.SSP_SAP.(behavField).(dataType).treatment,data.SSP_SAP.(behavField).(dataType).treatment,data.SSP_SAP.(behavField).(dataType).treatment);
%         Stats.(dataType).(behavField).Table.Frequency = cat(1,data.Blank_SAP.(behavField).(dataType).freqC001,data.Blank_SAP.(behavField).(dataType).freqC01,data.Blank_SAP.(behavField).(dataType).freqC05,...
%             data.SSP_SAP.(behavField).(dataType).freqC001,data.SSP_SAP.(behavField).(dataType).freqC01,data.SSP_SAP.(behavField).(dataType).freqC05);
%         Stats.(dataType).(behavField).Table.Coherence = cat(1,data.Blank_SAP.(behavField).(dataType).RH.C001,data.Blank_SAP.(behavField).(dataType).RH.C01,data.Blank_SAP.(behavField).(dataType).RH.C05,...
%             data.SSP_SAP.(behavField).(dataType).RH.C001,data.SSP_SAP.(behavField).(dataType).RH.C01,data.SSP_SAP.(behavField).(dataType).RH.C05);
%         Stats.(dataType).(behavField).FitFormula = 'Coherence ~ 1 + Treatment + Frequency + Frequency*Treatment + (1|Mouse)';
%         Stats.(dataType).(behavField).Stats = fitglme(Stats.(dataType).(behavField).Table,Stats.(dataType).(behavField).FitFormula);
%     end
% end
%% average HbT coherence
summaryFigure1 = figure;
sgtitle('Bilateral Coherence \DeltaHbT')
%% coherence^2 between bilateral HbT during rest
subplot(2,3,1);
p1 = semilogx(data.Naive.Rest.CBV_HbT.meanf,data.Naive.Rest.CBV_HbT.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Naive.Rest.CBV_HbT.meanf,data.Naive.Rest.CBV_HbT.meanC + data.Naive.Rest.CBV_HbT.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Naive.Rest.CBV_HbT.meanf,data.Naive.Rest.CBV_HbT.meanC - data.Naive.Rest.CBV_HbT.stdC,'color',colors('sapphire'),'LineWidth',0.5);
p2 = semilogx(data.Blank_SAP.Rest.CBV_HbT.meanf,data.Blank_SAP.Rest.CBV_HbT.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.Rest.CBV_HbT.meanf,data.Blank_SAP.Rest.CBV_HbT.meanC + data.Blank_SAP.Rest.CBV_HbT.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Rest.CBV_HbT.meanf,data.Blank_SAP.Rest.CBV_HbT.meanC - data.Blank_SAP.Rest.CBV_HbT.stdC,'color',colors('north texas green'),'LineWidth',0.5);
p3 = semilogx(data.SSP_SAP.Rest.CBV_HbT.meanf,data.SSP_SAP.Rest.CBV_HbT.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.Rest.CBV_HbT.meanf,data.SSP_SAP.Rest.CBV_HbT.meanC + data.SSP_SAP.Rest.CBV_HbT.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Rest.CBV_HbT.meanf,data.SSP_SAP.Rest.CBV_HbT.meanC - data.SSP_SAP.Rest.CBV_HbT.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Rest] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/10,0.5])
ylim([0,1])
set(gca,'box','off')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
%% coherence^2 between bilateral HbT during NREM
subplot(2,3,2);
semilogx(data.Naive.NREM.CBV_HbT.meanf,data.Naive.NREM.CBV_HbT.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Naive.NREM.CBV_HbT.meanf,data.Naive.NREM.CBV_HbT.meanC + data.Naive.NREM.CBV_HbT.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Naive.NREM.CBV_HbT.meanf,data.Naive.NREM.CBV_HbT.meanC - data.Naive.NREM.CBV_HbT.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.NREM.CBV_HbT.meanf,data.Blank_SAP.NREM.CBV_HbT.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.NREM.CBV_HbT.meanf,data.Blank_SAP.NREM.CBV_HbT.meanC + data.Blank_SAP.NREM.CBV_HbT.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.NREM.CBV_HbT.meanf,data.Blank_SAP.NREM.CBV_HbT.meanC - data.Blank_SAP.NREM.CBV_HbT.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.NREM.CBV_HbT.meanf,data.SSP_SAP.NREM.CBV_HbT.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.NREM.CBV_HbT.meanf,data.SSP_SAP.NREM.CBV_HbT.meanC + data.SSP_SAP.NREM.CBV_HbT.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.NREM.CBV_HbT.meanf,data.SSP_SAP.NREM.CBV_HbT.meanC - data.SSP_SAP.NREM.CBV_HbT.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[NREM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/30,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during REM
subplot(2,3,3);
semilogx(data.Naive.REM.CBV_HbT.meanf,data.Naive.REM.CBV_HbT.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Naive.REM.CBV_HbT.meanf,data.Naive.REM.CBV_HbT.meanC + data.Naive.REM.CBV_HbT.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Naive.REM.CBV_HbT.meanf,data.Naive.REM.CBV_HbT.meanC - data.Naive.REM.CBV_HbT.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.REM.CBV_HbT.meanf,data.Blank_SAP.REM.CBV_HbT.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.REM.CBV_HbT.meanf,data.Blank_SAP.REM.CBV_HbT.meanC + data.Blank_SAP.REM.CBV_HbT.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.REM.CBV_HbT.meanf,data.Blank_SAP.REM.CBV_HbT.meanC - data.Blank_SAP.REM.CBV_HbT.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.REM.CBV_HbT.meanf,data.SSP_SAP.REM.CBV_HbT.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.REM.CBV_HbT.meanf,data.SSP_SAP.REM.CBV_HbT.meanC + data.SSP_SAP.REM.CBV_HbT.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.REM.CBV_HbT.meanf,data.SSP_SAP.REM.CBV_HbT.meanC - data.SSP_SAP.REM.CBV_HbT.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[REM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/60,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Awake
subplot(2,3,4);
semilogx(data.Naive.Awake.CBV_HbT.meanf,data.Naive.Awake.CBV_HbT.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Naive.Awake.CBV_HbT.meanf,data.Naive.Awake.CBV_HbT.meanC + data.Naive.Awake.CBV_HbT.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Naive.Awake.CBV_HbT.meanf,data.Naive.Awake.CBV_HbT.meanC - data.Naive.Awake.CBV_HbT.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Awake.CBV_HbT.meanf,data.Blank_SAP.Awake.CBV_HbT.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.Awake.CBV_HbT.meanf,data.Blank_SAP.Awake.CBV_HbT.meanC + data.Blank_SAP.Awake.CBV_HbT.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Awake.CBV_HbT.meanf,data.Blank_SAP.Awake.CBV_HbT.meanC - data.Blank_SAP.Awake.CBV_HbT.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Awake.CBV_HbT.meanf,data.SSP_SAP.Awake.CBV_HbT.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.Awake.CBV_HbT.meanf,data.SSP_SAP.Awake.CBV_HbT.meanC + data.SSP_SAP.Awake.CBV_HbT.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Awake.CBV_HbT.meanf,data.SSP_SAP.Awake.CBV_HbT.meanC - data.SSP_SAP.Awake.CBV_HbT.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(2,3,5);
semilogx(data.Naive.Sleep.CBV_HbT.meanf,data.Naive.Sleep.CBV_HbT.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Naive.Sleep.CBV_HbT.meanf,data.Naive.Sleep.CBV_HbT.meanC + data.Naive.Sleep.CBV_HbT.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Naive.Sleep.CBV_HbT.meanf,data.Naive.Sleep.CBV_HbT.meanC - data.Naive.Sleep.CBV_HbT.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Sleep.CBV_HbT.meanf,data.Blank_SAP.Sleep.CBV_HbT.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.Sleep.CBV_HbT.meanf,data.Blank_SAP.Sleep.CBV_HbT.meanC + data.Blank_SAP.Sleep.CBV_HbT.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Sleep.CBV_HbT.meanf,data.Blank_SAP.Sleep.CBV_HbT.meanC - data.Blank_SAP.Sleep.CBV_HbT.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Sleep.CBV_HbT.meanf,data.SSP_SAP.Sleep.CBV_HbT.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.Sleep.CBV_HbT.meanf,data.SSP_SAP.Sleep.CBV_HbT.meanC + data.SSP_SAP.Sleep.CBV_HbT.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Sleep.CBV_HbT.meanf,data.SSP_SAP.Sleep.CBV_HbT.meanC - data.SSP_SAP.Sleep.CBV_HbT.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(2,3,6);
semilogx(data.Naive.All.CBV_HbT.meanf,data.Naive.All.CBV_HbT.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Naive.All.CBV_HbT.meanf,data.Naive.All.CBV_HbT.meanC + data.Naive.All.CBV_HbT.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Naive.All.CBV_HbT.meanf,data.Naive.All.CBV_HbT.meanC - data.Naive.All.CBV_HbT.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.All.CBV_HbT.meanf,data.Blank_SAP.All.CBV_HbT.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.All.CBV_HbT.meanf,data.Blank_SAP.All.CBV_HbT.meanC + data.Blank_SAP.All.CBV_HbT.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.All.CBV_HbT.meanf,data.Blank_SAP.All.CBV_HbT.meanC - data.Blank_SAP.All.CBV_HbT.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.All.CBV_HbT.meanf,data.SSP_SAP.All.CBV_HbT.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.All.CBV_HbT.meanf,data.SSP_SAP.All.CBV_HbT.meanC + data.SSP_SAP.All.CBV_HbT.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.All.CBV_HbT.meanf,data.SSP_SAP.All.CBV_HbT.meanC - data.SSP_SAP.All.CBV_HbT.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Bilateral Coherence - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'AverageBilateralCoherence_HbT']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'AverageBilateralCoherence_HbT'])
end
%% individual HbT coherence
summaryFigure2 = figure;
sgtitle('Bilateral Coherence \DeltaHbT - individual animals')
%% coherence^2 between bilateral HbT during rest
subplot(2,3,1);
% Naive
for aa = 1:size(data.Naive.Rest.CBV_HbT.C,2)
    semilogx(data.Naive.Rest.CBV_HbT.meanf,data.Naive.Rest.CBV_HbT.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Rest.CBV_HbT.C,2)
    semilogx(data.Blank_SAP.Rest.CBV_HbT.meanf,data.Blank_SAP.Rest.CBV_HbT.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Rest.CBV_HbT.C,2)
    semilogx(data.SSP_SAP.Rest.CBV_HbT.meanf,data.SSP_SAP.Rest.CBV_HbT.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Rest] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/10,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during NREM
subplot(2,3,2);
% Naive
for aa = 1:size(data.Naive.NREM.CBV_HbT.C,2)
    semilogx(data.Naive.NREM.CBV_HbT.meanf,data.Naive.NREM.CBV_HbT.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.NREM.CBV_HbT.C,2)
    semilogx(data.Blank_SAP.NREM.CBV_HbT.meanf,data.Blank_SAP.NREM.CBV_HbT.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.NREM.CBV_HbT.C,2)
    semilogx(data.SSP_SAP.NREM.CBV_HbT.meanf,data.SSP_SAP.NREM.CBV_HbT.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[NREM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/30,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during REM
subplot(2,3,3);
% Naive
for aa = 1:size(data.Naive.REM.CBV_HbT.C,2)
    semilogx(data.Naive.REM.CBV_HbT.meanf,data.Naive.REM.CBV_HbT.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.REM.CBV_HbT.C,2)
    semilogx(data.Blank_SAP.REM.CBV_HbT.meanf,data.Blank_SAP.REM.CBV_HbT.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.REM.CBV_HbT.C,2)
    semilogx(data.SSP_SAP.REM.CBV_HbT.meanf,data.SSP_SAP.REM.CBV_HbT.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[REM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.01,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Awake
subplot(2,3,4);
% Naive
for aa = 1:size(data.Naive.Awake.CBV_HbT.C,2)
    semilogx(data.Naive.Awake.CBV_HbT.meanf,data.Naive.Awake.CBV_HbT.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Awake.CBV_HbT.C,2)
    semilogx(data.Blank_SAP.Awake.CBV_HbT.meanf,data.Blank_SAP.Awake.CBV_HbT.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Awake.CBV_HbT.C,2)
    semilogx(data.SSP_SAP.Awake.CBV_HbT.meanf,data.SSP_SAP.Awake.CBV_HbT.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(2,3,5);
% Naive
for aa = 1:size(data.Naive.Sleep.CBV_HbT.C,2)
    semilogx(data.Naive.Sleep.CBV_HbT.meanf,data.Naive.Sleep.CBV_HbT.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Sleep.CBV_HbT.C,2)
    semilogx(data.Blank_SAP.Sleep.CBV_HbT.meanf,data.Blank_SAP.Sleep.CBV_HbT.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Sleep.CBV_HbT.C,2)
    semilogx(data.SSP_SAP.Sleep.CBV_HbT.meanf,data.SSP_SAP.Sleep.CBV_HbT.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(2,3,6);
% Naive
for aa = 1:size(data.Naive.All.CBV_HbT.C,2)
    semilogx(data.Naive.All.CBV_HbT.meanf,data.Naive.All.CBV_HbT.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.All.CBV_HbT.C,2)
    semilogx(data.Blank_SAP.All.CBV_HbT.meanf,data.Blank_SAP.All.CBV_HbT.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.All.CBV_HbT.C,2)
    semilogx(data.SSP_SAP.All.CBV_HbT.meanf,data.SSP_SAP.All.CBV_HbT.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Bilateral Coherence - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'IndividualCoherence_HbT']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'IndividualCoherence_HbT'])
end
%% average gamma-band coherence
summaryFigure3 = figure;
sgtitle('Bilateral Coherence Gamma-band [30-100 Hz] (envelope)')
%% coherence^2 between bilateral gamma-band during Rest
subplot(2,3,1);
p1 = semilogx(data.Naive.Rest.gammaBandPower.meanf,data.Naive.Rest.gammaBandPower.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Naive.Rest.gammaBandPower.meanf,data.Naive.Rest.gammaBandPower.meanC + data.Naive.Rest.gammaBandPower.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Naive.Rest.gammaBandPower.meanf,data.Naive.Rest.gammaBandPower.meanC - data.Naive.Rest.gammaBandPower.stdC,'color',colors('sapphire'),'LineWidth',0.5);
p2 = semilogx(data.Blank_SAP.Rest.gammaBandPower.meanf,data.Blank_SAP.Rest.gammaBandPower.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.Rest.gammaBandPower.meanf,data.Blank_SAP.Rest.gammaBandPower.meanC + data.Blank_SAP.Rest.gammaBandPower.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Rest.gammaBandPower.meanf,data.Blank_SAP.Rest.gammaBandPower.meanC - data.Blank_SAP.Rest.gammaBandPower.stdC,'color',colors('north texas green'),'LineWidth',0.5);
p3 = semilogx(data.SSP_SAP.Rest.gammaBandPower.meanf,data.SSP_SAP.Rest.gammaBandPower.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.Rest.gammaBandPower.meanf,data.SSP_SAP.Rest.gammaBandPower.meanC + data.SSP_SAP.Rest.gammaBandPower.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Rest.gammaBandPower.meanf,data.SSP_SAP.Rest.gammaBandPower.meanC - data.SSP_SAP.Rest.gammaBandPower.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Rest] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([1/10,0.5])
ylim([0,1])
set(gca,'box','off')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
%% coherence^2 between bilateral gamma-band during NREM
subplot(2,3,2);
semilogx(data.Naive.NREM.gammaBandPower.meanf,data.Naive.NREM.gammaBandPower.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Naive.NREM.gammaBandPower.meanf,data.Naive.NREM.gammaBandPower.meanC + data.Naive.NREM.gammaBandPower.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Naive.NREM.gammaBandPower.meanf,data.Naive.NREM.gammaBandPower.meanC - data.Naive.NREM.gammaBandPower.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.NREM.gammaBandPower.meanf,data.Blank_SAP.NREM.gammaBandPower.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.NREM.gammaBandPower.meanf,data.Blank_SAP.NREM.gammaBandPower.meanC + data.Blank_SAP.NREM.gammaBandPower.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.NREM.gammaBandPower.meanf,data.Blank_SAP.NREM.gammaBandPower.meanC - data.Blank_SAP.NREM.gammaBandPower.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.NREM.gammaBandPower.meanf,data.SSP_SAP.NREM.gammaBandPower.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.NREM.gammaBandPower.meanf,data.SSP_SAP.NREM.gammaBandPower.meanC + data.SSP_SAP.NREM.gammaBandPower.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.NREM.gammaBandPower.meanf,data.SSP_SAP.NREM.gammaBandPower.meanC - data.SSP_SAP.NREM.gammaBandPower.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[NREM] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([1/30,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during REM
subplot(2,3,3);
semilogx(data.Naive.REM.gammaBandPower.meanf,data.Naive.REM.gammaBandPower.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Naive.REM.gammaBandPower.meanf,data.Naive.REM.gammaBandPower.meanC + data.Naive.REM.gammaBandPower.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Naive.REM.gammaBandPower.meanf,data.Naive.REM.gammaBandPower.meanC - data.Naive.REM.gammaBandPower.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.REM.gammaBandPower.meanf,data.Blank_SAP.REM.gammaBandPower.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.REM.gammaBandPower.meanf,data.Blank_SAP.REM.gammaBandPower.meanC + data.Blank_SAP.REM.gammaBandPower.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.REM.gammaBandPower.meanf,data.Blank_SAP.REM.gammaBandPower.meanC - data.Blank_SAP.REM.gammaBandPower.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.REM.gammaBandPower.meanf,data.SSP_SAP.REM.gammaBandPower.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.REM.gammaBandPower.meanf,data.SSP_SAP.REM.gammaBandPower.meanC + data.SSP_SAP.REM.gammaBandPower.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.REM.gammaBandPower.meanf,data.SSP_SAP.REM.gammaBandPower.meanC - data.SSP_SAP.REM.gammaBandPower.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[REM] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([1/60,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during Awake
subplot(2,3,4);
semilogx(data.Naive.Awake.gammaBandPower.meanf,data.Naive.Awake.gammaBandPower.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Naive.Awake.gammaBandPower.meanf,data.Naive.Awake.gammaBandPower.meanC + data.Naive.Awake.gammaBandPower.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Naive.Awake.gammaBandPower.meanf,data.Naive.Awake.gammaBandPower.meanC - data.Naive.Awake.gammaBandPower.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Awake.gammaBandPower.meanf,data.Blank_SAP.Awake.gammaBandPower.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.Awake.gammaBandPower.meanf,data.Blank_SAP.Awake.gammaBandPower.meanC + data.Blank_SAP.Awake.gammaBandPower.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Awake.gammaBandPower.meanf,data.Blank_SAP.Awake.gammaBandPower.meanC - data.Blank_SAP.Awake.gammaBandPower.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Awake.gammaBandPower.meanf,data.SSP_SAP.Awake.gammaBandPower.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.Awake.gammaBandPower.meanf,data.SSP_SAP.Awake.gammaBandPower.meanC + data.SSP_SAP.Awake.gammaBandPower.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Awake.gammaBandPower.meanf,data.SSP_SAP.Awake.gammaBandPower.meanC - data.SSP_SAP.Awake.gammaBandPower.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during Sleep
subplot(2,3,5);
semilogx(data.Naive.Sleep.gammaBandPower.meanf,data.Naive.Sleep.gammaBandPower.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Naive.Sleep.gammaBandPower.meanf,data.Naive.Sleep.gammaBandPower.meanC + data.Naive.Sleep.gammaBandPower.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Naive.Sleep.gammaBandPower.meanf,data.Naive.Sleep.gammaBandPower.meanC - data.Naive.Sleep.gammaBandPower.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Sleep.gammaBandPower.meanf,data.Blank_SAP.Sleep.gammaBandPower.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.Sleep.gammaBandPower.meanf,data.Blank_SAP.Sleep.gammaBandPower.meanC + data.Blank_SAP.Sleep.gammaBandPower.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.Sleep.gammaBandPower.meanf,data.Blank_SAP.Sleep.gammaBandPower.meanC - data.Blank_SAP.Sleep.gammaBandPower.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Sleep.gammaBandPower.meanf,data.SSP_SAP.Sleep.gammaBandPower.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.Sleep.gammaBandPower.meanf,data.SSP_SAP.Sleep.gammaBandPower.meanC + data.SSP_SAP.Sleep.gammaBandPower.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.Sleep.gammaBandPower.meanf,data.SSP_SAP.Sleep.gammaBandPower.meanC - data.SSP_SAP.Sleep.gammaBandPower.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during All data
subplot(2,3,6);
semilogx(data.Naive.All.gammaBandPower.meanf,data.Naive.All.gammaBandPower.meanC,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Naive.All.gammaBandPower.meanf,data.Naive.All.gammaBandPower.meanC + data.Naive.All.gammaBandPower.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Naive.All.gammaBandPower.meanf,data.Naive.All.gammaBandPower.meanC - data.Naive.All.gammaBandPower.stdC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Blank_SAP.All.gammaBandPower.meanf,data.Blank_SAP.All.gammaBandPower.meanC,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.Blank_SAP.All.gammaBandPower.meanf,data.Blank_SAP.All.gammaBandPower.meanC + data.Blank_SAP.All.gammaBandPower.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.Blank_SAP.All.gammaBandPower.meanf,data.Blank_SAP.All.gammaBandPower.meanC - data.Blank_SAP.All.gammaBandPower.stdC,'color',colors('north texas green'),'LineWidth',0.5);
semilogx(data.SSP_SAP.All.gammaBandPower.meanf,data.SSP_SAP.All.gammaBandPower.meanC,'color',colors('electric purple'),'LineWidth',2);
semilogx(data.SSP_SAP.All.gammaBandPower.meanf,data.SSP_SAP.All.gammaBandPower.meanC + data.SSP_SAP.All.gammaBandPower.stdC,'color',colors('electric purple'),'LineWidth',0.5);
semilogx(data.SSP_SAP.All.gammaBandPower.meanf,data.SSP_SAP.All.gammaBandPower.meanC - data.SSP_SAP.All.gammaBandPower.stdC,'color',colors('electric purple'),'LineWidth',0.5);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Bilateral Coherence - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure3,[dirpath 'AverageBilateralCoherence_Gamma']);
    set(summaryFigure3,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'AverageBilateralCoherence_Gamma'])
end
%% individual gamma-band coherence
summaryFigure4 = figure;
sgtitle('Bilateral Coherence Gamma-band [30-100 Hz] - individual animals')
%% coherence^2 between bilateral gamma-band during rest
subplot(2,3,1);
% Naive
for aa = 1:size(data.Naive.Rest.gammaBandPower.C,2)
    semilogx(data.Naive.Rest.gammaBandPower.meanf,data.Naive.Rest.gammaBandPower.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Rest.gammaBandPower.C,2)
    semilogx(data.Blank_SAP.Rest.gammaBandPower.meanf,data.Blank_SAP.Rest.gammaBandPower.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Rest.gammaBandPower.C,2)
    semilogx(data.SSP_SAP.Rest.gammaBandPower.meanf,data.SSP_SAP.Rest.gammaBandPower.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Rest] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([1/10,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during NREM
subplot(2,3,2);
% Naive
for aa = 1:size(data.Naive.NREM.gammaBandPower.C,2)
    semilogx(data.Naive.NREM.gammaBandPower.meanf,data.Naive.NREM.gammaBandPower.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.NREM.gammaBandPower.C,2)
    semilogx(data.Blank_SAP.NREM.gammaBandPower.meanf,data.Blank_SAP.NREM.gammaBandPower.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.NREM.gammaBandPower.C,2)
    semilogx(data.SSP_SAP.NREM.gammaBandPower.meanf,data.SSP_SAP.NREM.gammaBandPower.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[NREM] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([1/30,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during REM
subplot(2,3,3);
% Naive
for aa = 1:size(data.Naive.REM.gammaBandPower.C,2)
    semilogx(data.Naive.REM.gammaBandPower.meanf,data.Naive.REM.gammaBandPower.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.REM.gammaBandPower.C,2)
    semilogx(data.Blank_SAP.REM.gammaBandPower.meanf,data.Blank_SAP.REM.gammaBandPower.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.REM.gammaBandPower.C,2)
    semilogx(data.SSP_SAP.REM.gammaBandPower.meanf,data.SSP_SAP.REM.gammaBandPower.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[REM] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([1/60,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during Awake
subplot(2,3,4);
% Naive
for aa = 1:size(data.Naive.Awake.gammaBandPower.C,2)
    semilogx(data.Naive.Awake.gammaBandPower.meanf,data.Naive.Awake.gammaBandPower.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Awake.gammaBandPower.C,2)
    semilogx(data.Blank_SAP.Awake.gammaBandPower.meanf,data.Blank_SAP.Awake.gammaBandPower.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Awake.gammaBandPower.C,2)
    semilogx(data.SSP_SAP.Awake.gammaBandPower.meanf,data.SSP_SAP.Awake.gammaBandPower.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during Sleep
subplot(2,3,5);
% Naive
for aa = 1:size(data.Naive.Sleep.gammaBandPower.C,2)
    semilogx(data.Naive.Sleep.gammaBandPower.meanf,data.Naive.Sleep.gammaBandPower.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Sleep.gammaBandPower.C,2)
    semilogx(data.Blank_SAP.Sleep.gammaBandPower.meanf,data.Blank_SAP.Sleep.gammaBandPower.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Sleep.gammaBandPower.C,2)
    semilogx(data.SSP_SAP.Sleep.gammaBandPower.meanf,data.SSP_SAP.Sleep.gammaBandPower.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during All data
subplot(2,3,6);
% Naive
for aa = 1:size(data.Naive.All.gammaBandPower.C,2)
    semilogx(data.Naive.All.gammaBandPower.meanf,data.Naive.All.gammaBandPower.C(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.All.gammaBandPower.C,2)
    semilogx(data.Blank_SAP.All.gammaBandPower.meanf,data.Blank_SAP.All.gammaBandPower.C(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.All.gammaBandPower.C,2)
    semilogx(data.SSP_SAP.All.gammaBandPower.meanf,data.SSP_SAP.All.gammaBandPower.C(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Bilateral Coherence - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure4,[dirpath 'IndividualBilateralCoherence_Gamma']);
    set(summaryFigure4,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'IndividualBilateralCoherence_Gamma'])
end
%% HbT and gamma-band coherence stats
summaryFigure5 = figure;
sgtitle('Bilateral Coherence Statics')
%% Alert HbT Stats
ax1 = subplot(2,3,1);
xInds = ones(1,length(animalIDs.Blank_SAP));
s1 = scatter(xInds*1,data.Blank_SAP.Awake.CBV_HbT.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(data.Blank_SAP.Awake.CBV_HbT.C001),std(data.Blank_SAP.Awake.CBV_HbT.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(xInds*2,data.SSP_SAP.Awake.CBV_HbT.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,mean(data.SSP_SAP.Awake.CBV_HbT.C001),std(data.Blank_SAP.Awake.CBV_HbT.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.Blank_SAP.Awake.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,mean(data.Blank_SAP.Awake.CBV_HbT.C01),std(data.Blank_SAP.Awake.CBV_HbT.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.SSP_SAP.Awake.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e4 = errorbar(4,mean(data.SSP_SAP.Awake.CBV_HbT.C01),std(data.Blank_SAP.Awake.CBV_HbT.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.Blank_SAP.Awake.CBV_HbT.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e5 = errorbar(5,mean(data.Blank_SAP.Awake.CBV_HbT.C05),std(data.Blank_SAP.Awake.CBV_HbT.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.SSP_SAP.Awake.CBV_HbT.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e6 = errorbar(6,mean(data.SSP_SAP.Awake.CBV_HbT.C05),std(data.Blank_SAP.Awake.CBV_HbT.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
% stat lines
plot([1,2],[1,1],'k');
text(1.5,1,'ns','FontSize',16)
plot([3,4],[1,1],'k');
text(3.5,1,'ns','FontSize',16)
plot([5,6],[1,1],'k');
text(5.5,1,'*','FontSize',16)
title({'Alert HbT',''})
ylabel('Average coherence^2')
legend([s1,s2],'Blank-SAP','SSP-SAP')
set(gca,'xtick',[1.5,3.5,5.5])
xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
axis square
axis tight
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% Asleep HbT stats
ax2 = subplot(2,3,2);
xInds2 = ones(1,length(data.Blank_SAP.Sleep.CBV_HbT.C001));
scatter(xInds2*1,data.Blank_SAP.Sleep.CBV_HbT.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(data.Blank_SAP.Sleep.CBV_HbT.C001),std(data.Blank_SAP.Sleep.CBV_HbT.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds2*2,data.SSP_SAP.Sleep.CBV_HbT.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,mean(data.SSP_SAP.Sleep.CBV_HbT.C001),std(data.Blank_SAP.Sleep.CBV_HbT.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds2*3,data.Blank_SAP.Sleep.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,mean(data.Blank_SAP.Sleep.CBV_HbT.C01),std(data.Blank_SAP.Sleep.CBV_HbT.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds2*4,data.SSP_SAP.Sleep.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e4 = errorbar(4,mean(data.SSP_SAP.Sleep.CBV_HbT.C01),std(data.Blank_SAP.Sleep.CBV_HbT.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds2*5,data.Blank_SAP.Sleep.CBV_HbT.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e5 = errorbar(5,mean(data.Blank_SAP.Sleep.CBV_HbT.C05),std(data.Blank_SAP.Sleep.CBV_HbT.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds2*6,data.SSP_SAP.Sleep.CBV_HbT.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e6 = errorbar(6,mean(data.SSP_SAP.Sleep.CBV_HbT.C05),std(data.Blank_SAP.Sleep.CBV_HbT.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
% stat lines
plot([1,2],[1,1],'k');
text(1.5,1,'ns','FontSize',16)
plot([3,4],[1,1],'k');
text(3.5,1,'*','FontSize',16)
plot([5,6],[1,1],'k');
text(5.5,1,'*','FontSize',16)
title({'Asleep HbT',''})
ylabel('Average coherence^2')
set(gca,'xtick',[1.5,3.5,5.5])
xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
axis square
axis tight
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% All HbT stats
ax3 = subplot(2,3,3);
xInds = ones(1,length(animalIDs.Blank_SAP));
scatter(xInds*1,data.Blank_SAP.All.CBV_HbT.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(data.Blank_SAP.All.CBV_HbT.C001),std(data.Blank_SAP.All.CBV_HbT.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.SSP_SAP.All.CBV_HbT.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,mean(data.SSP_SAP.All.CBV_HbT.C001),std(data.Blank_SAP.All.CBV_HbT.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.Blank_SAP.All.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,mean(data.Blank_SAP.All.CBV_HbT.C01),std(data.Blank_SAP.All.CBV_HbT.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.SSP_SAP.All.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e4 = errorbar(4,mean(data.SSP_SAP.All.CBV_HbT.C01),std(data.Blank_SAP.All.CBV_HbT.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.Blank_SAP.All.CBV_HbT.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e5 = errorbar(5,mean(data.Blank_SAP.All.CBV_HbT.C05),std(data.Blank_SAP.All.CBV_HbT.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.SSP_SAP.All.CBV_HbT.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e6 = errorbar(6,mean(data.SSP_SAP.All.CBV_HbT.C05),std(data.Blank_SAP.All.CBV_HbT.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
% stat lines
plot([1,2],[1,1],'k');
text(1.5,1,'ns','FontSize',16)
plot([3,4],[1,1],'k');
text(3.5,1,'ns','FontSize',16)
plot([5,6],[1,1],'k');
text(5.5,1,'***','FontSize',16)
title({'All HbT',''})
ylabel('Average coherence^2')
set(gca,'xtick',[1.5,3.5,5.5])
xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
axis square
axis tight
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% Alert Gamma-band power Stats
ax4 = subplot(2,3,4);
scatter(xInds*1,data.Blank_SAP.Awake.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(data.Blank_SAP.Awake.gammaBandPower.C001),std(data.Blank_SAP.Awake.gammaBandPower.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.SSP_SAP.Awake.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,mean(data.SSP_SAP.Awake.gammaBandPower.C001),std(data.Blank_SAP.Awake.gammaBandPower.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.Blank_SAP.Awake.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,mean(data.Blank_SAP.Awake.gammaBandPower.C01),std(data.Blank_SAP.Awake.gammaBandPower.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.SSP_SAP.Awake.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e4 = errorbar(4,mean(data.SSP_SAP.Awake.gammaBandPower.C01),std(data.Blank_SAP.Awake.gammaBandPower.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.Blank_SAP.Awake.gammaBandPower.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e5 = errorbar(5,mean(data.Blank_SAP.Awake.gammaBandPower.C05),std(data.Blank_SAP.Awake.gammaBandPower.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.SSP_SAP.Awake.gammaBandPower.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e6 = errorbar(6,mean(data.SSP_SAP.Awake.gammaBandPower.C05),std(data.Blank_SAP.Awake.gammaBandPower.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
% stat lines
plot([1,2],[1,1],'k');
text(1.5,1,'ns','FontSize',16)
plot([3,4],[1,1],'k');
text(3.5,1,'ns','FontSize',16)
plot([5,6],[1,1],'k');
text(5.5,1,'ns','FontSize',16)
title({'Alert Gamma-band power',''})
ylabel('Average coherence^2')
set(gca,'xtick',[1.5,3.5,5.5])
xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
axis square
axis tight
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% Asleep Gamma-band power stats
ax5 = subplot(2,3,5);
scatter(xInds2*1,data.Blank_SAP.Sleep.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(data.Blank_SAP.Sleep.gammaBandPower.C001),std(data.Blank_SAP.Sleep.gammaBandPower.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds2*2,data.SSP_SAP.Sleep.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,mean(data.SSP_SAP.Sleep.gammaBandPower.C001),std(data.Blank_SAP.Sleep.gammaBandPower.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds2*3,data.Blank_SAP.Sleep.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,mean(data.Blank_SAP.Sleep.gammaBandPower.C01),std(data.Blank_SAP.Sleep.gammaBandPower.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds2*4,data.SSP_SAP.Sleep.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e4 = errorbar(4,mean(data.SSP_SAP.Sleep.gammaBandPower.C01),std(data.Blank_SAP.Sleep.gammaBandPower.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds2*5,data.Blank_SAP.Sleep.gammaBandPower.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e5 = errorbar(5,mean(data.Blank_SAP.Sleep.gammaBandPower.C05),std(data.Blank_SAP.Sleep.gammaBandPower.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds2*6,data.SSP_SAP.Sleep.gammaBandPower.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e6 = errorbar(6,mean(data.SSP_SAP.Sleep.gammaBandPower.C05),std(data.Blank_SAP.Sleep.gammaBandPower.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
% stat lines
plot([1,2],[1,1],'k');
text(1.5,1,'ns','FontSize',16)
plot([3,4],[1,1],'k');
text(3.5,1,'*','FontSize',16)
plot([5,6],[1,1],'k');
text(5.5,1,'ns','FontSize',16)
title({'Asleep Gamma-band power',''})
ylabel('Average coherence^2')
set(gca,'xtick',[1.5,3.5,5.5])
xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
axis square
axis tight
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% All Gamma-band power stats
ax6 = subplot(2,3,6);
scatter(xInds*1,data.Blank_SAP.All.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(data.Blank_SAP.All.gammaBandPower.C001),std(data.Blank_SAP.All.gammaBandPower.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.SSP_SAP.All.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,mean(data.SSP_SAP.All.gammaBandPower.C001),std(data.Blank_SAP.All.gammaBandPower.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.Blank_SAP.All.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,mean(data.Blank_SAP.All.gammaBandPower.C01),std(data.Blank_SAP.All.gammaBandPower.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.SSP_SAP.All.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e4 = errorbar(4,mean(data.SSP_SAP.All.gammaBandPower.C01),std(data.Blank_SAP.All.gammaBandPower.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.Blank_SAP.All.gammaBandPower.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e5 = errorbar(5,mean(data.Blank_SAP.All.gammaBandPower.C05),std(data.Blank_SAP.All.gammaBandPower.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.SSP_SAP.All.gammaBandPower.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e6 = errorbar(6,mean(data.SSP_SAP.All.gammaBandPower.C05),std(data.Blank_SAP.All.gammaBandPower.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
% stat lines
plot([1,2],[1,1],'k');
text(1.5,1,'ns','FontSize',16)
plot([3,4],[1,1],'k');
text(3.5,1,'ns','FontSize',16)
plot([5,6],[1,1],'k');
text(5.5,1,'ns','FontSize',16)
title({'All Gamma-band power',''})
ylabel('Average coherence^2')
set(gca,'xtick',[1.5,3.5,5.5])
xticklabels({'0.03:0.01 Hz','0.01:0.1 Hz','0.1:0.5 Hz'})
axis square
axis tight
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Bilateral Coherence - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure5,[dirpath 'AverageBilateralCoherence_Statistics']);
    set(summaryFigure5,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'AverageBilateralCoherence_Statistics'])
    %% statistical diary
    diaryFile = [dirpath 'AverageBilateralCoherence_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % Awake stats
    disp('======================================================================================================================')
    disp('GLME statistics for gamma Coherence^2 for Awake data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.gammaBandPower.Awake.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.gammaBandPower.Awake.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.gammaBandPower.Awake.C05.Stats)
    % Sleep stats
    disp('======================================================================================================================')
    disp('GLME statistics for gamma Coherence^2 for Sleep data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.gammaBandPower.Sleep.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.gammaBandPower.Sleep.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.gammaBandPower.Sleep.C05.Stats)
    % All stats
    disp('======================================================================================================================')
    disp('GLME statistics for gamma Coherence^2 for All data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.gammaBandPower.All.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.gammaBandPower.All.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.gammaBandPower.All.C05.Stats)
    % Awake stats
    disp('======================================================================================================================')
    disp('GLME statistics for CBV_HbT Coherence^2 for Awake data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.CBV_HbT.Awake.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.CBV_HbT.Awake.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.CBV_HbT.Awake.C05.Stats)
    % Sleep stats
    disp('======================================================================================================================')
    disp('GLME statistics for CBV_HbT Coherence^2 for Sleep data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.CBV_HbT.Sleep.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.CBV_HbT.Sleep.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.CBV_HbT.Sleep.C05.Stats)
    % All stats
    disp('======================================================================================================================')
    disp('GLME statistics for CBV_HbT Coherence^2 for All data')
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.CBV_HbT.All.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.CBV_HbT.All.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.CBV_HbT.All.C05.Stats)
    diary off
end

end
