function [] = NeuralHemoCoherence_Bilateral_IOS(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% set-up and process data
resultsStruct = 'Results_NeuralHemoCoher';
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
behavFields = {'Awake','Sleep','All'};
dataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
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
                data.(treatment).(behavField).(dataType).LH.C = [];
                data.(treatment).(behavField).(dataType).LH.f = [];
                data.(treatment).(behavField).(dataType).LH.confC = [];
                data.(treatment).(behavField).(dataType).RH.C = [];
                data.(treatment).(behavField).(dataType).RH.f = [];
                data.(treatment).(behavField).(dataType).RH.confC = [];
                data.(treatment).(behavField).(dataType).animalID = {};
                data.(treatment).(behavField).(dataType).treatment = {};
                data.(treatment).(behavField).(dataType).freqC001 = {};
                data.(treatment).(behavField).(dataType).freqC01 = {};
                data.(treatment).(behavField).(dataType).freqC05 = {};
            end
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(Results_NeuralHemoCoher.(animalID).(behavField).(dataType).adjLH.C) == false
                % concatenate C/f for existing data - exclude any empty sets
                data.(treatment).(behavField).(dataType).LH.C = cat(2,data.(treatment).(behavField).(dataType).LH.C,Results_NeuralHemoCoher.(animalID).(behavField).(dataType).adjLH.C.^2);
                data.(treatment).(behavField).(dataType).LH.f = cat(1,data.(treatment).(behavField).(dataType).LH.f,Results_NeuralHemoCoher.(animalID).(behavField).(dataType).adjLH.f);
                data.(treatment).(behavField).(dataType).LH.confC = cat(1,data.(treatment).(behavField).(dataType).LH.confC,Results_NeuralHemoCoher.(animalID).(behavField).(dataType).adjLH.confC);
                data.(treatment).(behavField).(dataType).RH.C = cat(2,data.(treatment).(behavField).(dataType).RH.C,Results_NeuralHemoCoher.(animalID).(behavField).(dataType).adjRH.C.^2);
                data.(treatment).(behavField).(dataType).RH.f = cat(1,data.(treatment).(behavField).(dataType).RH.f,Results_NeuralHemoCoher.(animalID).(behavField).(dataType).adjRH.f);
                data.(treatment).(behavField).(dataType).RH.confC = cat(1,data.(treatment).(behavField).(dataType).RH.confC,Results_NeuralHemoCoher.(animalID).(behavField).(dataType).adjRH.confC);
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
freqs = {'C001','C01','C05'};
for aa = 1:length(freqs)
    freq = freqs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % statistics - generalized linear mixed effects model
            Stats.(dataType).(behavField).(freq).tableSize = cat(1,data.Blank_SAP.(behavField).(dataType).RH.(freq),data.SSP_SAP.(behavField).(dataType).RH.(freq));
            Stats.(dataType).(behavField).(freq).Table = table('Size',[size(Stats.(dataType).(behavField).(freq).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Treatment','Coherence'});
            Stats.(dataType).(behavField).(freq).Table.Mouse = cat(1,data.Blank_SAP.(behavField).(dataType).animalID,data.SSP_SAP.(behavField).(dataType).animalID);
            Stats.(dataType).(behavField).(freq).Table.Treatment = cat(1,data.Blank_SAP.(behavField).(dataType).treatment,data.SSP_SAP.(behavField).(dataType).treatment);
            Stats.(dataType).(behavField).(freq).Table.Coherence = cat(1,data.Blank_SAP.(behavField).(dataType).RH.(freq),data.SSP_SAP.(behavField).(dataType).RH.(freq));
            Stats.(dataType).(behavField).(freq).FitFormula = 'Coherence ~ 1 + Treatment + (1|Mouse)';
            Stats.(dataType).(behavField).(freq).Stats = fitglme(Stats.(dataType).(behavField).(freq).Table,Stats.(dataType).(behavField).(freq).FitFormula);
        end
    end
end
%% statistics - generalized linear mixed effects model
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
%% neural-hemo coherence figure
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    figName = ['summaryFigure' num2str(aa)]; %#ok<NASGU>
    figName = figure;
    sgtitle([dataType '- \DeltaHbT Coherence'])
    %% LH coherence^2 between neural-HbT during Awake
    ax1 = subplot(3,2,1);
    s1 = semilogx(data.Naive.Awake.(dataType).LH.meanf,data.Naive.Awake.(dataType).LH.meanC,'color',colors('sapphire'),'LineWidth',2);
    hold on
    semilogx(data.Naive.Awake.(dataType).LH.meanf,data.Naive.Awake.(dataType).LH.meanC + data.Naive.Awake.(dataType).LH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
    semilogx(data.Naive.Awake.(dataType).LH.meanf,data.Naive.Awake.(dataType).LH.meanC - data.Naive.Awake.(dataType).LH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
    s2 = semilogx(data.Blank_SAP.Awake.(dataType).LH.meanf,data.Blank_SAP.Awake.(dataType).LH.meanC,'color',colors('north texas green'),'LineWidth',2);
    semilogx(data.Blank_SAP.Awake.(dataType).LH.meanf,data.Blank_SAP.Awake.(dataType).LH.meanC + data.Blank_SAP.Awake.(dataType).LH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
    semilogx(data.Blank_SAP.Awake.(dataType).LH.meanf,data.Blank_SAP.Awake.(dataType).LH.meanC - data.Blank_SAP.Awake.(dataType).LH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
    s3 = semilogx(data.SSP_SAP.Awake.(dataType).LH.meanf,data.SSP_SAP.Awake.(dataType).LH.meanC,'color',colors('electric purple'),'LineWidth',2);
    semilogx(data.SSP_SAP.Awake.(dataType).LH.meanf,data.SSP_SAP.Awake.(dataType).LH.meanC + data.SSP_SAP.Awake.(dataType).LH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
    semilogx(data.SSP_SAP.Awake.(dataType).LH.meanf,data.SSP_SAP.Awake.(dataType).LH.meanC - data.SSP_SAP.Awake.(dataType).LH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
    ylabel('Coherence^2')
    xlabel('Freq (Hz)')
    title({'[Alert] LH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
    xlim([0.003,0.5])
    % ylim([0,1])
    legend([s1,s2,s3],'Naive','Blank-SAP','SSP-SAP')
    set(gca,'box','off')
    %% RH coherence^2 between neural-HbT during Awake
    ax2 = subplot(3,2,2);
    semilogx(data.Naive.Awake.(dataType).RH.meanf,data.Naive.Awake.(dataType).RH.meanC,'color',colors('sapphire'),'LineWidth',2);
    hold on
    semilogx(data.Naive.Awake.(dataType).RH.meanf,data.Naive.Awake.(dataType).RH.meanC + data.Naive.Awake.(dataType).RH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
    semilogx(data.Naive.Awake.gammaBandPower.RH.meanf,data.Naive.Awake.(dataType).RH.meanC - data.Naive.Awake.(dataType).RH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
    semilogx(data.Blank_SAP.Awake.(dataType).RH.meanf,data.Blank_SAP.Awake.(dataType).RH.meanC,'color',colors('north texas green'),'LineWidth',2);
    semilogx(data.Blank_SAP.Awake.(dataType).RH.meanf,data.Blank_SAP.Awake.(dataType).RH.meanC + data.Blank_SAP.Awake.(dataType).RH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
    semilogx(data.Blank_SAP.Awake.(dataType).RH.meanf,data.Blank_SAP.Awake.(dataType).RH.meanC - data.Blank_SAP.Awake.(dataType).RH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
    semilogx(data.SSP_SAP.Awake.(dataType).RH.meanf,data.SSP_SAP.Awake.(dataType).RH.meanC,'color',colors('electric purple'),'LineWidth',2);
    semilogx(data.SSP_SAP.Awake.(dataType).RH.meanf,data.SSP_SAP.Awake.(dataType).RH.meanC + data.SSP_SAP.Awake.(dataType).RH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
    semilogx(data.SSP_SAP.Awake.(dataType).RH.meanf,data.SSP_SAP.Awake.(dataType).RH.meanC - data.SSP_SAP.Awake.(dataType).RH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
    ylabel('Coherence^2')
    xlabel('Freq (Hz)')
    title({'[Alert] RH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
    xlim([0.003,0.5])
    % ylim([0,1])
    set(gca,'box','off')
    %% LH coherence^2 between neural-HbT during Sleep
    ax3 = subplot(3,2,3);
    semilogx(data.Naive.Sleep.(dataType).LH.meanf,data.Naive.Sleep.(dataType).LH.meanC,'color',colors('sapphire'),'LineWidth',2);
    hold on
    semilogx(data.Naive.Sleep.(dataType).LH.meanf,data.Naive.Sleep.(dataType).LH.meanC + data.Naive.Sleep.(dataType).LH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
    semilogx(data.Naive.Sleep.(dataType).LH.meanf,data.Naive.Sleep.(dataType).LH.meanC - data.Naive.Sleep.(dataType).LH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
    semilogx(data.Blank_SAP.Sleep.(dataType).LH.meanf,data.Blank_SAP.Sleep.(dataType).LH.meanC,'color',colors('north texas green'),'LineWidth',2);
    semilogx(data.Blank_SAP.Sleep.(dataType).LH.meanf,data.Blank_SAP.Sleep.(dataType).LH.meanC + data.Blank_SAP.Sleep.(dataType).LH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
    semilogx(data.Blank_SAP.Sleep.(dataType).LH.meanf,data.Blank_SAP.Sleep.(dataType).LH.meanC - data.Blank_SAP.Sleep.(dataType).LH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
    semilogx(data.SSP_SAP.Sleep.(dataType).LH.meanf,data.SSP_SAP.Sleep.(dataType).LH.meanC,'color',colors('electric purple'),'LineWidth',2);
    semilogx(data.SSP_SAP.Sleep.(dataType).LH.meanf,data.SSP_SAP.Sleep.(dataType).LH.meanC + data.SSP_SAP.Sleep.(dataType).LH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
    semilogx(data.SSP_SAP.Sleep.(dataType).LH.meanf,data.SSP_SAP.Sleep.(dataType).LH.meanC - data.SSP_SAP.Sleep.(dataType).LH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
    ylabel('Coherence^2')
    xlabel('Freq (Hz)')
    title({'[Asleep] LH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
    xlim([0.003,0.5])
    % ylim([0,1])
    set(gca,'box','off')
    %% RH coherence^2 between neural-HbT during Sleep
    ax4 = subplot(3,2,4);
    semilogx(data.Naive.Sleep.(dataType).RH.meanf,data.Naive.Sleep.(dataType).RH.meanC,'color',colors('sapphire'),'LineWidth',2);
    hold on
    semilogx(data.Naive.Sleep.(dataType).RH.meanf,data.Naive.Sleep.(dataType).RH.meanC + data.Naive.Sleep.(dataType).RH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
    semilogx(data.Naive.Sleep.(dataType).RH.meanf,data.Naive.Sleep.(dataType).RH.meanC - data.Naive.Sleep.(dataType).RH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
    semilogx(data.Blank_SAP.Sleep.(dataType).RH.meanf,data.Blank_SAP.Sleep.(dataType).RH.meanC,'color',colors('north texas green'),'LineWidth',2);
    semilogx(data.Blank_SAP.Sleep.(dataType).RH.meanf,data.Blank_SAP.Sleep.(dataType).RH.meanC + data.Blank_SAP.Sleep.(dataType).RH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
    semilogx(data.Blank_SAP.Sleep.(dataType).RH.meanf,data.Blank_SAP.Sleep.(dataType).RH.meanC - data.Blank_SAP.Sleep.(dataType).RH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
    semilogx(data.SSP_SAP.Sleep.(dataType).RH.meanf,data.SSP_SAP.Sleep.(dataType).RH.meanC,'color',colors('electric purple'),'LineWidth',2);
    semilogx(data.SSP_SAP.Sleep.(dataType).RH.meanf,data.SSP_SAP.Sleep.(dataType).RH.meanC + data.SSP_SAP.Sleep.(dataType).RH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
    semilogx(data.SSP_SAP.Sleep.(dataType).RH.meanf,data.SSP_SAP.Sleep.(dataType).RH.meanC - data.SSP_SAP.Sleep.(dataType).RH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
    ylabel('Coherence^2')
    xlabel('Freq (Hz)')
    title({'[Asleep] RH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
    xlim([0.003,0.5])
    % ylim([0,1])
    set(gca,'box','off')
    %% LH coherence^2 between neural-HbT during All
    ax5 = subplot(3,2,5);
    semilogx(data.Naive.All.(dataType).LH.meanf,data.Naive.All.(dataType).LH.meanC,'color',colors('sapphire'),'LineWidth',2);
    hold on
    semilogx(data.Naive.All.(dataType).LH.meanf,data.Naive.All.(dataType).LH.meanC + data.Naive.All.(dataType).LH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
    semilogx(data.Naive.All.(dataType).LH.meanf,data.Naive.All.(dataType).LH.meanC - data.Naive.All.(dataType).LH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
    semilogx(data.Blank_SAP.All.(dataType).LH.meanf,data.Blank_SAP.All.(dataType).LH.meanC,'color',colors('north texas green'),'LineWidth',2);
    semilogx(data.Blank_SAP.All.(dataType).LH.meanf,data.Blank_SAP.All.(dataType).LH.meanC + data.Blank_SAP.All.(dataType).LH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
    semilogx(data.Blank_SAP.All.(dataType).LH.meanf,data.Blank_SAP.All.(dataType).LH.meanC - data.Blank_SAP.All.(dataType).LH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
    semilogx(data.SSP_SAP.All.(dataType).LH.meanf,data.SSP_SAP.All.(dataType).LH.meanC,'color',colors('electric purple'),'LineWidth',2);
    semilogx(data.SSP_SAP.All.(dataType).LH.meanf,data.SSP_SAP.All.(dataType).LH.meanC + data.SSP_SAP.All.(dataType).LH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
    semilogx(data.SSP_SAP.All.(dataType).LH.meanf,data.SSP_SAP.All.(dataType).LH.meanC - data.SSP_SAP.All.(dataType).LH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
    ylabel('Coherence^2')
    xlabel('Freq (Hz)')
    title({'[All] LH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
    xlim([0.003,0.5])
    % ylim([0,1])
    set(gca,'box','off')
    %% RH coherence^2 between neural-HbT during All
    ax6 = subplot(3,2,6);
    semilogx(data.Naive.All.(dataType).RH.meanf,data.Naive.All.(dataType).RH.meanC,'color',colors('sapphire'),'LineWidth',2);
    hold on
    semilogx(data.Naive.All.(dataType).RH.meanf,data.Naive.All.(dataType).RH.meanC + data.Naive.All.(dataType).RH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
    semilogx(data.Naive.All.(dataType).RH.meanf,data.Naive.All.(dataType).RH.meanC - data.Naive.All.(dataType).RH.stdC,'color',colors('sapphire'),'LineWidth',0.5);
    semilogx(data.Blank_SAP.All.(dataType).RH.meanf,data.Blank_SAP.All.(dataType).RH.meanC,'color',colors('north texas green'),'LineWidth',2);
    semilogx(data.Blank_SAP.All.(dataType).RH.meanf,data.Blank_SAP.All.(dataType).RH.meanC + data.Blank_SAP.All.(dataType).RH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
    semilogx(data.Blank_SAP.All.(dataType).RH.meanf,data.Blank_SAP.All.(dataType).RH.meanC - data.Blank_SAP.All.(dataType).RH.stdC,'color',colors('north texas green'),'LineWidth',0.5);
    semilogx(data.SSP_SAP.All.(dataType).RH.meanf,data.SSP_SAP.All.(dataType).RH.meanC,'color',colors('electric purple'),'LineWidth',2);
    semilogx(data.SSP_SAP.All.(dataType).RH.meanf,data.SSP_SAP.All.(dataType).RH.meanC + data.SSP_SAP.All.(dataType).RH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
    semilogx(data.SSP_SAP.All.(dataType).RH.meanf,data.SSP_SAP.All.(dataType).RH.meanC - data.SSP_SAP.All.(dataType).RH.stdC,'color',colors('electric purple'),'LineWidth',0.5);
    ylabel('Coherence^2')
    xlabel('Freq (Hz)')
    title({'[All] RH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
    xlim([0.003,0.5])
    % ylim([0,1])
    set(gca,'box','off')
    %% figure characteristics
    linkaxes([ax1,ax2],'xy')
    linkaxes([ax3,ax4],'xy')
    linkaxes([ax5,ax6],'xy')
    %% save figure(s)
    if saveFigs == true
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Neural Hemo Coherence - Bilateral IOS' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(figName,[dirpath 'AverageNeuralHemoCoherence_' dataType]);
        set(figName,'PaperPositionMode','auto');
        print('-painters','-dpdf','-fillpage',[dirpath 'AverageNeuralHemoCoherence_' dataType])
    end
    %% statistical diary
    diaryFile = [dirpath 'AverageNeuralHemoCoherence_' dataType '_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % Awake stats
    disp('======================================================================================================================')
    disp(['GLME statistics for ' dataType '-HbT Coherence^2 for Awake data'])
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.(dataType).Awake.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.(dataType).Awake.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.(dataType).Awake.C05.Stats)
    % Sleep stats
    disp('======================================================================================================================')
    disp(['GLME statistics for ' dataType '-HbT Coherence^2 for Sleep data'])
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.(dataType).Sleep.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.(dataType).Sleep.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.(dataType).Sleep.C05.Stats)
    % All stats
    disp('======================================================================================================================')
    disp(['GLME statistics for ' dataType '-HbT Coherence^2 for All data'])
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.(dataType).All.C001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.(dataType).All.C01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.(dataType).All.C05.Stats)
    diary off
end
%% HbT and gamma-band coherence stats
summaryFigure6 = figure;
sgtitle('Gamma-HbT Coherence Statics')
%% Alert HbT Stats
ax1 = subplot(1,3,1);
xInds = ones(1,length(animalIDs.Blank_SAP));
s1 = scatter(xInds*1,data.Blank_SAP.Awake.gammaBandPower.RH.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(data.Blank_SAP.Awake.gammaBandPower.RH.C001),std(data.Blank_SAP.Awake.gammaBandPower.RH.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(xInds*2,data.SSP_SAP.Awake.gammaBandPower.RH.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,mean(data.SSP_SAP.Awake.gammaBandPower.RH.C001),std(data.Blank_SAP.Awake.gammaBandPower.RH.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.Blank_SAP.Awake.gammaBandPower.RH.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,mean(data.Blank_SAP.Awake.gammaBandPower.RH.C01),std(data.Blank_SAP.Awake.gammaBandPower.RH.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.SSP_SAP.Awake.gammaBandPower.RH.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e4 = errorbar(4,mean(data.SSP_SAP.Awake.gammaBandPower.RH.C01),std(data.Blank_SAP.Awake.gammaBandPower.RH.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.Blank_SAP.Awake.gammaBandPower.RH.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e5 = errorbar(5,mean(data.Blank_SAP.Awake.gammaBandPower.RH.C05),std(data.Blank_SAP.Awake.gammaBandPower.RH.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.SSP_SAP.Awake.gammaBandPower.RH.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e6 = errorbar(6,mean(data.SSP_SAP.Awake.gammaBandPower.RH.C05),std(data.Blank_SAP.Awake.gammaBandPower.RH.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
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
ax2 = subplot(1,3,2);
xInds2 = ones(1,length(data.Blank_SAP.Sleep.gammaBandPower.RH.C001));
scatter(xInds2*1,data.Blank_SAP.Sleep.gammaBandPower.RH.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(data.Blank_SAP.Sleep.gammaBandPower.RH.C001),std(data.Blank_SAP.Sleep.gammaBandPower.RH.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds2*2,data.SSP_SAP.Sleep.gammaBandPower.RH.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,mean(data.SSP_SAP.Sleep.gammaBandPower.RH.C001),std(data.Blank_SAP.Sleep.gammaBandPower.RH.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds2*3,data.Blank_SAP.Sleep.gammaBandPower.RH.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,mean(data.Blank_SAP.Sleep.gammaBandPower.RH.C01),std(data.Blank_SAP.Sleep.gammaBandPower.RH.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds2*4,data.SSP_SAP.Sleep.gammaBandPower.RH.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e4 = errorbar(4,mean(data.SSP_SAP.Sleep.gammaBandPower.RH.C01),std(data.Blank_SAP.Sleep.gammaBandPower.RH.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds2*5,data.Blank_SAP.Sleep.gammaBandPower.RH.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e5 = errorbar(5,mean(data.Blank_SAP.Sleep.gammaBandPower.RH.C05),std(data.Blank_SAP.Sleep.gammaBandPower.RH.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds2*6,data.SSP_SAP.Sleep.gammaBandPower.RH.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e6 = errorbar(6,mean(data.SSP_SAP.Sleep.gammaBandPower.RH.C05),std(data.Blank_SAP.Sleep.gammaBandPower.RH.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
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
ax3 = subplot(1,3,3);
xInds = ones(1,length(animalIDs.Blank_SAP));
scatter(xInds*1,data.Blank_SAP.All.gammaBandPower.RH.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(data.Blank_SAP.All.gammaBandPower.RH.C001),std(data.Blank_SAP.All.gammaBandPower.RH.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.SSP_SAP.All.gammaBandPower.RH.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,mean(data.SSP_SAP.All.gammaBandPower.RH.C001),std(data.Blank_SAP.All.gammaBandPower.RH.C001),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.Blank_SAP.All.gammaBandPower.RH.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,mean(data.Blank_SAP.All.gammaBandPower.RH.C01),std(data.Blank_SAP.All.gammaBandPower.RH.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.SSP_SAP.All.gammaBandPower.RH.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e4 = errorbar(4,mean(data.SSP_SAP.All.gammaBandPower.RH.C01),std(data.Blank_SAP.All.gammaBandPower.RH.C01),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.Blank_SAP.All.gammaBandPower.RH.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e5 = errorbar(5,mean(data.Blank_SAP.All.gammaBandPower.RH.C05),std(data.Blank_SAP.All.gammaBandPower.RH.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.SSP_SAP.All.gammaBandPower.RH.C05,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e6 = errorbar(6,mean(data.SSP_SAP.All.gammaBandPower.RH.C05),std(data.Blank_SAP.All.gammaBandPower.RH.C05),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
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
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Neural Hemo Coherence - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure6,[dirpath 'AverageNeuralHemoCoherence_gammaBandPower_Statistics']);
    set(summaryFigure6,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'AverageNeuralHemoCoherence_gammaBandPower_Statistics'])
end
