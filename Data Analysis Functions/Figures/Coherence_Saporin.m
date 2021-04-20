function [AnalysisResults] = Coherence_Saporin(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: 
%________________________________________________________________________________________________________________________

%% set-up and process data
animalIDs = {'T141','T155','T156','T157','T142','T144','T159','T172','T150','T165','T166','T177','T179','T186','T187','T188','T189'};
C57BL6J_IDs = {'T141','T155','T156','T157','T186','T187','T188','T189'};
SSP_SAP_IDs = {'T142','T144','T159','T172'};
Blank_SAP_IDs = {'T150','T165','T166','T177','T179'};
treatments = {'C57BL6J','SSP_SAP','Blank_SAP'};
behavFields = {'Rest','NREM','REM','Awake','Sleep','All'};
dataTypes = {'CBV_HbT','gammaBandPower'};
%% average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    % recognize treatment based on animal group
    if ismember(animalID,C57BL6J_IDs) == true
        treatment = 'C57BL6J';
    elseif ismember(animalIDs{1,aa},SSP_SAP_IDs) == true
        treatment = 'SSP_SAP';
    elseif ismember(animalIDs{1,aa},Blank_SAP_IDs) == true
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
            end
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).Coherence.(behavField).(dataType).C) == false
                % concatenate C/f for existing data - exclude any empty sets
                data.(treatment).(behavField).(dataType).C = cat(2,data.(treatment).(behavField).(dataType).C,AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                data.(treatment).(behavField).(dataType).f = cat(1,data.(treatment).(behavField).(dataType).f,AnalysisResults.(animalID).Coherence.(behavField).(dataType).f);
                data.(treatment).(behavField).(dataType).confC = cat(1,data.(treatment).(behavField).(dataType).confC,AnalysisResults.(animalID).Coherence.(behavField).(dataType).confC);
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
%% average HbT coherence
summaryFigure1 = figure;
%% coherence^2 between bilateral HbT during rest
subplot(2,3,1);
p1 = semilogx(data.C57BL6J.Rest.CBV_HbT.meanf,data.C57BL6J.Rest.CBV_HbT.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
p2 = semilogx(data.Blank_SAP.Rest.CBV_HbT.meanf,data.Blank_SAP.Rest.CBV_HbT.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
p3 = semilogx(data.SSP_SAP.Rest.CBV_HbT.meanf,data.SSP_SAP.Rest.CBV_HbT.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Rest] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/10,0.5])
ylim([0,1])
set(gca,'box','off')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
%% coherence^2 between bilateral HbT during NREM
subplot(2,3,2);
semilogx(data.C57BL6J.NREM.CBV_HbT.meanf,data.C57BL6J.NREM.CBV_HbT.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.NREM.CBV_HbT.meanf,data.Blank_SAP.NREM.CBV_HbT.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.NREM.CBV_HbT.meanf,data.SSP_SAP.NREM.CBV_HbT.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[NREM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/30,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during REM
subplot(2,3,3);
semilogx(data.C57BL6J.REM.CBV_HbT.meanf,data.C57BL6J.REM.CBV_HbT.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.REM.CBV_HbT.meanf,data.Blank_SAP.REM.CBV_HbT.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.REM.CBV_HbT.meanf,data.SSP_SAP.REM.CBV_HbT.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[REM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/60,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Awake
subplot(2,3,4);
semilogx(data.C57BL6J.Awake.CBV_HbT.meanf,data.C57BL6J.Awake.CBV_HbT.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.Awake.CBV_HbT.meanf,data.Blank_SAP.Awake.CBV_HbT.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.Awake.CBV_HbT.meanf,data.SSP_SAP.Awake.CBV_HbT.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(2,3,5);
semilogx(data.C57BL6J.Sleep.CBV_HbT.meanf,data.C57BL6J.Sleep.CBV_HbT.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.Sleep.CBV_HbT.meanf,data.Blank_SAP.Sleep.CBV_HbT.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.Sleep.CBV_HbT.meanf,data.SSP_SAP.Sleep.CBV_HbT.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(2,3,6);
semilogx(data.C57BL6J.All.CBV_HbT.meanf,data.C57BL6J.All.CBV_HbT.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.All.CBV_HbT.meanf,data.Blank_SAP.All.CBV_HbT.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.All.CBV_HbT.meanf,data.SSP_SAP.All.CBV_HbT.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'Coherence_HbT']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Coherence_HbT'])
end
%% individual HbT coherence
summaryFigure2 = figure;
%% coherence^2 between bilateral HbT during rest
subplot(2,3,1);
% C57BL6J
for aa = 1:size(data.C57BL6J.Rest.CBV_HbT.C,2)
    semilogx(data.C57BL6J.Rest.CBV_HbT.meanf,data.C57BL6J.Rest.CBV_HbT.C(:,aa).^2,'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Rest.CBV_HbT.C,2)
    semilogx(data.Blank_SAP.Rest.CBV_HbT.meanf,data.Blank_SAP.Rest.CBV_HbT.C(:,aa).^2,'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Rest.CBV_HbT.C,2)
    semilogx(data.SSP_SAP.Rest.CBV_HbT.meanf,data.SSP_SAP.Rest.CBV_HbT.C(:,aa).^2,'color',colors('electric purple'),'LineWidth',0.5);
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
% C57BL6J
for aa = 1:size(data.C57BL6J.NREM.CBV_HbT.C,2)
    semilogx(data.C57BL6J.NREM.CBV_HbT.meanf,data.C57BL6J.NREM.CBV_HbT.C(:,aa).^2,'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.NREM.CBV_HbT.C,2)
    semilogx(data.Blank_SAP.NREM.CBV_HbT.meanf,data.Blank_SAP.NREM.CBV_HbT.C(:,aa).^2,'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.NREM.CBV_HbT.C,2)
    semilogx(data.SSP_SAP.NREM.CBV_HbT.meanf,data.SSP_SAP.NREM.CBV_HbT.C(:,aa).^2,'color',colors('electric purple'),'LineWidth',0.5);
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
% C57BL6J
for aa = 1:size(data.C57BL6J.REM.CBV_HbT.C,2)
    semilogx(data.C57BL6J.REM.CBV_HbT.meanf,data.C57BL6J.REM.CBV_HbT.C(:,aa).^2,'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.REM.CBV_HbT.C,2)
    semilogx(data.Blank_SAP.REM.CBV_HbT.meanf,data.Blank_SAP.REM.CBV_HbT.C(:,aa).^2,'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.REM.CBV_HbT.C,2)
    semilogx(data.SSP_SAP.REM.CBV_HbT.meanf,data.SSP_SAP.REM.CBV_HbT.C(:,aa).^2,'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[REM] Bilateral coherence^2','\Delta[HbT] (\muM)'})
xlim([1/60,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Awake
subplot(2,3,4);
% C57BL6J
for aa = 1:size(data.C57BL6J.Awake.CBV_HbT.C,2)
    semilogx(data.C57BL6J.Awake.CBV_HbT.meanf,data.C57BL6J.Awake.CBV_HbT.C(:,aa).^2,'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Awake.CBV_HbT.C,2)
    semilogx(data.Blank_SAP.Awake.CBV_HbT.meanf,data.Blank_SAP.Awake.CBV_HbT.C(:,aa).^2,'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Awake.CBV_HbT.C,2)
    semilogx(data.SSP_SAP.Awake.CBV_HbT.meanf,data.SSP_SAP.Awake.CBV_HbT.C(:,aa).^2,'color',colors('electric purple'),'LineWidth',0.5);
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
% C57BL6J
for aa = 1:size(data.C57BL6J.Sleep.CBV_HbT.C,2)
    semilogx(data.C57BL6J.Sleep.CBV_HbT.meanf,data.C57BL6J.Sleep.CBV_HbT.C(:,aa).^2,'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Sleep.CBV_HbT.C,2)
    semilogx(data.Blank_SAP.Sleep.CBV_HbT.meanf,data.Blank_SAP.Sleep.CBV_HbT.C(:,aa).^2,'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Sleep.CBV_HbT.C,2)
    semilogx(data.SSP_SAP.Sleep.CBV_HbT.meanf,data.SSP_SAP.Sleep.CBV_HbT.C(:,aa).^2,'color',colors('electric purple'),'LineWidth',0.5);
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
% C57BL6J
for aa = 1:size(data.C57BL6J.All.CBV_HbT.C,2)
    semilogx(data.C57BL6J.All.CBV_HbT.meanf,data.C57BL6J.All.CBV_HbT.C(:,aa).^2,'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.All.CBV_HbT.C,2)
    semilogx(data.Blank_SAP.All.CBV_HbT.meanf,data.Blank_SAP.All.CBV_HbT.C(:,aa).^2,'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.All.CBV_HbT.C,2)
    semilogx(data.SSP_SAP.All.CBV_HbT.meanf,data.SSP_SAP.All.CBV_HbT.C(:,aa).^2,'color',colors('electric purple'),'LineWidth',0.5);
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
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'indCoherence_HbT']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'indCoherence_HbT'])
end
%% average gamma-band coherence
summaryFigure3 = figure;
%% coherence^2 between bilateral gamma-band during Rest
subplot(2,3,1);
p1 = semilogx(data.C57BL6J.Rest.gammaBandPower.meanf,data.C57BL6J.Rest.gammaBandPower.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
p2 = semilogx(data.Blank_SAP.Rest.gammaBandPower.meanf,data.Blank_SAP.Rest.gammaBandPower.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
p3 = semilogx(data.SSP_SAP.Rest.gammaBandPower.meanf,data.SSP_SAP.Rest.gammaBandPower.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Rest] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([1/10,0.5])
ylim([0,1])
set(gca,'box','off')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
%% coherence^2 between bilateral gamma-band during NREM
subplot(2,3,2);
semilogx(data.C57BL6J.NREM.gammaBandPower.meanf,data.C57BL6J.NREM.gammaBandPower.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.NREM.gammaBandPower.meanf,data.Blank_SAP.NREM.gammaBandPower.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.NREM.gammaBandPower.meanf,data.SSP_SAP.NREM.gammaBandPower.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[NREM] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([1/30,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during REM
subplot(2,3,3);
semilogx(data.C57BL6J.REM.gammaBandPower.meanf,data.C57BL6J.REM.gammaBandPower.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.REM.gammaBandPower.meanf,data.Blank_SAP.REM.gammaBandPower.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.REM.gammaBandPower.meanf,data.SSP_SAP.REM.gammaBandPower.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[REM] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([1/60,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during Awake
subplot(2,3,4);
semilogx(data.C57BL6J.Awake.gammaBandPower.meanf,data.C57BL6J.Awake.gammaBandPower.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.Awake.gammaBandPower.meanf,data.Blank_SAP.Awake.gammaBandPower.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.Awake.gammaBandPower.meanf,data.SSP_SAP.Awake.gammaBandPower.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during Sleep
subplot(2,3,5);
semilogx(data.C57BL6J.Sleep.gammaBandPower.meanf,data.C57BL6J.Sleep.gammaBandPower.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.Sleep.gammaBandPower.meanf,data.Blank_SAP.Sleep.gammaBandPower.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.Sleep.gammaBandPower.meanf,data.SSP_SAP.Sleep.gammaBandPower.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral gamma-band during All data
subplot(2,3,6);
semilogx(data.C57BL6J.All.gammaBandPower.meanf,data.C57BL6J.All.gammaBandPower.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.All.gammaBandPower.meanf,data.Blank_SAP.All.gammaBandPower.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.All.gammaBandPower.meanf,data.SSP_SAP.All.gammaBandPower.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] Bilateral coherence^2','Gamma-band power [30-100 Hz] (%)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure3,[dirpath 'Coherence_gamma']);
    set(summaryFigure3,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Coherence_gamma'])
end
%% individual gamma-band coherence
summaryFigure4 = figure;
%% coherence^2 between bilateral gamma-band during rest
subplot(2,3,1);
% C57BL6J
for aa = 1:size(data.C57BL6J.Rest.gammaBandPower.C,2)
    semilogx(data.C57BL6J.Rest.gammaBandPower.meanf,data.C57BL6J.Rest.gammaBandPower.C(:,aa).^2,'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Rest.gammaBandPower.C,2)
    semilogx(data.Blank_SAP.Rest.gammaBandPower.meanf,data.Blank_SAP.Rest.gammaBandPower.C(:,aa).^2,'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Rest.gammaBandPower.C,2)
    semilogx(data.SSP_SAP.Rest.gammaBandPower.meanf,data.SSP_SAP.Rest.gammaBandPower.C(:,aa).^2,'color',colors('electric purple'),'LineWidth',0.5);
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
% C57BL6J
for aa = 1:size(data.C57BL6J.NREM.gammaBandPower.C,2)
    semilogx(data.C57BL6J.NREM.gammaBandPower.meanf,data.C57BL6J.NREM.gammaBandPower.C(:,aa).^2,'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.NREM.gammaBandPower.C,2)
    semilogx(data.Blank_SAP.NREM.gammaBandPower.meanf,data.Blank_SAP.NREM.gammaBandPower.C(:,aa).^2,'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.NREM.gammaBandPower.C,2)
    semilogx(data.SSP_SAP.NREM.gammaBandPower.meanf,data.SSP_SAP.NREM.gammaBandPower.C(:,aa).^2,'color',colors('electric purple'),'LineWidth',0.5);
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
% C57BL6J
for aa = 1:size(data.C57BL6J.REM.gammaBandPower.C,2)
    semilogx(data.C57BL6J.REM.gammaBandPower.meanf,data.C57BL6J.REM.gammaBandPower.C(:,aa).^2,'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.REM.gammaBandPower.C,2)
    semilogx(data.Blank_SAP.REM.gammaBandPower.meanf,data.Blank_SAP.REM.gammaBandPower.C(:,aa).^2,'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.REM.gammaBandPower.C,2)
    semilogx(data.SSP_SAP.REM.gammaBandPower.meanf,data.SSP_SAP.REM.gammaBandPower.C(:,aa).^2,'color',colors('electric purple'),'LineWidth',0.5);
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
% C57BL6J
for aa = 1:size(data.C57BL6J.Awake.gammaBandPower.C,2)
    semilogx(data.C57BL6J.Awake.gammaBandPower.meanf,data.C57BL6J.Awake.gammaBandPower.C(:,aa).^2,'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Awake.gammaBandPower.C,2)
    semilogx(data.Blank_SAP.Awake.gammaBandPower.meanf,data.Blank_SAP.Awake.gammaBandPower.C(:,aa).^2,'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Awake.gammaBandPower.C,2)
    semilogx(data.SSP_SAP.Awake.gammaBandPower.meanf,data.SSP_SAP.Awake.gammaBandPower.C(:,aa).^2,'color',colors('electric purple'),'LineWidth',0.5);
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
% C57BL6J
for aa = 1:size(data.C57BL6J.Sleep.gammaBandPower.C,2)
    semilogx(data.C57BL6J.Sleep.gammaBandPower.meanf,data.C57BL6J.Sleep.gammaBandPower.C(:,aa).^2,'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Sleep.gammaBandPower.C,2)
    semilogx(data.Blank_SAP.Sleep.gammaBandPower.meanf,data.Blank_SAP.Sleep.gammaBandPower.C(:,aa).^2,'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Sleep.gammaBandPower.C,2)
    semilogx(data.SSP_SAP.Sleep.gammaBandPower.meanf,data.SSP_SAP.Sleep.gammaBandPower.C(:,aa).^2,'color',colors('electric purple'),'LineWidth',0.5);
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
% C57BL6J
for aa = 1:size(data.C57BL6J.All.gammaBandPower.C,2)
    semilogx(data.C57BL6J.All.gammaBandPower.meanf,data.C57BL6J.All.gammaBandPower.C(:,aa).^2,'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.All.gammaBandPower.C,2)
    semilogx(data.Blank_SAP.All.gammaBandPower.meanf,data.Blank_SAP.All.gammaBandPower.C(:,aa).^2,'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.All.gammaBandPower.C,2)
    semilogx(data.SSP_SAP.All.gammaBandPower.meanf,data.SSP_SAP.All.gammaBandPower.C(:,aa).^2,'color',colors('electric purple'),'LineWidth',0.5);
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
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure4,[dirpath 'indCoherence_gamma']);
    set(summaryFigure4,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'indCoherence_gamma'])
end

end
