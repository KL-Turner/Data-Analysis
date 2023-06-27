function [AnalysisResults] = WhiskHemoCoherence_Bilateral_IOS(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% set-up and process data
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
        % pre-allocate necessary variable fields
        data.(treatment).(behavField).dummCheck = 1;
        if isfield(data.(treatment).(behavField),'LH') == false
            data.(treatment).(behavField).LH.C = [];
            data.(treatment).(behavField).LH.f = [];
            data.(treatment).(behavField).LH.confC = [];
            data.(treatment).(behavField).RH.C = [];
            data.(treatment).(behavField).RH.f = [];
            data.(treatment).(behavField).RH.confC = [];
        end
        % don't concatenate empty arrays where there was no data for this behavior
        if isempty(AnalysisResults.(animalID).WhiskHemoCoherence.(behavField).adjLH.C) == false
            % concatenate C/f for existing data - exclude any empty sets
            data.(treatment).(behavField).LH.C = cat(2,data.(treatment).(behavField).LH.C,AnalysisResults.(animalID).WhiskHemoCoherence.(behavField).adjLH.C);
            data.(treatment).(behavField).LH.f = cat(1,data.(treatment).(behavField).LH.f,AnalysisResults.(animalID).WhiskHemoCoherence.(behavField).adjLH.f);
            data.(treatment).(behavField).LH.confC = cat(1,data.(treatment).(behavField).LH.confC,AnalysisResults.(animalID).WhiskHemoCoherence.(behavField).adjLH.confC);
            data.(treatment).(behavField).RH.C = cat(2,data.(treatment).(behavField).RH.C,AnalysisResults.(animalID).WhiskHemoCoherence.(behavField).adjRH.C);
            data.(treatment).(behavField).RH.f = cat(1,data.(treatment).(behavField).RH.f,AnalysisResults.(animalID).WhiskHemoCoherence.(behavField).adjRH.f);
            data.(treatment).(behavField).RH.confC = cat(1,data.(treatment).(behavField).RH.confC,AnalysisResults.(animalID).WhiskHemoCoherence.(behavField).adjRH.confC);
        end
    end
end
%% take mean/StD of C/f and determine confC line
for qq = 1:length(treatments)
    treatment = treatments{1,qq};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        data.(treatment).(behavField).LH.meanC = mean(data.(treatment).(behavField).LH.C,2);
        data.(treatment).(behavField).LH.stdC = std(data.(treatment).(behavField).LH.C,0,2);
        data.(treatment).(behavField).LH.meanf = mean(data.(treatment).(behavField).LH.f,1);
        data.(treatment).(behavField).LH.maxConfC = geomean(data.(treatment).(behavField).LH.confC);
        data.(treatment).(behavField).LH.maxConfC_Y = ones(length(data.(treatment).(behavField).LH.meanf),1)*data.(treatment).(behavField).LH.maxConfC;
        data.(treatment).(behavField).RH.meanC = mean(data.(treatment).(behavField).RH.C,2);
        data.(treatment).(behavField).RH.stdC = std(data.(treatment).(behavField).RH.C,0,2);
        data.(treatment).(behavField).RH.meanf = mean(data.(treatment).(behavField).RH.f,1);
        data.(treatment).(behavField).RH.maxConfC = geomean(data.(treatment).(behavField).RH.confC);
        data.(treatment).(behavField).RH.maxConfC_Y = ones(length(data.(treatment).(behavField).RH.meanf),1)*data.(treatment).(behavField).RH.maxConfC;
    end
end
%% average HbT coherence
summaryFigure1 = figure;
sgtitle('Binarized whisking - \DeltaHbT Coherence')
%% coherence^2 between bilateral HbT during Awake
subplot(3,2,1);
s1 = semilogx(data.Naive.Awake.LH.meanf,data.Naive.Awake.LH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
s2 = semilogx(data.Blank_SAP.Awake.LH.meanf,data.Blank_SAP.Awake.LH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
s3 = semilogx(data.SSP_SAP.Awake.LH.meanf,data.SSP_SAP.Awake.LH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title('[Alert] LH Whisk-hemo coherence^2')
xlim([0.003,0.5])
% ylim([0,1])
legend([s1,s2,s3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Awake
subplot(3,2,2);
semilogx(data.Naive.Awake.RH.meanf,data.Naive.Awake.RH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.Awake.RH.meanf,data.Blank_SAP.Awake.RH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.Awake.RH.meanf,data.SSP_SAP.Awake.RH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title('[Alert] RH Whisk-hemo coherence^2')
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(3,2,3);
semilogx(data.Naive.Sleep.LH.meanf,data.Naive.Sleep.LH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.Sleep.LH.meanf,data.Blank_SAP.Sleep.LH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.Sleep.LH.meanf,data.SSP_SAP.Sleep.LH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title('[Asleep] LH Whisk-hemo coherence^2')
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(3,2,4);
semilogx(data.Naive.Sleep.RH.meanf,data.Naive.Sleep.RH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.Sleep.RH.meanf,data.Blank_SAP.Sleep.RH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.Sleep.RH.meanf,data.SSP_SAP.Sleep.RH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title('[Asleep] RH Whisk-hemo coherence^2')
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(3,2,5);
semilogx(data.Naive.All.LH.meanf,data.Naive.All.LH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.All.LH.meanf,data.Blank_SAP.All.LH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.All.LH.meanf,data.SSP_SAP.All.LH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title('[All] LH Whisk-hemo coherence^2')
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(3,2,6);
semilogx(data.Naive.All.RH.meanf,data.Naive.All.RH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.All.RH.meanf,data.Blank_SAP.All.RH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.All.RH.meanf,data.SSP_SAP.All.RH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title('[All] RH Whisk-hemo coherence^2')
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Whisk Hemo Coherence' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'WhiskHemoCoherence_Gamma']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'WhiskHemoCoherence_Gamma'])
end

end
