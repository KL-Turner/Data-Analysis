function [] = Coherence_Pupil(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% variables for loops
resultsStruct = 'Results_Coherence';
load(resultsStruct);
animalIDs = fieldnames(Results_Coherence);
behavFields = {'Rest','NREM','REM','Awake','Asleep','All'};
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
colorRest = [(0/256),(166/256),(81/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
colorAlert = [(255/256),(191/256),(0/256)];
colorAsleep = [(0/256),(128/256),(255/256)];
colorAll = [(183/256),(115/256),(51/256)];
%% pre-allocate data structure
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.(behavField).(dataType).C = [];
        data.(behavField).(dataType).f = [];
    end
end
%% power spectra during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(Results_Coherence.(animalID).(behavField).(dataType).LH_HbT.C) == false
                data.(behavField).(dataType).C = cat(2,data.(behavField).(dataType).C,Results_Coherence.(animalID).(behavField).(dataType).LH_HbT.C,Results_Coherence.(animalID).(behavField).(dataType).RH_HbT.C);
                data.(behavField).(dataType).f = cat(1,data.(behavField).(dataType).f,Results_Coherence.(animalID).(behavField).(dataType).LH_HbT.f,Results_Coherence.(animalID).(behavField).(dataType).RH_HbT.f);
            end
        end
    end
end
%% take mean/StD of S/f
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.(behavField).(dataType).meanC = mean(data.(behavField).(dataType).C,2);
        data.(behavField).(dataType).stdC = std(data.(behavField).(dataType).C,0,2);
        data.(behavField).(dataType).meanf = mean(data.(behavField).(dataType).f,1);
    end
end
%% figures
summaryFigure = figure;
sgtitle('Coherence between HbT and Pupil data')
subplot(2,2,1);
L1 = semilogx(data.Rest.mmArea.meanf,data.Rest.mmArea.meanC,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.1,0.1 - 0.005,1],'FaceColor','w','EdgeColor','w')
L2 = semilogx(data.NREM.mmArea.meanf,data.NREM.mmArea.meanC,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.1,1/30 - 0.005,1],'FaceColor','w','EdgeColor','w')
L3 = semilogx(data.REM.mmArea.meanf,data.REM.mmArea.meanC,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.1,1/60 - 0.005,1],'FaceColor','w','EdgeColor','w')
L4 = semilogx(data.Awake.mmArea.meanf,data.Awake.mmArea.meanC,'color',colorAlert,'LineWidth',2);
L5 = semilogx(data.Asleep.mmArea.meanf,data.Asleep.mmArea.meanC,'color',colorAsleep,'LineWidth',2);
L6 = semilogx(data.All.mmArea.meanf,data.All.mmArea.meanC,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('mmArea')
ylabel('Coherence')
xlabel('Freq (Hz)')
legend([L1,L2,L3,L4,L5,L6],'Rest','NREM','REM','Alert','Asleep','All','Location','SouthEast')
% axis square
xlim([0.003,1])
ylim([0,1])
set(gca,'box','off')
subplot(2,2,2);
semilogx(data.Rest.mmDiameter.meanf,data.Rest.mmDiameter.meanC,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.1,0.1 - 0.005,1],'FaceColor','w','EdgeColor','w')
semilogx(data.NREM.mmDiameter.meanf,data.NREM.mmDiameter.meanC,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.1,1/30 - 0.005,1],'FaceColor','w','EdgeColor','w')
semilogx(data.REM.mmDiameter.meanf,data.REM.mmDiameter.meanC,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.1,1/60 - 0.005,1],'FaceColor','w','EdgeColor','w')
semilogx(data.Awake.mmDiameter.meanf,data.Awake.mmDiameter.meanC,'color',colorAlert,'LineWidth',2);
semilogx(data.Asleep.mmDiameter.meanf,data.Asleep.mmDiameter.meanC,'color',colorAsleep,'LineWidth',2);
semilogx(data.All.mmDiameter.meanf,data.All.mmDiameter.meanC,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('mmDiameter')
ylabel('Coherence')
xlabel('Freq (Hz)')
% axis square
xlim([0.003,1])
ylim([0,1])
set(gca,'box','off')
subplot(2,2,3);
semilogx(data.Rest.zArea.meanf,data.Rest.zArea.meanC,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.1,0.1 - 0.005,1],'FaceColor','w','EdgeColor','w')
semilogx(data.NREM.zArea.meanf,data.NREM.zArea.meanC,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.1,1/30 - 0.005,1],'FaceColor','w','EdgeColor','w')
semilogx(data.REM.zArea.meanf,data.REM.zArea.meanC,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.1,1/60 - 0.005,1],'FaceColor','w','EdgeColor','w')
semilogx(data.Awake.zArea.meanf,data.Awake.zArea.meanC,'color',colorAlert,'LineWidth',2);
semilogx(data.Asleep.zArea.meanf,data.Asleep.zArea.meanC,'color',colorAsleep,'LineWidth',2);
semilogx(data.All.zArea.meanf,data.All.zArea.meanC,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('zArea')
ylabel('Coherence')
xlabel('Freq (Hz)')
% axis square
xlim([0.003,1])
ylim([0,1])
set(gca,'box','off')
subplot(2,2,4);
semilogx(data.Rest.zDiameter.meanf,data.Rest.zDiameter.meanC,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.1,0.1 - 0.005,1],'FaceColor','w','EdgeColor','w')
semilogx(data.NREM.zDiameter.meanf,data.NREM.zDiameter.meanC,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.1,1/30 - 0.005,1],'FaceColor','w','EdgeColor','w')
semilogx(data.REM.zDiameter.meanf,data.REM.zDiameter.meanC,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.1,1/60 - 0.005,1],'FaceColor','w','EdgeColor','w')
semilogx(data.Awake.zDiameter.meanf,data.Awake.zDiameter.meanC,'color',colorAlert,'LineWidth',2);
semilogx(data.Asleep.zDiameter.meanf,data.Asleep.zDiameter.meanC,'color',colorAsleep,'LineWidth',2);
semilogx(data.All.zDiameter.meanf,data.All.zDiameter.meanC,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('zDiameter')
ylabel('Coherence')
xlabel('Freq (Hz)')
% axis square
xlim([0.003,1])
ylim([0,1])
set(gca,'box','off')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Pupil-HbT Coherence' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Pupil_Coherence']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Pupil_Coherence'])
end

end
