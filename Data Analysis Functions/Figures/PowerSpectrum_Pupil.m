function [] = PowerSpectrum_Pupil(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% variables for loops
resultsStruct = 'Results_PowerSpectrum';
load(resultsStruct);
animalIDs = fieldnames(Results_PowerSpectrum);
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
        data.(behavField).(dataType).S = [];
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
            if isempty(Results_PowerSpectrum.(animalID).(behavField).(dataType).S) == false
                data.(behavField).(dataType).S = cat(2,data.(behavField).(dataType).S,Results_PowerSpectrum.(animalID).(behavField).(dataType).S);
                data.(behavField).(dataType).f = cat(1,data.(behavField).(dataType).f,Results_PowerSpectrum.(animalID).(behavField).(dataType).f);
            end
        end
    end
end
%% take mean/StD of S/f
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.(behavField).(dataType).meanS = mean(data.(behavField).(dataType).S,2);
        data.(behavField).(dataType).stdS = std(data.(behavField).(dataType).S,0,2);
        data.(behavField).(dataType).meanf = mean(data.(behavField).(dataType).f,1);
    end
end
%% figures
summaryFigure = figure;
sgtitle('Power Spectra of Pupil data')
subplot(2,2,1);
L1 = loglog(data.Rest.mmArea.meanf,data.Rest.mmArea.meanS,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,min(data.Rest.mmArea.meanS)*100,0.1 - 0.005,max(data.Rest.mmArea.meanS)*10],'FaceColor','w','EdgeColor','w')
L2 = loglog(data.NREM.mmArea.meanf,data.NREM.mmArea.meanS,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,min(data.NREM.mmArea.meanS)*100,1/30 - 0.005,max(data.NREM.mmArea.meanS)*10],'FaceColor','w','EdgeColor','w')
L3 = loglog(data.REM.mmArea.meanf,data.REM.mmArea.meanS,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,min(data.REM.mmArea.meanS)*100,1/60 - 0.005,max(data.REM.mmArea.meanS)*10],'FaceColor','w','EdgeColor','w')
L4 = loglog(data.Awake.mmArea.meanf,data.Awake.mmArea.meanS,'color',colorAlert,'LineWidth',2);
L5 = loglog(data.Asleep.mmArea.meanf,data.Asleep.mmArea.meanS,'color',colorAsleep,'LineWidth',2);
L6 = loglog(data.All.mmArea.meanf,data.All.mmArea.meanS,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('mmArea')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([L1,L2,L3,L4,L5,L6],'Rest','NREM','REM','Alert','Asleep','All','Location','NorthEast')
% axis square
axis tight
xlim([0.003,1])
set(gca,'box','off')
subplot(2,2,2);
loglog(data.Rest.mmDiameter.meanf,data.Rest.mmDiameter.meanS,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,min(data.Rest.mmDiameter.meanS)*100,0.1 - 0.005,max(data.Rest.mmDiameter.meanS)*10],'FaceColor','w','EdgeColor','w')
loglog(data.NREM.mmDiameter.meanf,data.NREM.mmDiameter.meanS,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,min(data.NREM.mmDiameter.meanS)*100,1/30 - 0.005,max(data.NREM.mmDiameter.meanS)*10],'FaceColor','w','EdgeColor','w')
loglog(data.REM.mmDiameter.meanf,data.REM.mmDiameter.meanS,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,min(data.REM.mmDiameter.meanS)*100,1/60 - 0.005,max(data.REM.mmDiameter.meanS)*10],'FaceColor','w','EdgeColor','w')
loglog(data.Awake.mmDiameter.meanf,data.Awake.mmDiameter.meanS,'color',colorAlert,'LineWidth',2);
loglog(data.Asleep.mmDiameter.meanf,data.Asleep.mmDiameter.meanS,'color',colorAsleep,'LineWidth',2);
loglog(data.All.mmDiameter.meanf,data.All.mmDiameter.meanS,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('mmDiameter')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
% axis square
axis tight
xlim([0.003,1])
set(gca,'box','off')
subplot(2,2,3);
loglog(data.Rest.zArea.meanf,data.Rest.zArea.meanS,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,min(data.Rest.zArea.meanS)*100,0.1 - 0.005,max(data.Rest.zArea.meanS)*10],'FaceColor','w','EdgeColor','w')
loglog(data.NREM.zArea.meanf,data.NREM.zArea.meanS,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,min(data.NREM.zArea.meanS)*100,1/30 - 0.005,max(data.NREM.zArea.meanS)*10],'FaceColor','w','EdgeColor','w')
loglog(data.REM.zArea.meanf,data.REM.zArea.meanS,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,min(data.REM.zArea.meanS)*100,1/60 - 0.005,max(data.REM.zArea.meanS)*10],'FaceColor','w','EdgeColor','w')
loglog(data.Awake.zArea.meanf,data.Awake.zArea.meanS,'color',colorAlert,'LineWidth',2);
loglog(data.Asleep.zArea.meanf,data.Asleep.zArea.meanS,'color',colorAsleep,'LineWidth',2);
loglog(data.All.zArea.meanf,data.All.zArea.meanS,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('zArea')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
% axis square
axis tight
xlim([0.003,1])
set(gca,'box','off')
subplot(2,2,4);
loglog(data.Rest.zDiameter.meanf,data.Rest.zDiameter.meanS,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,min(data.Rest.zDiameter.meanS)*100,0.1 - 0.005,max(data.Rest.zDiameter.meanS)*10],'FaceColor','w','EdgeColor','w')
loglog(data.NREM.zDiameter.meanf,data.NREM.zDiameter.meanS,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,min(data.NREM.zDiameter.meanS)*100,1/30 - 0.005,max(data.NREM.zDiameter.meanS)*10],'FaceColor','w','EdgeColor','w')
loglog(data.REM.zDiameter.meanf,data.REM.zDiameter.meanS,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,min(data.REM.zDiameter.meanS)*100,1/60 - 0.005,max(data.REM.zDiameter.meanS)*10],'FaceColor','w','EdgeColor','w')
loglog(data.Awake.zDiameter.meanf,data.Awake.zDiameter.meanS,'color',colorAlert,'LineWidth',2);
loglog(data.Asleep.zDiameter.meanf,data.Asleep.zDiameter.meanS,'color',colorAsleep,'LineWidth',2);
loglog(data.All.zDiameter.meanf,data.All.zDiameter.meanS,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('zDiameter')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
% axis square
axis tight
xlim([0.003,1])
set(gca,'box','off')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Pupil Power Spectrum' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Pupil_PowerSpectrum']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Pupil_PowerSpectrum'])
end

end
