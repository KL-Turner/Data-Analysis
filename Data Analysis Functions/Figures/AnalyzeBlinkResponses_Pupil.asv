function [] = AnalyzeBlinkResponses_Pupil(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

resultsStruct = 'Results_Evoked';
load(resultsStruct);
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
solenoidNames = {'LPadSol','RPadSol','AudSol'};
timeVector = (0:12*30)/30 - 2;
animalIDs = fieldnames(Results_Evoked);
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    % whisk lengths
    for bb = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,bb};
        data.(animalID).(whiskDataType).area = Results_Evoked.(animalID).Whisk.pupilArea.(whiskDataType).pupilArea.mean;
    end
    % solenoids
    for bb = 1:length(solenoidNames)
        solenoidName = solenoidNames{1,bb};
        data.(animalID).(solenoidName).area = Results_Evoked.(animalID).Stim.(solenoidName).pupilArea.mean;
    end
end
%
procData.controlSolenoid.area = []; procData.stimSolenoid.area = [];
procData.briefWhisk.area = []; procData.interWhisk.area = []; procData.extendWhisk.area = [];
for cc = 1:length(animalIDs)
    animalID = animalIDs{cc,1};
    % whisking
    procData.briefWhisk.area = cat(1,procData.briefWhisk.area,data.(animalID).ShortWhisks.area);
    procData.interWhisk.area = cat(1,procData.interWhisk.area,data.(animalID).IntermediateWhisks.area);
    procData.extendWhisk.area = cat(1,procData.extendWhisk.area,data.(animalID).LongWhisks.area);
    % solenoids
    procData.stimSolenoid.area = cat(1,procData.stimSolenoid.area,data.(animalID).LPadSol.area,data.(animalID).RPadSol.area);
    procData.controlSolenoid.area = cat(1,procData.controlSolenoid.area,data.(animalID).AudSol.area);
end
% 
procData.briefWhisk.meanArea = mean(procData.briefWhisk.area,1);
procData.briefWhisk.stdArea = std(procData.briefWhisk.area,0,1)./sqrt(size(procData.briefWhisk.area,1));
procData.interWhisk.meanArea = mean(procData.interWhisk.area,1);
procData.interWhisk.stdArea = std(procData.interWhisk.area,0,1)./sqrt(size(procData.interWhisk.area,1));
procData.extendWhisk.meanArea = mean(procData.extendWhisk.area,1);
procData.extendWhisk.stdArea = std(procData.extendWhisk.area,0,1)./sqrt(size(procData.extendWhisk.area,1));
procData.stimSolenoid.meanArea = mean(procData.stimSolenoid.area,1);
procData.stimSolenoid.stdArea = std(procData.stimSolenoid.area,0,1)./sqrt(size(procData.stimSolenoid.area,1));
procData.controlSolenoid.meanArea = mean(procData.controlSolenoid.area,1);
procData.controlSolenoid.stdArea = std(procData.controlSolenoid.area,0,1)./sqrt(size(procData.controlSolenoid.area,1));
% average whisk-evoked figures
summaryFigure1 = figure;
% brief whisks
ax1 = subplot(1,3,1);
% Naives
plot(timeVector,procData.briefWhisk.meanArea,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.meanArea + procData.briefWhisk.stdArea,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.meanArea - procData.briefWhisk.stdArea,'color',colors('smoky black'),'LineWidth',0.5)
title('Brief Whisks')
ylabel('\DeltaArea (pixels)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% intermediate whisks
ax2 = subplot(1,3,2);
plot(timeVector,procData.interWhisk.meanArea,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.meanArea + procData.interWhisk.stdArea,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.meanArea - procData.interWhisk.stdArea,'color',colors('smoky black'),'LineWidth',0.5)
title('Intermediate Whisks')
ylabel('\DeltaArea (pixels)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% extended whisks
ax3 = subplot(1,3,3);
plot(timeVector,procData.extendWhisk.meanArea,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.meanArea + procData.extendWhisk.stdArea,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.meanArea - procData.extendWhisk.stdArea,'color',colors('smoky black'),'LineWidth',0.5)
title('Extended Whisks')
ylabel('\DeltaArea (pixels)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% figure characteristics
linkaxes([ax1,ax2,ax3],'xy')
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Whisking-evoked Pupil Area' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'Whisking_PupilArea']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Whisking_PupilArea'])
end
% average whisk-evoked figures
summaryFigure2 = figure;
% brief whisks
ax4 = subplot(1,2,1);
% Naives
plot(timeVector,procData.stimSolenoid.meanArea,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.meanArea + procData.stimSolenoid.stdArea,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.meanArea - procData.stimSolenoid.stdArea,'color',colors('smoky black'),'LineWidth',0.5)
title('Whisker Stimulus')
ylabel('\DeltaArea (pixels)')
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% intermediate whisks
ax5 = subplot(1,2,2);
plot(timeVector,procData.controlSolenoid.meanArea,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.meanArea + procData.controlSolenoid.stdArea,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.meanArea - procData.controlSolenoid.stdArea,'color',colors('smoky black'),'LineWidth',0.5)
title('Auditory Stimulus')
ylabel('\DeltaArea (pixels)')
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% figure characteristics
linkaxes([ax4,ax5],'xy')
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Stimulus-evoked Pupil Area' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'Stimulus_PupilArea']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Stimulus_PupilArea'])
end

end
