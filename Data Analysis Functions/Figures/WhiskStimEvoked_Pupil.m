function [] = WhiskStimEvoked_Pupil(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

resultsStruct = 'Results_Evoked';
load(resultsStruct);
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter','LH_HbT','RH_HbT'};
timeVector = (0:12*30)/30 - 2;
animalIDs = fieldnames(Results_Evoked);
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    % pre-allocate data types
        data.blinkWhisk.(dataType) = [];
    data.briefWhisk.(dataType) = [];
    data.interWhisk.(dataType) = [];
    data.extendWhisk.(dataType) = [];
        data.controlSolenoid.(dataType) = [];
    data.stimSolenoid.(dataType) = [];
end
%%
for cc = 1:length(animalIDs)
    animalID = animalIDs{cc,1};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        % whisking
        data.blinkWhisk.(dataType) = cat(1,data.blinkWhisk.(dataType),Results_Evoked.(animalID).Whisk.(dataType).BlinkWhisks.mean);
        data.briefWhisk.(dataType) = cat(1,data.briefWhisk.(dataType),Results_Evoked.(animalID).Whisk.(dataType).ShortWhisks.mean);
        data.interWhisk.(dataType) = cat(1,data.interWhisk.(dataType),Results_Evoked.(animalID).Whisk.(dataType).IntermediateWhisks.mean);
        data.extendWhisk.(dataType) = cat(1,data.extendWhisk.(dataType),Results_Evoked.(animalID).Whisk.(dataType).LongWhisks.mean);
        % solenoids
        if strcmp(dataType,'LH_HbT') == true
            data.stimSolenoid.(dataType) = cat(1,data.stimSolenoid.(dataType),Results_Evoked.(animalID).Stim.(dataType).RPadSol.mean);
        elseif strcmp(dataType,'RH_HbT') == true
            data.stimSolenoid.(dataType) = cat(1,data.stimSolenoid.(dataType),Results_Evoked.(animalID).Stim.(dataType).LPadSol.mean);
        else
            data.stimSolenoid.(dataType) = cat(1,data.stimSolenoid.(dataType),Results_Evoked.(animalID).Stim.(dataType).LPadSol.mean,Results_Evoked.(animalID).Stim.(dataType).RPadSol.mean);
        end
        data.controlSolenoid.(dataType) = cat(1,data.controlSolenoid.(dataType),Results_Evoked.(animalID).Stim.(dataType).AudSol.mean);
    end
end
data.blinkWhisk.HbT = cat(1,data.blinkWhisk.LH_HbT,data.blinkWhisk.RH_HbT);
data.briefWhisk.HbT = cat(1,data.briefWhisk.LH_HbT,data.blinkWhisk.RH_HbT);
data.interWhisk.HbT = cat(1,data.interWhisk.LH_HbT,data.blinkWhisk.RH_HbT);
data.extendWhisk.HbT = cat(1,data.interWhisk.LH_HbT,data.blinkWhisk.RH_HbT);
data.stimSolenoid.HbT = cat(1,data.stimSolenoid.LH_HbT,data.stimSolenoid.RH_HbT);
data.controlSolenoid.HbT = cat(1,data.controlSolenoid.LH_HbT,data.controlSolenoid.RH_HbT);
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter','HbT'};
for dd = 1:length(dataTypes)
    dataType = dataTypes{1,dd};
    procData.blinkWhisk.(dataType).mean = mean(data.blinkWhisk.(dataType),1);
    procData.blinkWhisk.(dataType).std = std(data.blinkWhisk.(dataType),0,1)./sqrt(size(data.blinkWhisk.(dataType),1));
    procData.briefWhisk.(dataType).mean = mean(data.briefWhisk.(dataType),1);
    procData.briefWhisk.(dataType).std = std(data.briefWhisk.(dataType),0,1)./sqrt(size(data.briefWhisk.(dataType),1));
    procData.interWhisk.(dataType).mean = mean(data.interWhisk.(dataType),1);
    procData.interWhisk.(dataType).std = std(data.interWhisk.(dataType),0,1)./sqrt(size(data.interWhisk.(dataType),1));
    procData.extendWhisk.(dataType).mean = mean(data.extendWhisk.(dataType),1);
    procData.extendWhisk.(dataType).std = std(data.extendWhisk.(dataType),0,1)./sqrt(size(data.extendWhisk.(dataType),1));   
    procData.stimSolenoid.(dataType).mean = mean(data.stimSolenoid.(dataType),1);
    procData.stimSolenoid.(dataType).std = std(data.stimSolenoid.(dataType),0,1)./sqrt(size(data.stimSolenoid.(dataType),1));
    procData.controlSolenoid.(dataType).mean = mean(data.controlSolenoid.(dataType),1);
    procData.controlSolenoid.(dataType).std = std(data.controlSolenoid.(dataType),0,1)./sqrt(size(data.controlSolenoid.(dataType),1));
end
%% average whisk-evoked figures
summaryFigure1 = figure;
sgtitle('Whisking evoked pupil response')
%% brief whisks
ax1 = subplot(4,4,1);
plot(timeVector,procData.blinkWhisk.mmArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.blinkWhisk.mmArea.mean + procData.blinkWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.blinkWhisk.mmArea.mean - procData.blinkWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Short Whisks')
ylabel('\DeltaArea (mm^2)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
ax2 = subplot(4,4,2);
plot(timeVector,procData.briefWhisk.mmArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.mmArea.mean + procData.briefWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.mmArea.mean - procData.briefWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Brief Whisks')
ylabel('\DeltaArea (mm^2)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% intermediate whisks
ax3 = subplot(4,4,3);
plot(timeVector,procData.interWhisk.mmArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.mmArea.mean + procData.interWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.mmArea.mean - procData.interWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Intermediate Whisks')
ylabel('\DeltaArea (mm^2)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% extended whisks
ax4 = subplot(4,4,4);
plot(timeVector,procData.extendWhisk.mmArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.mmArea.mean + procData.extendWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.mmArea.mean - procData.extendWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Extended Whisks')
ylabel('\DeltaArea (mm^2)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% figure characteristics
ax5 = subplot(4,4,5);
plot(timeVector,procData.blinkWhisk.zArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.blinkWhisk.zArea.mean + procData.blinkWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.blinkWhisk.zArea.mean - procData.blinkWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Brief Whisks')
ylabel('\DeltaZ Units')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
ax6 = subplot(4,4,6);
plot(timeVector,procData.briefWhisk.zArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.zArea.mean + procData.briefWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.zArea.mean - procData.briefWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Brief Whisks')
ylabel('\DeltaZ Units')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% intermediate whisks
ax7 = subplot(4,4,7);
plot(timeVector,procData.interWhisk.zArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.zArea.mean + procData.interWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.zArea.mean - procData.interWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Intermediate Whisks')
ylabel('\DeltaZ Units')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% extended whisks
ax8 = subplot(4,4,8);
plot(timeVector,procData.extendWhisk.zArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.zArea.mean + procData.extendWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.zArea.mean - procData.extendWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Extended Whisks')
ylabel('\DeltaZ Units')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% figure characteristics
ax9 = subplot(4,4,9);
plot(timeVector,procData.blinkWhisk.mmDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.blinkWhisk.mmDiameter.mean + procData.blinkWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.blinkWhisk.mmDiameter.mean - procData.blinkWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Brief Whisks')
ylabel('\DeltaDiameter (mm)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
ax10 = subplot(4,4,10);
plot(timeVector,procData.briefWhisk.mmDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.mmDiameter.mean + procData.briefWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.mmDiameter.mean - procData.briefWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Brief Whisks')
ylabel('\DeltaDiameter (mm)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% intermediate whisks
ax11 = subplot(4,4,11);
plot(timeVector,procData.interWhisk.mmDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.mmDiameter.mean + procData.interWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.mmDiameter.mean - procData.interWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Intermediate Whisks')
ylabel('\DeltaDiameter (mm)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% extended whisks
ax12 = subplot(4,4,12);
plot(timeVector,procData.extendWhisk.mmDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.mmDiameter.mean + procData.extendWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.mmDiameter.mean - procData.extendWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Extended Whisks')
ylabel('\DeltaDiameter (mm)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% figure characteristics
ax13 = subplot(4,4,13);
plot(timeVector,procData.blinkWhisk.zDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.blinkWhisk.zDiameter.mean + procData.blinkWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.blinkWhisk.zDiameter.mean - procData.blinkWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Brief Whisks')
ylabel('\DeltaZ Units')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
ax14 = subplot(4,4,14);
plot(timeVector,procData.briefWhisk.zDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.zDiameter.mean + procData.briefWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.zDiameter.mean - procData.briefWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Brief Whisks')
ylabel('\DeltaZ Units')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% intermediate whisks
ax15 = subplot(4,4,15);
plot(timeVector,procData.interWhisk.zDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.zDiameter.mean + procData.interWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.zDiameter.mean - procData.interWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Intermediate Whisks')
ylabel('\DeltaZ Units')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% extended whisks
ax16 = subplot(4,4,16);
plot(timeVector,procData.extendWhisk.zDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.zDiameter.mean + procData.extendWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.zDiameter.mean - procData.extendWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Extended Whisks')
ylabel('\DeltaZ Units')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% figure characteristics
linkaxes([ax1,ax2,ax3,ax4],'xy')
linkaxes([ax5,ax6,ax7,ax8],'xy')
linkaxes([ax9,ax10,ax11,ax12],'xy')
linkaxes([ax13,ax14,ax15,ax16],'xy')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Whisking Evoked Pupil Data' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'Pupil_WhiskingEvoked']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Pupil_WhiskingEvoked'])
end


%%

%% average whisk-evoked figures
summaryFigure2 = figure;
% whisker stim
ax1 = subplot(4,2,1);
plot(timeVector,procData.stimSolenoid.mmArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.mmArea.mean + procData.stimSolenoid.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.mmArea.mean - procData.stimSolenoid.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Whisker Stimulus')
ylabel('\DeltaArea (mm^2)')
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% auditory stim
ax2 = subplot(4,2,2);
plot(timeVector,procData.controlSolenoid.mmArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.mmArea.mean + procData.controlSolenoid.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.mmArea.mean - procData.controlSolenoid.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Auditory Stimulus')
ylabel('\DeltaArea (mm^2)')
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% whisker stim
ax3 = subplot(4,2,3);
plot(timeVector,procData.stimSolenoid.zArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.zArea.mean + procData.stimSolenoid.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.zArea.mean - procData.stimSolenoid.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Whisker Stimulus')
ylabel('\DeltaZ Units')
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% auditory stim
ax4 = subplot(4,2,4);
plot(timeVector,procData.controlSolenoid.zArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.zArea.mean + procData.controlSolenoid.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.zArea.mean - procData.controlSolenoid.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Auditory Stimulus')
ylabel('\DeltaZ Units')
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% whisker stim
ax5 = subplot(4,2,5);
plot(timeVector,procData.stimSolenoid.mmDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.mmDiameter.mean + procData.stimSolenoid.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.mmDiameter.mean - procData.stimSolenoid.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Whisker Stimulus')
ylabel('\DeltaDiameter (mm)')
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% auditory stim
ax6 = subplot(4,2,6);
plot(timeVector,procData.controlSolenoid.mmDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.mmDiameter.mean + procData.controlSolenoid.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.mmDiameter.mean - procData.controlSolenoid.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Auditory Stimulus')
ylabel('\DeltaDiameter (mm)')
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% whisker stim
ax7 = subplot(4,2,7);
plot(timeVector,procData.stimSolenoid.zDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.zDiameter.mean + procData.stimSolenoid.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.zDiameter.mean - procData.stimSolenoid.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Whisker Stimulus')
ylabel('\DeltaZ Units')
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% auditory stim
ax8 = subplot(4,2,8);
plot(timeVector,procData.controlSolenoid.zDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.zDiameter.mean + procData.controlSolenoid.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.zDiameter.mean - procData.controlSolenoid.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
title('Auditory Stimulus')
ylabel('\DeltaZ Units')
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% figure characteristics
linkaxes([ax1,ax2],'xy')
linkaxes([ax3,ax4],'xy')
linkaxes([ax5,ax6],'xy')
linkaxes([ax7,ax8],'xy')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Stimulus Evoked Pupil Data' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'Pupil_StimulusEvoked']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Pupil_StimulusEvoked'])
end

end
