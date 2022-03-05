function [] = FigS2_TBD(rootFolder,saveFigs,delim)
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
%% average stimulus-evoked figures
summaryFigure2 = figure;
% whisker stim
subplot(4,2,1);
plot(timeVector,procData.stimSolenoid.mmArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.mmArea.mean + procData.stimSolenoid.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.mmArea.mean - procData.stimSolenoid.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaArea (mm^2)')
ylim([-0.05,0.15])
yyaxis right
plot(timeVector,procData.stimSolenoid.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.HbT.mean + procData.stimSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.HbT.mean - procData.stimSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-7,20])
title('Whisker Stimulus')
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% auditory stim
subplot(4,2,2);
plot(timeVector,procData.controlSolenoid.mmArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.mmArea.mean + procData.controlSolenoid.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.mmArea.mean - procData.controlSolenoid.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaArea (mm^2)')
ylim([-0.05,0.15])
yyaxis right
plot(timeVector,procData.controlSolenoid.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.HbT.mean + procData.controlSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.HbT.mean - procData.controlSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-7,20])
title('Auditory Stimulus')
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% whisker stim
subplot(4,2,3);
plot(timeVector,procData.stimSolenoid.zArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.zArea.mean + procData.stimSolenoid.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.zArea.mean - procData.stimSolenoid.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaArea (Z Units)')
ylim([-0.7,2.2])
yyaxis right
plot(timeVector,procData.stimSolenoid.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.HbT.mean + procData.stimSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.HbT.mean - procData.stimSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-7,20])
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% auditory stim
subplot(4,2,4);
plot(timeVector,procData.controlSolenoid.zArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.zArea.mean + procData.controlSolenoid.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.zArea.mean - procData.controlSolenoid.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaArea (Z Units)')
ylim([-0.7,2.2])
yyaxis right
plot(timeVector,procData.controlSolenoid.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.HbT.mean + procData.controlSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.HbT.mean - procData.controlSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-7,20])
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% whisker stim
subplot(4,2,5);
plot(timeVector,procData.stimSolenoid.mmDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.mmDiameter.mean + procData.stimSolenoid.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.mmDiameter.mean - procData.stimSolenoid.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaDiameter (mm)')
ylim([-0.05,0.15])
yyaxis right
plot(timeVector,procData.stimSolenoid.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.HbT.mean + procData.stimSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.HbT.mean - procData.stimSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-7,20])
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% auditory stim
subplot(4,2,6);
plot(timeVector,procData.controlSolenoid.mmDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.mmDiameter.mean + procData.controlSolenoid.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.mmDiameter.mean - procData.controlSolenoid.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaDiameter (mm)')
ylim([-0.05,0.15])
yyaxis right
plot(timeVector,procData.controlSolenoid.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.HbT.mean + procData.controlSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.HbT.mean - procData.controlSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-7,20])
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% whisker stim
subplot(4,2,7);
plot(timeVector,procData.stimSolenoid.zDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.zDiameter.mean + procData.stimSolenoid.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.zDiameter.mean - procData.stimSolenoid.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaDiameter (Z Units)')
ylim([-0.7,2.2])
yyaxis right
plot(timeVector,procData.stimSolenoid.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.stimSolenoid.HbT.mean + procData.stimSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.HbT.mean - procData.stimSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-7,20])
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% auditory stim
subplot(4,2,8);
plot(timeVector,procData.controlSolenoid.zDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.zDiameter.mean + procData.controlSolenoid.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.zDiameter.mean - procData.controlSolenoid.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaDiameter (Z Units)')
ylim([-0.7,2.2])
yyaxis right
plot(timeVector,procData.controlSolenoid.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.HbT.mean + procData.controlSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.HbT.mean - procData.controlSolenoid.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-7,20])
xlabel('Peri-stim time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
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
