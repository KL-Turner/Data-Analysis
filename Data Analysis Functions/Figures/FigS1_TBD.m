function [] = FigS1_TBD(rootFolder,saveFigs,delim)
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
data.briefWhisk.HbT = cat(1,data.briefWhisk.LH_HbT,data.briefWhisk.RH_HbT);
data.interWhisk.HbT = cat(1,data.interWhisk.LH_HbT,data.interWhisk.RH_HbT);
data.extendWhisk.HbT = cat(1,data.extendWhisk.LH_HbT,data.extendWhisk.RH_HbT);
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
subplot(4,4,1);
plot(timeVector,procData.blinkWhisk.mmArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.blinkWhisk.mmArea.mean + procData.blinkWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.blinkWhisk.mmArea.mean - procData.blinkWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaArea (mm^2)')
ylim([-0.1,0.3])
yyaxis right
plot(timeVector,procData.blinkWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.blinkWhisk.HbT.mean + procData.blinkWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.blinkWhisk.HbT.mean - procData.blinkWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
title('Short Whisks')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
subplot(4,4,2);
plot(timeVector,procData.briefWhisk.mmArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.mmArea.mean + procData.briefWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.mmArea.mean - procData.briefWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaArea (mm^2)')
ylim([-0.1,0.3])
yyaxis right
plot(timeVector,procData.briefWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.HbT.mean + procData.briefWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.HbT.mean - procData.briefWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom')
ylim([-5,15])
title('Brief Whisks')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% intermediate whisks
subplot(4,4,3);
plot(timeVector,procData.interWhisk.mmArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.mmArea.mean + procData.interWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.mmArea.mean - procData.interWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaArea (mm^2)')
ylim([-0.1,0.3])
yyaxis right
plot(timeVector,procData.interWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.HbT.mean + procData.interWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.HbT.mean - procData.interWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
title('Intermediate Whisks')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% extended whisks
subplot(4,4,4);
plot(timeVector,procData.extendWhisk.mmArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.mmArea.mean + procData.extendWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.mmArea.mean - procData.extendWhisk.mmArea.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaArea (mm^2)')
ylim([-0.1,0.3])
yyaxis right
plot(timeVector,procData.extendWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.HbT.mean + procData.extendWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.HbT.mean - procData.extendWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
title('Extended Whisks')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% figure characteristics
subplot(4,4,5);
plot(timeVector,procData.blinkWhisk.zArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.blinkWhisk.zArea.mean + procData.blinkWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.blinkWhisk.zArea.mean - procData.blinkWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaArea (Z Units)')
ylim([-1.5,5])
yyaxis right
plot(timeVector,procData.blinkWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.blinkWhisk.HbT.mean + procData.blinkWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.blinkWhisk.HbT.mean - procData.blinkWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
subplot(4,4,6);
plot(timeVector,procData.briefWhisk.zArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.zArea.mean + procData.briefWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.zArea.mean - procData.briefWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaArea (Z Units)')
ylim([-1.5,5])
yyaxis right
plot(timeVector,procData.briefWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.HbT.mean + procData.briefWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.HbT.mean - procData.briefWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% intermediate whisks
subplot(4,4,7);
plot(timeVector,procData.interWhisk.zArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.zArea.mean + procData.interWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.zArea.mean - procData.interWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaArea (Z Units)')
ylim([-1.5,5])
yyaxis right
plot(timeVector,procData.interWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.HbT.mean + procData.interWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.HbT.mean - procData.interWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% extended whisks
subplot(4,4,8);
plot(timeVector,procData.extendWhisk.zArea.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.zArea.mean + procData.extendWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.zArea.mean - procData.extendWhisk.zArea.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaArea (Z Units)')
ylim([-1.5,5])
yyaxis right
plot(timeVector,procData.extendWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.HbT.mean + procData.extendWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.HbT.mean - procData.extendWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% figure characteristics
subplot(4,4,9);
plot(timeVector,procData.blinkWhisk.mmDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.blinkWhisk.mmDiameter.mean + procData.blinkWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.blinkWhisk.mmDiameter.mean - procData.blinkWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaDiameter (mm)')
ylim([-0.1,0.3])
yyaxis right
plot(timeVector,procData.blinkWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.blinkWhisk.HbT.mean + procData.blinkWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.blinkWhisk.HbT.mean - procData.blinkWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
subplot(4,4,10);
plot(timeVector,procData.briefWhisk.mmDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.mmDiameter.mean + procData.briefWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.mmDiameter.mean - procData.briefWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaDiameter (mm)')
ylim([-0.1,0.3])
yyaxis right
plot(timeVector,procData.briefWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.HbT.mean + procData.briefWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.HbT.mean - procData.briefWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% intermediate whisks
subplot(4,4,11);
plot(timeVector,procData.interWhisk.mmDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.mmDiameter.mean + procData.interWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.mmDiameter.mean - procData.interWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaDiameter (mm)')
ylim([-0.1,0.3])
yyaxis right
plot(timeVector,procData.interWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.HbT.mean + procData.interWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.HbT.mean - procData.interWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% extended whisks
subplot(4,4,12);
plot(timeVector,procData.extendWhisk.mmDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.mmDiameter.mean + procData.extendWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.mmDiameter.mean - procData.extendWhisk.mmDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaDiameter (mm)')
ylim([-0.1,0.3])
yyaxis right
plot(timeVector,procData.extendWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.HbT.mean + procData.extendWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.HbT.mean - procData.extendWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% figure characteristics
subplot(4,4,13);
plot(timeVector,procData.blinkWhisk.zDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.blinkWhisk.zDiameter.mean + procData.blinkWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.blinkWhisk.zDiameter.mean - procData.blinkWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaDiameter (Z Units)')
ylim([-1.5,5])
yyaxis right
plot(timeVector,procData.blinkWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.blinkWhisk.HbT.mean + procData.blinkWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.blinkWhisk.HbT.mean - procData.blinkWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
subplot(4,4,14);
plot(timeVector,procData.briefWhisk.zDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.zDiameter.mean + procData.briefWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.zDiameter.mean - procData.briefWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaDiameter (Z Units)')
ylim([-1.5,5])
yyaxis right
plot(timeVector,procData.briefWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.briefWhisk.HbT.mean + procData.briefWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.briefWhisk.HbT.mean - procData.briefWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% intermediate whisks
subplot(4,4,15);
plot(timeVector,procData.interWhisk.zDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.zDiameter.mean + procData.interWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.zDiameter.mean - procData.interWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaDiameter (Z Units)')
ylim([-1.5,5])
yyaxis right
plot(timeVector,procData.interWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.HbT.mean + procData.interWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.HbT.mean - procData.interWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
% extended whisks
subplot(4,4,16);
plot(timeVector,procData.extendWhisk.zDiameter.mean,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.zDiameter.mean + procData.extendWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.zDiameter.mean - procData.extendWhisk.zDiameter.std,'color',colors('smoky black'),'LineWidth',0.5)
ylabel('\DeltaDiameter (Z Units)')
ylim([-1.5,5])
yyaxis right
plot(timeVector,procData.extendWhisk.HbT.mean,'color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(timeVector,procData.extendWhisk.HbT.mean + procData.extendWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.extendWhisk.HbT.mean - procData.extendWhisk.HbT.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
ylabel('\DeltaHbT (\muM)','rotation',-90,'VerticalAlignment','bottom') 
ylim([-5,15])
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
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

end
