function [AnalysisResults] = PowerSpec_Saporin(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: 
%________________________________________________________________________________________________________________________

%% set-up and process data
expGroups = {'C57BL6J','SSP-SAP','Blank-SAP'};
setName = 'IOS Set A';
animalIDs.all = {};
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs.all = horzcat(animalIDs.all,{folderList.name});
    animalIDs.(strrep(expGroups{1,aa},'-','_')) = {folderList.name};
end
treatments = {'C57BL6J','SSP_SAP','Blank_SAP'};
behavFields = {'Rest','NREM','REM','Awake','Sleep','All'};
dataTypes = {'CBV_HbT','gammaBandPower'};
%% power spectra during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,aa};
    % recognize treatment based on animal group
    if ismember(animalID,animalIDs.C57BL6J) == true
        treatment = 'C57BL6J';
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
                data.(treatment).(behavField).(dataType).adjLH.S = [];
                data.(treatment).(behavField).(dataType).adjLH.f = [];
                data.(treatment).(behavField).(dataType).adjRH.S = [];
                data.(treatment).(behavField).(dataType).adjRH.f = [];
            end
            data.(treatment).(behavField).(dataType).adjLH.S = cat(2,data.(treatment).(behavField).(dataType).adjLH.S,AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.S);
            data.(treatment).(behavField).(dataType).adjLH.f = cat(1,data.(treatment).(behavField).(dataType).adjLH.f,AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.f);
            data.(treatment).(behavField).(dataType).adjRH.S = cat(2,data.(treatment).(behavField).(dataType).adjRH.S,AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.S);
            data.(treatment).(behavField).(dataType).adjRH.f = cat(1,data.(treatment).(behavField).(dataType).adjRH.f,AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.f);
        end
    end
end
%% find the peak of the resting PSD for each animal/hemisphere
for aa = 1:length(treatments)
    treatment = treatments{1,aa};
    for cc = 1:length(dataTypes)
        dataType = dataTypes{1,cc};
        for ee = 1:size(data.(treatment).Rest.(dataType).adjLH.S,2)
            data.(treatment).baseline.(dataType).LH(ee,1) = max(data.(treatment).Rest.(dataType).adjLH.S(:,ee));
            data.(treatment).baseline.(dataType).RH(ee,1) = max(data.(treatment).Rest.(dataType).adjRH.S(:,ee));
        end
    end
end
%% DC-shift each animal/hemisphere/behavior PSD with respect to the resting peak
for aa = 1:length(treatments)
    treatment = treatments{1,aa};
    for dd = 1:length(behavFields)
        behavField = behavFields{1,dd};
        for jj = 1:length(dataTypes)
            dataType = dataTypes{1,jj};
            for ee = 1:size(data.(treatment).(behavField).(dataType).adjLH.S,2)
                data.(treatment).(behavField).(dataType).adjLH.normS(:,ee) = (data.(treatment).(behavField).(dataType).adjLH.S(:,ee))*(1/(data.(treatment).baseline.(dataType).LH(ee,1)));
                data.(treatment).(behavField).(dataType).adjRH.normS(:,ee) = (data.(treatment).(behavField).(dataType).adjRH.S(:,ee))*(1/(data.(treatment).baseline.(dataType).RH(ee,1)));
            end
        end
    end
end
%% take mean/StD of S/f
for aa = 1:length(treatments)
    treatment = treatments{1,aa};
    for h = 1:length(behavFields)
        behavField = behavFields{1,h};
        for jj = 1:length(dataTypes)
            dataType = dataTypes{1,jj};
            data.(treatment).(behavField).(dataType).adjLH.meanCortS = mean(data.(treatment).(behavField).(dataType).adjLH.normS,2);
            data.(treatment).(behavField).(dataType).adjLH.stdCortS = std(data.(treatment).(behavField).(dataType).adjLH.normS,0,2);
            data.(treatment).(behavField).(dataType).adjLH.meanCortf = mean(data.(treatment).(behavField).(dataType).adjLH.f,1);
            data.(treatment).(behavField).(dataType).adjRH.meanCortS = mean(data.(treatment).(behavField).(dataType).adjRH.normS,2);
            data.(treatment).(behavField).(dataType).adjRH.stdCortS = std(data.(treatment).(behavField).(dataType).adjRH.normS,0,2);
            data.(treatment).(behavField).(dataType).adjRH.meanCortf = mean(data.(treatment).(behavField).(dataType).adjRH.f,1);
        end
    end
end
%% average HbT power
summaryFigure1 = figure;
sgtitle('\DeltaHbT (\muM) cortical power spectra')
%% LH power spectra of HbT power during Rest
ax1 = subplot(3,4,1);
p1 = loglog(data.C57BL6J.Rest.CBV_HbT.adjLH.meanCortf,data.C57BL6J.Rest.CBV_HbT.adjLH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
p2 = loglog(data.Blank_SAP.Rest.CBV_HbT.adjLH.meanCortf,data.Blank_SAP.Rest.CBV_HbT.adjLH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
p3 = loglog(data.SSP_SAP.Rest.CBV_HbT.adjLH.meanCortf,data.SSP_SAP.Rest.CBV_HbT.adjLH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Rest] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
xlim([1/10,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during Rest
ax2 = subplot(3,4,2);
loglog(data.C57BL6J.Rest.CBV_HbT.adjRH.meanCortf,data.C57BL6J.Rest.CBV_HbT.adjRH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Rest.CBV_HbT.adjRH.meanCortf,data.Blank_SAP.Rest.CBV_HbT.adjRH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Rest.CBV_HbT.adjRH.meanCortf,data.SSP_SAP.Rest.CBV_HbT.adjRH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Rest] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/10,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during NREM
ax3 = subplot(3,4,3);
loglog(data.C57BL6J.NREM.CBV_HbT.adjLH.meanCortf,data.C57BL6J.NREM.CBV_HbT.adjLH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.NREM.CBV_HbT.adjLH.meanCortf,data.Blank_SAP.NREM.CBV_HbT.adjLH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.NREM.CBV_HbT.adjLH.meanCortf,data.SSP_SAP.NREM.CBV_HbT.adjLH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[NREM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
xlim([1/30,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during NREM
ax4 = subplot(3,4,4);
loglog(data.C57BL6J.NREM.CBV_HbT.adjRH.meanCortf,data.C57BL6J.NREM.CBV_HbT.adjRH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.NREM.CBV_HbT.adjRH.meanCortf,data.Blank_SAP.NREM.CBV_HbT.adjRH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.NREM.CBV_HbT.adjRH.meanCortf,data.SSP_SAP.NREM.CBV_HbT.adjRH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[NREM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/30,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during REM
ax5 = subplot(3,4,5);
loglog(data.C57BL6J.REM.CBV_HbT.adjLH.meanCortf,data.C57BL6J.REM.CBV_HbT.adjLH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.REM.CBV_HbT.adjLH.meanCortf,data.Blank_SAP.REM.CBV_HbT.adjLH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.REM.CBV_HbT.adjLH.meanCortf,data.SSP_SAP.REM.CBV_HbT.adjLH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[REM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
xlim([1/60,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during REM
ax6 = subplot(3,4,6);
loglog(data.C57BL6J.REM.CBV_HbT.adjRH.meanCortf,data.C57BL6J.REM.CBV_HbT.adjRH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.REM.CBV_HbT.adjRH.meanCortf,data.Blank_SAP.REM.CBV_HbT.adjRH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.REM.CBV_HbT.adjRH.meanCortf,data.SSP_SAP.REM.CBV_HbT.adjRH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[REM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/60,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during Awake
ax7 = subplot(3,4,7);
loglog(data.C57BL6J.Awake.CBV_HbT.adjLH.meanCortf,data.C57BL6J.Awake.CBV_HbT.adjLH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Awake.CBV_HbT.adjLH.meanCortf,data.Blank_SAP.Awake.CBV_HbT.adjLH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Awake.CBV_HbT.adjLH.meanCortf,data.SSP_SAP.Awake.CBV_HbT.adjLH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Alert] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during Awake
ax8 = subplot(3,4,8);
loglog(data.C57BL6J.Awake.CBV_HbT.adjRH.meanCortf,data.C57BL6J.Awake.CBV_HbT.adjRH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Awake.CBV_HbT.adjRH.meanCortf,data.Blank_SAP.Awake.CBV_HbT.adjRH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Awake.CBV_HbT.adjRH.meanCortf,data.SSP_SAP.Awake.CBV_HbT.adjRH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Alert] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during Sleep
ax9 = subplot(3,4,9);
loglog(data.C57BL6J.Sleep.CBV_HbT.adjLH.meanCortf,data.C57BL6J.Sleep.CBV_HbT.adjLH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Sleep.CBV_HbT.adjLH.meanCortf,data.Blank_SAP.Sleep.CBV_HbT.adjLH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Sleep.CBV_HbT.adjLH.meanCortf,data.SSP_SAP.Sleep.CBV_HbT.adjLH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Asleep] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during Sleep
ax10 = subplot(3,4,10);
loglog(data.C57BL6J.Sleep.CBV_HbT.adjRH.meanCortf,data.C57BL6J.Sleep.CBV_HbT.adjRH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Sleep.CBV_HbT.adjRH.meanCortf,data.Blank_SAP.Sleep.CBV_HbT.adjRH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Sleep.CBV_HbT.adjRH.meanCortf,data.SSP_SAP.Sleep.CBV_HbT.adjRH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Asleep] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during All data
ax11 = subplot(3,4,11);
loglog(data.C57BL6J.All.CBV_HbT.adjLH.meanCortf,data.C57BL6J.All.CBV_HbT.adjLH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.All.CBV_HbT.adjLH.meanCortf,data.Blank_SAP.All.CBV_HbT.adjLH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.All.CBV_HbT.adjLH.meanCortf,data.SSP_SAP.All.CBV_HbT.adjLH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[All] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during All data
ax12 = subplot(3,4,12);
loglog(data.C57BL6J.All.CBV_HbT.adjRH.meanCortf,data.C57BL6J.All.CBV_HbT.adjRH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.All.CBV_HbT.adjRH.meanCortf,data.Blank_SAP.All.CBV_HbT.adjRH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.All.CBV_HbT.adjRH.meanCortf,data.SSP_SAP.All.CBV_HbT.adjRH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[All] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% figure characteristics
linkaxes([ax1,ax2],'y')
linkaxes([ax3,ax4],'y')
linkaxes([ax5,ax6],'y')
linkaxes([ax7,ax8],'y')
linkaxes([ax9,ax10],'y')
linkaxes([ax11,ax12],'y')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'PowerSpec_HbT']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'PowerSpec_HbT'])
end
%% individual HbT power
summaryFigure2 = figure;
sgtitle('\DeltaHbT (\muM) cortical power spectra - individual animals')
%% LH power spectra of HbT power during Rest
ax1 = subplot(3,4,1);
% C57BL6J
for aa = 1:size(data.C57BL6J.Rest.CBV_HbT.adjLH.normS,2)
    loglog(data.C57BL6J.Rest.CBV_HbT.adjLH.meanCortf,data.C57BL6J.Rest.CBV_HbT.adjLH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Rest.CBV_HbT.adjLH.normS,2)
    loglog(data.Blank_SAP.Rest.CBV_HbT.adjLH.meanCortf,data.Blank_SAP.Rest.CBV_HbT.adjLH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Rest.CBV_HbT.adjLH.normS,2)
    loglog(data.SSP_SAP.Rest.CBV_HbT.adjLH.meanCortf,data.SSP_SAP.Rest.CBV_HbT.adjLH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Rest] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/10,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during Rest
ax2 = subplot(3,4,2);
% C57BL6J
for aa = 1:size(data.C57BL6J.Rest.CBV_HbT.adjRH.normS,2)
    loglog(data.C57BL6J.Rest.CBV_HbT.adjRH.meanCortf,data.C57BL6J.Rest.CBV_HbT.adjRH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Rest.CBV_HbT.adjRH.normS,2)
    loglog(data.Blank_SAP.Rest.CBV_HbT.adjRH.meanCortf,data.Blank_SAP.Rest.CBV_HbT.adjRH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Rest.CBV_HbT.adjRH.normS,2)
    loglog(data.SSP_SAP.Rest.CBV_HbT.adjRH.meanCortf,data.SSP_SAP.Rest.CBV_HbT.adjRH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Rest] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/10,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during NREM
ax3 = subplot(3,4,3);
% C57BL6J
for aa = 1:size(data.C57BL6J.NREM.CBV_HbT.adjLH.normS,2)
    loglog(data.C57BL6J.NREM.CBV_HbT.adjLH.meanCortf,data.C57BL6J.NREM.CBV_HbT.adjLH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.NREM.CBV_HbT.adjLH.normS,2)
    loglog(data.Blank_SAP.NREM.CBV_HbT.adjLH.meanCortf,data.Blank_SAP.NREM.CBV_HbT.adjLH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.NREM.CBV_HbT.adjLH.normS,2)
    loglog(data.SSP_SAP.NREM.CBV_HbT.adjLH.meanCortf,data.SSP_SAP.NREM.CBV_HbT.adjLH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[NREM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/30,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during NREM
ax4 = subplot(3,4,4);
% C57BL6J
for aa = 1:size(data.C57BL6J.NREM.CBV_HbT.adjRH.normS,2)
    loglog(data.C57BL6J.NREM.CBV_HbT.adjRH.meanCortf,data.C57BL6J.NREM.CBV_HbT.adjRH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.NREM.CBV_HbT.adjRH.normS,2)
    loglog(data.Blank_SAP.NREM.CBV_HbT.adjRH.meanCortf,data.Blank_SAP.NREM.CBV_HbT.adjRH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.NREM.CBV_HbT.adjRH.normS,2)
    loglog(data.SSP_SAP.NREM.CBV_HbT.adjRH.meanCortf,data.SSP_SAP.NREM.CBV_HbT.adjRH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[NREM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/30,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during REM
ax5 = subplot(3,4,5);
% C57BL6J
for aa = 1:size(data.C57BL6J.REM.CBV_HbT.adjLH.normS,2)
    loglog(data.C57BL6J.REM.CBV_HbT.adjLH.meanCortf,data.C57BL6J.REM.CBV_HbT.adjLH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.REM.CBV_HbT.adjLH.normS,2)
    loglog(data.Blank_SAP.REM.CBV_HbT.adjLH.meanCortf,data.Blank_SAP.REM.CBV_HbT.adjLH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.REM.CBV_HbT.adjLH.normS,2)
    loglog(data.SSP_SAP.REM.CBV_HbT.adjLH.meanCortf,data.SSP_SAP.REM.CBV_HbT.adjLH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[REM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/60,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during REM
ax6 = subplot(3,4,6);
% C57BL6J
for aa = 1:size(data.C57BL6J.REM.CBV_HbT.adjRH.normS,2)
    loglog(data.C57BL6J.REM.CBV_HbT.adjRH.meanCortf,data.C57BL6J.REM.CBV_HbT.adjRH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.REM.CBV_HbT.adjRH.normS,2)
    loglog(data.Blank_SAP.REM.CBV_HbT.adjRH.meanCortf,data.Blank_SAP.REM.CBV_HbT.adjRH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.REM.CBV_HbT.adjRH.normS,2)
    loglog(data.SSP_SAP.REM.CBV_HbT.adjRH.meanCortf,data.SSP_SAP.REM.CBV_HbT.adjRH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[REM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/60,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during Awake
ax7 = subplot(3,4,7);
% C57BL6J
for aa = 1:size(data.C57BL6J.Awake.CBV_HbT.adjLH.normS,2)
    loglog(data.C57BL6J.Awake.CBV_HbT.adjLH.meanCortf,data.C57BL6J.Awake.CBV_HbT.adjLH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Awake.CBV_HbT.adjLH.normS,2)
    loglog(data.Blank_SAP.Awake.CBV_HbT.adjLH.meanCortf,data.Blank_SAP.Awake.CBV_HbT.adjLH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Awake.CBV_HbT.adjLH.normS,2)
    loglog(data.SSP_SAP.Awake.CBV_HbT.adjLH.meanCortf,data.SSP_SAP.Awake.CBV_HbT.adjLH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Alert] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during Awake
ax8 = subplot(3,4,8);
% C57BL6J
for aa = 1:size(data.C57BL6J.Awake.CBV_HbT.adjRH.normS,2)
    loglog(data.C57BL6J.Awake.CBV_HbT.adjRH.meanCortf,data.C57BL6J.Awake.CBV_HbT.adjRH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Awake.CBV_HbT.adjRH.normS,2)
    loglog(data.Blank_SAP.Awake.CBV_HbT.adjRH.meanCortf,data.Blank_SAP.Awake.CBV_HbT.adjRH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Awake.CBV_HbT.adjRH.normS,2)
    loglog(data.SSP_SAP.Awake.CBV_HbT.adjRH.meanCortf,data.SSP_SAP.Awake.CBV_HbT.adjRH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Alert] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during Sleep
ax9 = subplot(3,4,9);
% C57BL6J
for aa = 1:size(data.C57BL6J.Sleep.CBV_HbT.adjLH.normS,2)
    loglog(data.C57BL6J.Sleep.CBV_HbT.adjLH.meanCortf,data.C57BL6J.Sleep.CBV_HbT.adjLH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Sleep.CBV_HbT.adjLH.normS,2)
    loglog(data.Blank_SAP.Sleep.CBV_HbT.adjLH.meanCortf,data.Blank_SAP.Sleep.CBV_HbT.adjLH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Sleep.CBV_HbT.adjLH.normS,2)
    loglog(data.SSP_SAP.Sleep.CBV_HbT.adjLH.meanCortf,data.SSP_SAP.Sleep.CBV_HbT.adjLH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Asleep] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during Sleep
ax10 = subplot(3,4,10);
% C57BL6J
for aa = 1:size(data.C57BL6J.Sleep.CBV_HbT.adjRH.normS,2)
    loglog(data.C57BL6J.Sleep.CBV_HbT.adjRH.meanCortf,data.C57BL6J.Sleep.CBV_HbT.adjRH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Sleep.CBV_HbT.adjRH.normS,2)
    loglog(data.Blank_SAP.Sleep.CBV_HbT.adjRH.meanCortf,data.Blank_SAP.Sleep.CBV_HbT.adjRH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Sleep.CBV_HbT.adjRH.normS,2)
    loglog(data.SSP_SAP.Sleep.CBV_HbT.adjRH.meanCortf,data.SSP_SAP.Sleep.CBV_HbT.adjRH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Asleep] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during All data
ax11 = subplot(3,4,11);
% C57BL6J
for aa = 1:size(data.C57BL6J.All.CBV_HbT.adjLH.normS,2)
    loglog(data.C57BL6J.All.CBV_HbT.adjLH.meanCortf,data.C57BL6J.All.CBV_HbT.adjLH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.All.CBV_HbT.adjLH.normS,2)
    loglog(data.Blank_SAP.All.CBV_HbT.adjLH.meanCortf,data.Blank_SAP.All.CBV_HbT.adjLH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.All.CBV_HbT.adjLH.normS,2)
    loglog(data.SSP_SAP.All.CBV_HbT.adjLH.meanCortf,data.SSP_SAP.All.CBV_HbT.adjLH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[All] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during All data
ax12 = subplot(3,4,12);
% C57BL6J
for aa = 1:size(data.C57BL6J.All.CBV_HbT.adjRH.normS,2)
    loglog(data.C57BL6J.All.CBV_HbT.adjRH.meanCortf,data.C57BL6J.All.CBV_HbT.adjRH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.All.CBV_HbT.adjRH.normS,2)
    loglog(data.Blank_SAP.All.CBV_HbT.adjRH.meanCortf,data.Blank_SAP.All.CBV_HbT.adjRH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.All.CBV_HbT.adjRH.normS,2)
    loglog(data.SSP_SAP.All.CBV_HbT.adjRH.meanCortf,data.SSP_SAP.All.CBV_HbT.adjRH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[All] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% figure characteristics
linkaxes([ax1,ax2],'y')
linkaxes([ax3,ax4],'y')
linkaxes([ax5,ax6],'y')
linkaxes([ax7,ax8],'y')
linkaxes([ax9,ax10],'y')
linkaxes([ax11,ax12],'y')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'indPowerSpec_HbT']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'indPowerSpec_HbT'])
end
%% average gamma-band power
summaryFigure3 = figure;
sgtitle('Gamma-band [30-100] Hz (envelope) cortical power spectra')
%% LH power spectra of gamma-band power during Rest
ax1 = subplot(3,4,1);
p1 = loglog(data.C57BL6J.Rest.gammaBandPower.adjLH.meanCortf,data.C57BL6J.Rest.gammaBandPower.adjLH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
p2 = loglog(data.Blank_SAP.Rest.gammaBandPower.adjLH.meanCortf,data.Blank_SAP.Rest.gammaBandPower.adjLH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
p3 = loglog(data.SSP_SAP.Rest.gammaBandPower.adjLH.meanCortf,data.SSP_SAP.Rest.gammaBandPower.adjLH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Rest] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
xlim([1/10,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Rest
ax2 = subplot(3,4,2);
loglog(data.C57BL6J.Rest.gammaBandPower.adjRH.meanCortf,data.C57BL6J.Rest.gammaBandPower.adjRH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Rest.gammaBandPower.adjRH.meanCortf,data.Blank_SAP.Rest.gammaBandPower.adjRH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Rest.gammaBandPower.adjRH.meanCortf,data.SSP_SAP.Rest.gammaBandPower.adjRH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Rest] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/10,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during NREM
ax3 = subplot(3,4,3);
loglog(data.C57BL6J.NREM.gammaBandPower.adjLH.meanCortf,data.C57BL6J.NREM.gammaBandPower.adjLH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.NREM.gammaBandPower.adjLH.meanCortf,data.Blank_SAP.NREM.gammaBandPower.adjLH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.NREM.gammaBandPower.adjLH.meanCortf,data.SSP_SAP.NREM.gammaBandPower.adjLH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[NREM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
xlim([1/30,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during NREM
ax4 = subplot(3,4,4);
loglog(data.C57BL6J.NREM.gammaBandPower.adjRH.meanCortf,data.C57BL6J.NREM.gammaBandPower.adjRH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.NREM.gammaBandPower.adjRH.meanCortf,data.Blank_SAP.NREM.gammaBandPower.adjRH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.NREM.gammaBandPower.adjRH.meanCortf,data.SSP_SAP.NREM.gammaBandPower.adjRH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[NREM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/30,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during REM
ax5 = subplot(3,4,5);
loglog(data.C57BL6J.REM.gammaBandPower.adjLH.meanCortf,data.C57BL6J.REM.gammaBandPower.adjLH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.REM.gammaBandPower.adjLH.meanCortf,data.Blank_SAP.REM.gammaBandPower.adjLH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.REM.gammaBandPower.adjLH.meanCortf,data.SSP_SAP.REM.gammaBandPower.adjLH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[REM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
xlim([1/60,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during REM
ax6 = subplot(3,4,6);
loglog(data.C57BL6J.REM.gammaBandPower.adjRH.meanCortf,data.C57BL6J.REM.gammaBandPower.adjRH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.REM.gammaBandPower.adjRH.meanCortf,data.Blank_SAP.REM.gammaBandPower.adjRH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.REM.gammaBandPower.adjRH.meanCortf,data.SSP_SAP.REM.gammaBandPower.adjRH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[REM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/60,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during Awake
ax7 = subplot(3,4,7);
loglog(data.C57BL6J.Awake.gammaBandPower.adjLH.meanCortf,data.C57BL6J.Awake.gammaBandPower.adjLH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Awake.gammaBandPower.adjLH.meanCortf,data.Blank_SAP.Awake.gammaBandPower.adjLH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Awake.gammaBandPower.adjLH.meanCortf,data.SSP_SAP.Awake.gammaBandPower.adjLH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Alert] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Awake
ax8 = subplot(3,4,8);
loglog(data.C57BL6J.Awake.gammaBandPower.adjRH.meanCortf,data.C57BL6J.Awake.gammaBandPower.adjRH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Awake.gammaBandPower.adjRH.meanCortf,data.Blank_SAP.Awake.gammaBandPower.adjRH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Awake.gammaBandPower.adjRH.meanCortf,data.SSP_SAP.Awake.gammaBandPower.adjRH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Alert] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during Sleep
ax9 = subplot(3,4,9);
loglog(data.C57BL6J.Sleep.gammaBandPower.adjLH.meanCortf,data.C57BL6J.Sleep.gammaBandPower.adjLH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Sleep.gammaBandPower.adjLH.meanCortf,data.Blank_SAP.Sleep.gammaBandPower.adjLH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Sleep.gammaBandPower.adjLH.meanCortf,data.SSP_SAP.Sleep.gammaBandPower.adjLH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Asleep] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Sleep
ax10 = subplot(3,4,10);
loglog(data.C57BL6J.Sleep.gammaBandPower.adjRH.meanCortf,data.C57BL6J.Sleep.gammaBandPower.adjRH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Sleep.gammaBandPower.adjRH.meanCortf,data.Blank_SAP.Sleep.gammaBandPower.adjRH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Sleep.gammaBandPower.adjRH.meanCortf,data.SSP_SAP.Sleep.gammaBandPower.adjRH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Asleep] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during All data
ax11 = subplot(3,4,11);
loglog(data.C57BL6J.All.gammaBandPower.adjLH.meanCortf,data.C57BL6J.All.gammaBandPower.adjLH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.All.gammaBandPower.adjLH.meanCortf,data.Blank_SAP.All.gammaBandPower.adjLH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.All.gammaBandPower.adjLH.meanCortf,data.SSP_SAP.All.gammaBandPower.adjLH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[All] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during All data
ax12 = subplot(3,4,12);
loglog(data.C57BL6J.All.gammaBandPower.adjRH.meanCortf,data.C57BL6J.All.gammaBandPower.adjRH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.All.gammaBandPower.adjRH.meanCortf,data.Blank_SAP.All.gammaBandPower.adjRH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.All.gammaBandPower.adjRH.meanCortf,data.SSP_SAP.All.gammaBandPower.adjRH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[All] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% figure characteristics
linkaxes([ax1,ax2],'y')
linkaxes([ax3,ax4],'y')
linkaxes([ax5,ax6],'y')
linkaxes([ax7,ax8],'y')
linkaxes([ax9,ax10],'y')
linkaxes([ax11,ax12],'y')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure3,[dirpath 'PowerSpec_Gamma']);
    set(summaryFigure3,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'PowerSpec_Gamma'])
end
%% individual gamma-band power
summaryFigure4 = figure;
sgtitle('Gamma-band [30-100] Hz (envelope) cortical power spectra')
%% LH power spectra of gamma-band power during Rest
ax1 = subplot(3,4,1);
% C57BL6J
for aa = 1:size(data.C57BL6J.Rest.gammaBandPower.adjLH.normS,2)
    loglog(data.C57BL6J.Rest.gammaBandPower.adjLH.meanCortf,data.C57BL6J.Rest.gammaBandPower.adjLH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Rest.gammaBandPower.adjLH.normS,2)
    loglog(data.Blank_SAP.Rest.gammaBandPower.adjLH.meanCortf,data.Blank_SAP.Rest.gammaBandPower.adjLH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Rest.gammaBandPower.adjLH.normS,2)
    loglog(data.SSP_SAP.Rest.gammaBandPower.adjLH.meanCortf,data.SSP_SAP.Rest.gammaBandPower.adjLH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Rest] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/10,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Rest
ax2 = subplot(3,4,2);
% C57BL6J
for aa = 1:size(data.C57BL6J.Rest.gammaBandPower.adjRH.normS,2)
    loglog(data.C57BL6J.Rest.gammaBandPower.adjRH.meanCortf,data.C57BL6J.Rest.gammaBandPower.adjRH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Rest.gammaBandPower.adjRH.normS,2)
    loglog(data.Blank_SAP.Rest.gammaBandPower.adjRH.meanCortf,data.Blank_SAP.Rest.gammaBandPower.adjRH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Rest.gammaBandPower.adjRH.normS,2)
    loglog(data.SSP_SAP.Rest.gammaBandPower.adjRH.meanCortf,data.SSP_SAP.Rest.gammaBandPower.adjRH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Rest] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/10,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during NREM
ax3 = subplot(3,4,3);
% C57BL6J
for aa = 1:size(data.C57BL6J.NREM.gammaBandPower.adjLH.normS,2)
    loglog(data.C57BL6J.NREM.gammaBandPower.adjLH.meanCortf,data.C57BL6J.NREM.gammaBandPower.adjLH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.NREM.gammaBandPower.adjLH.normS,2)
    loglog(data.Blank_SAP.NREM.gammaBandPower.adjLH.meanCortf,data.Blank_SAP.NREM.gammaBandPower.adjLH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.NREM.gammaBandPower.adjLH.normS,2)
    loglog(data.SSP_SAP.NREM.gammaBandPower.adjLH.meanCortf,data.SSP_SAP.NREM.gammaBandPower.adjLH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[NREM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/30,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during NREM
ax4 = subplot(3,4,4);
% C57BL6J
for aa = 1:size(data.C57BL6J.NREM.gammaBandPower.adjRH.normS,2)
    loglog(data.C57BL6J.NREM.gammaBandPower.adjRH.meanCortf,data.C57BL6J.NREM.gammaBandPower.adjRH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.NREM.gammaBandPower.adjRH.normS,2)
    loglog(data.Blank_SAP.NREM.gammaBandPower.adjRH.meanCortf,data.Blank_SAP.NREM.gammaBandPower.adjRH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.NREM.gammaBandPower.adjRH.normS,2)
    loglog(data.SSP_SAP.NREM.gammaBandPower.adjRH.meanCortf,data.SSP_SAP.NREM.gammaBandPower.adjRH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[NREM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/30,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during REM
ax5 = subplot(3,4,5);
% C57BL6J
for aa = 1:size(data.C57BL6J.REM.gammaBandPower.adjLH.normS,2)
    loglog(data.C57BL6J.REM.gammaBandPower.adjLH.meanCortf,data.C57BL6J.REM.gammaBandPower.adjLH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.REM.gammaBandPower.adjLH.normS,2)
    loglog(data.Blank_SAP.REM.gammaBandPower.adjLH.meanCortf,data.Blank_SAP.REM.gammaBandPower.adjLH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.REM.gammaBandPower.adjLH.normS,2)
    loglog(data.SSP_SAP.REM.gammaBandPower.adjLH.meanCortf,data.SSP_SAP.REM.gammaBandPower.adjLH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[REM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/60,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during REM
ax6 = subplot(3,4,6);
% C57BL6J
for aa = 1:size(data.C57BL6J.REM.gammaBandPower.adjRH.normS,2)
    loglog(data.C57BL6J.REM.gammaBandPower.adjRH.meanCortf,data.C57BL6J.REM.gammaBandPower.adjRH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.REM.gammaBandPower.adjRH.normS,2)
    loglog(data.Blank_SAP.REM.gammaBandPower.adjRH.meanCortf,data.Blank_SAP.REM.gammaBandPower.adjRH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.REM.gammaBandPower.adjRH.normS,2)
    loglog(data.SSP_SAP.REM.gammaBandPower.adjRH.meanCortf,data.SSP_SAP.REM.gammaBandPower.adjRH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[REM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/60,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during Awake
ax7 = subplot(3,4,7);
% C57BL6J
for aa = 1:size(data.C57BL6J.Awake.gammaBandPower.adjLH.normS,2)
    loglog(data.C57BL6J.Awake.gammaBandPower.adjLH.meanCortf,data.C57BL6J.Awake.gammaBandPower.adjLH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Awake.gammaBandPower.adjLH.normS,2)
    loglog(data.Blank_SAP.Awake.gammaBandPower.adjLH.meanCortf,data.Blank_SAP.Awake.gammaBandPower.adjLH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Awake.gammaBandPower.adjLH.normS,2)
    loglog(data.SSP_SAP.Awake.gammaBandPower.adjLH.meanCortf,data.SSP_SAP.Awake.gammaBandPower.adjLH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Alert] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Awake
ax8 = subplot(3,4,8);
% C57BL6J
for aa = 1:size(data.C57BL6J.Awake.gammaBandPower.adjRH.normS,2)
    loglog(data.C57BL6J.Awake.gammaBandPower.adjRH.meanCortf,data.C57BL6J.Awake.gammaBandPower.adjRH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Awake.gammaBandPower.adjRH.normS,2)
    loglog(data.Blank_SAP.Awake.gammaBandPower.adjRH.meanCortf,data.Blank_SAP.Awake.gammaBandPower.adjRH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Awake.gammaBandPower.adjRH.normS,2)
    loglog(data.SSP_SAP.Awake.gammaBandPower.adjRH.meanCortf,data.SSP_SAP.Awake.gammaBandPower.adjRH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Alert] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during Sleep
ax9 = subplot(3,4,9);
% C57BL6J
for aa = 1:size(data.C57BL6J.Sleep.gammaBandPower.adjLH.normS,2)
    loglog(data.C57BL6J.Sleep.gammaBandPower.adjLH.meanCortf,data.C57BL6J.Sleep.gammaBandPower.adjLH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Sleep.gammaBandPower.adjLH.normS,2)
    loglog(data.Blank_SAP.Sleep.gammaBandPower.adjLH.meanCortf,data.Blank_SAP.Sleep.gammaBandPower.adjLH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Sleep.gammaBandPower.adjLH.normS,2)
    loglog(data.SSP_SAP.Sleep.gammaBandPower.adjLH.meanCortf,data.SSP_SAP.Sleep.gammaBandPower.adjLH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Asleep] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Sleep
ax10 = subplot(3,4,10);
% C57BL6J
for aa = 1:size(data.C57BL6J.Sleep.gammaBandPower.adjRH.normS,2)
    loglog(data.C57BL6J.Sleep.gammaBandPower.adjRH.meanCortf,data.C57BL6J.Sleep.gammaBandPower.adjRH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Sleep.gammaBandPower.adjRH.normS,2)
    loglog(data.Blank_SAP.Sleep.gammaBandPower.adjRH.meanCortf,data.Blank_SAP.Sleep.gammaBandPower.adjRH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Sleep.gammaBandPower.adjRH.normS,2)
    loglog(data.SSP_SAP.Sleep.gammaBandPower.adjRH.meanCortf,data.SSP_SAP.Sleep.gammaBandPower.adjRH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Asleep] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during All data
ax11 = subplot(3,4,11);
% C57BL6J
for aa = 1:size(data.C57BL6J.All.gammaBandPower.adjLH.normS,2)
    loglog(data.C57BL6J.All.gammaBandPower.adjLH.meanCortf,data.C57BL6J.All.gammaBandPower.adjLH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.All.gammaBandPower.adjLH.normS,2)
    loglog(data.Blank_SAP.All.gammaBandPower.adjLH.meanCortf,data.Blank_SAP.All.gammaBandPower.adjLH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.All.gammaBandPower.adjLH.normS,2)
    loglog(data.SSP_SAP.All.gammaBandPower.adjLH.meanCortf,data.SSP_SAP.All.gammaBandPower.adjLH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[All] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during All data
ax12 = subplot(3,4,12);
% C57BL6J
for aa = 1:size(data.C57BL6J.All.gammaBandPower.adjRH.normS,2)
    loglog(data.C57BL6J.All.gammaBandPower.adjRH.meanCortf,data.C57BL6J.All.gammaBandPower.adjRH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.All.gammaBandPower.adjRH.normS,2)
    loglog(data.Blank_SAP.All.gammaBandPower.adjRH.meanCortf,data.Blank_SAP.All.gammaBandPower.adjRH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.All.gammaBandPower.adjRH.normS,2)
    loglog(data.SSP_SAP.All.gammaBandPower.adjRH.meanCortf,data.SSP_SAP.All.gammaBandPower.adjRH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[All] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% figure characteristics
linkaxes([ax1,ax2],'y')
linkaxes([ax3,ax4],'y')
linkaxes([ax5,ax6],'y')
linkaxes([ax7,ax8],'y')
linkaxes([ax9,ax10],'y')
linkaxes([ax11,ax12],'y')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure4,[dirpath 'indPowerSpec_Gamma']);
    set(summaryFigure4,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'indPowerSpec_Gamma'])
end

end
