function [] = PowerSpectrumLFP_Bilateral_IOS(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% set-up and process data
resultsStruct = 'Results_PowerSpecLFP';
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
behavFields = {'Alert','Asleep','All'};
%% power spectra during different behaviors
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
            data.(treatment).(behavField).LH.S = [];
            data.(treatment).(behavField).LH.deltaS = [];
            data.(treatment).(behavField).LH.f = [];
            data.(treatment).(behavField).RH.S = [];
            data.(treatment).(behavField).RH.deltaS = [];
            data.(treatment).(behavField).RH.f = [];
            data.(treatment).(behavField).animalID = {};
            data.(treatment).(behavField).treatment = {};
        end
        data.(treatment).(behavField).LH.f = cat(1,data.(treatment).(behavField).LH.f,Results_PowerSpecLFP.(animalID).(behavField).LH.f);
        data.(treatment).(behavField).LH.S = cat(2,data.(treatment).(behavField).LH.S,Results_PowerSpecLFP.(animalID).(behavField).LH.S);
        data.(treatment).(behavField).RH.f = cat(1,data.(treatment).(behavField).RH.f,Results_PowerSpecLFP.(animalID).(behavField).RH.f);
        data.(treatment).(behavField).RH.S = cat(2,data.(treatment).(behavField).RH.S,Results_PowerSpecLFP.(animalID).(behavField).RH.S);
        index = find(round(Results_PowerSpecLFP.(animalID).(behavField).LH.f,2) == 4);
        if isempty(index) == false
            deltaIndex = index(end);
            data.(treatment).(behavField).LH.deltaS = cat(1,data.(treatment).(behavField).LH.deltaS,mean(Results_PowerSpecLFP.(animalID).(behavField).LH.S(1:deltaIndex)));
            data.(treatment).(behavField).RH.deltaS = cat(1,data.(treatment).(behavField).RH.deltaS,mean(Results_PowerSpecLFP.(animalID).(behavField).RH.S(1:deltaIndex)));
            data.(treatment).(behavField).animalID = cat(1,data.(treatment).(behavField).animalID,animalID);
            data.(treatment).(behavField).treatment = cat(1,data.(treatment).(behavField).treatment,treatment);
        end
    end
end
%% take mean/StD of S/f
for aa = 1:length(treatments)
    treatment = treatments{1,aa};
    for h = 1:length(behavFields)
        behavField = behavFields{1,h};
        data.(treatment).(behavField).LH.meanCortS = mean(data.(treatment).(behavField).LH.S,2);
        data.(treatment).(behavField).LH.stdCortS = std(data.(treatment).(behavField).LH.S,0,2);
        data.(treatment).(behavField).LH.meanCortf = mean(data.(treatment).(behavField).LH.f,1);
%         data.(treatment).(behavField).LH.meanDeltaS = mean(data.(treatment).(behavField).LH.deltaS,1);
%         data.(treatment).(behavField).LH.stdDeltaS = std(data.(treatment).(behavField).LH.deltaS,0,1);
        data.(treatment).(behavField).RH.meanCortS = mean(data.(treatment).(behavField).RH.S,2);
        data.(treatment).(behavField).RH.stdCortS = std(data.(treatment).(behavField).RH.S,0,2);
        data.(treatment).(behavField).RH.meanCortf = mean(data.(treatment).(behavField).RH.f,1);
%         data.(treatment).(behavField).RH.meanDeltaS = mean(data.(treatment).(behavField).RH.deltaS,1);
%         data.(treatment).(behavField).RH.stdDeltaS = std(data.(treatment).(behavField).RH.deltaS,0,1);
    end
end
%% statistics - generalized linear mixed effects model
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    Stats.(behavField).tableSize = cat(1,data.Blank_SAP.(behavField).RH.deltaS,data.SSP_SAP.(behavField).RH.deltaS,data.Naive.(behavField).RH.deltaS);
    Stats.(behavField).Table = table('Size',[size(Stats.(behavField).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Treatment','deltaS'});
    Stats.(behavField).Table.Mouse = cat(1,data.Blank_SAP.(behavField).animalID,data.SSP_SAP.(behavField).animalID,data.Naive.(behavField).animalID);
    Stats.(behavField).Table.Treatment = cat(1,data.Blank_SAP.(behavField).treatment,data.SSP_SAP.(behavField).treatment,data.Naive.(behavField).treatment);
    Stats.(behavField).Table.deltaS = cat(1,data.Blank_SAP.(behavField).RH.deltaS,data.SSP_SAP.(behavField).RH.deltaS,data.Naive.(behavField).RH.deltaS);
    Stats.(behavField).FitFormula = 'deltaS ~ 1 + Treatment + (1|Mouse)';
    Stats.(behavField).Stats = fitglme(Stats.(behavField).Table,Stats.(behavField).FitFormula);
end
%% average gamma-band power
summaryFigure1 = figure;
sgtitle('LFP Power Spectra [1-100 Hz]')
%% LH power spectra of gamma-band power during Alert
ax1 = subplot(3,2,1);
L1 = loglog(data.Naive.Alert.LH.meanCortf,data.Naive.Alert.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
% loglog(data.Naive.Alert.LH.meanCortf,data.Naive.Alert.LH.stdCortS + data.Naive.Alert.LH.meanCortS,'color',colors('sapphire'),'LineWidth',0.5);
% loglog(data.Naive.Alert.LH.meanCortf,data.Naive.Alert.LH.stdCortS - data.Naive.Alert.LH.meanCortS,'color',colors('sapphire'),'LineWidth',0.5);
hold on
L2 = loglog(data.Blank_SAP.Alert.LH.meanCortf,data.Blank_SAP.Alert.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
% loglog(data.Blank_SAP.Alert.LH.meanCortf,data.Blank_SAP.Alert.LH.stdCortS + data.Blank_SAP.Alert.LH.meanCortS,'color',colors('north texas green'),'LineWidth',0.5);
% loglog(data.Blank_SAP.Alert.LH.meanCortf,data.Blank_SAP.Alert.LH.stdCortS - data.Blank_SAP.Alert.LH.meanCortS,'color',colors('north texas green'),'LineWidth',0.5);
L3 = loglog(data.SSP_SAP.Alert.LH.meanCortf,data.SSP_SAP.Alert.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
% loglog(data.SSP_SAP.Alert.LH.meanCortf,data.SSP_SAP.Alert.LH.stdCortS + data.SSP_SAP.Alert.LH.meanCortS,'color',colors('electric purple'),'LineWidth',0.5);
% loglog(data.SSP_SAP.Alert.LH.meanCortf,data.SSP_SAP.Alert.LH.stdCortS - data.SSP_SAP.Alert.LH.meanCortS,'color',colors('electric purple'),'LineWidth',0.5);
title('[Alert] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1,100])
legend([L1,L2,L3],'Naive','SSP-SAP','Blank-SAP')
set(gca,'box','off')
%% RH power spectra of gamma-band power during Alert
ax2 = subplot(3,2,2);
loglog(data.Naive.Alert.RH.meanCortf,data.Naive.Alert.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
% loglog(data.Naive.Alert.RH.meanCortf,data.Naive.Alert.RH.stdCortS + data.Naive.Alert.RH.meanCortS,'color',colors('sapphire'),'LineWidth',0.5);
% loglog(data.Naive.Alert.RH.meanCortf,data.Naive.Alert.RH.stdCortS - data.Naive.Alert.RH.meanCortS,'color',colors('sapphire'),'LineWidth',0.5);
hold on
loglog(data.Blank_SAP.Alert.RH.meanCortf,data.Blank_SAP.Alert.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
% loglog(data.Blank_SAP.Alert.RH.meanCortf,data.Blank_SAP.Alert.RH.stdCortS + data.Blank_SAP.Alert.RH.meanCortS,'color',colors('north texas green'),'LineWidth',0.5);
% loglog(data.Blank_SAP.Alert.RH.meanCortf,data.Blank_SAP.Alert.RH.stdCortS - data.Blank_SAP.Alert.RH.meanCortS,'color',colors('north texas green'),'LineWidth',0.5);
loglog(data.SSP_SAP.Alert.RH.meanCortf,data.SSP_SAP.Alert.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
% loglog(data.SSP_SAP.Alert.RH.meanCortf,data.SSP_SAP.Alert.RH.stdCortS + data.SSP_SAP.Alert.RH.meanCortS,'color',colors('electric purple'),'LineWidth',0.5);
% loglog(data.SSP_SAP.Alert.RH.meanCortf,data.SSP_SAP.Alert.RH.stdCortS - data.SSP_SAP.Alert.RH.meanCortS,'color',colors('electric purple'),'LineWidth',0.5);
title('[Alert] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1,100])
set(gca,'box','off')
%% LH power spectra of gamma-band power during Asleep
ax3 = subplot(3,2,3);
loglog(data.Naive.Asleep.LH.meanCortf,data.Naive.Asleep.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
% loglog(data.Naive.Asleep.LH.meanCortf,data.Naive.Asleep.LH.stdCortS + data.Naive.Asleep.LH.meanCortS,'color',colors('sapphire'),'LineWidth',0.5);
% loglog(data.Naive.Asleep.LH.meanCortf,data.Naive.Asleep.LH.stdCortS - data.Naive.Asleep.LH.meanCortS,'color',colors('sapphire'),'LineWidth',0.5);
hold on
loglog(data.Blank_SAP.Asleep.LH.meanCortf,data.Blank_SAP.Asleep.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
% loglog(data.Blank_SAP.Asleep.LH.meanCortf,data.Blank_SAP.Asleep.LH.stdCortS + data.Blank_SAP.Asleep.LH.meanCortS,'color',colors('north texas green'),'LineWidth',0.5);
% loglog(data.Blank_SAP.Asleep.LH.meanCortf,data.Blank_SAP.Asleep.LH.stdCortS - data.Blank_SAP.Asleep.LH.meanCortS,'color',colors('north texas green'),'LineWidth',0.5);
loglog(data.SSP_SAP.Asleep.LH.meanCortf,data.SSP_SAP.Asleep.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
% loglog(data.SSP_SAP.Asleep.LH.meanCortf,data.SSP_SAP.Asleep.LH.stdCortS + data.SSP_SAP.Asleep.LH.meanCortS,'color',colors('electric purple'),'LineWidth',0.5);
% loglog(data.SSP_SAP.Asleep.LH.meanCortf,data.SSP_SAP.Asleep.LH.stdCortS - data.SSP_SAP.Asleep.LH.meanCortS,'color',colors('electric purple'),'LineWidth',0.5);
title('[Asleep] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1,100])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Asleep
ax4 = subplot(3,2,4);
loglog(data.Naive.Asleep.RH.meanCortf,data.Naive.Asleep.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
% loglog(data.Naive.Asleep.RH.meanCortf,data.Naive.Asleep.RH.stdCortS + data.Naive.Asleep.RH.meanCortS,'color',colors('sapphire'),'LineWidth',0.5);
% loglog(data.Naive.Asleep.RH.meanCortf,data.Naive.Asleep.RH.stdCortS - data.Naive.Asleep.RH.meanCortS,'color',colors('sapphire'),'LineWidth',0.5);
hold on
loglog(data.Blank_SAP.Asleep.RH.meanCortf,data.Blank_SAP.Asleep.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
% loglog(data.Blank_SAP.Asleep.RH.meanCortf,data.Blank_SAP.Asleep.RH.stdCortS + data.Blank_SAP.Asleep.RH.meanCortS,'color',colors('north texas green'),'LineWidth',0.5);
% loglog(data.Blank_SAP.Asleep.RH.meanCortf,data.Blank_SAP.Asleep.RH.stdCortS - data.Blank_SAP.Asleep.RH.meanCortS,'color',colors('north texas green'),'LineWidth',0.5);
loglog(data.SSP_SAP.Asleep.RH.meanCortf,data.SSP_SAP.Asleep.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
% loglog(data.SSP_SAP.Asleep.RH.meanCortf,data.SSP_SAP.Asleep.RH.stdCortS + data.SSP_SAP.Asleep.RH.meanCortS,'color',colors('electric purple'),'LineWidth',0.5);
% loglog(data.SSP_SAP.Asleep.RH.meanCortf,data.SSP_SAP.Asleep.RH.stdCortS - data.SSP_SAP.Asleep.RH.meanCortS,'color',colors('electric purple'),'LineWidth',0.5);
title('[Asleep] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1,100])
set(gca,'box','off')
%% LH power spectra of gamma-band power during All data
ax5 = subplot(3,2,5);
loglog(data.Naive.All.LH.meanCortf,data.Naive.All.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
% loglog(data.Naive.All.LH.meanCortf,data.Naive.All.LH.stdCortS + data.Naive.All.LH.meanCortS,'color',colors('sapphire'),'LineWidth',0.5);
% loglog(data.Naive.All.LH.meanCortf,data.Naive.All.LH.stdCortS - data.Naive.All.LH.meanCortS,'color',colors('sapphire'),'LineWidth',0.5);
hold on
loglog(data.Blank_SAP.All.LH.meanCortf,data.Blank_SAP.All.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
% loglog(data.Blank_SAP.All.LH.meanCortf,data.Blank_SAP.All.LH.stdCortS + data.Blank_SAP.All.LH.meanCortS,'color',colors('north texas green'),'LineWidth',0.5);
% loglog(data.Blank_SAP.All.LH.meanCortf,data.Blank_SAP.All.LH.stdCortS - data.Blank_SAP.All.LH.meanCortS,'color',colors('north texas green'),'LineWidth',0.5);
loglog(data.SSP_SAP.All.LH.meanCortf,data.SSP_SAP.All.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
% loglog(data.SSP_SAP.All.LH.meanCortf,data.SSP_SAP.All.LH.stdCortS + data.SSP_SAP.All.LH.meanCortS,'color',colors('electric purple'),'LineWidth',0.5);
% loglog(data.SSP_SAP.All.LH.meanCortf,data.SSP_SAP.All.LH.stdCortS - data.SSP_SAP.All.LH.meanCortS,'color',colors('electric purple'),'LineWidth',0.5);
title('[All] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1,100])
set(gca,'box','off')
%% RH power spectra of gamma-band power during All data
ax6 = subplot(3,2,6);
loglog(data.Naive.All.RH.meanCortf,data.Naive.All.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
% loglog(data.Naive.All.RH.meanCortf,data.Naive.All.RH.stdCortS + data.Naive.All.RH.meanCortS,'color',colors('sapphire'),'LineWidth',0.5);
% loglog(data.Naive.All.RH.meanCortf,data.Naive.All.RH.stdCortS - data.Naive.All.RH.meanCortS,'color',colors('sapphire'),'LineWidth',0.5);
hold on
loglog(data.Blank_SAP.All.RH.meanCortf,data.Blank_SAP.All.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
% loglog(data.Blank_SAP.All.RH.meanCortf,data.Blank_SAP.All.RH.stdCortS + data.Blank_SAP.All.RH.meanCortS,'color',colors('north texas green'),'LineWidth',0.5);
% loglog(data.Blank_SAP.All.RH.meanCortf,data.Blank_SAP.All.RH.stdCortS - data.Blank_SAP.All.RH.meanCortS,'color',colors('north texas green'),'LineWidth',0.5);
loglog(data.SSP_SAP.All.RH.meanCortf,data.SSP_SAP.All.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
% loglog(data.SSP_SAP.All.RH.meanCortf,data.SSP_SAP.All.RH.stdCortS + data.SSP_SAP.All.RH.meanCortS,'color',colors('electric purple'),'LineWidth',0.5);
% loglog(data.SSP_SAP.All.RH.meanCortf,data.SSP_SAP.All.RH.stdCortS - data.SSP_SAP.All.RH.meanCortS,'color',colors('electric purple'),'LineWidth',0.5);
title('[All] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1,100])
set(gca,'box','off')
%% figure characteristics
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'xy')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'LFP Power Spectrum - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'PowerSpec_LFP']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'PowerSpec_LFP'])
end
%% power spec stats
summaryFigure2 = figure;
ax1 = subplot(1,3,1);
s1 = scatter(ones(1,length(data.Naive.Alert.RH.deltaS))*1,data.Naive.Alert.RH.deltaS,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(data.Naive.Alert.RH.deltaS),std(data.Naive.Alert.RH.deltaS),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(data.SSP_SAP.Alert.RH.deltaS))*2,data.SSP_SAP.Alert.RH.deltaS,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,mean(data.SSP_SAP.Alert.RH.deltaS),std(data.SSP_SAP.Alert.RH.deltaS),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(data.Blank_SAP.Alert.RH.deltaS))*3,data.Blank_SAP.Alert.RH.deltaS,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,mean(data.Blank_SAP.Alert.RH.deltaS),std(data.Blank_SAP.Alert.RH.deltaS),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
% stat lines
% plot([1,3],[1,1],'k');
% text(2,1,'ns','FontSize',16)
% plot([2,3],[1,1],'k');
% text(2.5,1,'***','FontSize',16)
ylabel('Power (a.u.)')
legend([s1,s2,s3],'Naive','SSP-SAP','Blank-SAP')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,4])
set(gca,'box','off')
ax2 = subplot(1,3,2);
s1 = scatter(ones(1,length(data.Naive.Asleep.RH.deltaS))*1,data.Naive.Asleep.RH.deltaS,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(data.Naive.Asleep.RH.deltaS),std(data.Naive.Asleep.RH.deltaS),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(data.SSP_SAP.Asleep.RH.deltaS))*2,data.SSP_SAP.Asleep.RH.deltaS,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,mean(data.SSP_SAP.Asleep.RH.deltaS),std(data.SSP_SAP.Asleep.RH.deltaS),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(data.Blank_SAP.Asleep.RH.deltaS))*3,data.Blank_SAP.Asleep.RH.deltaS,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,mean(data.Blank_SAP.Asleep.RH.deltaS),std(data.Blank_SAP.Asleep.RH.deltaS),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
% stat lines
% plot([1,3],[1,1],'k');
% text(2,1,'ns','FontSize',16)
% plot([2,3],[1,1],'k');
% text(2.5,1,'***','FontSize',16)
ylabel('Power (a.u.)')
legend([s1,s2,s3],'Naive','SSP-SAP','Blank-SAP')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,4])
set(gca,'box','off')
ax3 = subplot(1,3,3);
s1 = scatter(ones(1,length(data.Naive.All.RH.deltaS))*1,data.Naive.All.RH.deltaS,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(data.Naive.All.RH.deltaS),std(data.Naive.All.RH.deltaS),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(data.SSP_SAP.All.RH.deltaS))*2,data.SSP_SAP.All.RH.deltaS,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,mean(data.SSP_SAP.All.RH.deltaS),std(data.SSP_SAP.All.RH.deltaS),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(data.Blank_SAP.All.RH.deltaS))*3,data.Blank_SAP.All.RH.deltaS,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e3 = errorbar(3,mean(data.Blank_SAP.All.RH.deltaS),std(data.Blank_SAP.All.RH.deltaS),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
% stat lines
% plot([1,3],[1,1],'k');
% text(2,1,'ns','FontSize',16)
% plot([2,3],[1,1],'k');
% text(2.5,1,'***','FontSize',16)
ylabel('Power (a.u.)')
legend([s1,s2,s3],'Naive','SSP-SAP','Blank-SAP')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,4])
set(gca,'box','off')
linkaxes([ax1,ax2,ax3],'xy')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'LFP Power Spectrum - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'LFP_Statistics']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'LFP_Statistics'])
    %% statistical diary
    diaryFile = [dirpath 'LFP_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % Awake stats
    disp('======================================================================================================================')
    disp('GLME statistics for Alert')
    disp('======================================================================================================================')
    disp(Stats.Alert.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('======================================================================================================================')
    disp('GLME statistics for Asleep')
    disp('======================================================================================================================')
    disp(Stats.Asleep.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('======================================================================================================================')
    disp('GLME statistics for All')
    disp('======================================================================================================================')
    disp(Stats.All.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
end

end