function [] = PowerSpectrum_Bilateral_IOS(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% set-up and process data
resultsStruct = 'Results_PowerSpec';
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
% behavFields = {'Rest','NREM','REM','Awake','Asleep','All'};
behavFields = {'Awake','Asleep','All'};
dataTypes = {'CBV_HbT','gammaBandPower'};
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
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % pre-allocate necessary variable fields
            data.(treatment).(behavField).dummCheck = 1;
            if isfield(data.(treatment).(behavField),dataType) == false
                data.(treatment).(behavField).(dataType).LH.S = [];
                data.(treatment).(behavField).(dataType).LH.f = [];
                data.(treatment).(behavField).(dataType).RH.S = [];
                data.(treatment).(behavField).(dataType).RH.f = [];
                data.(treatment).(behavField).(dataType).animalID = {};
                data.(treatment).(behavField).(dataType).treatment = {};
                data.(treatment).(behavField).(dataType).freqS001 = {};
                data.(treatment).(behavField).(dataType).freqS01 = {};
                data.(treatment).(behavField).(dataType).freqS05 = {};
            end
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(Results_PowerSpec.(animalID).(behavField).(dataType).LH.S) == false
                data.(treatment).(behavField).(dataType).LH.S = cat(2,data.(treatment).(behavField).(dataType).LH.S,Results_PowerSpec.(animalID).(behavField).(dataType).adjLH.S);
                data.(treatment).(behavField).(dataType).LH.f = cat(1,data.(treatment).(behavField).(dataType).LH.f,Results_PowerSpec.(animalID).(behavField).(dataType).adjLH.f);
                data.(treatment).(behavField).(dataType).RH.S = cat(2,data.(treatment).(behavField).(dataType).RH.S,Results_PowerSpec.(animalID).(behavField).(dataType).adjRH.S);
                data.(treatment).(behavField).(dataType).RH.f = cat(1,data.(treatment).(behavField).(dataType).RH.f,Results_PowerSpec.(animalID).(behavField).(dataType).adjRH.f);
                data.(treatment).(behavField).(dataType).animalID = cat(1,data.(treatment).(behavField).(dataType).animalID,animalID);
                data.(treatment).(behavField).(dataType).treatment = cat(1,data.(treatment).(behavField).(dataType).treatment,treatment);
                data.(treatment).(behavField).(dataType).freqS001 = cat(1,data.(treatment).(behavField).(dataType).freqS001,'S001');
                data.(treatment).(behavField).(dataType).freqS01 = cat(1,data.(treatment).(behavField).(dataType).freqS01,'S01');
                data.(treatment).(behavField).(dataType).freqS05 = cat(1,data.(treatment).(behavField).(dataType).freqS05,'S05');
            end
        end
    end
end
%% find the peak of the resting PSD for each animal/hemisphere
for aa = 1:length(treatments)
    treatment = treatments{1,aa};
    for cc = 1:length(dataTypes)
        dataType = dataTypes{1,cc};
        for ee = 1:size(data.(treatment).Rest.(dataType).LH.S,2)
            data.(treatment).baseline.(dataType).LH(ee,1) = max(data.(treatment).Rest.(dataType).LH.S(:,ee));
            data.(treatment).baseline.(dataType).RH(ee,1) = max(data.(treatment).Rest.(dataType).RH.S(:,ee));
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
            for ee = 1:size(data.(treatment).(behavField).(dataType).LH.S,2)
                data.(treatment).(behavField).(dataType).LH.normS(:,ee) = (data.(treatment).(behavField).(dataType).LH.S(:,ee))*(1/(data.(treatment).baseline.(dataType).LH(ee,1)));
                data.(treatment).(behavField).(dataType).RH.normS(:,ee) = (data.(treatment).(behavField).(dataType).RH.S(:,ee))*(1/(data.(treatment).baseline.(dataType).RH(ee,1)));
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
            data.(treatment).(behavField).(dataType).LH.meanCortS = mean(data.(treatment).(behavField).(dataType).LH.normS,2);
            data.(treatment).(behavField).(dataType).LH.stdCortS = std(data.(treatment).(behavField).(dataType).LH.normS,0,2);
            data.(treatment).(behavField).(dataType).LH.meanCortf = mean(data.(treatment).(behavField).(dataType).LH.f,1);
            data.(treatment).(behavField).(dataType).RH.meanCortS = mean(data.(treatment).(behavField).(dataType).RH.normS,2);
            data.(treatment).(behavField).(dataType).RH.stdCortS = std(data.(treatment).(behavField).(dataType).RH.normS,0,2);
            data.(treatment).(behavField).(dataType).RH.meanCortf = mean(data.(treatment).(behavField).(dataType).RH.f,1);
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
            data.(treatment).(behavField).(dataType).LH.meanC = mean(data.(treatment).(behavField).(dataType).LH.C,2);
            data.(treatment).(behavField).(dataType).LH.stdC = std(data.(treatment).(behavField).(dataType).LH.C,0,2);
            data.(treatment).(behavField).(dataType).LH.meanf = mean(data.(treatment).(behavField).(dataType).LH.f,1);
            data.(treatment).(behavField).(dataType).LH.maxConfC = geomean(data.(treatment).(behavField).(dataType).LH.confC);
            data.(treatment).(behavField).(dataType).LH.maxConfC_Y = ones(length(data.(treatment).(behavField).(dataType).LH.meanf),1)*data.(treatment).(behavField).(dataType).LH.maxConfC;
            data.(treatment).(behavField).(dataType).RH.meanC = mean(data.(treatment).(behavField).(dataType).RH.C,2);
            data.(treatment).(behavField).(dataType).RH.stdC = std(data.(treatment).(behavField).(dataType).RH.C,0,2);
            data.(treatment).(behavField).(dataType).RH.meanf = mean(data.(treatment).(behavField).(dataType).RH.f,1);
            data.(treatment).(behavField).(dataType).RH.maxConfC = geomean(data.(treatment).(behavField).(dataType).RH.confC);
            data.(treatment).(behavField).(dataType).RH.maxConfC_Y = ones(length(data.(treatment).(behavField).(dataType).RH.meanf),1)*data.(treatment).(behavField).(dataType).RH.maxConfC;
        end
    end
end
%% find Hz peaks in power
for qq = 1:length(treatments)
    treatment = treatments{1,qq};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        for ff = 1:length(dataTypes)
            dataType = dataTypes{1,ff};
            for gg = 1:size(data.(treatment).(behavField).(dataType).RH.normS,2)
                F = round(data.(treatment).(behavField).(dataType).RH.f(gg,:),3);
                RH_S = data.(treatment).(behavField).(dataType).RH.normS(:,gg);
                index001 = find(F == 0.01);
                index01 = find(F == 0.1);
                index05 = find(F == 0.5);
                data.(treatment).(behavField).(dataType).RH.S001(gg,1) = mean(RH_S(1:index001(1)));
                data.(treatment).(behavField).(dataType).RH.S01(gg,1) = mean(RH_S(index001(1) + 1:index01(1)));
                data.(treatment).(behavField).(dataType).RH.S05(gg,1) = mean(RH_S(index01(1) + 1:index05(1)));
            end
        end
    end
end
%% statistics - generalized linear mixed effects model
freqs = {'S001','S01','S05'};
for aa = 1:length(freqs)
    freq = freqs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % statistics - generalized linear mixed effects model
            Stats.(dataType).(behavField).(freq).tableSize = cat(1,data.Blank_SAP.(behavField).(dataType).RH.(freq),data.SSP_SAP.(behavField).(dataType).RH.(freq));
            Stats.(dataType).(behavField).(freq).Table = table('Size',[size(Stats.(dataType).(behavField).(freq).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Treatment','Power'});
            Stats.(dataType).(behavField).(freq).Table.Mouse = cat(1,data.Blank_SAP.(behavField).(dataType).animalID,data.SSP_SAP.(behavField).(dataType).animalID);
            Stats.(dataType).(behavField).(freq).Table.Treatment = cat(1,data.Blank_SAP.(behavField).(dataType).treatment,data.SSP_SAP.(behavField).(dataType).treatment);
            Stats.(dataType).(behavField).(freq).Table.Power = cat(1,data.Blank_SAP.(behavField).(dataType).RH.(freq),data.SSP_SAP.(behavField).(dataType).RH.(freq));
            Stats.(dataType).(behavField).(freq).FitFormula = 'Power ~ 1 + Treatment + (1|Mouse)';
            Stats.(dataType).(behavField).(freq).Stats = fitglme(Stats.(dataType).(behavField).(freq).Table,Stats.(dataType).(behavField).(freq).FitFormula);
        end
    end
end
%% stats
% for bb = 1:length(behavFields)
%     behavField = behavFields{1,bb};
%     for cc = 1:length(dataTypes)
%         dataType = dataTypes{1,cc};
%         % statistics - generalized linear mixed effects model
%         Stats.(dataType).(behavField).tableSize = cat(1,data.Blank_SAP.(behavField).(dataType).RH.S001,data.Blank_SAP.(behavField).(dataType).RH.S01,data.Blank_SAP.(behavField).(dataType).RH.S05,...
%             data.SSP_SAP.(behavField).(dataType).RH.S001,data.SSP_SAP.(behavField).(dataType).RH.S01,data.SSP_SAP.(behavField).(dataType).RH.S05);
%         Stats.(dataType).(behavField).Table = table('Size',[size(Stats.(dataType).(behavField).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Treatment','Frequency','Power'});
%         Stats.(dataType).(behavField).Table.Mouse = cat(1,data.Blank_SAP.(behavField).(dataType).animalID,data.Blank_SAP.(behavField).(dataType).animalID,data.Blank_SAP.(behavField).(dataType).animalID,...
%             data.SSP_SAP.(behavField).(dataType).animalID,data.SSP_SAP.(behavField).(dataType).animalID,data.SSP_SAP.(behavField).(dataType).animalID);
%         Stats.(dataType).(behavField).Table.Treatment = cat(1,data.Blank_SAP.(behavField).(dataType).treatment,data.Blank_SAP.(behavField).(dataType).treatment,data.Blank_SAP.(behavField).(dataType).treatment,...,
%             data.SSP_SAP.(behavField).(dataType).treatment,data.SSP_SAP.(behavField).(dataType).treatment,data.SSP_SAP.(behavField).(dataType).treatment);
%         Stats.(dataType).(behavField).Table.Frequency = cat(1,data.Blank_SAP.(behavField).(dataType).freqS001,data.Blank_SAP.(behavField).(dataType).freqS01,data.Blank_SAP.(behavField).(dataType).freqS05,...
%             data.SSP_SAP.(behavField).(dataType).freqS001,data.SSP_SAP.(behavField).(dataType).freqS01,data.SSP_SAP.(behavField).(dataType).freqS05);
%         Stats.(dataType).(behavField).Table.Power = cat(1,data.Blank_SAP.(behavField).(dataType).RH.S001,data.Blank_SAP.(behavField).(dataType).RH.S01,data.Blank_SAP.(behavField).(dataType).RH.S05,...
%             data.SSP_SAP.(behavField).(dataType).RH.S001,data.SSP_SAP.(behavField).(dataType).RH.S01,data.SSP_SAP.(behavField).(dataType).RH.S05);
%         Stats.(dataType).(behavField).FitFormula = 'Power ~ 1 + Treatment + Frequency + Frequency*Treatment + (1|Mouse)';
%         Stats.(dataType).(behavField).Stats = fitglme(Stats.(dataType).(behavField).Table,Stats.(dataType).(behavField).FitFormula);
%     end
% end
%% average HbT power
summaryFigure1 = figure;
sgtitle('\DeltaHbT (\muM) cortical power spectra')
%% LH power spectra of HbT power during Rest
ax1 = subplot(3,4,1);
p1 = loglog(data.Naive.Rest.CBV_HbT.LH.meanCortf,data.Naive.Rest.CBV_HbT.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
p2 = loglog(data.Blank_SAP.Rest.CBV_HbT.LH.meanCortf,data.Blank_SAP.Rest.CBV_HbT.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
p3 = loglog(data.SSP_SAP.Rest.CBV_HbT.LH.meanCortf,data.SSP_SAP.Rest.CBV_HbT.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Rest] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
xlim([1/10,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during Rest
ax2 = subplot(3,4,2);
loglog(data.Naive.Rest.CBV_HbT.RH.meanCortf,data.Naive.Rest.CBV_HbT.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Rest.CBV_HbT.RH.meanCortf,data.Blank_SAP.Rest.CBV_HbT.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Rest.CBV_HbT.RH.meanCortf,data.SSP_SAP.Rest.CBV_HbT.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Rest] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/10,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during NREM
ax3 = subplot(3,4,3);
loglog(data.Naive.NREM.CBV_HbT.LH.meanCortf,data.Naive.NREM.CBV_HbT.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.NREM.CBV_HbT.LH.meanCortf,data.Blank_SAP.NREM.CBV_HbT.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.NREM.CBV_HbT.LH.meanCortf,data.SSP_SAP.NREM.CBV_HbT.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[NREM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
xlim([1/30,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during NREM
ax4 = subplot(3,4,4);
loglog(data.Naive.NREM.CBV_HbT.RH.meanCortf,data.Naive.NREM.CBV_HbT.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.NREM.CBV_HbT.RH.meanCortf,data.Blank_SAP.NREM.CBV_HbT.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.NREM.CBV_HbT.RH.meanCortf,data.SSP_SAP.NREM.CBV_HbT.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[NREM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/30,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during REM
ax5 = subplot(3,4,5);
loglog(data.Naive.REM.CBV_HbT.LH.meanCortf,data.Naive.REM.CBV_HbT.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.REM.CBV_HbT.LH.meanCortf,data.Blank_SAP.REM.CBV_HbT.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.REM.CBV_HbT.LH.meanCortf,data.SSP_SAP.REM.CBV_HbT.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[REM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
xlim([1/60,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during REM
ax6 = subplot(3,4,6);
loglog(data.Naive.REM.CBV_HbT.RH.meanCortf,data.Naive.REM.CBV_HbT.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.REM.CBV_HbT.RH.meanCortf,data.Blank_SAP.REM.CBV_HbT.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.REM.CBV_HbT.RH.meanCortf,data.SSP_SAP.REM.CBV_HbT.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[REM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/60,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during Awake
ax7 = subplot(3,4,7);
loglog(data.Naive.Awake.CBV_HbT.LH.meanCortf,data.Naive.Awake.CBV_HbT.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Awake.CBV_HbT.LH.meanCortf,data.Blank_SAP.Awake.CBV_HbT.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Awake.CBV_HbT.LH.meanCortf,data.SSP_SAP.Awake.CBV_HbT.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Alert] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during Awake
ax8 = subplot(3,4,8);
loglog(data.Naive.Awake.CBV_HbT.RH.meanCortf,data.Naive.Awake.CBV_HbT.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Awake.CBV_HbT.RH.meanCortf,data.Blank_SAP.Awake.CBV_HbT.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Awake.CBV_HbT.RH.meanCortf,data.SSP_SAP.Awake.CBV_HbT.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Alert] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during Asleep
ax9 = subplot(3,4,9);
loglog(data.Naive.Asleep.CBV_HbT.LH.meanCortf,data.Naive.Asleep.CBV_HbT.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Asleep.CBV_HbT.LH.meanCortf,data.Blank_SAP.Asleep.CBV_HbT.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Asleep.CBV_HbT.LH.meanCortf,data.SSP_SAP.Asleep.CBV_HbT.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[AAsleep] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during Asleep
ax10 = subplot(3,4,10);
loglog(data.Naive.Asleep.CBV_HbT.RH.meanCortf,data.Naive.Asleep.CBV_HbT.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Asleep.CBV_HbT.RH.meanCortf,data.Blank_SAP.Asleep.CBV_HbT.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Asleep.CBV_HbT.RH.meanCortf,data.SSP_SAP.Asleep.CBV_HbT.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[AAsleep] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during All data
ax11 = subplot(3,4,11);
loglog(data.Naive.All.CBV_HbT.LH.meanCortf,data.Naive.All.CBV_HbT.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.All.CBV_HbT.LH.meanCortf,data.Blank_SAP.All.CBV_HbT.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.All.CBV_HbT.LH.meanCortf,data.SSP_SAP.All.CBV_HbT.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[All] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during All data
ax12 = subplot(3,4,12);
loglog(data.Naive.All.CBV_HbT.RH.meanCortf,data.Naive.All.CBV_HbT.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.All.CBV_HbT.RH.meanCortf,data.Blank_SAP.All.CBV_HbT.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.All.CBV_HbT.RH.meanCortf,data.SSP_SAP.All.CBV_HbT.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
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
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Power Spectrum - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'AveragePowerSpec_HbT']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'AveragePowerSpec_HbT'])
end
%% individual HbT power
summaryFigure2 = figure;
sgtitle('\DeltaHbT (\muM) cortical power spectra - individual animals')
%% LH power spectra of HbT power during Rest
ax1 = subplot(3,4,1);
% Naive
for aa = 1:size(data.Naive.Rest.CBV_HbT.LH.normS,2)
    loglog(data.Naive.Rest.CBV_HbT.LH.meanCortf,data.Naive.Rest.CBV_HbT.LH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Rest.CBV_HbT.LH.normS,2)
    loglog(data.Blank_SAP.Rest.CBV_HbT.LH.meanCortf,data.Blank_SAP.Rest.CBV_HbT.LH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Rest.CBV_HbT.LH.normS,2)
    loglog(data.SSP_SAP.Rest.CBV_HbT.LH.meanCortf,data.SSP_SAP.Rest.CBV_HbT.LH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Rest] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/10,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during Rest
ax2 = subplot(3,4,2);
% Naive
for aa = 1:size(data.Naive.Rest.CBV_HbT.RH.normS,2)
    loglog(data.Naive.Rest.CBV_HbT.RH.meanCortf,data.Naive.Rest.CBV_HbT.RH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Rest.CBV_HbT.RH.normS,2)
    loglog(data.Blank_SAP.Rest.CBV_HbT.RH.meanCortf,data.Blank_SAP.Rest.CBV_HbT.RH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Rest.CBV_HbT.RH.normS,2)
    loglog(data.SSP_SAP.Rest.CBV_HbT.RH.meanCortf,data.SSP_SAP.Rest.CBV_HbT.RH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Rest] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/10,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during NREM
ax3 = subplot(3,4,3);
% Naive
for aa = 1:size(data.Naive.NREM.CBV_HbT.LH.normS,2)
    loglog(data.Naive.NREM.CBV_HbT.LH.meanCortf,data.Naive.NREM.CBV_HbT.LH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.NREM.CBV_HbT.LH.normS,2)
    loglog(data.Blank_SAP.NREM.CBV_HbT.LH.meanCortf,data.Blank_SAP.NREM.CBV_HbT.LH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.NREM.CBV_HbT.LH.normS,2)
    loglog(data.SSP_SAP.NREM.CBV_HbT.LH.meanCortf,data.SSP_SAP.NREM.CBV_HbT.LH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[NREM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/30,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during NREM
ax4 = subplot(3,4,4);
% Naive
for aa = 1:size(data.Naive.NREM.CBV_HbT.RH.normS,2)
    loglog(data.Naive.NREM.CBV_HbT.RH.meanCortf,data.Naive.NREM.CBV_HbT.RH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.NREM.CBV_HbT.RH.normS,2)
    loglog(data.Blank_SAP.NREM.CBV_HbT.RH.meanCortf,data.Blank_SAP.NREM.CBV_HbT.RH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.NREM.CBV_HbT.RH.normS,2)
    loglog(data.SSP_SAP.NREM.CBV_HbT.RH.meanCortf,data.SSP_SAP.NREM.CBV_HbT.RH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[NREM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/30,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during REM
ax5 = subplot(3,4,5);
% Naive
for aa = 1:size(data.Naive.REM.CBV_HbT.LH.normS,2)
    loglog(data.Naive.REM.CBV_HbT.LH.meanCortf,data.Naive.REM.CBV_HbT.LH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.REM.CBV_HbT.LH.normS,2)
    loglog(data.Blank_SAP.REM.CBV_HbT.LH.meanCortf,data.Blank_SAP.REM.CBV_HbT.LH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.REM.CBV_HbT.LH.normS,2)
    loglog(data.SSP_SAP.REM.CBV_HbT.LH.meanCortf,data.SSP_SAP.REM.CBV_HbT.LH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[REM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/60,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during REM
ax6 = subplot(3,4,6);
% Naive
for aa = 1:size(data.Naive.REM.CBV_HbT.RH.normS,2)
    loglog(data.Naive.REM.CBV_HbT.RH.meanCortf,data.Naive.REM.CBV_HbT.RH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.REM.CBV_HbT.RH.normS,2)
    loglog(data.Blank_SAP.REM.CBV_HbT.RH.meanCortf,data.Blank_SAP.REM.CBV_HbT.RH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.REM.CBV_HbT.RH.normS,2)
    loglog(data.SSP_SAP.REM.CBV_HbT.RH.meanCortf,data.SSP_SAP.REM.CBV_HbT.RH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[REM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/60,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during Awake
ax7 = subplot(3,4,7);
% Naive
for aa = 1:size(data.Naive.Awake.CBV_HbT.LH.normS,2)
    loglog(data.Naive.Awake.CBV_HbT.LH.meanCortf,data.Naive.Awake.CBV_HbT.LH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Awake.CBV_HbT.LH.normS,2)
    loglog(data.Blank_SAP.Awake.CBV_HbT.LH.meanCortf,data.Blank_SAP.Awake.CBV_HbT.LH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Awake.CBV_HbT.LH.normS,2)
    loglog(data.SSP_SAP.Awake.CBV_HbT.LH.meanCortf,data.SSP_SAP.Awake.CBV_HbT.LH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Alert] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during Awake
ax8 = subplot(3,4,8);
% Naive
for aa = 1:size(data.Naive.Awake.CBV_HbT.RH.normS,2)
    loglog(data.Naive.Awake.CBV_HbT.RH.meanCortf,data.Naive.Awake.CBV_HbT.RH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Awake.CBV_HbT.RH.normS,2)
    loglog(data.Blank_SAP.Awake.CBV_HbT.RH.meanCortf,data.Blank_SAP.Awake.CBV_HbT.RH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Awake.CBV_HbT.RH.normS,2)
    loglog(data.SSP_SAP.Awake.CBV_HbT.RH.meanCortf,data.SSP_SAP.Awake.CBV_HbT.RH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Alert] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during Asleep
ax9 = subplot(3,4,9);
% Naive
for aa = 1:size(data.Naive.Asleep.CBV_HbT.LH.normS,2)
    loglog(data.Naive.Asleep.CBV_HbT.LH.meanCortf,data.Naive.Asleep.CBV_HbT.LH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Asleep.CBV_HbT.LH.normS,2)
    loglog(data.Blank_SAP.Asleep.CBV_HbT.LH.meanCortf,data.Blank_SAP.Asleep.CBV_HbT.LH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Asleep.CBV_HbT.LH.normS,2)
    loglog(data.SSP_SAP.Asleep.CBV_HbT.LH.meanCortf,data.SSP_SAP.Asleep.CBV_HbT.LH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[AAsleep] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during Asleep
ax10 = subplot(3,4,10);
% Naive
for aa = 1:size(data.Naive.Asleep.CBV_HbT.RH.normS,2)
    loglog(data.Naive.Asleep.CBV_HbT.RH.meanCortf,data.Naive.Asleep.CBV_HbT.RH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Asleep.CBV_HbT.RH.normS,2)
    loglog(data.Blank_SAP.Asleep.CBV_HbT.RH.meanCortf,data.Blank_SAP.Asleep.CBV_HbT.RH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Asleep.CBV_HbT.RH.normS,2)
    loglog(data.SSP_SAP.Asleep.CBV_HbT.RH.meanCortf,data.SSP_SAP.Asleep.CBV_HbT.RH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[AAsleep] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of HbT power during All data
ax11 = subplot(3,4,11);
% Naive
for aa = 1:size(data.Naive.All.CBV_HbT.LH.normS,2)
    loglog(data.Naive.All.CBV_HbT.LH.meanCortf,data.Naive.All.CBV_HbT.LH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.All.CBV_HbT.LH.normS,2)
    loglog(data.Blank_SAP.All.CBV_HbT.LH.meanCortf,data.Blank_SAP.All.CBV_HbT.LH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.All.CBV_HbT.LH.normS,2)
    loglog(data.SSP_SAP.All.CBV_HbT.LH.meanCortf,data.SSP_SAP.All.CBV_HbT.LH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[All] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of HbT power during All data
ax12 = subplot(3,4,12);
% Naive
for aa = 1:size(data.Naive.All.CBV_HbT.RH.normS,2)
    loglog(data.Naive.All.CBV_HbT.RH.meanCortf,data.Naive.All.CBV_HbT.RH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.All.CBV_HbT.RH.normS,2)
    loglog(data.Blank_SAP.All.CBV_HbT.RH.meanCortf,data.Blank_SAP.All.CBV_HbT.RH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.All.CBV_HbT.RH.normS,2)
    loglog(data.SSP_SAP.All.CBV_HbT.RH.meanCortf,data.SSP_SAP.All.CBV_HbT.RH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
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
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Power Spectrum - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'IndividualPowerSpec_HbT']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'IndividualPowerSpec_HbT'])
end
%% average gamma-band power
summaryFigure3 = figure;
sgtitle('Gamma-band [30-100] Hz (envelope) cortical power spectra')
%% LH power spectra of gamma-band power during Rest
ax1 = subplot(3,4,1);
p1 = loglog(data.Naive.Rest.gammaBandPower.LH.meanCortf,data.Naive.Rest.gammaBandPower.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
p2 = loglog(data.Blank_SAP.Rest.gammaBandPower.LH.meanCortf,data.Blank_SAP.Rest.gammaBandPower.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
p3 = loglog(data.SSP_SAP.Rest.gammaBandPower.LH.meanCortf,data.SSP_SAP.Rest.gammaBandPower.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Rest] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
xlim([1/10,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Rest
ax2 = subplot(3,4,2);
loglog(data.Naive.Rest.gammaBandPower.RH.meanCortf,data.Naive.Rest.gammaBandPower.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Rest.gammaBandPower.RH.meanCortf,data.Blank_SAP.Rest.gammaBandPower.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Rest.gammaBandPower.RH.meanCortf,data.SSP_SAP.Rest.gammaBandPower.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Rest] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/10,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during NREM
ax3 = subplot(3,4,3);
loglog(data.Naive.NREM.gammaBandPower.LH.meanCortf,data.Naive.NREM.gammaBandPower.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.NREM.gammaBandPower.LH.meanCortf,data.Blank_SAP.NREM.gammaBandPower.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.NREM.gammaBandPower.LH.meanCortf,data.SSP_SAP.NREM.gammaBandPower.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[NREM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
xlim([1/30,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during NREM
ax4 = subplot(3,4,4);
loglog(data.Naive.NREM.gammaBandPower.RH.meanCortf,data.Naive.NREM.gammaBandPower.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.NREM.gammaBandPower.RH.meanCortf,data.Blank_SAP.NREM.gammaBandPower.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.NREM.gammaBandPower.RH.meanCortf,data.SSP_SAP.NREM.gammaBandPower.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[NREM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/30,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during REM
ax5 = subplot(3,4,5);
loglog(data.Naive.REM.gammaBandPower.LH.meanCortf,data.Naive.REM.gammaBandPower.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.REM.gammaBandPower.LH.meanCortf,data.Blank_SAP.REM.gammaBandPower.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.REM.gammaBandPower.LH.meanCortf,data.SSP_SAP.REM.gammaBandPower.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[REM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
xlim([1/60,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during REM
ax6 = subplot(3,4,6);
loglog(data.Naive.REM.gammaBandPower.RH.meanCortf,data.Naive.REM.gammaBandPower.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.REM.gammaBandPower.RH.meanCortf,data.Blank_SAP.REM.gammaBandPower.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.REM.gammaBandPower.RH.meanCortf,data.SSP_SAP.REM.gammaBandPower.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[REM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/60,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during Awake
ax7 = subplot(3,4,7);
loglog(data.Naive.Awake.gammaBandPower.LH.meanCortf,data.Naive.Awake.gammaBandPower.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Awake.gammaBandPower.LH.meanCortf,data.Blank_SAP.Awake.gammaBandPower.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Awake.gammaBandPower.LH.meanCortf,data.SSP_SAP.Awake.gammaBandPower.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Alert] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Awake
ax8 = subplot(3,4,8);
loglog(data.Naive.Awake.gammaBandPower.RH.meanCortf,data.Naive.Awake.gammaBandPower.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Awake.gammaBandPower.RH.meanCortf,data.Blank_SAP.Awake.gammaBandPower.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Awake.gammaBandPower.RH.meanCortf,data.SSP_SAP.Awake.gammaBandPower.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Alert] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during Asleep
ax9 = subplot(3,4,9);
loglog(data.Naive.Asleep.gammaBandPower.LH.meanCortf,data.Naive.Asleep.gammaBandPower.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Asleep.gammaBandPower.LH.meanCortf,data.Blank_SAP.Asleep.gammaBandPower.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Asleep.gammaBandPower.LH.meanCortf,data.SSP_SAP.Asleep.gammaBandPower.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[AAsleep] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Asleep
ax10 = subplot(3,4,10);
loglog(data.Naive.Asleep.gammaBandPower.RH.meanCortf,data.Naive.Asleep.gammaBandPower.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Asleep.gammaBandPower.RH.meanCortf,data.Blank_SAP.Asleep.gammaBandPower.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Asleep.gammaBandPower.RH.meanCortf,data.SSP_SAP.Asleep.gammaBandPower.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[AAsleep] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during All data
ax11 = subplot(3,4,11);
loglog(data.Naive.All.gammaBandPower.LH.meanCortf,data.Naive.All.gammaBandPower.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.All.gammaBandPower.LH.meanCortf,data.Blank_SAP.All.gammaBandPower.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.All.gammaBandPower.LH.meanCortf,data.SSP_SAP.All.gammaBandPower.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[All] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during All data
ax12 = subplot(3,4,12);
loglog(data.Naive.All.gammaBandPower.RH.meanCortf,data.Naive.All.gammaBandPower.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.All.gammaBandPower.RH.meanCortf,data.Blank_SAP.All.gammaBandPower.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.All.gammaBandPower.RH.meanCortf,data.SSP_SAP.All.gammaBandPower.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
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
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Power Spectrum - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure3,[dirpath 'AveragePowerSpec_Gamma']);
    set(summaryFigure3,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'AveragePowerSpec_Gamma'])
end
%% individual gamma-band power
summaryFigure4 = figure;
sgtitle('Gamma-band [30-100] Hz (envelope) cortical power spectra')
%% LH power spectra of gamma-band power during Rest
ax1 = subplot(3,4,1);
% Naive
for aa = 1:size(data.Naive.Rest.gammaBandPower.LH.normS,2)
    loglog(data.Naive.Rest.gammaBandPower.LH.meanCortf,data.Naive.Rest.gammaBandPower.LH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Rest.gammaBandPower.LH.normS,2)
    loglog(data.Blank_SAP.Rest.gammaBandPower.LH.meanCortf,data.Blank_SAP.Rest.gammaBandPower.LH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Rest.gammaBandPower.LH.normS,2)
    loglog(data.SSP_SAP.Rest.gammaBandPower.LH.meanCortf,data.SSP_SAP.Rest.gammaBandPower.LH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Rest] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/10,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Rest
ax2 = subplot(3,4,2);
% Naive
for aa = 1:size(data.Naive.Rest.gammaBandPower.RH.normS,2)
    loglog(data.Naive.Rest.gammaBandPower.RH.meanCortf,data.Naive.Rest.gammaBandPower.RH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Rest.gammaBandPower.RH.normS,2)
    loglog(data.Blank_SAP.Rest.gammaBandPower.RH.meanCortf,data.Blank_SAP.Rest.gammaBandPower.RH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Rest.gammaBandPower.RH.normS,2)
    loglog(data.SSP_SAP.Rest.gammaBandPower.RH.meanCortf,data.SSP_SAP.Rest.gammaBandPower.RH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Rest] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/10,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during NREM
ax3 = subplot(3,4,3);
% Naive
for aa = 1:size(data.Naive.NREM.gammaBandPower.LH.normS,2)
    loglog(data.Naive.NREM.gammaBandPower.LH.meanCortf,data.Naive.NREM.gammaBandPower.LH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.NREM.gammaBandPower.LH.normS,2)
    loglog(data.Blank_SAP.NREM.gammaBandPower.LH.meanCortf,data.Blank_SAP.NREM.gammaBandPower.LH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.NREM.gammaBandPower.LH.normS,2)
    loglog(data.SSP_SAP.NREM.gammaBandPower.LH.meanCortf,data.SSP_SAP.NREM.gammaBandPower.LH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[NREM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/30,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during NREM
ax4 = subplot(3,4,4);
% Naive
for aa = 1:size(data.Naive.NREM.gammaBandPower.RH.normS,2)
    loglog(data.Naive.NREM.gammaBandPower.RH.meanCortf,data.Naive.NREM.gammaBandPower.RH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.NREM.gammaBandPower.RH.normS,2)
    loglog(data.Blank_SAP.NREM.gammaBandPower.RH.meanCortf,data.Blank_SAP.NREM.gammaBandPower.RH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.NREM.gammaBandPower.RH.normS,2)
    loglog(data.SSP_SAP.NREM.gammaBandPower.RH.meanCortf,data.SSP_SAP.NREM.gammaBandPower.RH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[NREM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/30,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during REM
ax5 = subplot(3,4,5);
% Naive
for aa = 1:size(data.Naive.REM.gammaBandPower.LH.normS,2)
    loglog(data.Naive.REM.gammaBandPower.LH.meanCortf,data.Naive.REM.gammaBandPower.LH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.REM.gammaBandPower.LH.normS,2)
    loglog(data.Blank_SAP.REM.gammaBandPower.LH.meanCortf,data.Blank_SAP.REM.gammaBandPower.LH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.REM.gammaBandPower.LH.normS,2)
    loglog(data.SSP_SAP.REM.gammaBandPower.LH.meanCortf,data.SSP_SAP.REM.gammaBandPower.LH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[REM] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/60,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during REM
ax6 = subplot(3,4,6);
% Naive
for aa = 1:size(data.Naive.REM.gammaBandPower.RH.normS,2)
    loglog(data.Naive.REM.gammaBandPower.RH.meanCortf,data.Naive.REM.gammaBandPower.RH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.REM.gammaBandPower.RH.normS,2)
    loglog(data.Blank_SAP.REM.gammaBandPower.RH.meanCortf,data.Blank_SAP.REM.gammaBandPower.RH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.REM.gammaBandPower.RH.normS,2)
    loglog(data.SSP_SAP.REM.gammaBandPower.RH.meanCortf,data.SSP_SAP.REM.gammaBandPower.RH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[REM] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1/60,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during Awake
ax7 = subplot(3,4,7);
% Naive
for aa = 1:size(data.Naive.Awake.gammaBandPower.LH.normS,2)
    loglog(data.Naive.Awake.gammaBandPower.LH.meanCortf,data.Naive.Awake.gammaBandPower.LH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Awake.gammaBandPower.LH.normS,2)
    loglog(data.Blank_SAP.Awake.gammaBandPower.LH.meanCortf,data.Blank_SAP.Awake.gammaBandPower.LH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Awake.gammaBandPower.LH.normS,2)
    loglog(data.SSP_SAP.Awake.gammaBandPower.LH.meanCortf,data.SSP_SAP.Awake.gammaBandPower.LH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Alert] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Awake
ax8 = subplot(3,4,8);
% Naive
for aa = 1:size(data.Naive.Awake.gammaBandPower.RH.normS,2)
    loglog(data.Naive.Awake.gammaBandPower.RH.meanCortf,data.Naive.Awake.gammaBandPower.RH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Awake.gammaBandPower.RH.normS,2)
    loglog(data.Blank_SAP.Awake.gammaBandPower.RH.meanCortf,data.Blank_SAP.Awake.gammaBandPower.RH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Awake.gammaBandPower.RH.normS,2)
    loglog(data.SSP_SAP.Awake.gammaBandPower.RH.meanCortf,data.SSP_SAP.Awake.gammaBandPower.RH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[Alert] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during Asleep
ax9 = subplot(3,4,9);
% Naive
for aa = 1:size(data.Naive.Asleep.gammaBandPower.LH.normS,2)
    loglog(data.Naive.Asleep.gammaBandPower.LH.meanCortf,data.Naive.Asleep.gammaBandPower.LH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Asleep.gammaBandPower.LH.normS,2)
    loglog(data.Blank_SAP.Asleep.gammaBandPower.LH.meanCortf,data.Blank_SAP.Asleep.gammaBandPower.LH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Asleep.gammaBandPower.LH.normS,2)
    loglog(data.SSP_SAP.Asleep.gammaBandPower.LH.meanCortf,data.SSP_SAP.Asleep.gammaBandPower.LH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[AAsleep] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Asleep
ax10 = subplot(3,4,10);
% Naive
for aa = 1:size(data.Naive.Asleep.gammaBandPower.RH.normS,2)
    loglog(data.Naive.Asleep.gammaBandPower.RH.meanCortf,data.Naive.Asleep.gammaBandPower.RH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Asleep.gammaBandPower.RH.normS,2)
    loglog(data.Blank_SAP.Asleep.gammaBandPower.RH.meanCortf,data.Blank_SAP.Asleep.gammaBandPower.RH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Asleep.gammaBandPower.RH.normS,2)
    loglog(data.SSP_SAP.Asleep.gammaBandPower.RH.meanCortf,data.SSP_SAP.Asleep.gammaBandPower.RH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[AAsleep] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% LH power spectra of gamma-band power during All data
ax11 = subplot(3,4,11);
% Naive
for aa = 1:size(data.Naive.All.gammaBandPower.LH.normS,2)
    loglog(data.Naive.All.gammaBandPower.LH.meanCortf,data.Naive.All.gammaBandPower.LH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.All.gammaBandPower.LH.normS,2)
    loglog(data.Blank_SAP.All.gammaBandPower.LH.meanCortf,data.Blank_SAP.All.gammaBandPower.LH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.All.gammaBandPower.LH.normS,2)
    loglog(data.SSP_SAP.All.gammaBandPower.LH.meanCortf,data.SSP_SAP.All.gammaBandPower.LH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('[All] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5])
set(gca,'box','off')
%% RH power spectra of gamma-band power during All data
ax12 = subplot(3,4,12);
% Naive
for aa = 1:size(data.Naive.All.gammaBandPower.RH.normS,2)
    loglog(data.Naive.All.gammaBandPower.RH.meanCortf,data.Naive.All.gammaBandPower.RH.normS(:,aa),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.All.gammaBandPower.RH.normS,2)
    loglog(data.Blank_SAP.All.gammaBandPower.RH.meanCortf,data.Blank_SAP.All.gammaBandPower.RH.normS(:,aa),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.All.gammaBandPower.RH.normS,2)
    loglog(data.SSP_SAP.All.gammaBandPower.RH.meanCortf,data.SSP_SAP.All.gammaBandPower.RH.normS(:,aa),'color',colors('electric purple'),'LineWidth',0.5);
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
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Power Spectrum - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure4,[dirpath 'IndividualPowerSpec_Gamma']);
    set(summaryFigure4,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'IndividualPowerSpec_Gamma'])
end

%%
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    %% statistical diary
    diaryFile = [dirpath 'PowerSpec_' dataType '_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % Awake stats
    disp('======================================================================================================================')
    disp(['GLME statistics for ' dataType 'Power spectrum for Awake data'])
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.(dataType).Awake.S001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.(dataType).Awake.S01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.(dataType).Awake.S05.Stats)
    % Sleep stats
    disp('======================================================================================================================')
    disp(['GLME statistics for ' dataType 'Power spectrum for Sleep data'])
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.(dataType).Sleep.S001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.(dataType).Sleep.S01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.(dataType).Sleep.S05.Stats)
    % All stats
    disp('======================================================================================================================')
    disp(['GLME statistics for ' dataType 'Power spectrum for All data'])
    disp('======================================================================================================================')
    disp('0 -> 0.01 Hz')
    disp(Stats.(dataType).All.S001.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.01 -> 0.1 Hz')
    disp(Stats.(dataType).All.S01.Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp('0.1 -> 0.5 Hz')
    disp(Stats.(dataType).All.S05.Stats)
    diary off
end

end
