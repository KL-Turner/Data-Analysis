function [] = AwakeEvokedResponses_GarborgTBD(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figures and supporting information for Figure Panel 2
%________________________________________________________________________________________________________________________

dataLocation = [rootFolder delim 'Analysis Structures'];
cd(dataLocation)
%% stimulus and whisking evoked pupil changes
resultsStruct = 'Results_Evoked.mat';
load(resultsStruct);
dataTypes = {'SSS','lSSS','rSSS'};
timeVector = (0:12*30)/30 - 2;
animalIDs = fieldnames(Results_Evoked);
% pre-allocate
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    data.Evoked.controlSolenoid.(dataType).data = [];
    data.Evoked.stimSolenoid.(dataType).data = [];
    data.Evoked.briefWhisk.(dataType).data = [];
    data.Evoked.interWhisk.(dataType).data = [];
    data.Evoked.extendWhisk.(dataType).data = [];
end
% concatenate each stimuli type from each animal
for cc = 1:length(animalIDs)
    animalID = animalIDs{cc,1};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        % whisking
        data.Evoked.briefWhisk.(dataType).data = cat(1,data.Evoked.briefWhisk.(dataType).data,Results_Evoked.(animalID).Whisk.(dataType).ShortWhisks.mean);
        data.Evoked.interWhisk.(dataType).data = cat(1,data.Evoked.interWhisk.(dataType).data,Results_Evoked.(animalID).Whisk.(dataType).IntermediateWhisks.mean);
        data.Evoked.extendWhisk.(dataType).data = cat(1,data.Evoked.extendWhisk.(dataType).data,Results_Evoked.(animalID).Whisk.(dataType).LongWhisks.mean);
        % solenoids
        data.Evoked.stimSolenoid.(dataType).data = cat(1,data.Evoked.stimSolenoid.(dataType).data,Results_Evoked.(animalID).Stim.(dataType).LPadSol.mean,Results_Evoked.(animalID).Stim.(dataType).RPadSol.mean);
        data.Evoked.controlSolenoid.(dataType).data = cat(1,data.Evoked.controlSolenoid.(dataType).data,Results_Evoked.(animalID).Stim.(dataType).AudSol.mean);
    end
end
% mean and standard error for each stimulation
for dd = 1:length(dataTypes)
    dataType = dataTypes{1,dd};
    data.Evoked.briefWhisk.(dataType).mean = mean(data.Evoked.briefWhisk.(dataType).data,1);
    data.Evoked.briefWhisk.(dataType).sem = std(data.Evoked.briefWhisk.(dataType).data,0,1)./sqrt(size(data.Evoked.briefWhisk.(dataType).data,1));
    data.Evoked.briefWhisk.(dataType).peak = max(data.Evoked.briefWhisk.(dataType).data,[],2);
    data.Evoked.briefWhisk.(dataType).meanPeak = mean(data.Evoked.briefWhisk.(dataType).peak,1);
    data.Evoked.briefWhisk.(dataType).stdPeak = std(data.Evoked.briefWhisk.(dataType).peak,0,1);
    data.Evoked.interWhisk.(dataType).mean = mean(data.Evoked.interWhisk.(dataType).data,1);
    data.Evoked.interWhisk.(dataType).sem = std(data.Evoked.interWhisk.(dataType).data,0,1)./sqrt(size(data.Evoked.interWhisk.(dataType).data,1));
    data.Evoked.interWhisk.(dataType).peak = max(data.Evoked.interWhisk.(dataType).data,[],2);
    data.Evoked.interWhisk.(dataType).meanPeak = mean(data.Evoked.interWhisk.(dataType).peak,1);
    data.Evoked.interWhisk.(dataType).stdPeak = std(data.Evoked.interWhisk.(dataType).peak,0,1);
    data.Evoked.extendWhisk.(dataType).mean = mean(data.Evoked.extendWhisk.(dataType).data,1);
    data.Evoked.extendWhisk.(dataType).sem = std(data.Evoked.extendWhisk.(dataType).data,0,1)./sqrt(size(data.Evoked.extendWhisk.(dataType).data,1));
    data.Evoked.extendWhisk.(dataType).peak = max(data.Evoked.extendWhisk.(dataType).data,[],2);
    data.Evoked.extendWhisk.(dataType).meanPeak = mean(data.Evoked.extendWhisk.(dataType).peak,1);
    data.Evoked.extendWhisk.(dataType).stdPeak = std(data.Evoked.extendWhisk.(dataType).peak,0,1);
    data.Evoked.stimSolenoid.(dataType).mean = mean(data.Evoked.stimSolenoid.(dataType).data,1);
    data.Evoked.stimSolenoid.(dataType).sem = std(data.Evoked.stimSolenoid.(dataType).data,0,1)./sqrt(size(data.Evoked.stimSolenoid.(dataType).data,1));
    data.Evoked.stimSolenoid.(dataType).peak = max(data.Evoked.stimSolenoid.(dataType).data,[],2);
    data.Evoked.stimSolenoid.(dataType).meanPeak = mean(data.Evoked.stimSolenoid.(dataType).peak,1);
    data.Evoked.stimSolenoid.(dataType).stdPeak = std(data.Evoked.stimSolenoid.(dataType).peak,0,1);
    data.Evoked.controlSolenoid.(dataType).mean = mean(data.Evoked.controlSolenoid.(dataType).data,1);
    data.Evoked.controlSolenoid.(dataType).sem = std(data.Evoked.controlSolenoid.(dataType).data,0,1)./sqrt(size(data.Evoked.controlSolenoid.(dataType).data,1));
    data.Evoked.controlSolenoid.(dataType).peak = max(data.Evoked.controlSolenoid.(dataType).data,[],2);
    data.Evoked.controlSolenoid.(dataType).meanPeak = mean(data.Evoked.controlSolenoid.(dataType).peak,1);
    data.Evoked.controlSolenoid.(dataType).stdPeak = std(data.Evoked.controlSolenoid.(dataType).peak,0,1);
end
%% figures
FigA = figure('Name','Figure Panel A - Garborg et al. TBD','Units','Normalized','OuterPosition',[0,0,1,1]);
%% stimulus evoked changes
ax1 = subplot(2,2,1);

p1 = plot(timeVector,data.Evoked.stimSolenoid.lSSS.mean,'color',colors('red'),'LineWidth',2);
hold on
plot(timeVector,data.Evoked.stimSolenoid.lSSS.mean + data.Evoked.stimSolenoid.lSSS.sem,'color',colors('red'),'LineWidth',0.5)
plot(timeVector,data.Evoked.stimSolenoid.lSSS.mean - data.Evoked.stimSolenoid.lSSS.sem,'color',colors('red'),'LineWidth',0.5)

p2 = plot(timeVector,data.Evoked.stimSolenoid.SSS.mean,'color',colors('green'),'LineWidth',2);
plot(timeVector,data.Evoked.stimSolenoid.SSS.mean + data.Evoked.stimSolenoid.SSS.sem,'color',colors('green'),'LineWidth',0.5)
plot(timeVector,data.Evoked.stimSolenoid.SSS.mean - data.Evoked.stimSolenoid.SSS.sem,'color',colors('green'),'LineWidth',0.5)

p3 = plot(timeVector,data.Evoked.stimSolenoid.rSSS.mean,'color',colors('blue'),'LineWidth',2);
plot(timeVector,data.Evoked.stimSolenoid.rSSS.mean + data.Evoked.stimSolenoid.rSSS.sem,'color',colors('blue'),'LineWidth',0.5)
plot(timeVector,data.Evoked.stimSolenoid.rSSS.mean - data.Evoked.stimSolenoid.rSSS.sem,'color',colors('blue'),'LineWidth',0.5)

title('Evoked diameter (z-units)')
ylabel('\Deltaz-units')
xlabel('Time (s)')
legend([p1,p2,p3],'lSSS','SSS','rSSS')
set(gca,'box','off')
xlim([-2,10])
axis square
ax1.TickLength = [0.03,0.03];

%% stimulus evoked changes
ax2 = subplot(2,2,2);

plot(timeVector,data.Evoked.controlSolenoid.lSSS.mean,'color',colors('red'),'LineWidth',2);
hold on
plot(timeVector,data.Evoked.controlSolenoid.lSSS.mean + data.Evoked.controlSolenoid.lSSS.sem,'color',colors('red'),'LineWidth',0.5)
plot(timeVector,data.Evoked.controlSolenoid.lSSS.mean - data.Evoked.controlSolenoid.lSSS.sem,'color',colors('red'),'LineWidth',0.5)

plot(timeVector,data.Evoked.controlSolenoid.SSS.mean,'color',colors('green'),'LineWidth',2);
plot(timeVector,data.Evoked.controlSolenoid.SSS.mean + data.Evoked.controlSolenoid.SSS.sem,'color',colors('green'),'LineWidth',0.5)
plot(timeVector,data.Evoked.controlSolenoid.SSS.mean - data.Evoked.controlSolenoid.SSS.sem,'color',colors('green'),'LineWidth',0.5)

plot(timeVector,data.Evoked.controlSolenoid.rSSS.mean,'color',colors('blue'),'LineWidth',2);
plot(timeVector,data.Evoked.controlSolenoid.rSSS.mean + data.Evoked.controlSolenoid.rSSS.sem,'color',colors('blue'),'LineWidth',0.5)
plot(timeVector,data.Evoked.controlSolenoid.rSSS.mean - data.Evoked.controlSolenoid.rSSS.sem,'color',colors('blue'),'LineWidth',0.5)

title('control solenoid (z-units)')
ylabel('\Deltaz-units')
xlabel('Time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
ax2.TickLength = [0.03,0.03];

%%
ax4 = subplot(2,3,4);
plot(timeVector,data.Evoked.briefWhisk.lSSS.mean,'color',colors('red'),'LineWidth',2);
hold on
plot(timeVector,data.Evoked.briefWhisk.lSSS.mean + data.Evoked.briefWhisk.lSSS.sem,'color',colors('red'),'LineWidth',0.5)
plot(timeVector,data.Evoked.briefWhisk.lSSS.mean - data.Evoked.briefWhisk.lSSS.sem,'color',colors('red'),'LineWidth',0.5)

plot(timeVector,data.Evoked.briefWhisk.SSS.mean,'color',colors('green'),'LineWidth',2);
plot(timeVector,data.Evoked.briefWhisk.SSS.mean + data.Evoked.briefWhisk.SSS.sem,'color',colors('green'),'LineWidth',0.5)
plot(timeVector,data.Evoked.briefWhisk.SSS.mean - data.Evoked.briefWhisk.SSS.sem,'color',colors('green'),'LineWidth',0.5)

plot(timeVector,data.Evoked.briefWhisk.rSSS.mean,'color',colors('blue'),'LineWidth',2);
plot(timeVector,data.Evoked.briefWhisk.rSSS.mean + data.Evoked.briefWhisk.lSSS.sem,'color',colors('blue'),'LineWidth',0.5)
plot(timeVector,data.Evoked.briefWhisk.rSSS.mean - data.Evoked.briefWhisk.lSSS.sem,'color',colors('blue'),'LineWidth',0.5)

ylabel('Auditory sol')
ylabel('\DeltaF/F (%)')
xlabel('Time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
ax4.TickLength = [0.03,0.03];

%% stimulus evoked changes
ax5 = subplot(2,3,5);
plot(timeVector,data.Evoked.interWhisk.lSSS.mean,'color',colors('red'),'LineWidth',2);
hold on
plot(timeVector,data.Evoked.interWhisk.lSSS.mean + data.Evoked.interWhisk.lSSS.sem,'color',colors('red'),'LineWidth',0.5)
plot(timeVector,data.Evoked.interWhisk.lSSS.mean - data.Evoked.interWhisk.lSSS.sem,'color',colors('red'),'LineWidth',0.5)

plot(timeVector,data.Evoked.interWhisk.SSS.mean,'color',colors('green'),'LineWidth',2);
plot(timeVector,data.Evoked.interWhisk.SSS.mean + data.Evoked.interWhisk.SSS.sem,'color',colors('green'),'LineWidth',0.5)
plot(timeVector,data.Evoked.interWhisk.SSS.mean - data.Evoked.interWhisk.SSS.sem,'color',colors('green'),'LineWidth',0.5)

plot(timeVector,data.Evoked.interWhisk.rSSS.mean,'color',colors('blue'),'LineWidth',2);
plot(timeVector,data.Evoked.interWhisk.rSSS.mean + data.Evoked.interWhisk.lSSS.sem,'color',colors('blue'),'LineWidth',0.5)
plot(timeVector,data.Evoked.interWhisk.rSSS.mean - data.Evoked.interWhisk.lSSS.sem,'color',colors('blue'),'LineWidth',0.5)

ylabel('Intermed whisk')
ylabel('\DeltaF/F (%)')
xlabel('Time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
ax5.TickLength = [0.03,0.03];
%%
ax6 = subplot(2,3,6);
plot(timeVector,data.Evoked.extendWhisk.lSSS.mean,'color',colors('red'),'LineWidth',2);
hold on
plot(timeVector,data.Evoked.extendWhisk.lSSS.mean + data.Evoked.extendWhisk.lSSS.sem,'color',colors('red'),'LineWidth',0.5)
plot(timeVector,data.Evoked.extendWhisk.lSSS.mean - data.Evoked.extendWhisk.lSSS.sem,'color',colors('red'),'LineWidth',0.5)

plot(timeVector,data.Evoked.extendWhisk.SSS.mean,'color',colors('green'),'LineWidth',2);
plot(timeVector,data.Evoked.extendWhisk.SSS.mean + data.Evoked.extendWhisk.SSS.sem,'color',colors('green'),'LineWidth',0.5)
plot(timeVector,data.Evoked.extendWhisk.SSS.mean - data.Evoked.extendWhisk.SSS.sem,'color',colors('green'),'LineWidth',0.5)

plot(timeVector,data.Evoked.extendWhisk.rSSS.mean,'color',colors('blue'),'LineWidth',2);
plot(timeVector,data.Evoked.extendWhisk.rSSS.mean + data.Evoked.extendWhisk.rSSS.sem,'color',colors('blue'),'LineWidth',0.5)
plot(timeVector,data.Evoked.extendWhisk.rSSS.mean - data.Evoked.extendWhisk.rSSS.sem,'color',colors('blue'),'LineWidth',0.5)

ylabel('Extended whisk')
ylabel('\DeltaF/F (%)')
xlabel('Time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
ax6.TickLength = [0.03,0.03];

figure;
%% stimulus evoked changes
for aa = 1:size(data.Evoked.stimSolenoid.SSS.data,1)
    plot(timeVector,data.Evoked.stimSolenoid.SSS.data(aa,:),'LineWidth',1);
    hold on
end
title('Evoked diameter (z-units)')
ylabel('\DeltaF/F')
xlabel('Time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
ax1.TickLength = [0.03,0.03];
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(FigA,[dirpath 'Fig2_JNeurosci2022']);
    set(FigA,'PaperPositionMode','auto');
    print('-vector','-dpdf','-bestfit',[dirpath 'Fig2_JNeurosci2022'])
    % text diary
    diaryFile = [dirpath 'Fig2_Text.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % mm diameter statistical diary
    disp('======================================================================================================================')
    disp('GLME stats for mm diameter during Rest, Whisk, Stim, NREM, and REM')
    disp('======================================================================================================================')
    disp(mmDiameterStatsA)
    disp(mmDiameterStatsB)
    disp(mmDiameterStatsC)
    disp(mmDiameterStatsD)
    disp(mmDiameterStatsE)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  diameter (mm): ' num2str(round(data.Diameter.Rest.meanDiameter,2)) ' ± ' num2str(round(data.Diameter.Whisk.stdDiameter,2)) ' (n = ' num2str(length(data.Diameter.Rest.mmDiameter)) ') mice']); disp(' ')
    disp(['Whisk diameter (mm): ' num2str(round(data.Diameter.Whisk.meanDiameter,2)) ' ± ' num2str(round(data.Diameter.Whisk.stdDiameter,2)) ' (n = ' num2str(length(data.Diameter.Whisk.mmDiameter)) ') mice']); disp(' ')
    disp(['Stim  diameter (mm): ' num2str(round(data.Diameter.Stim.meanDiameter,2)) ' ± ' num2str(round(data.Diameter.Stim.stdDiameter,2)) ' (n = ' num2str(length(data.Diameter.Stim.mmDiameter)) ') mice']); disp(' ')
    disp(['NREM  diameter (mm): ' num2str(round(data.Diameter.NREM.meanDiameter,2)) ' ± ' num2str(round(data.Diameter.NREM.stdDiameter,2)) ' (n = ' num2str(length(data.Diameter.NREM.mmDiameter)) ') mice']); disp(' ')
    disp(['REM   diameter (mm): ' num2str(round(data.Diameter.REM.meanDiameter,2)) ' ± ' num2str(round(data.Diameter.REM.stdDiameter,2)) ' (n = ' num2str(length(data.Diameter.REM.mmDiameter)) ') mice']); disp(' ')
    disp(['*p < ' num2str(alpha10A) ' **p < ' num2str(alpha10B) ' ***p < ' num2str(alpha10C)]);
    disp('----------------------------------------------------------------------------------------------------------------------')
    % z-unit Diameter statistical diary
    disp('======================================================================================================================')
    disp('GLME stats for z-unit diameter during Rest, Whisk, Stim, NREM, and REM')
    disp('======================================================================================================================')
    disp(lSSSStatsA)
    disp(lSSSStatsB)
    disp(lSSSStatsC)
    disp(lSSSStatsD)
    disp(lSSSStatsE)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  diameter (z-unit): ' num2str(round(data.Diameter.Rest.meanlSSS,2)) ' ± ' num2str(round(data.Diameter.Whisk.stdlSSS,2)) ' (n = ' num2str(length(data.Diameter.Rest.lSSS)) ') mice']); disp(' ')
    disp(['Whisk diameter (z-unit): ' num2str(round(data.Diameter.Whisk.meanlSSS,2)) ' ± ' num2str(round(data.Diameter.Whisk.stdlSSS,2)) ' (n = ' num2str(length(data.Diameter.Whisk.lSSS)) ') mice']); disp(' ')
    disp(['Stim  diameter (z-unit): ' num2str(round(data.Diameter.Stim.meanlSSS,2)) ' ± ' num2str(round(data.Diameter.Stim.stdlSSS,2)) ' (n = ' num2str(length(data.Diameter.Stim.lSSS)) ') mice']); disp(' ')
    disp(['NREM  diameter (z-unit): ' num2str(round(data.Diameter.NREM.meanlSSS,2)) ' ± ' num2str(round(data.Diameter.NREM.stdlSSS,2)) ' (n = ' num2str(length(data.Diameter.NREM.lSSS)) ') mice']); disp(' ')
    disp(['REM   diameter (z-unit): ' num2str(round(data.Diameter.REM.meanlSSS,2)) ' ± ' num2str(round(data.Diameter.REM.stdlSSS,2)) ' (n = ' num2str(length(data.Diameter.REM.lSSS)) ') mice']); disp(' ')
    disp(['*p < ' num2str(alpha10A) ' **p < ' num2str(alpha10B) ' ***p < ' num2str(alpha10C)]);
    disp('----------------------------------------------------------------------------------------------------------------------')
    % peak whisk/stim
    disp('======================================================================================================================')
    disp('Peak change in z-unit post-stimulation/whisking')
    disp('======================================================================================================================')
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Whisker stimulation z-unit increase ' num2str(data.Evoked.stimSolenoid.lSSS.meanPeak) ' ± ' num2str(data.Evoked.stimSolenoid.lSSS.stdPeak) ' (n = ' num2str(length(data.Evoked.stimSolenoid.lSSS.peak)/2) ') mice']); disp(' ')
    disp(['Auditory stimulation z-unit increase ' num2str(data.Evoked.controlSolenoid.lSSS.meanPeak) ' ± ' num2str(data.Evoked.controlSolenoid.lSSS.stdPeak) ' (n = ' num2str(length(data.Evoked.controlSolenoid.lSSS.peak)) ') mice']); disp(' ')
    disp(['Volitional whisking z-unit increase ' num2str(data.Evoked.interWhisk.lSSS.meanPeak) ' ± ' num2str(data.Evoked.interWhisk.lSSS.stdPeak) ' (n = ' num2str(length(data.Evoked.interWhisk.lSSS.peak)) ') mice']); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
end
cd(rootFolder)
end
