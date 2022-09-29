function [] = StimEvoked_PulseTrain_2PLSM(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% set-up
resultsStruct = 'Results_VesselEvoked';
load(resultsStruct);
expGroups = {'SSP-SAP','Blank-SAP'};
setName = '2PLSM Set B';
animalIDs.all = {};
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs.all = horzcat(animalIDs.all,{folderList.name});
    animalIDs.(strrep(expGroups{1,aa},'-','_')) = {folderList.name};
end
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
treatments = {'SSP_SAP','Blank_SAP'};
data = [];
cortVariables = {'diameter','baseline','timeVector','count'};
%% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs.all)
    % recognize treatment based on animal group
    if ismember(animalIDs.all{1,aa},animalIDs.SSP_SAP) == true
        treatment = 'SSP_SAP';
    elseif ismember(animalIDs.all{1,aa},animalIDs.Blank_SAP) == true
        treatment = 'Blank_SAP';
    end
    for bb = 1:length(solenoidNames)
        % left, right hemishpere hemo & neural data
        vesselIDs = fieldnames(Results_VesselEvoked.(animalIDs.all{1,aa}).Stim.(solenoidNames{1,bb}));
        for cc = 1:length(vesselIDs)
            % pre-allocate necessary variable fields
            data.(treatment).(solenoidNames{1,bb}).dummCheck = 1;
            for dd = 1:length(cortVariables)
                if isfield(data.(treatment).(solenoidNames{1,bb}),(cortVariables{1,dd})) == false
                    data.(treatment).(solenoidNames{1,bb}).(cortVariables{1,dd}) = [];
                end
            end
            data.(treatment).(solenoidNames{1,bb}).diameter = cat(1,data.(treatment).(solenoidNames{1,bb}).diameter,Results_VesselEvoked.(animalIDs.all{1,aa}).Stim.(solenoidNames{1,bb}).(vesselIDs{cc,1}).mean);
            data.(treatment).(solenoidNames{1,bb}).timeVector = cat(1,data.(treatment).(solenoidNames{1,bb}).timeVector,Results_VesselEvoked.(animalIDs.all{1,aa}).Stim.(solenoidNames{1,bb}).(vesselIDs{cc,1}).timeVector);
            data.(treatment).(solenoidNames{1,bb}).count = cat(1,data.(treatment).(solenoidNames{1,bb}).count,Results_VesselEvoked.(animalIDs.all{1,aa}).Stim.(solenoidNames{1,bb}).(vesselIDs{cc,1}).count);
            data.(treatment).(solenoidNames{1,bb}).baseline = cat(1,data.(treatment).(solenoidNames{1,bb}).baseline,Results_VesselEvoked.(animalIDs.all{1,aa}).Stim.(solenoidNames{1,bb}).(vesselIDs{cc,1}).baseline);
        end
    end
end
%% concatenate the data from the contra and ipsi data
% contra
for ff = 1:length(treatments)
    for gg = 1:length(cortVariables)
        data.(treatments{1,ff}).Contra.(cortVariables{1,gg}) = data.(treatments{1,ff}).LPadSol.(cortVariables{1,gg});
    end
end
% Ipsi
for ff = 1:length(treatments)
    for gg = 1:length(cortVariables)
        data.(treatments{1,ff}).Ipsi.(cortVariables{1,gg}) = data.(treatments{1,ff}).RPadSol.(cortVariables{1,gg});
    end
end
% auditory
for ff = 1:length(treatments)
    for gg = 1:length(cortVariables)
        data.(treatments{1,ff}).Auditory.(cortVariables{1,gg}) = data.(treatments{1,ff}).AudSol.(cortVariables{1,gg});
    end
end
%% take the averages of each field through the proper dimension
for ee = 1:length(treatments)
    treatment = treatments{1,ee};
    for ff = 1:length(compDataTypes)
        compDataType = compDataTypes{1,ff};
        data.(treatment).(compDataType).meanDiameter = mean(data.(treatment).(compDataType).diameter,1);
        data.(treatment).(compDataType).stdErrDiameter = std(data.(treatment).(compDataType).diameter,0,1)./sqrt(size(data.(treatment).(compDataType).diameter,1));
        data.(treatment).(compDataType).meanTimeVector = mean(data.(treatment).(compDataType).timeVector,1);
        data.(treatment).(compDataType).meanCount = mean(data.(treatment).(compDataType).count,1);
        data.(treatment).(compDataType).stdCount = std(data.(treatment).(compDataType).count,0,1);
        data.(treatment).baselines = data.(treatment).LPadSol.baseline;   % should be the same for each data type
    end
end
%% time to peak and peak
for gg = 1:length(treatments)
    treatment = treatments{1,gg};
    compData.(treatment).baselines = data.(treatment).baselines;
    offset = 2;
    vesselSamplingRate = 5;
    for hh = 1:length(data.(treatment).baselines)
        [compData.(treatment).peakDiameter(hh,1),timeToPeak] = max(data.(treatment).Contra.diameter(hh,offset*vesselSamplingRate:end));
        compData.(treatment).timeToPeak(hh,1) = timeToPeak/vesselSamplingRate + offset;
    end
end
%% average stim-evoked figures
summaryFigure1 = figure;
sgtitle('Stimulus-evoked \DeltaD/D repsonses')
%% RH contra stim
ax1 = subplot(1,3,1);
% Blank-SAP
p1 = plot(data.Blank_SAP.Contra.meanTimeVector,data.Blank_SAP.Contra.meanDiameter,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(data.Blank_SAP.Contra.meanTimeVector,data.Blank_SAP.Contra.meanDiameter + data.Blank_SAP.Contra.stdErrDiameter,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Contra.meanTimeVector,data.Blank_SAP.Contra.meanDiameter - data.Blank_SAP.Contra.stdErrDiameter,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
p2 = plot(data.SSP_SAP.Contra.meanTimeVector,data.SSP_SAP.Contra.meanDiameter,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Contra.meanTimeVector,data.SSP_SAP.Contra.meanDiameter + data.SSP_SAP.Contra.stdErrDiameter,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Contra.meanTimeVector,data.SSP_SAP.Contra.meanDiameter - data.SSP_SAP.Contra.stdErrDiameter,'color',colors('electric purple'),'LineWidth',0.5)
title('Contra Stim')
ylabel('\DeltaD/D (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
axis square
%% RH ipsi stim
ax2 = subplot(1,3,2);
% Blank-SAP
plot(data.Blank_SAP.Ipsi.meanTimeVector,data.Blank_SAP.Ipsi.meanDiameter,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(data.Blank_SAP.Ipsi.meanTimeVector,data.Blank_SAP.Ipsi.meanDiameter + data.Blank_SAP.Ipsi.stdErrDiameter,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Ipsi.meanTimeVector,data.Blank_SAP.Ipsi.meanDiameter - data.Blank_SAP.Ipsi.stdErrDiameter,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Ipsi.meanTimeVector,data.SSP_SAP.Ipsi.meanDiameter,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Ipsi.meanTimeVector,data.SSP_SAP.Ipsi.meanDiameter + data.SSP_SAP.Ipsi.stdErrDiameter,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Ipsi.meanTimeVector,data.SSP_SAP.Ipsi.meanDiameter - data.SSP_SAP.Ipsi.stdErrDiameter,'color',colors('electric purple'),'LineWidth',0.5)
title('Ipsi Stim')
ylabel('\DeltaD/D (%)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
%% RH auditory stim
ax3 = subplot(1,3,3);
% Blank-SAP
plot(data.Blank_SAP.Auditory.meanTimeVector,data.Blank_SAP.Auditory.meanDiameter,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(data.Blank_SAP.Auditory.meanTimeVector,data.Blank_SAP.Auditory.meanDiameter + data.Blank_SAP.Auditory.stdErrDiameter,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Auditory.meanTimeVector,data.Blank_SAP.Auditory.meanDiameter - data.Blank_SAP.Auditory.stdErrDiameter,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Auditory.meanTimeVector,data.SSP_SAP.Auditory.meanDiameter,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Auditory.meanTimeVector,data.SSP_SAP.Auditory.meanDiameter + data.SSP_SAP.Auditory.stdErrDiameter,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Auditory.meanTimeVector,data.SSP_SAP.Auditory.meanDiameter - data.SSP_SAP.Auditory.stdErrDiameter,'color',colors('electric purple'),'LineWidth',0.5)
title('Auditory Stim')
ylabel('\DeltaD/D (%)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
%% figure characteristics
linkaxes([ax1,ax2,ax3],'xy')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Stimulus Evoked - Pulse Train 2PLSM' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'AverageStimEvoked_2PLSM']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'AverageStimEvoked_2PLSM'])
end
%% average stim-evoked figures
summaryFigure2 = figure;
sgtitle('Stimulus-evoked \DeltaD/D repsonses - individual animals')
%% RH contra stim
ax1 = subplot(1,3,1);
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Contra.diameter,1)
    p1 = plot(data.Blank_SAP.Contra.meanTimeVector,data.Blank_SAP.Contra.diameter(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Contra.diameter,1)
    p2 = plot(data.SSP_SAP.Contra.meanTimeVector,data.SSP_SAP.Contra.diameter(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('Contra Stim')
ylabel('\DeltaD/D (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
axis square
%% RH ipsi stim
ax2 = subplot(1,3,2);
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Contra.diameter,1)
    plot(data.Blank_SAP.Ipsi.meanTimeVector,data.Blank_SAP.Ipsi.diameter(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Contra.diameter,1)
    plot(data.SSP_SAP.Ipsi.meanTimeVector,data.SSP_SAP.Ipsi.diameter(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('Ipsi Stim')
ylabel('\DeltaD/D (%)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
%% RH auditory stim
ax3 = subplot(1,3,3);
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Contra.diameter,1)
    plot(data.Blank_SAP.Auditory.meanTimeVector,data.Blank_SAP.Auditory.diameter(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Contra.diameter,1)
    plot(data.SSP_SAP.Auditory.meanTimeVector,data.SSP_SAP.Auditory.diameter(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('Auditory Stim')
ylabel('\DeltaD/D (%)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
%% figure characteristics
linkaxes([ax1,ax2,ax3],'xy')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Stimulus Evoked - Pulse Train 2PLSM' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'IndividualStimEvoked_2PLSM']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'IndividualStimEvoked_2PLSM'])
end
%% average stim-evoked figures
summaryFigure3 = figure;
sgtitle('Stimulus-evoked \DeltaD/D repsonses')
%% RH contra stim
subplot(1,2,1);
% Blank-SAP
s1 = scatter(compData.Blank_SAP.baselines,compData.Blank_SAP.peakDiameter,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'));
hold on
% SSP-SAP
s2 = scatter(compData.SSP_SAP.baselines,compData.SSP_SAP.peakDiameter,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'));
title('Baseline vs. Peak Diameter')
xlabel('Baseline diameter (\mum)')
ylabel('Peak dilation \DeltaD/D (%)')
legend([s1,s2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
axis square
%% RH contra stim
subplot(1,2,2);
% Blank-SAP
scatter(compData.Blank_SAP.baselines,compData.Blank_SAP.timeToPeak,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'));
hold on
% SSP-SAP
scatter(compData.SSP_SAP.baselines,compData.SSP_SAP.timeToPeak,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'));
title('Baseline vs. Time-to-Peak')
xlabel('Baseline diameter (\mum)')
ylabel('Time-to-Peak (sec)')
set(gca,'box','off')
axis square
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Stimulus Evoked - Pulse Train 2PLSM' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure3,[dirpath 'DiameterVsDilation_2PLSM']);
    set(summaryFigure3,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'DiameterVsDilation_2PLSM'])
end

end
