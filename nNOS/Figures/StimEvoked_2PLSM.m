function [] = StimEvoked_2PLSM(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% set-up
resultsStruct = 'Results_Evoked';
load(resultsStruct);
expGroups = {'SSP-SAP','Blank-SAP'};
setName = 'IOS Set B';
animalIDs.all = {};
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs.all = horzcat(animalIDs.all,{folderList.name});
    animalIDs.(strrep(expGroups{1,aa},'-','_')) = {folderList.name};
end
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
hemispheres = {'adjBarrels'};
treatments = {'SSP_SAP','Blank_SAP'};
data = [];
cortVariables = {'HbT','timeVector','count'};
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
        for cc = 1:length(hemispheres)
            % pre-allocate necessary variable fields
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).dummCheck = 1;
            for dd = 1:length(cortVariables)
                if isfield(data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}),cortVariables{1,dd}) == false
                    data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).(cortVariables{1,dd}) = [];
                end
            end
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).HbT = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).HbT,Results_Evoked.(animalIDs.all{1,aa}).Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).CBV_HbT.HbT);
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).timeVector = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).timeVector,Results_Evoked.(animalIDs.all{1,aa}).Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).timeVector);
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).count = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).count,Results_Evoked.(animalIDs.all{1,aa}).Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).count);
        end
    end
end
%% concatenate the data from the contra and ipsi data
% contra
for ff = 1:length(treatments)
    for gg = 1:length(cortVariables)
        data.(treatments{1,ff}).Contra.adjBarrels.(cortVariables{1,gg}) = data.(treatments{1,ff}).LPadSol.adjBarrels.(cortVariables{1,gg});
    end
end
% Ipsi
for ff = 1:length(treatments)
    for gg = 1:length(cortVariables)
        data.(treatments{1,ff}).Ipsi.adjBarrels.(cortVariables{1,gg}) = data.(treatments{1,ff}).RPadSol.adjBarrels.(cortVariables{1,gg});
    end
end
% auditory
for ff = 1:length(treatments)
    for gg = 1:length(cortVariables)
        data.(treatments{1,ff}).Auditory.adjBarrels.(cortVariables{1,gg}) = data.(treatments{1,ff}).AudSol.adjBarrels.(cortVariables{1,gg});
    end
end
%% take the averages of each field through the proper dimension
for ee = 1:length(treatments)
    treatment = treatments{1,ee};
    for gg = 1:length(hemispheres)
        hemisphere = hemispheres{1,gg};
        for ff = 1:length(compDataTypes)
            compDataType = compDataTypes{1,ff};
            data.(treatment).(compDataType).(hemisphere).meanHbT = mean(data.(treatment).(compDataType).(hemisphere).HbT,1);
            data.(treatment).(compDataType).(hemisphere).stdHbT = std(data.(treatment).(compDataType).(hemisphere).HbT,0,1);
            data.(treatment).(compDataType).(hemisphere).meanTimeVector = mean(data.(treatment).(compDataType).(hemisphere).timeVector,1);
            data.(treatment).(compDataType).(hemisphere).meanCount = mean(data.(treatment).(compDataType).(hemisphere).count,1);
            data.(treatment).(compDataType).(hemisphere).stdCount = std(data.(treatment).(compDataType).(hemisphere).count,0,1);
        end
    end
end
%% average stim-evoked figures
summaryFigure1 = figure;
sgtitle('Stimulus-evoked \DeltaHbT repsonses')
%% RH contra stim
ax1 = subplot(1,3,1);
% Blank-SAP
p1 = plot(data.Blank_SAP.Contra.adjBarrels.meanTimeVector,data.Blank_SAP.Contra.adjBarrels.meanHbT,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(data.Blank_SAP.Contra.adjBarrels.meanTimeVector,data.Blank_SAP.Contra.adjBarrels.meanHbT + data.Blank_SAP.Contra.adjBarrels.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Contra.adjBarrels.meanTimeVector,data.Blank_SAP.Contra.adjBarrels.meanHbT - data.Blank_SAP.Contra.adjBarrels.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
p2 = plot(data.SSP_SAP.Contra.adjBarrels.meanTimeVector,data.SSP_SAP.Contra.adjBarrels.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Contra.adjBarrels.meanTimeVector,data.SSP_SAP.Contra.adjBarrels.meanHbT + data.SSP_SAP.Contra.adjBarrels.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Contra.adjBarrels.meanTimeVector,data.SSP_SAP.Contra.adjBarrels.meanHbT - data.SSP_SAP.Contra.adjBarrels.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('Contra Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
axis square
%% RH ipsi stim
ax2 = subplot(1,3,2);
% Blank-SAP
plot(data.Blank_SAP.Ipsi.adjBarrels.meanTimeVector,data.Blank_SAP.Ipsi.adjBarrels.meanHbT,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(data.Blank_SAP.Ipsi.adjBarrels.meanTimeVector,data.Blank_SAP.Ipsi.adjBarrels.meanHbT + data.Blank_SAP.Ipsi.adjBarrels.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Ipsi.adjBarrels.meanTimeVector,data.Blank_SAP.Ipsi.adjBarrels.meanHbT - data.Blank_SAP.Ipsi.adjBarrels.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Ipsi.adjBarrels.meanTimeVector,data.SSP_SAP.Ipsi.adjBarrels.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Ipsi.adjBarrels.meanTimeVector,data.SSP_SAP.Ipsi.adjBarrels.meanHbT + data.SSP_SAP.Ipsi.adjBarrels.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Ipsi.adjBarrels.meanTimeVector,data.SSP_SAP.Ipsi.adjBarrels.meanHbT - data.SSP_SAP.Ipsi.adjBarrels.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('Ipsi Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
%% RH auditory stim
ax3 = subplot(1,3,3);
% Blank-SAP
plot(data.Blank_SAP.Auditory.adjBarrels.meanTimeVector,data.Blank_SAP.Auditory.adjBarrels.meanHbT,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(data.Blank_SAP.Auditory.adjBarrels.meanTimeVector,data.Blank_SAP.Auditory.adjBarrels.meanHbT + data.Blank_SAP.Auditory.adjBarrels.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Auditory.adjBarrels.meanTimeVector,data.Blank_SAP.Auditory.adjBarrels.meanHbT - data.Blank_SAP.Auditory.adjBarrels.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Auditory.adjBarrels.meanTimeVector,data.SSP_SAP.Auditory.adjBarrels.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Auditory.adjBarrels.meanTimeVector,data.SSP_SAP.Auditory.adjBarrels.meanHbT + data.SSP_SAP.Auditory.adjBarrels.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Auditory.adjBarrels.meanTimeVector,data.SSP_SAP.Auditory.adjBarrels.meanHbT - data.SSP_SAP.Auditory.adjBarrels.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('Auditory Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
%% figure characteristics
linkaxes([ax1,ax2,ax3],'xy')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Stim Evoked - IOS Pulse Train' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'AverageStimEvoked_HbT']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'AverageStimEvoked_HbT'])
end
%% individual stim-evoked figures
summaryFigure2 = figure;
sgtitle('Stimulus-evoked \DeltaHbT repsonses - individual animals')
%% contra stim
ax1 = subplot(1,3,1);
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Contra.adjBarrels.HbT,1)
    p1 = plot(data.Blank_SAP.Contra.adjBarrels.meanTimeVector,data.Blank_SAP.Contra.adjBarrels.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Contra.adjBarrels.HbT,1)
    p2 = plot(data.SSP_SAP.Contra.adjBarrels.meanTimeVector,data.SSP_SAP.Contra.adjBarrels.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('Contra Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
axis square
%% LH ipsi stim
ax2 = subplot(1,3,2);
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Ipsi.adjBarrels.HbT,1)
    plot(data.Blank_SAP.Ipsi.adjBarrels.meanTimeVector,data.Blank_SAP.Ipsi.adjBarrels.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Ipsi.adjBarrels.HbT,1)
    plot(data.SSP_SAP.Ipsi.adjBarrels.meanTimeVector,data.SSP_SAP.Ipsi.adjBarrels.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('Ipsi Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
%% auditory stim
ax3 = subplot(1,3,3);
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Auditory.adjBarrels.HbT,1)
    plot(data.Blank_SAP.Auditory.adjBarrels.meanTimeVector,data.Blank_SAP.Auditory.adjBarrels.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Auditory.adjBarrels.HbT,1)
    plot(data.SSP_SAP.Auditory.adjBarrels.meanTimeVector,data.SSP_SAP.Auditory.adjBarrels.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('Auditory Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
%% figure characteristics
linkaxes([ax1,ax2,ax3],'xy')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Stim Evoked - IOS Pulse Train' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'IndividualStimEvoked_HbT']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'IndividualStimEvoked_HbT'])
end

end
