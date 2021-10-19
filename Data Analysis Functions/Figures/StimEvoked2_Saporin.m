function [AnalysisResults] = StimEvoked2_Saporin(rootFolder,saveFigs,delim,AnalysisResults)
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
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).HbT = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).HbT,AnalysisResults.(animalIDs.all{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).CBV_HbT.HbT);
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).timeVector = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).timeVector,AnalysisResults.(animalIDs.all{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).timeVector);
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).count = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).count,AnalysisResults.(animalIDs.all{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).count);
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
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'Stim_Evoked_HbT']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Stim_Evoked_HbT'])
end
% %% individual stim-evoked figures
% summaryFigure2 = figure;
% sgtitle('Stimulus-evoked \DeltaHbT repsonses - individual animals')
% %% LH contra stim
% ax1 = subplot(3,2,1);
% % C57BL6Js
% for aa = 1:size(data.C57BL6J.Contra.adjLH.HbT,1)
%     plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.Blank_SAP.Contra.adjLH.HbT,1)
%     plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.SSP_SAP.Contra.adjLH.HbT,1)
%     plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% title('LH (UnRx) \DeltaHbT - Contra Stim')
% ylabel('\DeltaHbT (\muM)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% RH contra stim
% ax2 = subplot(3,2,2);
% % C57BL6Js
% for aa = 1:size(data.C57BL6J.Contra.adjRH.HbT,1)
%     plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.Blank_SAP.Contra.adjRH.HbT,1)
%     plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.SSP_SAP.Contra.adjRH.HbT,1)
%     plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% title('RH (Rx) \DeltaHbT - Contra Stim')
% ylabel('\DeltaHbT (\muM)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% LH ipsi stim
% ax3 = subplot(3,2,3);
% % C57BL6Js
% for aa = 1:size(data.C57BL6J.Ipsi.adjLH.HbT,1)
%     plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.Blank_SAP.Ipsi.adjLH.HbT,1)
%     plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.SSP_SAP.Ipsi.adjLH.HbT,1)
%     plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% title('LH (UnRx) \DeltaHbT - Ipsi Stim')
% ylabel('\DeltaHbT (\muM)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% RH ipsi stim
% ax4 = subplot(3,2,4);
% % C57BL6Js
% for aa = 1:size(data.C57BL6J.Ipsi.adjRH.HbT,1)
%     plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.Blank_SAP.Ipsi.adjRH.HbT,1)
%     plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.SSP_SAP.Ipsi.adjRH.HbT,1)
%     plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% title('RH (Rx) \DeltaHbT - Ipsi Stim')
% ylabel('\DeltaHbT (\muM)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% LH auditory stim
% ax5 = subplot(3,2,5);
% % C57BL6Js
% for aa = 1:size(data.C57BL6J.Auditory.adjLH.HbT,1)
%     plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.Blank_SAP.Auditory.adjLH.HbT,1)
%     plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.SSP_SAP.Auditory.adjLH.HbT,1)
%     plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% title('LH (UnRx) \DeltaHbT - Auditory Stim')
% ylabel('\DeltaHbT (\muM)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% RH auditory stim
% ax6 = subplot(3,2,6);
% % C57BL6Js
% for aa = 1:size(data.C57BL6J.Auditory.adjRH.HbT,1)
%     plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.Blank_SAP.Auditory.adjRH.HbT,1)
%     plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.SSP_SAP.Auditory.adjRH.HbT,1)
%     plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% title('RH (Rx) \DeltaHbT - Auditory Stim')
% ylabel('\DeltaHbT (\muM)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% figure characteristics
% linkaxes([ax1,ax2],'xy')
% linkaxes([ax3,ax4],'xy')
% linkaxes([ax5,ax6],'xy')
% %% save figure(s)
% if strcmp(saveFigs,'y') == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure2,[dirpath 'indStim_Evoked_HbT']);
%     set(summaryFigure2,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'indStim_Evoked_HbT'])
% end
% %% average whisk-evoked figures
% summaryFigure3 = figure;
% sgtitle('Stimulus-evoked cortical MUA [300-3000 Hz]  repsonses')
% %% LH short whisks
% ax1 = subplot(3,2,1);
% % C57BL6Js
% p1 = plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
% hold on
% plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.meanCortMUA + data.C57BL6J.Contra.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.meanCortMUA - data.C57BL6J.Contra.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% % Blank-SAP
% p2 = plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
% plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.meanCortMUA + data.Blank_SAP.Contra.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.meanCortMUA - data.Blank_SAP.Contra.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% % SSP-SAP
% p3 = plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
% plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.meanCortMUA + data.SSP_SAP.Contra.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
% plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.meanCortMUA - data.SSP_SAP.Contra.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
% title('LH (UnRx) MUA [300-3000 Hz] - Contra Stim')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% RH short whisks
% ax2 = subplot(3,2,2);
% % C57BL6Js
% plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
% hold on
% plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.meanCortMUA + data.C57BL6J.Contra.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.meanCortMUA - data.C57BL6J.Contra.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% % Blank-SAP
% plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
% plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.meanCortMUA + data.Blank_SAP.Contra.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.meanCortMUA - data.Blank_SAP.Contra.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% % SSP-SAP
% plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
% plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.meanCortMUA + data.SSP_SAP.Contra.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
% plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.meanCortMUA - data.SSP_SAP.Contra.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
% title('RH (Rx) MUA [300-3000 Hz] - Contra Stim')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% LH intermediate whisks
% ax3 = subplot(3,2,3);
% % C57BL6Js
% plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
% hold on
% plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.meanCortMUA + data.C57BL6J.Ipsi.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.meanCortMUA - data.C57BL6J.Ipsi.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% % Blank-SAP
% plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
% plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.meanCortMUA + data.Blank_SAP.Ipsi.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.meanCortMUA - data.Blank_SAP.Ipsi.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% % SSP-SAP
% plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
% plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.meanCortMUA + data.SSP_SAP.Ipsi.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
% plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.meanCortMUA - data.SSP_SAP.Ipsi.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
% title('LH (UnRx) MUA [300-3000 Hz] - Ipsi Stim')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% RH intermediate whisks
% ax4 = subplot(3,2,4);
% % C57BL6Js
% plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
% hold on
% plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.meanCortMUA + data.C57BL6J.Ipsi.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.meanCortMUA - data.C57BL6J.Ipsi.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% % Blank-SAP
% plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
% plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.meanCortMUA + data.Blank_SAP.Ipsi.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.meanCortMUA - data.Blank_SAP.Ipsi.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% % SSP-SAP
% plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
% plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.meanCortMUA + data.SSP_SAP.Ipsi.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
% plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.meanCortMUA - data.SSP_SAP.Ipsi.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
% title('RH (Rx) MUA [300-3000 Hz] - Ipsi Stim')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% LH long whisks
% ax5 = subplot(3,2,5);
% % C57BL6Js
% plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
% hold on
% plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.meanCortMUA + data.C57BL6J.Auditory.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.meanCortMUA - data.C57BL6J.Auditory.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% % Blank-SAP
% plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
% plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.meanCortMUA + data.Blank_SAP.Auditory.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.meanCortMUA - data.Blank_SAP.Auditory.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% % SSP-SAP
% plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
% plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.meanCortMUA + data.SSP_SAP.Auditory.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
% plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.meanCortMUA - data.SSP_SAP.Auditory.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
% title('LH (UnRx) MUA [300-3000 Hz] - Auditory Stim')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% RH long whisks
% ax6 = subplot(3,2,6);
% % C57BL6Js
% plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
% hold on
% plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.meanCortMUA + data.C57BL6J.Auditory.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.meanCortMUA - data.C57BL6J.Auditory.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% % Blank-SAP
% plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
% plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.meanCortMUA + data.Blank_SAP.Auditory.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.meanCortMUA - data.Blank_SAP.Auditory.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% % SSP-SAP
% plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
% plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.meanCortMUA + data.SSP_SAP.Auditory.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
% plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.meanCortMUA - data.SSP_SAP.Auditory.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
% title('RH (Rx) MUA [300-3000 Hz] - Auditory Stim')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% figure characteristics
% linkaxes([ax1,ax2],'xy')
% linkaxes([ax3,ax4],'xy')
% linkaxes([ax5,ax6],'xy')
% %% save figure(s)
% if strcmp(saveFigs,'y') == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure3,[dirpath 'Stim_Evoked_MUA']);
%     set(summaryFigure3,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'Stim_Evoked_MUA'])
% end
% %% individual whisk-evoked figures
% summaryFigure4 = figure;
% sgtitle('Stimulus-evoked cortical MUA [300-3000 Hz] repsonses - individual animals')
% %% LH short whisks
% ax1 = subplot(3,2,1);
% % C57BL6Js
% for aa = 1:size(data.C57BL6J.Contra.adjLH.cortMUA,1)
%     plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.Blank_SAP.Contra.adjLH.cortMUA,1)
%     plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.SSP_SAP.Contra.adjLH.cortMUA,1)
%     plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% title('LH (UnRx) MUA [300-3000 Hz] - Contra Stim')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% set(gca,'box','off')
% xlim([-2,10])
% %% RH short whisks
% ax2 = subplot(3,2,2);
% % C57BL6Js
% for aa = 1:size(data.C57BL6J.Contra.adjRH.cortMUA,1)
%     plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.Blank_SAP.Contra.adjRH.cortMUA,1)
%     plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.SSP_SAP.Contra.adjRH.cortMUA,1)
%     plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% title('RH (Rx) MUA [300-3000 Hz] - Contra Stim')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% LH intermediate whisks
% ax3 = subplot(3,2,3);
% % C57BL6Js
% for aa = 1:size(data.C57BL6J.Ipsi.adjLH.cortMUA,1)
%     plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.Blank_SAP.Ipsi.adjLH.cortMUA,1)
%     plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.SSP_SAP.Ipsi.adjLH.cortMUA,1)
%     plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% title('LH (UnRx) MUA [300-3000 Hz] - Ipsi Stim')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% RH intermediate whisks
% ax4 = subplot(3,2,4);
% % C57BL6Js
% for aa = 1:size(data.C57BL6J.Ipsi.adjRH.cortMUA,1)
%     plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.Blank_SAP.Ipsi.adjRH.cortMUA,1)
%     plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.SSP_SAP.Ipsi.adjRH.cortMUA,1)
%     plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% title('RH (Rx) MUA [300-3000 Hz] - Ipsi Stim')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% LH long whisks
% ax5 = subplot(3,2,5);
% % C57BL6Js
% for aa = 1:size(data.C57BL6J.Auditory.adjLH.cortMUA,1)
%     plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.Blank_SAP.Auditory.adjLH.cortMUA,1)
%     plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.SSP_SAP.Auditory.adjLH.cortMUA,1)
%     plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% title('LH (UnRx) MUA [300-3000 Hz] - Auditory Stim')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% RH long whisks
% ax6 = subplot(3,2,6);
% % C57BL6Js
% for aa = 1:size(data.C57BL6J.Auditory.adjRH.cortMUA,1)
%     plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
%     hold on
% end
% % Blank-SAP
% for aa = 1:size(data.Blank_SAP.Auditory.adjRH.cortMUA,1)
%     plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
%     hold on
% end
% % SSP-SAP
% for aa = 1:size(data.SSP_SAP.Auditory.adjRH.cortMUA,1)
%     plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
%     hold on
% end
% title('RH (Rx) MUA [300-3000 Hz] - Auditory Stim')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% set(gca,'box','off')
% xlim([-2,10])
% %% figure characteristics
% linkaxes([ax1,ax2],'xy')
% linkaxes([ax3,ax4],'xy')
% linkaxes([ax5,ax6],'xy')
% %% save figure(s)
% if strcmp(saveFigs,'y') == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure4,[dirpath 'indStim_Evoked_MUA']);
%     set(summaryFigure4,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'indStim_Evoked_MUA'])
% end
% %% average neural responses
% summaryFigure5 = figure;
% sgtitle('Stimulus-evoked cortical neural (LFP) repsonses')
% %% stim cortical LFP
% subplot(3,6,1);
% imagesc(data.C57BL6J.Contra.adjLH.meanCortT,data.C57BL6J.Contra.adjLH.meanCortF,data.C57BL6J.Contra.adjLH.meanCortS)
% title('C57BL6J Contra Stim LH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% stim cortical LFP
% subplot(3,6,2);
% imagesc(data.C57BL6J.Contra.adjRH.meanCortT,data.C57BL6J.Contra.adjRH.meanCortF,data.C57BL6J.Contra.adjRH.meanCortS)
% title('C57BL6J Contra Stim RH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% stim cortical LFP
% subplot(3,6,3);
% imagesc(data.SSP_SAP.Contra.adjLH.meanCortT,data.SSP_SAP.Contra.adjLH.meanCortF,data.SSP_SAP.Contra.adjLH.meanCortS)
% title('SSP-SAP Contra Stim LH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% stim cortical LFP
% subplot(3,6,4);
% imagesc(data.SSP_SAP.Contra.adjRH.meanCortT,data.SSP_SAP.Contra.adjRH.meanCortF,data.SSP_SAP.Contra.adjRH.meanCortS)
% title('SSP-SAP Contra Stim RH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% stim cortical LFP
% subplot(3,6,5);
% imagesc(data.Blank_SAP.Contra.adjLH.meanCortT,data.Blank_SAP.Contra.adjLH.meanCortF,data.Blank_SAP.Contra.adjLH.meanCortS)
% title('Blank-SAP Contra Stim LH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% stim cortical LFP
% subplot(3,6,6);
% imagesc(data.Blank_SAP.Contra.adjRH.meanCortT,data.Blank_SAP.Contra.adjRH.meanCortF,data.Blank_SAP.Contra.adjRH.meanCortS)
% title('Blank-SAP Contra Stim RH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% Ipsi cortical LFP
% subplot(3,6,7);
% imagesc(data.C57BL6J.Ipsi.adjLH.meanCortT,data.C57BL6J.Ipsi.adjLH.meanCortF,data.C57BL6J.Ipsi.adjLH.meanCortS)
% title('C57BL6J Ipsi Stim LH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% Ipsi cortical LFP
% subplot(3,6,8);
% imagesc(data.C57BL6J.Ipsi.adjRH.meanCortT,data.C57BL6J.Ipsi.adjRH.meanCortF,data.C57BL6J.Ipsi.adjRH.meanCortS)
% title('C57BL6J Ipsi Stim RH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% Ipsi cortical LFP
% subplot(3,6,9);
% imagesc(data.SSP_SAP.Ipsi.adjLH.meanCortT,data.SSP_SAP.Ipsi.adjLH.meanCortF,data.SSP_SAP.Ipsi.adjLH.meanCortS)
% title('SSP-SAP Ipsi Stim LH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% Ipsi cortical LFP
% subplot(3,6,10);
% imagesc(data.SSP_SAP.Ipsi.adjRH.meanCortT,data.SSP_SAP.Ipsi.adjRH.meanCortF,data.SSP_SAP.Ipsi.adjRH.meanCortS)
% title('SSP-SAP Ipsi Stim RH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% Ipsi cortical LFP
% subplot(3,6,11);
% imagesc(data.Blank_SAP.Ipsi.adjLH.meanCortT,data.Blank_SAP.Ipsi.adjLH.meanCortF,data.Blank_SAP.Ipsi.adjLH.meanCortS)
% title('Blank-SAP Ipsi Stim LH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% Ipsi cortical LFP
% subplot(3,6,12);
% imagesc(data.Blank_SAP.Ipsi.adjRH.meanCortT,data.Blank_SAP.Ipsi.adjRH.meanCortF,data.Blank_SAP.Ipsi.adjRH.meanCortS)
% title('Blank-SAP Ipsi Stim RH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% Auditory cortical LFP
% subplot(3,6,13);
% imagesc(data.C57BL6J.Auditory.adjLH.meanCortT,data.C57BL6J.Auditory.adjLH.meanCortF,data.C57BL6J.Auditory.adjLH.meanCortS)
% title('C57BL6J Auditory Stim LH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% stim cortical LFP
% subplot(3,6,14);
% imagesc(data.C57BL6J.Auditory.adjRH.meanCortT,data.C57BL6J.Auditory.adjRH.meanCortF,data.C57BL6J.Auditory.adjRH.meanCortS)
% title('C57BL6J Auditory Stim RH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% stim cortical LFP
% subplot(3,6,15);
% imagesc(data.SSP_SAP.Auditory.adjLH.meanCortT,data.SSP_SAP.Auditory.adjLH.meanCortF,data.SSP_SAP.Auditory.adjLH.meanCortS)
% title('SSP-SAP Auditory Stim LH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% stim cortical LFP
% subplot(3,6,16);
% imagesc(data.SSP_SAP.Auditory.adjRH.meanCortT,data.SSP_SAP.Auditory.adjRH.meanCortF,data.SSP_SAP.Auditory.adjRH.meanCortS)
% title('SSP-SAP Auditory Stim RH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% stim cortical LFP
% subplot(3,6,17);
% imagesc(data.Blank_SAP.Auditory.adjLH.meanCortT,data.Blank_SAP.Auditory.adjLH.meanCortF,data.Blank_SAP.Auditory.adjLH.meanCortS)
% title('Blank-SAP Auditory Stim LH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% stim cortical LFP
% subplot(3,6,18);
% imagesc(data.Blank_SAP.Auditory.adjRH.meanCortT,data.Blank_SAP.Auditory.adjRH.meanCortF,data.Blank_SAP.Auditory.adjRH.meanCortS)
% title('Blank-SAP Auditory Stim RH')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% % caxis([-50,100])
% set(gca,'Ticklength',[0,0])
% axis xy
% set(gca,'box','off')
% %% save figure(s)
% if strcmp(saveFigs,'y') == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure5,[dirpath 'Stim_Evoked_LFP']);
%     set(summaryFigure5,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'Stim_Evoked_LFP'])
% end

end
