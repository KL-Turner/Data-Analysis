function [AnalysisResults] = NeuralHemoCoherence_Saporin(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: 
%________________________________________________________________________________________________________________________

%% set-up and process data
animalIDs = {'T141','T155','T156','T157','T142','T144','T159','T172','T150','T165','T166','T177','T179','T186','T187','T188','T189'};
C57BL6J_IDs = {'T141','T155','T156','T157','T186','T187','T188','T189'};
SSP_SAP_IDs = {'T142','T144','T159','T172'};
Blank_SAP_IDs = {'T150','T165','T166','T177','T179'};
treatments = {'C57BL6J','SSP_SAP','Blank_SAP'};
behavFields = {'Rest','NREM','REM','Awake','Sleep','All'};
dataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
%% average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    % recognize treatment based on animal group
    if ismember(animalID,C57BL6J_IDs) == true
        treatment = 'C57BL6J';
    elseif ismember(animalIDs{1,aa},SSP_SAP_IDs) == true
        treatment = 'SSP_SAP';
    elseif ismember(animalIDs{1,aa},Blank_SAP_IDs) == true
        treatment = 'Blank_SAP';
    end
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % pre-allocate necessary variable fields
            data.(treatment).(behavField).dummCheck = 1;
            if isfield(data.(treatment).(behavField),dataType) == false
                data.(treatment).(behavField).(dataType).LH.C = [];
                data.(treatment).(behavField).(dataType).LH.f = [];
                data.(treatment).(behavField).(dataType).LH.confC = [];
                data.(treatment).(behavField).(dataType).RH.C = [];
                data.(treatment).(behavField).(dataType).RH.f = [];
                data.(treatment).(behavField).(dataType).RH.confC = [];
            end
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C) == false
                % concatenate C/f for existing data - exclude any empty sets
                data.(treatment).(behavField).(dataType).LH.C = cat(2,data.(treatment).(behavField).(dataType).LH.C,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C);
                data.(treatment).(behavField).(dataType).LH.f = cat(1,data.(treatment).(behavField).(dataType).LH.f,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.f);
                data.(treatment).(behavField).(dataType).LH.confC = cat(1,data.(treatment).(behavField).(dataType).LH.confC,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.confC);
                data.(treatment).(behavField).(dataType).RH.C = cat(2,data.(treatment).(behavField).(dataType).RH.C,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.C);
                data.(treatment).(behavField).(dataType).RH.f = cat(1,data.(treatment).(behavField).(dataType).RH.f,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.f);
                data.(treatment).(behavField).(dataType).RH.confC = cat(1,data.(treatment).(behavField).(dataType).RH.confC,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.confC);
            end
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
%% average HbT coherence
summaryFigure1 = figure;
sgtitle('Gamma-band - \DeltaHbT Coherence')
% %% coherence^2 between bilateral HbT during rest
% subplot(3,4,1);
% p1 = semilogx(data.C57BL6J.Rest.gammaBandPower.LH.meanf,data.C57BL6J.Rest.gammaBandPower.LH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
% hold on
% p2 = semilogx(data.Blank_SAP.Rest.gammaBandPower.LH.meanf,data.Blank_SAP.Rest.gammaBandPower.LH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
% p3 = semilogx(data.SSP_SAP.Rest.gammaBandPower.LH.meanf,data.SSP_SAP.Rest.gammaBandPower.LH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Rest] Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
% xlim([1/10,0.5])
% % ylim([0,1])
% set(gca,'box','off')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% %% coherence^2 between bilateral HbT during rest
% subplot(3,4,2);
% p1 = semilogx(data.C57BL6J.Rest.gammaBandPower.RH.meanf,data.C57BL6J.Rest.gammaBandPower.RH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
% hold on
% p2 = semilogx(data.Blank_SAP.Rest.gammaBandPower.RH.meanf,data.Blank_SAP.Rest.gammaBandPower.RH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
% p3 = semilogx(data.SSP_SAP.Rest.gammaBandPower.RH.meanf,data.SSP_SAP.Rest.gammaBandPower.RH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[Rest] Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
% xlim([1/10,0.5])
% % ylim([0,1])
% set(gca,'box','off')
% legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
% %% coherence^2 between bilateral HbT during NREM
% subplot(3,4,3);
% semilogx(data.C57BL6J.NREM.gammaBandPower.LH.meanf,data.C57BL6J.NREM.gammaBandPower.LH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Blank_SAP.NREM.gammaBandPower.LH.meanf,data.Blank_SAP.NREM.gammaBandPower.LH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.SSP_SAP.NREM.gammaBandPower.LH.meanf,data.SSP_SAP.NREM.gammaBandPower.LH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[NREM] Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
% xlim([1/30,0.5])
% % ylim([0,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral HbT during NREM
% subplot(3,4,4);
% semilogx(data.C57BL6J.NREM.gammaBandPower.RH.meanf,data.C57BL6J.NREM.gammaBandPower.RH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Blank_SAP.NREM.gammaBandPower.RH.meanf,data.Blank_SAP.NREM.gammaBandPower.RH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.SSP_SAP.NREM.gammaBandPower.RH.meanf,data.SSP_SAP.NREM.gammaBandPower.RH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[NREM] Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
% xlim([1/30,0.5])
% % ylim([0,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral HbT during REM
% subplot(3,4,5);
% semilogx(data.C57BL6J.REM.gammaBandPower.LH.meanf,data.C57BL6J.REM.gammaBandPower.LH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Blank_SAP.REM.gammaBandPower.LH.meanf,data.Blank_SAP.REM.gammaBandPower.LH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.SSP_SAP.REM.gammaBandPower.LH.meanf,data.SSP_SAP.REM.gammaBandPower.LH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[REM] Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
% xlim([1/60,0.5])
% % ylim([0,1])
% set(gca,'box','off')
% %% coherence^2 between bilateral HbT during REM
% subplot(3,4,6);
% semilogx(data.C57BL6J.REM.gammaBandPower.RH.meanf,data.C57BL6J.REM.gammaBandPower.RH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
% hold on
% semilogx(data.Blank_SAP.REM.gammaBandPower.RH.meanf,data.Blank_SAP.REM.gammaBandPower.RH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
% semilogx(data.SSP_SAP.REM.gammaBandPower.RH.meanf,data.SSP_SAP.REM.gammaBandPower.RH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[REM] Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
% xlim([1/60,0.5])
% % ylim([0,1])
% set(gca,'box','off')
%% coherence^2 between bilateral HbT during Awake
subplot(3,2,1);
s1 = semilogx(data.C57BL6J.Awake.gammaBandPower.LH.meanf,data.C57BL6J.Awake.gammaBandPower.LH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
s2 = semilogx(data.Blank_SAP.Awake.gammaBandPower.LH.meanf,data.Blank_SAP.Awake.gammaBandPower.LH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
s3 = semilogx(data.SSP_SAP.Awake.gammaBandPower.LH.meanf,data.SSP_SAP.Awake.gammaBandPower.LH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] LH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
% ylim([0,1])
legend([s1,s2,s3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Awake
subplot(3,2,2);
semilogx(data.C57BL6J.Awake.gammaBandPower.RH.meanf,data.C57BL6J.Awake.gammaBandPower.RH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.Awake.gammaBandPower.RH.meanf,data.Blank_SAP.Awake.gammaBandPower.RH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.Awake.gammaBandPower.RH.meanf,data.SSP_SAP.Awake.gammaBandPower.RH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Alert] RH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(3,2,3);
semilogx(data.C57BL6J.Sleep.gammaBandPower.LH.meanf,data.C57BL6J.Sleep.gammaBandPower.LH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.Sleep.gammaBandPower.LH.meanf,data.Blank_SAP.Sleep.gammaBandPower.LH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.Sleep.gammaBandPower.LH.meanf,data.SSP_SAP.Sleep.gammaBandPower.LH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] LH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during Sleep
subplot(3,2,4);
semilogx(data.C57BL6J.Sleep.gammaBandPower.RH.meanf,data.C57BL6J.Sleep.gammaBandPower.RH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.Sleep.gammaBandPower.RH.meanf,data.Blank_SAP.Sleep.gammaBandPower.RH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.Sleep.gammaBandPower.RH.meanf,data.SSP_SAP.Sleep.gammaBandPower.RH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[Asleep] RH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(3,2,5);
semilogx(data.C57BL6J.All.gammaBandPower.LH.meanf,data.C57BL6J.All.gammaBandPower.LH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.All.gammaBandPower.LH.meanf,data.Blank_SAP.All.gammaBandPower.LH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.All.gammaBandPower.LH.meanf,data.SSP_SAP.All.gammaBandPower.LH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] LH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% coherence^2 between bilateral HbT during All data
subplot(3,2,6);
semilogx(data.C57BL6J.All.gammaBandPower.RH.meanf,data.C57BL6J.All.gammaBandPower.RH.meanC.^2,'color',colors('sapphire'),'LineWidth',2);
hold on
semilogx(data.Blank_SAP.All.gammaBandPower.RH.meanf,data.Blank_SAP.All.gammaBandPower.RH.meanC.^2,'color',colors('north texas green'),'LineWidth',2);
semilogx(data.SSP_SAP.All.gammaBandPower.RH.meanf,data.SSP_SAP.All.gammaBandPower.RH.meanC.^2,'color',colors('electric purple'),'LineWidth',2);
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[All] RH Neural-hemo coherence^2','Gamma-band \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
% ylim([0,1])
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'NeuralHemoCoherence_Gamma']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'NeuralHemoCoherence_Gamma'])
end

end
