function [AnalysisResults] = PowerSpec2_Saporin(rootFolder,saveFigs,delim,AnalysisResults)
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
behavFields = {'Alert','Asleep','All'};
dataTypes = {'gammaBandPower'};
%% power spectra during different behaviors
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
                data.(treatment).(behavField).(dataType).LH.S = [];
                data.(treatment).(behavField).(dataType).LH.f = [];
                data.(treatment).(behavField).(dataType).RH.S = [];
                data.(treatment).(behavField).(dataType).RH.f = [];
            end
            data.(treatment).(behavField).(dataType).LH.S = cat(2,data.(treatment).(behavField).(dataType).LH.S,AnalysisResults.(animalID).PowerSpectra2.(behavField).(dataType).LH.S);
            data.(treatment).(behavField).(dataType).LH.f = cat(1,data.(treatment).(behavField).(dataType).LH.f,AnalysisResults.(animalID).PowerSpectra2.(behavField).(dataType).LH.f);
            data.(treatment).(behavField).(dataType).RH.S = cat(2,data.(treatment).(behavField).(dataType).RH.S,AnalysisResults.(animalID).PowerSpectra2.(behavField).(dataType).RH.S);
            data.(treatment).(behavField).(dataType).RH.f = cat(1,data.(treatment).(behavField).(dataType).RH.f,AnalysisResults.(animalID).PowerSpectra2.(behavField).(dataType).RH.f);
        end
    end
end
%% find the peak of the resting PSD for each animal/hemisphere
for aa = 1:length(treatments)
    treatment = treatments{1,aa};
    for cc = 1:length(dataTypes)
        dataType = dataTypes{1,cc};
        for ee = 1:size(data.(treatment).Alert.(dataType).LH.S,2)
            data.(treatment).baseline.(dataType).LH(ee,1) = max(data.(treatment).Alert.(dataType).LH.S(:,ee));
            data.(treatment).baseline.(dataType).RH(ee,1) = max(data.(treatment).Alert.(dataType).RH.S(:,ee));
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
                data.(treatment).(behavField).(dataType).LH.normS(:,ee) = (data.(treatment).(behavField).(dataType).LH.S(:,ee));%*(1/(data.(treatment).baseline.(dataType).LH(ee,1)));
                data.(treatment).(behavField).(dataType).RH.normS(:,ee) = (data.(treatment).(behavField).(dataType).RH.S(:,ee));%*(1/(data.(treatment).baseline.(dataType).RH(ee,1)));
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
%% average gamma-band power
summaryFigure1 = figure;
sgtitle('LFP Power Spectra [1-100 Hz]')
%% LH power spectra of gamma-band power during Alert
ax1 = subplot(3,2,1);
L1 = loglog(data.C57BL6J.Alert.gammaBandPower.LH.meanCortf,data.C57BL6J.Alert.gammaBandPower.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
L2 = loglog(data.Blank_SAP.Alert.gammaBandPower.LH.meanCortf,data.Blank_SAP.Alert.gammaBandPower.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
L3 = loglog(data.SSP_SAP.Alert.gammaBandPower.LH.meanCortf,data.SSP_SAP.Alert.gammaBandPower.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Alert] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1,100])
legend([L1,L2,L3],'C57BL6J','SSP-SAP','Blank-SAP')
set(gca,'box','off')
%% RH power spectra of gamma-band power during Alert
ax2 = subplot(3,2,2);
loglog(data.C57BL6J.Alert.gammaBandPower.RH.meanCortf,data.C57BL6J.Alert.gammaBandPower.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Alert.gammaBandPower.RH.meanCortf,data.Blank_SAP.Alert.gammaBandPower.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Alert.gammaBandPower.RH.meanCortf,data.SSP_SAP.Alert.gammaBandPower.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Alert] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1,100])
set(gca,'box','off')
%% LH power spectra of gamma-band power during Asleep
ax3 = subplot(3,2,3);
loglog(data.C57BL6J.Asleep.gammaBandPower.LH.meanCortf,data.C57BL6J.Asleep.gammaBandPower.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Asleep.gammaBandPower.LH.meanCortf,data.Blank_SAP.Asleep.gammaBandPower.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Asleep.gammaBandPower.LH.meanCortf,data.SSP_SAP.Asleep.gammaBandPower.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Asleep] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1,100])
set(gca,'box','off')
%% RH power spectra of gamma-band power during Asleep
ax4 = subplot(3,2,4);
loglog(data.C57BL6J.Asleep.gammaBandPower.RH.meanCortf,data.C57BL6J.Asleep.gammaBandPower.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.Asleep.gammaBandPower.RH.meanCortf,data.Blank_SAP.Asleep.gammaBandPower.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.Asleep.gammaBandPower.RH.meanCortf,data.SSP_SAP.Asleep.gammaBandPower.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[Asleep] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1,100])
set(gca,'box','off')
%% LH power spectra of gamma-band power during All data
ax5 = subplot(3,2,5);
loglog(data.C57BL6J.All.gammaBandPower.LH.meanCortf,data.C57BL6J.All.gammaBandPower.LH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.All.gammaBandPower.LH.meanCortf,data.Blank_SAP.All.gammaBandPower.LH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.All.gammaBandPower.LH.meanCortf,data.SSP_SAP.All.gammaBandPower.LH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[All] LH (UnRx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1,100])
set(gca,'box','off')
%% RH power spectra of gamma-band power during All data
ax6 = subplot(3,2,6);
loglog(data.C57BL6J.All.gammaBandPower.RH.meanCortf,data.C57BL6J.All.gammaBandPower.RH.meanCortS,'color',colors('sapphire'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.All.gammaBandPower.RH.meanCortf,data.Blank_SAP.All.gammaBandPower.RH.meanCortS,'color',colors('north texas green'),'LineWidth',2);
loglog(data.SSP_SAP.All.gammaBandPower.RH.meanCortf,data.SSP_SAP.All.gammaBandPower.RH.meanCortS,'color',colors('electric purple'),'LineWidth',2);
title('[All] RH (Rx)')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([1,100])
set(gca,'box','off')
%% figure characteristics
linkaxes([ax1,ax2],'y')
linkaxes([ax3,ax4],'y')
linkaxes([ax5,ax6],'y')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'PowerSpec_LFP']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'PowerSpec_LFP'])
end

end
