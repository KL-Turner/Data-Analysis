function [AnalysisResults] = FigS15_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel S15 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
IOSanimalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
%% mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,aa};
    for bb = 1:length(transitions)
        transition = transitions{1,bb};
        data.(animalID).(transition) = [];
        uniqueFileDates = AnalysisResults.(animalID).Transitions.(transition).fileDates;
        for cc = 1:length(uniqueFileDates)
            data.(animalID).(transition).(uniqueFileDates{cc,1}).LH_HbT = [];
            data.(animalID).(transition).(uniqueFileDates{cc,1}).RH_HbT = [];
        end
        for dd = 1:length(AnalysisResults.(animalID).Transitions.(transition).indFileDate)
            strDay = AnalysisResults.(animalID).Transitions.(transition).indFileDate{dd,1};
            data.(animalID).(transition).(strDay).LH_HbT = cat(1,data.(animalID).(transition).(strDay).LH_HbT,AnalysisResults.(animalID).Transitions.(transition).LH_HbT(dd,:));
            data.(animalID).(transition).(strDay).RH_HbT = cat(1,data.(animalID).(transition).(strDay).RH_HbT,AnalysisResults.(animalID).Transitions.(transition).RH_HbT(dd,:));
        end
    end
end
% put together L/R for each day
for ee = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,ee};
    for ff = 1:length(transitions)
        transition = transitions{1,ff};
        uniqueFileDates = fieldnames(data.(animalID).(transition));
        for gg = 1:length(uniqueFileDates)
            strDay = uniqueFileDates{gg,1};
            procData.(animalID).(transition).(strDay).HbT = cat(1,data.(animalID).(transition).(strDay).LH_HbT,data.(animalID).(transition).(strDay).RH_HbT);
        end
    end
end
% take average for each animal's behavioral transition per day
for ee = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,ee};
    for ff = 1:length(transitions)
        transition = transitions{1,ff};
        uniqueFileDates = fieldnames(data.(animalID).(transition));
        for gg = 1:length(uniqueFileDates)
            strDay = uniqueFileDates{gg,1};
            if isempty(procData.(animalID).(transition).(strDay).HbT) == false
                finData.(transition){gg,1}(ee,:) = mean(procData.(animalID).(transition).(strDay).HbT,1);
            else
                finData.(transition){gg,1}(ee,:) = NaN(1,1800);
            end
        end
    end
end
% patch Animal T110 who only had 5 days
for hh = 1:length(transitions)
    transition = transitions{1,hh};
    finData.(transition){6,1}(8,:) = NaN(1,1800);
end
% take average across animals
for ii = 1:length(transitions)
    transition = transitions{1,ii};
    for jj = 1:6
        meanData.(transition){jj,1} = nanmean(finData.(transition){jj,1},1);
        % check nans
        ll = 1;
        for kk = 1:size(finData.(transition){jj,1},1)
            if isnan(finData.(transition){jj,1}(kk,1)) == false
                nCount.(transition){jj,1} = ll; %#ok<STRNU>
                ll = ll + 1;
            end
        end
    end
end
T1 = -30:(1/30):30;
T1 = T1(1:end - 1);
%% Fig. S15
summaryFigure = figure('Name','FigS15 (a-d)');
sgtitle('Figure S15 - Turner et al. 2020')
%% [S15a] Awake to NREM
ax1 = subplot(2,2,1);
for kk = 1:6
    p(kk) = plot(T1,meanData.AWAKEtoNREM{kk,1},'-','LineWidth',2); %#ok<*AGROW>
    hold on
end
title('[S15a] Awake to NREM transition')
xlabel('Time (s)')
ylabel('\Delta[HbT] (\muM)')
legend([p(1),p(2),p(3),p(4),p(5),p(6)],'Day 1','Day 2','Day 3','Day 4','Day 5','Day 6')
xlim([-30,30])
ylim([-5,45])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [S15b] NREM to Awake
ax2 = subplot(2,2,2);
for kk = 1:6
    plot(T1,meanData.NREMtoAWAKE{kk,1},'-','LineWidth',2);
    hold on
end
title('[S15b] NREM to Awake transition')
xlabel('Time (s)')
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
ylim([-5,45])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [S15c] NREM to REM
ax3 = subplot(2,2,3);
for kk = 1:6
    plot(T1,meanData.NREMtoREM{kk,1},'-','LineWidth',2);
    hold on
end
title('[S15c] NREM to REM transition')
xlabel('Time (s)')
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
ylim([35,80])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [S15d] REM to Awake
ax4 = subplot(2,2,4);
for kk = 1:6
    plot(T1,meanData.REMtoAWAKE{kk,1},'-','LineWidth',2);
    hold on
end
title('[S15d] REM to Awake transition')
xlabel('Time (s)')
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
ylim([0,100])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'FigS15']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'FigS15'])
end

end
