function [AnalysisResults] = FigS9_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel S9 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
colorRfcAwake = [(0/256),(64/256),(64/256)];
colorRfcNREM = [(0/256),(174/256),(239/256)];
colorRfcREM = [(190/256),(30/256),(45/256)];
% colorRest = [(0/256),(166/256),(81/256)];
% colorWhisk = [(31/256),(120/256),(179/256)];
% colorStim = [(255/256),(28/256),(206/256)];
% colorNREM = [(191/256),(0/256),(255/256)];
% colorREM = [(254/256),(139/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
% mean HbT and heart rate comparison between behaviors
behavFields = {'Awake','NREM','REM'};
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        data.BehavioralDistributions.(behavField).EMG(aa,1) = mean(AnalysisResults.(animalID).BehaviorDistributions.(behavField).EMG);
        data.BehavioralDistributions.(behavField).Whisk(aa,1) = mean(AnalysisResults.(animalID).BehaviorDistributions.(behavField).Whisk);
        data.BehavioralDistributions.(behavField).Heart(aa,1) = mean(AnalysisResults.(animalID).BehaviorDistributions.(behavField).HR);
        data.BehavioralDistributions.(behavField).animalIDs{aa,1} = animalID;
        data.BehavioralDistributions.(behavField).behaviors{aa,1} = behavField;
    end
end
% take the mean and standard deviation of each set of signals
for cc = 1:length(behavFields)
    behavField = behavFields{1,cc};
    data.BehavioralDistributions.(behavField).meanEMG = mean(data.BehavioralDistributions.(behavField).EMG,1);
    data.BehavioralDistributions.(behavField).stdEMG = std(data.BehavioralDistributions.(behavField).EMG,0,1);
    data.BehavioralDistributions.(behavField).meanWhisk = mean(data.BehavioralDistributions.(behavField).Whisk,1);
    data.BehavioralDistributions.(behavField).stdWhisk = std(data.BehavioralDistributions.(behavField).Whisk,0,1);
    data.BehavioralDistributions.(behavField).meanHeart = mean(data.BehavioralDistributions.(behavField).Heart,1);
    data.BehavioralDistributions.(behavField).stdHeart = std(data.BehavioralDistributions.(behavField).Heart,0,1);
end
%% statistics - generalized linear mixed-effects model
% EMG
EMGtableSize = cat(1,data.BehavioralDistributions.Awake.EMG,data.BehavioralDistributions.NREM.EMG,data.BehavioralDistributions.REM.EMG);
EMGTable = table('Size',[size(EMGtableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','EMG','Behavior'});
EMGTable.Mouse = cat(1,data.BehavioralDistributions.Awake.animalIDs,data.BehavioralDistributions.NREM.animalIDs,data.BehavioralDistributions.REM.animalIDs);
EMGTable.EMG = cat(1,data.BehavioralDistributions.Awake.EMG,data.BehavioralDistributions.NREM.EMG,data.BehavioralDistributions.REM.EMG);
EMGTable.Behavior = cat(1,data.BehavioralDistributions.Awake.behaviors,data.BehavioralDistributions.NREM.behaviors,data.BehavioralDistributions.REM.behaviors);
EMGFitFormula = 'EMG ~ 1 + Behavior + (1|Mouse)';
EMGStats = fitglme(EMGTable,EMGFitFormula);
% whisker
WhisktableSize = cat(1,data.BehavioralDistributions.Awake.Whisk,data.BehavioralDistributions.NREM.Whisk,data.BehavioralDistributions.REM.Whisk);
WhiskTable = table('Size',[size(WhisktableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','Whisk','Behavior'});
WhiskTable.Mouse = cat(1,data.BehavioralDistributions.Awake.animalIDs,data.BehavioralDistributions.NREM.animalIDs,data.BehavioralDistributions.REM.animalIDs);
WhiskTable.Whisk = cat(1,data.BehavioralDistributions.Awake.Whisk,data.BehavioralDistributions.NREM.Whisk,data.BehavioralDistributions.REM.Whisk);
WhiskTable.Behavior = cat(1,data.BehavioralDistributions.Awake.behaviors,data.BehavioralDistributions.NREM.behaviors,data.BehavioralDistributions.REM.behaviors);
WhiskFitFormula = 'Whisk ~ 1 + Behavior + (1|Mouse)';
WhiskStats = fitglme(WhiskTable,WhiskFitFormula);
% heart rate
HearttableSize = cat(1,data.BehavioralDistributions.Awake.Heart,data.BehavioralDistributions.NREM.Heart,data.BehavioralDistributions.REM.Heart);
HeartTable = table('Size',[size(HearttableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','Heart','Behavior'});
HeartTable.Mouse = cat(1,data.BehavioralDistributions.Awake.animalIDs,data.BehavioralDistributions.NREM.animalIDs,data.BehavioralDistributions.REM.animalIDs);
HeartTable.Heart = cat(1,data.BehavioralDistributions.Awake.Heart,data.BehavioralDistributions.NREM.Heart,data.BehavioralDistributions.REM.Heart);
HeartTable.Behavior = cat(1,data.BehavioralDistributions.Awake.behaviors,data.BehavioralDistributions.NREM.behaviors,data.BehavioralDistributions.REM.behaviors);
HeartFitFormula = 'Heart ~ 1 + Behavior + (1|Mouse)';
HeartStats = fitglme(HeartTable,HeartFitFormula);
%% Fig. S9
summaryFigure = figure('Name','FigS9 (a-c)');
sgtitle('Figure S9 - Turner et al. 2020')
%% [S9a] mean EMG power during different behaviors
ax1 = subplot(1,3,1);
xInds = ones(1,length(animalIDs));
s1 = scatter(xInds*1,data.BehavioralDistributions.Awake.EMG,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcAwake,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.BehavioralDistributions.Awake.meanEMG,data.BehavioralDistributions.Awake.stdEMG,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(xInds*2,data.BehavioralDistributions.NREM.EMG,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcNREM,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.BehavioralDistributions.NREM.meanEMG,data.BehavioralDistributions.NREM.stdEMG,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(xInds*3,data.BehavioralDistributions.REM.EMG,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.BehavioralDistributions.REM.meanEMG,data.BehavioralDistributions.REM.stdEMG,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title('[S9a] Mean EMG power')
ylabel('EMG power (a.u.)')
legend([s1,s2,s3],'rfc-Awake','rfc-NREM','rfc-REM')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
ylim([-2,0.5])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [S9b] mean Whisker variance during different behaviors
ax2 = subplot(1,3,2);
xInds = ones(1,length(animalIDs));
scatter(xInds*1,data.BehavioralDistributions.Awake.Whisk,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcAwake,'jitter','on','jitterAmount',0.25);
hold on
e4 = errorbar(1,data.BehavioralDistributions.Awake.meanWhisk,data.BehavioralDistributions.Awake.stdWhisk,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*2,data.BehavioralDistributions.NREM.Whisk,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcNREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(2,data.BehavioralDistributions.NREM.meanWhisk,data.BehavioralDistributions.NREM.stdWhisk,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*3,data.BehavioralDistributions.REM.Whisk,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcREM,'jitter','on','jitterAmount',0.25);
e6 = errorbar(3,data.BehavioralDistributions.REM.meanWhisk,data.BehavioralDistributions.REM.stdWhisk,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title('[S9b] Mean whisker variance')
ylabel('Whisker variance (deg^2)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
ylim([0,45])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [S9c] mean heart rate during different behaviors
ax3 = subplot(1,3,3);
xInds = ones(1,length(animalIDs));
scatter(xInds*1,data.BehavioralDistributions.Awake.Heart,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcAwake,'jitter','on','jitterAmount',0.25);
hold on
e7 = errorbar(1,data.BehavioralDistributions.Awake.meanHeart,data.BehavioralDistributions.Awake.stdHeart,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*2,data.BehavioralDistributions.NREM.Heart,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(2,data.BehavioralDistributions.NREM.meanHeart,data.BehavioralDistributions.NREM.stdHeart,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*3,data.BehavioralDistributions.REM.Heart,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(3,data.BehavioralDistributions.REM.meanHeart,data.BehavioralDistributions.REM.stdHeart,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
title('[S9c] Mean heart rate')
ylabel('Heart rate (Hz)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
ylim([5,9])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'FigS9']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'FigS9'])
    %% statistical diary
    diaryFile = [dirpath 'FigS9_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % EMG statistical diary
    disp('======================================================================================================================')
    disp('[S9a] Generalized linear mixed-effects model statistics for mean EMG during Not Asleep, NREM, and REM')
    disp('======================================================================================================================')
    disp(EMGStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake EMG pwr: ' num2str(round(data.BehavioralDistributions.Awake.meanEMG,1)) ' +/- ' num2str(round(data.BehavioralDistributions.Awake.stdEMG,1))]); disp(' ')
    disp(['NREM EMG pwr: ' num2str(round(data.BehavioralDistributions.NREM.meanEMG,1)) ' +/- ' num2str(round(data.BehavioralDistributions.NREM.stdEMG,1))]); disp(' ')
    disp(['REM EMG pwr: ' num2str(round(data.BehavioralDistributions.REM.meanEMG,1)) ' +/- ' num2str(round(data.BehavioralDistributions.REM.stdEMG,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % whisker variance statistical diary
    disp('======================================================================================================================')
    disp('[S9b] Generalized linear mixed-effects model statistics for whisker angle variance during Not Asleep, NREM, and REM')
    disp('======================================================================================================================')
    disp(WhiskStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake whisk var: ' num2str(round(data.BehavioralDistributions.Awake.meanWhisk,1)) ' +/- ' num2str(round(data.BehavioralDistributions.Awake.stdWhisk,1))]); disp(' ')
    disp(['NREM whisk var: ' num2str(round(data.BehavioralDistributions.NREM.meanWhisk,1)) ' +/- ' num2str(round(data.BehavioralDistributions.NREM.stdWhisk,1))]); disp(' ')
    disp(['REM whisk var: ' num2str(round(data.BehavioralDistributions.REM.meanWhisk,1)) ' +/- ' num2str(round(data.BehavioralDistributions.REM.stdWhisk,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % heart rate statistical diary
    disp('======================================================================================================================')
    disp('[S9c] Generalized linear mixed-effects model statistics for mean heart rate during Not Asleep, NREM, and REM')
    disp('======================================================================================================================')
    disp(HeartStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake heart rate (Hz): ' num2str(round(data.BehavioralDistributions.Awake.meanHeart,1)) ' +/- ' num2str(round(data.BehavioralDistributions.Awake.stdHeart,1))]); disp(' ')
    disp(['NREM heart rate (Hz): ' num2str(round(data.BehavioralDistributions.NREM.meanHeart,1)) ' +/- ' num2str(round(data.BehavioralDistributions.NREM.stdHeart,1))]); disp(' ')
    disp(['REM heart rate (Hz): ' num2str(round(data.BehavioralDistributions.REM.meanHeart,1)) ' +/- ' num2str(round(data.BehavioralDistributions.REM.stdHeart,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
end

end
