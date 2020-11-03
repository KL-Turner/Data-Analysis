function [AnalysisResults] = Fig1_S4_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1_S4 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
% colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
colorRest = [(0/256),(166/256),(81/256)];
% colorWhisk = [(31/256),(120/256),(179/256)];
% colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
TwoP_animalIDs = {'T115','T116','T117','T118','T125','T126'};
behavFields = {'Rest','NREM','REM'};
%% HbT comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH)
            data.HbT.(animalID).(behavField).maxLH(cc,1) = max(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{cc,1});
            data.HbT.(animalID).(behavField).maxRH(cc,1) = max(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{cc,1});
            data.HbT.(animalID).(behavField).p2pLH(cc,1) = abs(max(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{cc,1})) + abs(min(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{cc,1}));
            data.HbT.(animalID).(behavField).p2pRH(cc,1) = abs(max(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{cc,1})) + abs(min(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{cc,1}));
        end
    end
end
% put data into arrays and prep for stats
for dd = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,dd};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.HbT,behavField) == false
            data.HbT.(behavField).catMax = [];
            data.HbT.(behavField).catP2P = [];
            data.HbT.(behavField).animalID = {};
            data.HbT.(behavField).behavior = {};
            data.HbT.(behavField).hemisphere = {};
        end
        data.HbT.(behavField).catMax = cat(1,data.HbT.(behavField).catMax,mean(data.HbT.(animalID).(behavField).maxLH),mean(data.HbT.(animalID).(behavField).maxRH));
        data.HbT.(behavField).catP2P = cat(1,data.HbT.(behavField).catP2P,mean(data.HbT.(animalID).(behavField).p2pLH),mean(data.HbT.(animalID).(behavField).p2pRH));
        data.HbT.(behavField).animalID = cat(1,data.HbT.(behavField).animalID,animalID,animalID);
        data.HbT.(behavField).behavior = cat(1,data.HbT.(behavField).behavior,behavField,behavField);
        data.HbT.(behavField).hemisphere = cat(1,data.HbT.(behavField).hemisphere,'LH','RH');
    end
end
% take mean/StD
for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.HbT.(behavField).meanMax = mean(data.HbT.(behavField).catMax,1);
    data.HbT.(behavField).stdMax = std(data.HbT.(behavField).catMax,0,1);
    data.HbT.(behavField).meanP2P = mean(data.HbT.(behavField).catP2P,1);
    data.HbT.(behavField).stdP2P = std(data.HbT.(behavField).catP2P,0,1);
end
%% statistics - generalized linear mixed-effects model
% peak-to-peak
tableSize = cat(1,data.HbT.Rest.animalID,data.HbT.NREM.animalID,data.HbT.REM.animalID);
HbT_p2pTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','P2Pamplitude'});
HbT_p2pTable.Mouse = cat(1,data.HbT.Rest.animalID,data.HbT.NREM.animalID,data.HbT.REM.animalID);
HbT_p2pTable.Hemisphere = cat(1,data.HbT.Rest.hemisphere,data.HbT.NREM.hemisphere,data.HbT.REM.hemisphere);
HbT_p2pTable.Behavior = cat(1,data.HbT.Rest.behavior,data.HbT.NREM.behavior,data.HbT.REM.behavior);
HbT_p2pTable.P2Pamplitude = cat(1,data.HbT.Rest.catP2P,data.HbT.NREM.catP2P,data.HbT.REM.catP2P);
HbT_p2pFitFormula = 'P2Pamplitude ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
HbT_p2pStats = fitglme(HbT_p2pTable,HbT_p2pFitFormula);
% max
tableSize = cat(1,data.HbT.Rest.animalID,data.HbT.NREM.animalID,data.HbT.REM.animalID);
HbT_maxTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Peak'});
HbT_maxTable.Mouse = cat(1,data.HbT.Rest.animalID,data.HbT.NREM.animalID,data.HbT.REM.animalID);
HbT_maxTable.Hemisphere = cat(1,data.HbT.Rest.hemisphere,data.HbT.NREM.hemisphere,data.HbT.REM.hemisphere);
HbT_maxTable.Behavior = cat(1,data.HbT.Rest.behavior,data.HbT.NREM.behavior,data.HbT.REM.behavior);
HbT_maxTable.Peak = cat(1,data.HbT.Rest.catMax,data.HbT.NREM.catMax,data.HbT.REM.catMax);
HbT_maxFitFormula = 'Peak ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
HbT_maxStats = fitglme(HbT_maxTable,HbT_maxFitFormula);
%% vessel diameter comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    for bb =  1:length(behavFields)
        behavField = behavFields{1,bb};
        if isfield(AnalysisResults.(animalID).MeanVesselDiameter,behavField) == true
            vesselIDs = fieldnames(AnalysisResults.(animalID).MeanVesselDiameter.(behavField));
            for cc = 1:length(vesselIDs)
                vID = vesselIDs{cc,1};
                if strcmp(vID(1),'V') == false
                    for dd = 1:length(AnalysisResults.(animalID).MeanVesselDiameter.(behavField).(vID).indEvents)
                        data.TwoP.(animalID).(behavField).(vID).max(dd,1) = max(AnalysisResults.(animalID).MeanVesselDiameter.(behavField).(vID).indEvents{dd,1});
                        data.TwoP.(animalID).(behavField).(vID).p2p(dd,1) = abs(max(AnalysisResults.(animalID).MeanVesselDiameter.(behavField).(vID).indEvents{dd,1})) + abs(min(AnalysisResults.(animalID).MeanVesselDiameter.(behavField).(vID).indEvents{dd,1}));
                    end
                end
            end
        end
    end
end
% put data into arrays and prep for stats
for ee = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,ee};
    for ff = 1:length(behavFields)
        behavField = behavFields{1,ff};
        if isfield(data.TwoP,behavField) == false
            data.TwoP.(behavField).catMax = [];
            data.TwoP.(behavField).catP2P = [];
            data.TwoP.(behavField).animalID = {};
            data.TwoP.(behavField).vessel = {};
            data.TwoP.(behavField).behavior = {};
        end
        if isfield(data.TwoP.(animalID),behavField) == true
            vesselIDs = fieldnames(data.TwoP.(animalID).(behavField));
            for gg = 1:length(vesselIDs)
                vID = vesselIDs{gg,1};
                data.TwoP.(behavField).catMax = cat(1,data.TwoP.(behavField).catMax,mean(data.TwoP.(animalID).(behavField).(vID).max));
                data.TwoP.(behavField).catP2P = cat(1,data.TwoP.(behavField).catP2P,mean(data.TwoP.(animalID).(behavField).(vID).p2p));
                data.TwoP.(behavField).animalID = cat(1,data.TwoP.(behavField).animalID,animalID);
                data.TwoP.(behavField).vessel = cat(1,data.TwoP.(behavField).vessel,vID);
                data.TwoP.(behavField).behavior = cat(1,data.TwoP.(behavField).behavior,behavField);
            end
        end
    end
end
% take mean/StD
for hh = 1:length(behavFields)
    behavField = behavFields{1,hh};
    data.TwoP.(behavField).meanMax = mean(data.TwoP.(behavField).catMax,1);
    data.TwoP.(behavField).stdMax = std(data.TwoP.(behavField).catMax,0,1);
    data.TwoP.(behavField).meanP2P = mean(data.TwoP.(behavField).catP2P,1);
    data.TwoP.(behavField).stdP2P = std(data.TwoP.(behavField).catP2P,0,1);
end
%% statistics - generalized linear mixed-effects model
% peak-to-peak
tableSize = cat(1,data.TwoP.Rest.animalID,data.TwoP.NREM.animalID,data.TwoP.REM.animalID);
TwoP_p2pTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','P2Pamplitude'});
TwoP_p2pTable.Mouse = cat(1,data.TwoP.Rest.animalID,data.TwoP.NREM.animalID,data.TwoP.REM.animalID);
TwoP_p2pTable.Vessel = cat(1,data.TwoP.Rest.vessel,data.TwoP.NREM.vessel,data.TwoP.REM.vessel);
TwoP_p2pTable.Behavior = cat(1,data.TwoP.Rest.behavior,data.TwoP.NREM.behavior,data.TwoP.REM.behavior);
TwoP_p2pTable.P2Pamplitude = cat(1,data.TwoP.Rest.catP2P,data.TwoP.NREM.catP2P,data.TwoP.REM.catP2P);
TwoP_p2pFitFormula = 'P2Pamplitude ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Vessel)';
TwoP_p2pStats = fitglme(TwoP_p2pTable,TwoP_p2pFitFormula);
% max
tableSize = cat(1,data.TwoP.Rest.animalID,data.TwoP.NREM.animalID,data.TwoP.REM.animalID);
TwoP_maxTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Peak'});
TwoP_maxTable.Mouse = cat(1,data.TwoP.Rest.animalID,data.TwoP.NREM.animalID,data.TwoP.REM.animalID);
TwoP_maxTable.Vessel = cat(1,data.TwoP.Rest.vessel,data.TwoP.NREM.vessel,data.TwoP.REM.vessel);
TwoP_maxTable.Behavior = cat(1,data.TwoP.Rest.behavior,data.TwoP.NREM.behavior,data.TwoP.REM.behavior);
TwoP_maxTable.Peak = cat(1,data.TwoP.Rest.catMax,data.TwoP.NREM.catMax,data.TwoP.REM.catMax);
TwoP_maxFitFormula = 'Peak ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Vessel)';
TwoP_maxStats = fitglme(TwoP_maxTable,TwoP_maxFitFormula);
%% Fig. 1_S4
summaryFigure = figure('Name','Fig1_S4 (a-d)');
sgtitle('Figure 1_S4 - Turner et al. 2020')
%% [1_S4a] peak-to-peak HbT
ax1 = subplot(2,2,1);
xInds = ones(1,length(IOS_animalIDs)*2);
s1= scatter(xInds*1,data.HbT.Rest.catP2P,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.HbT.Rest.meanP2P,data.HbT.Rest.stdP2P,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(xInds*2,data.HbT.NREM.catP2P,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(2,data.HbT.NREM.meanP2P,data.HbT.NREM.stdP2P,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s3 = scatter(xInds*3,data.HbT.REM.catP2P,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(3,data.HbT.REM.meanP2P,data.HbT.REM.stdP2P,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[1_S4a] Mean Peak-to-Peak \Delta[HbT] (\muM)','during arousal-states'})
ylabel('Peak-to-peak \Delta[HbT] (\muM)')
legend([s1,s2,s3],'Rest','NREM','REM')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
ylim([0,180])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [1_S4b] peak HbT
ax2 = subplot(2,2,2);
xInds = ones(1,length(IOS_animalIDs)*2);
scatter(xInds*1,data.HbT.Rest.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.HbT.Rest.meanMax,data.HbT.Rest.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.HbT.NREM.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(2,data.HbT.NREM.meanMax,data.HbT.NREM.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*3,data.HbT.REM.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(3,data.HbT.REM.meanMax,data.HbT.REM.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[1_S4b] Mean Peak \Delta[HbT] (\muM)','during arousal-states'})
ylabel('Peak \Delta[HbT] (\muM)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
ylim([0,140])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [1_S4c] peak-to-peak D/D
ax3 = subplot(2,2,3);
scatter(ones(1,length(data.TwoP.Rest.catP2P))*1,data.TwoP.Rest.catP2P,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.TwoP.Rest.meanP2P,data.TwoP.Rest.stdP2P,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.TwoP.NREM.catP2P))*2,data.TwoP.NREM.catP2P,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(2,data.TwoP.NREM.meanP2P,data.TwoP.NREM.stdP2P,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.TwoP.REM.catP2P))*3,data.TwoP.REM.catP2P,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(3,data.TwoP.REM.meanP2P,data.TwoP.REM.stdP2P,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[1_S4c] Mean Peak-to-Peak \DeltaD/D (%)','during arousal-states'})
ylabel('Peak-to-peak \DeltaD/D (%)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
ylim([0,90])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [1_S4d] peak D/D
ax4 = subplot(2,2,4);
scatter(ones(1,length(data.TwoP.Rest.catP2P))*1,data.TwoP.Rest.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.TwoP.Rest.meanMax,data.TwoP.Rest.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.TwoP.NREM.catP2P))*2,data.TwoP.NREM.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(2,data.TwoP.NREM.meanMax,data.TwoP.NREM.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.TwoP.REM.catP2P))*3,data.TwoP.REM.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(3,data.TwoP.REM.meanMax,data.TwoP.REM.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[1_S4d] Mean Peak \DeltaD/D (%)','during arousal-states'})
ylabel('Peak \DeltaD/D (%)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
ylim([0,70])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    set(summaryFigure,'PaperPositionMode','auto');
    savefig(summaryFigure,[dirpath 'Fig1_S4']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig1_S4'])
    %% Text diary
    diaryFile = [dirpath 'Fig1_S4_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % peak-to-peak HbT statistical diary
    disp('======================================================================================================================')
    disp('[1_S4a] Generalized linear mixed-effects model statistics for peak-to-peak HbT during Rest, NREM, and REM')
    disp('======================================================================================================================')
    disp(HbT_p2pStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest P2P [HbT] (uM): ' num2str(round(data.HbT.Rest.meanP2P,1)) ' +/- ' num2str(round(data.HbT.Rest.stdP2P,1))]); disp(' ')
    disp(['NREM P2P [HbT] (uM): ' num2str(round(data.HbT.NREM.meanP2P,1)) ' +/- ' num2str(round(data.HbT.NREM.stdP2P,1))]); disp(' ')
    disp(['REM P2P [HbT] (uM): ' num2str(round(data.HbT.REM.meanP2P,1)) ' +/- ' num2str(round(data.HbT.REM.stdP2P,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % peak HbT statistical diary
    disp('======================================================================================================================')
    disp('[1_S4b] Generalized linear mixed-effects model statistics for peak HbT during Rest, NREM, and REM')
    disp('======================================================================================================================')
    disp(HbT_maxStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest Peak [HbT] (uM): ' num2str(round(data.HbT.Rest.meanMax,1)) ' +/- ' num2str(round(data.HbT.Rest.stdMax,1))]); disp(' ')
    disp(['NREM  Peak [HbT] (uM): ' num2str(round(data.HbT.NREM.meanMax,1)) ' +/- ' num2str(round(data.HbT.NREM.stdMax,1))]); disp(' ')
    disp(['REM  Peak [HbT] (uM): ' num2str(round(data.HbT.REM.meanMax,1)) ' +/- ' num2str(round(data.HbT.REM.stdMax,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % peak-to-peak D/D statistical diary
    disp('======================================================================================================================')
    disp('[1_S4c] Generalized linear mixed-effects model statistics for peak-to-peak D/D during Rest, NREM, and REM')
    disp('======================================================================================================================')
    disp(TwoP_p2pStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest P2P D/D (%): ' num2str(round(data.TwoP.Rest.meanP2P,1)) ' +/- ' num2str(round(data.TwoP.Rest.stdP2P,1))]); disp(' ')
    disp(['NREM P2P D/D (%): ' num2str(round(data.TwoP.NREM.meanP2P,1)) ' +/- ' num2str(round(data.TwoP.NREM.stdP2P,1))]); disp(' ')
    disp(['REM P2P D/D (%): ' num2str(round(data.TwoP.REM.meanP2P,1)) ' +/- ' num2str(round(data.TwoP.REM.stdP2P,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % peak D/D statistical diary
    disp('======================================================================================================================')
    disp('[1_S4d] Generalized linear mixed-effects model statistics for peak D/D during Rest, NREM, and REM')
    disp('======================================================================================================================')
    disp(TwoP_maxStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest Peak D/D (%): ' num2str(round(data.TwoP.Rest.meanMax,1)) ' +/- ' num2str(round(data.TwoP.Rest.stdMax,1))]); disp(' ')
    disp(['NREM  Peak D/D (%): ' num2str(round(data.TwoP.NREM.meanMax,1)) ' +/- ' num2str(round(data.TwoP.NREM.stdMax,1))]); disp(' ')
    disp(['REM  Peak D/D (%): ' num2str(round(data.TwoP.REM.meanMax,1)) ' +/- ' num2str(round(data.TwoP.REM.stdMax,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
end

end
