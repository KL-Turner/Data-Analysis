function [AnalysisResults] = Fig6_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generate figure panel 6-S1 for Turner_Gheres_Proctor_Drew_eLife2020
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
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Rest','NREM','REM'};
%% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        data.(behavField).adjLH.HbTvLFPxcVals(:,:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvLFPxcVals;
        data.(behavField).adjLH.LFP_lags(:,:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.LFP_lags;
        data.(behavField).adjLH.F(:,:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.F;
        data.(behavField).adjRH.HbTvLFPxcVals(:,:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvLFPxcVals;
        data.(behavField).adjRH.LFP_lags(:,:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.LFP_lags;
        data.(behavField).adjRH.F(:,:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.F;
        data.(behavField).adjLH.HbTvMUAxcVals(:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvMUAxcVals;
        data.(behavField).adjLH.HbTvMUAxcVals_std(:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvMUAxcVals_std;
        data.(behavField).adjLH.MUA_lags(:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjLH.LFP_lags;
        data.(behavField).adjRH.HbTvMUAxcVals(:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvMUAxcVals;
        data.(behavField).adjRH.HbTvMUAxcVals_std(:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvMUAxcVals_std;
        data.(behavField).adjRH.MUA_lags(:,aa) = AnalysisResults.(animalID).XCorr.(behavField).adjRH.LFP_lags;
        data.(behavField).animalID{aa,1} = animalID;
        data.(behavField).behavior{aa,1} = behavField;
        data.(behavField).LH{aa,1} = 'LH';
        data.(behavField).RH{aa,1} = 'RH';
    end
end
% concatenate the data from the left and right hemispheres
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    data.(behavField).cat_HbTvLFPxcVals = cat(3,data.(behavField).adjLH.HbTvLFPxcVals, data.(behavField).adjRH.HbTvLFPxcVals);
    data.(behavField).cat_LFP_lags = cat(3,data.(behavField).adjLH.LFP_lags, data.(behavField).adjRH.LFP_lags);
    data.(behavField).cat_LFP_F = cat(3,data.(behavField).adjLH.F, data.(behavField).adjRH.F);
    data.(behavField).cat_HbTvMUAxcVals = cat(2,data.(behavField).adjLH.HbTvMUAxcVals, data.(behavField).adjRH.HbTvMUAxcVals);
    data.(behavField).cat_MUA_lags = cat(2,data.(behavField).adjLH.MUA_lags, data.(behavField).adjRH.MUA_lags);
end
% take the averages of each field through the proper dimension
for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.(behavField).meanHbTvLFPxcVals = mean(data.(behavField).cat_HbTvLFPxcVals,3);
    data.(behavField).meanLFP_lags = mean(data.(behavField).cat_LFP_lags,3);
    data.(behavField).meanLFP_F = mean(data.(behavField).cat_LFP_F,3);
    data.(behavField).meanHbTvMUAxcVals = mean(data.(behavField).cat_HbTvMUAxcVals,2);
    data.(behavField).stdHbTvMUAxcVals = std(data.(behavField).cat_HbTvMUAxcVals,0,2);
    data.(behavField).meanMUA_lags = mean(data.(behavField).cat_MUA_lags,2);
end
%% find max/time to peak for MUA/Gamma-band power
for gg = 1:length(behavFields)
    behavField = behavFields{1,gg};
    for hh = 1:size(data.(behavField).cat_HbTvLFPxcVals,3)
        % gamma-band
        LFPmat = data.(behavField).cat_HbTvLFPxcVals(:,:,hh);
        gammaArray = mean(LFPmat(49:end,:),1);
        [gammaMax,gammaIndex] = max(gammaArray);
        data.(behavField).gammaPeak(hh,1) = gammaMax;
        data.(behavField).gammaTTP(hh,1) = data.(behavField).meanLFP_lags(gammaIndex)/30;
        % mua
        muaArray = data.(behavField).cat_HbTvMUAxcVals(:,hh)';
        [muaMax,muaIndex] = max(muaArray);
        data.(behavField).muaPeak(hh,1) = muaMax;
        data.(behavField).muaTTP(hh,1) = data.(behavField).meanMUA_lags(muaIndex)/30;
    end
end
% mean/std
for ii = 1:length(behavFields)
    behavField = behavFields{1,ii};
    data.(behavField).meanMuaPeak = mean(data.(behavField).muaPeak,1);
    data.(behavField).stdMuaPeak = std(data.(behavField).muaPeak,0,1);
    data.(behavField).meanMuaTTP = mean(data.(behavField).muaTTP,1);
    data.(behavField).stdMuaTTP = std(data.(behavField).muaTTP,0,1);
    data.(behavField).meanGammaPeak = mean(data.(behavField).gammaPeak,1);
    data.(behavField).stdGammaPeak = std(data.(behavField).gammaPeak,0,1);
    data.(behavField).meanGammaTTP = mean(data.(behavField).gammaTTP,1);
    data.(behavField).stdGammaTTP = std(data.(behavField).gammaTTP,0,1);
end
%% statistics - generalized linear mixed effects model
muaPeakTableSize = cat(1,data.Rest.muaPeak,data.NREM.muaPeak,data.REM.muaPeak);
muaPeakTable = table('Size',[size(muaPeakTableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Peak','Behavior','Hemisphere'});
muaPeakTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
muaPeakTable.Peak = cat(1,data.Rest.muaPeak,data.NREM.muaPeak,data.REM.muaPeak);
muaPeakTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
muaPeakTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
muaPeakFitFormula = 'Peak ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
muaPeakStats = fitglme(muaPeakTable,muaPeakFitFormula); %#ok<*NASGU>
%% statistics - generalized linear mixed effects model
muaTTPTableSize = cat(1,data.Rest.muaTTP,data.NREM.muaTTP,data.REM.muaTTP);
muaTTPTable = table('Size',[size(muaTTPTableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Peak','Behavior','Hemisphere'});
muaTTPTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
muaTTPTable.TTP = cat(1,data.Rest.muaTTP,data.NREM.muaTTP,data.REM.muaTTP);
muaTTPTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
muaTTPTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
muaTTPFitFormula = 'TTP ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
muaTTPStats = fitglme(muaTTPTable,muaTTPFitFormula); %#ok<*NASGU>
%% statistics - generalized linear mixed effects model
gammaPeakTableSize = cat(1,data.Rest.gammaPeak,data.NREM.gammaPeak,data.REM.gammaPeak);
gammaPeakTable = table('Size',[size(gammaPeakTableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Peak','Behavior','Hemisphere'});
gammaPeakTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
gammaPeakTable.Peak = cat(1,data.Rest.gammaPeak,data.NREM.gammaPeak,data.REM.gammaPeak);
gammaPeakTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
gammaPeakTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
gammaPeakFitFormula = 'Peak ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
gammaPeakStats = fitglme(gammaPeakTable,gammaPeakFitFormula); %#ok<*NASGU>
%% statistics - generalized linear mixed effects model
gammaTTPTableSize = cat(1,data.Rest.gammaTTP,data.NREM.gammaTTP,data.REM.gammaTTP);
gammaTTPTable = table('Size',[size(gammaTTPTableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Peak','Behavior','Hemisphere'});
gammaTTPTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
gammaTTPTable.TTP = cat(1,data.Rest.gammaTTP,data.NREM.gammaTTP,data.REM.gammaTTP);
gammaTTPTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
gammaTTPTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
gammaTTPFitFormula = 'TTP ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
gammaTTPStats = fitglme(gammaTTPTable,gammaTTPFitFormula); %#ok<*NASGU>
%% Fig. 6-S1
summaryFigure = figure('Name','Fig6-S1 (a-d)');
sgtitle('Figure 6-S1 - Turner et al. 2020')
%% [6-S1a] Peak cross-corr MUA
ax1 = subplot(2,2,1);
xInds = ones(1,length(animalIDs)*2);
s1 = scatter(xInds*1,data.Rest.muaPeak,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.meanMuaPeak,data.Rest.stdMuaPeak,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(xInds*2,data.NREM.muaPeak,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.NREM.meanMuaPeak,data.NREM.stdMuaPeak,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(xInds*3,data.REM.muaPeak,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.REM.meanMuaPeak,data.REM.stdMuaPeak,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[6-S1a] Peak cross-correlation MUA vs. \Delta[HbT] (\muM)','during arousal-states'})
ylabel('Peak corr. coef.')
legend([s1,s2,s3],'Rest','NREM','REM')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
ylim([0,0.6])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [6-S1b] Time-to-peak MUA
ax2 = subplot(2,2,2);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.Rest.muaTTP,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.meanMuaTTP,data.Rest.stdMuaTTP,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.NREM.muaTTP,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.NREM.meanMuaTTP,data.NREM.stdMuaTTP,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.REM.muaTTP,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.REM.meanMuaTTP,data.REM.stdMuaTTP,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[6-S1b] Time-to-peak MUA vs. \Delta[HbT] (\muM)','during arousal-states'})
ylabel('Time-to-peak (s)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
ylim([0,3.5])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [6-S1c] Peak cross-corr gamma-band
ax3 = subplot(2,2,3);
scatter(xInds*1,data.Rest.gammaPeak,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.meanGammaPeak,data.Rest.stdGammaPeak,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.NREM.gammaPeak,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.NREM.meanGammaPeak,data.NREM.stdGammaPeak,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.REM.gammaPeak,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.REM.meanGammaPeak,data.REM.stdGammaPeak,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[6-S1c] Peak cross-correlation gamma-band vs. \Delta[HbT] (\muM)','during arousal-states'})
ylabel('Peak corr. coef.')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
ylim([0,0.4])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [6-S1d] Time-to-peak gamma-band
ax4 = subplot(2,2,4);
scatter(xInds*1,data.Rest.gammaTTP,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.meanGammaTTP,data.Rest.stdGammaTTP,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.NREM.gammaTTP,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.NREM.meanGammaTTP,data.NREM.stdGammaTTP,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.REM.gammaTTP,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.REM.meanGammaTTP,data.REM.stdGammaTTP,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[6-S1d] Time-to-peak gamma-band vs. \Delta[HbT] (\muM)','during arousal-states'})
ylabel('Time-to-peak (s)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
ylim([0,3.5])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    set(summaryFigure,'PaperPositionMode','auto');
    savefig(summaryFigure,[dirpath 'Fig6-S1']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig6-S1'])
    %% Text diary
    diaryFile = [dirpath 'Fig6-S1_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % peak-to-peak HbT statistical diary
    disp('======================================================================================================================')
    disp('[6-S1a] Generalized linear mixed-effects model statistics for peak MUA cross-corr during Rest, NREM, and REM')
    disp('======================================================================================================================')
    disp(muaPeakStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest peak [MUA] (corrCoef): ' num2str(round(data.Rest.meanMuaPeak,2)) ' +/- ' num2str(round(data.Rest.stdMuaPeak,2))]); disp(' ')
    disp(['NREM peak [MUA] (corrCoef): ' num2str(round(data.NREM.meanMuaPeak,2)) ' +/- ' num2str(round(data.NREM.stdMuaPeak,2))]); disp(' ')
    disp(['REM peak [MUA] (corrCoef): ' num2str(round(data.REM.meanMuaPeak,2)) ' +/- ' num2str(round(data.REM.stdMuaPeak,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % peak HbT statistical diary
    disp('======================================================================================================================')
    disp('[6-S1b] Generalized linear mixed-effects model statistics for MUA time-to-peak during Rest, NREM, and REM')
    disp('======================================================================================================================')
    disp(muaTTPStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest TTP [MUA] (s): ' num2str(round(data.Rest.meanMuaTTP,2)) ' +/- ' num2str(round(data.Rest.stdMuaTTP,2))]); disp(' ')
    disp(['NREM  TTP [MUA] (s): ' num2str(round(data.NREM.meanMuaTTP,2)) ' +/- ' num2str(round(data.NREM.stdMuaTTP,2))]); disp(' ')
    disp(['REM  TTP [MUA] (s): ' num2str(round(data.REM.meanMuaTTP,2)) ' +/- ' num2str(round(data.REM.stdMuaTTP,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % peak-to-peak D/D statistical diary
    disp('======================================================================================================================')
    disp('[6-S1c] Generalized linear mixed-effects model statistics for peak gamma-band cross-corr during Rest, NREM, and REM')
    disp('======================================================================================================================')
    disp(gammaPeakStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest peak [gamma] (corrCoef): ' num2str(round(data.Rest.meanGammaPeak,2)) ' +/- ' num2str(round(data.Rest.stdGammaPeak,2))]); disp(' ')
    disp(['NREM peak [gamma] (corrCoef): ' num2str(round(data.NREM.meanGammaPeak,2)) ' +/- ' num2str(round(data.NREM.stdGammaPeak,2))]); disp(' ')
    disp(['REM peak [gamma] (corrCoef): ' num2str(round(data.REM.meanGammaPeak,2)) ' +/- ' num2str(round(data.REM.stdGammaPeak,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % peak D/D statistical diary
    disp('======================================================================================================================')
    disp('[6-S1d] Generalized linear mixed-effects model statistics for gamma-band time-to-peak during Rest, NREM, and REM')
    disp('======================================================================================================================')
    disp(gammaTTPStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest TTP [gamma] (s): ' num2str(round(data.Rest.meanGammaTTP,2)) ' +/- ' num2str(round(data.Rest.stdGammaTTP,2))]); disp(' ')
    disp(['NREM  TTP [gamma] (s): ' num2str(round(data.NREM.meanGammaTTP,2)) ' +/- ' num2str(round(data.NREM.stdGammaTTP,2))]); disp(' ')
    disp(['REM  TTP [gamma] (s): ' num2str(round(data.REM.meanGammaTTP,2)) ' +/- ' num2str(round(data.REM.stdGammaTTP,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
end

end
