function [AnalysisResults] = XCorr_Saporin(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose:
%________________________________________________________________________________________________________________________

%% set-up and process data
expGroups = {'C57BL6J','SSP-SAP','Blank-SAP'};
setName = 'IOS Set A';
animalIDs.all = {};
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs.all = horzcat(animalIDs.all,{folderList.name});
    animalIDs.(strrep(expGroups{1,aa},'-','_')) = {folderList.name};
end
treatments = {'C57BL6J','SSP_SAP','Blank_SAP'};
behavFields = {'Rest','NREM','REM'};
variables = {'HbTvLFPxcVals','LFP_lags','F','HbTvMUAxcVals','MUA_lags'};
%% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,aa};
    % recognize treatment based on animal group
    if ismember(animalID,animalIDs.C57BL6J) == true
        treatment = 'C57BL6J';
    elseif ismember(animalID,animalIDs.SSP_SAP) == true
        treatment = 'SSP_SAP';
    elseif ismember(animalID,animalIDs.Blank_SAP) == true
        treatment = 'Blank_SAP';
    end
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        % pre-allocate necessary variable fields
        data.(treatment).(behavField).adjLH.dummCheck = 1;
        for dd = 1:length(variables)
            if isfield(data.(treatment).(behavField).adjLH,variables{1,dd}) == false
                data.(treatment).(behavField).adjLH.(variables{1,dd}) = [];
                data.(treatment).(behavField).adjRH.(variables{1,dd}) = [];
            end
        end
        data.(treatment).(behavField).adjLH.HbTvLFPxcVals = cat(3,data.(treatment).(behavField).adjLH.HbTvLFPxcVals,AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvLFPxcVals);
        data.(treatment).(behavField).adjLH.LFP_lags = cat(3,data.(treatment).(behavField).adjLH.LFP_lags,AnalysisResults.(animalID).XCorr.(behavField).adjLH.LFP_lags);
        data.(treatment).(behavField).adjLH.F = cat(3,data.(treatment).(behavField).adjLH.F,AnalysisResults.(animalID).XCorr.(behavField).adjLH.F);
        data.(treatment).(behavField).adjRH.HbTvLFPxcVals = cat(3,data.(treatment).(behavField).adjRH.HbTvLFPxcVals,AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvLFPxcVals);
        data.(treatment).(behavField).adjRH.LFP_lags = cat(3,data.(treatment).(behavField).adjRH.LFP_lags,AnalysisResults.(animalID).XCorr.(behavField).adjRH.LFP_lags);
        data.(treatment).(behavField).adjRH.F = cat(3,data.(treatment).(behavField).adjRH.F,AnalysisResults.(animalID).XCorr.(behavField).adjRH.F);
        data.(treatment).(behavField).adjLH.HbTvMUAxcVals = cat(1,data.(treatment).(behavField).adjLH.HbTvMUAxcVals,AnalysisResults.(animalID).XCorr.(behavField).adjLH.HbTvMUAxcVals);
        data.(treatment).(behavField).adjLH.MUA_lags = cat(1,data.(treatment).(behavField).adjLH.MUA_lags,AnalysisResults.(animalID).XCorr.(behavField).adjLH.LFP_lags);
        data.(treatment).(behavField).adjRH.HbTvMUAxcVals = cat(1,data.(treatment).(behavField).adjRH.HbTvMUAxcVals,AnalysisResults.(animalID).XCorr.(behavField).adjRH.HbTvMUAxcVals);
        data.(treatment).(behavField).adjRH.MUA_lags = cat(1,data.(treatment).(behavField).adjRH.MUA_lags,AnalysisResults.(animalID).XCorr.(behavField).adjRH.LFP_lags);
    end
end
%% take the averages of each field through the proper dimension
for dd = 1:length(treatments)
    treatment = treatments{1,dd};
    for ff = 1:length(behavFields)
        behavField = behavFields{1,ff};
        % LH
        data.(treatment).(behavField).adjLH.meanHbTvLFPxcVals = mean(data.(treatment).(behavField).adjLH.HbTvLFPxcVals,3);
        data.(treatment).(behavField).adjLH.meanLFP_lags = mean(data.(treatment).(behavField).adjLH.LFP_lags,3);
        data.(treatment).(behavField).adjLH.meanLFP_F = mean(data.(treatment).(behavField).adjLH.F,3);
        data.(treatment).(behavField).adjLH.meanHbTvMUAxcVals = mean(data.(treatment).(behavField).adjLH.HbTvMUAxcVals,2);
        data.(treatment).(behavField).adjLH.stdHbTvMUAxcVals = std(data.(treatment).(behavField).adjLH.HbTvMUAxcVals,0,2);
        data.(treatment).(behavField).adjLH.meanMUA_lags = mean(data.(treatment).(behavField).adjLH.MUA_lags,2);
        % RH
        data.(treatment).(behavField).adjRH.meanHbTvLFPxcVals = mean(data.(treatment).(behavField).adjRH.HbTvLFPxcVals,3);
        data.(treatment).(behavField).adjRH.meanLFP_lags = mean(data.(treatment).(behavField).adjRH.LFP_lags,3);
        data.(treatment).(behavField).adjRH.meanLFP_F = mean(data.(treatment).(behavField).adjRH.F,3);
        data.(treatment).(behavField).adjRH.meanHbTvMUAxcVals = mean(data.(treatment).(behavField).adjRH.HbTvMUAxcVals,2);
        data.(treatment).(behavField).adjRH.stdHbTvMUAxcVals = std(data.(treatment).(behavField).adjRH.HbTvMUAxcVals,0,2);
        data.(treatment).(behavField).adjRH.meanMUA_lags = mean(data.(treatment).(behavField).adjRH.MUA_lags,2);
    end
end
%% Awake Rest XCorr
freq = 30;
restLag = 5;
sleepLag = 5;
summaryFigure1 = figure;
sgtitle('Awake Rest LFP-HbT Cross-correlation')
%% LH LFP-HbT XCorr [C57BL6J] rest
subplot(3,2,1);
imagesc(data.C57BL6J.Rest.adjLH.meanLFP_lags,data.C57BL6J.Rest.adjLH.meanLFP_F,data.C57BL6J.Rest.adjLH.meanHbTvLFPxcVals)
title({'C57BL6J LH Awake Rest','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c1 = colorbar;
caxis([-0.1,0.1])
ylabel(c1,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% RH LFP-HbT XCorr [C57BL6J] rest
subplot(3,2,2);
imagesc(data.C57BL6J.Rest.adjRH.meanLFP_lags,data.C57BL6J.Rest.adjRH.meanLFP_F,data.C57BL6J.Rest.adjRH.meanHbTvLFPxcVals)
title({'C57BL6J RH Awake Rest','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c2 = colorbar;
caxis([-0.1,0.1])
ylabel(c2,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% LH LFP-HbT XCorr [SSP-SAP] rest
subplot(3,2,3);
imagesc(data.SSP_SAP.Rest.adjLH.meanLFP_lags,data.SSP_SAP.Rest.adjLH.meanLFP_F,data.SSP_SAP.Rest.adjLH.meanHbTvLFPxcVals)
title({'SSP-SAP UnRx LH Awake Rest','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c5 = colorbar;
caxis([-0.1,0.1])
ylabel(c5,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% RH LFP-HbT XCorr [SSP-SAP] rest
subplot(3,2,4);
imagesc(data.SSP_SAP.Rest.adjRH.meanLFP_lags,data.SSP_SAP.Rest.adjRH.meanLFP_F,data.SSP_SAP.Rest.adjRH.meanHbTvLFPxcVals)
title({'SSP-SAP treated RH Awake Rest','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c4 = colorbar;
caxis([-0.1,0.1])
ylabel(c4,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% LH LFP-HbT XCorr [Blank-SAP] rest
subplot(3,2,5);
imagesc(data.Blank_SAP.Rest.adjLH.meanLFP_lags,data.Blank_SAP.Rest.adjLH.meanLFP_F,data.Blank_SAP.Rest.adjLH.meanHbTvLFPxcVals)
title({'Blank-SAP UnRx LH Awake Rest','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c5 = colorbar;
caxis([-0.1,0.1])
ylabel(c5,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% RH LFP-HbT XCorr [Blank-SAP] rest
subplot(3,2,6);
imagesc(data.Blank_SAP.Rest.adjRH.meanLFP_lags,data.Blank_SAP.Rest.adjRH.meanLFP_F,data.Blank_SAP.Rest.adjRH.meanHbTvLFPxcVals)
title({'Blank-SAP treated RH Awake Rest','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c6 = colorbar;
caxis([-0.1,0.1])
ylabel(c6,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'XCorr_Rest']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'XCorr_Rest'])
end
%% NREM XCorr
summaryFigure2 = figure;
sgtitle('NREM LFP-HbT Cross-correlation')
%% LH LFP-HbT XCorr [C57BL6J] NREM
subplot(3,2,1);
imagesc(data.C57BL6J.NREM.adjLH.meanLFP_lags,data.C57BL6J.NREM.adjLH.meanLFP_F,data.C57BL6J.NREM.adjLH.meanHbTvLFPxcVals)
title({'C57BL6J LH NREM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c1 = colorbar;
caxis([-0.4,0.4])
ylabel(c1,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% RH LFP-HbT XCorr [C57BL6J] NREM
subplot(3,2,2);
imagesc(data.C57BL6J.NREM.adjRH.meanLFP_lags,data.C57BL6J.NREM.adjRH.meanLFP_F,data.C57BL6J.NREM.adjRH.meanHbTvLFPxcVals)
title({'C57BL6J RH NREM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c2 = colorbar;
caxis([-0.4,0.4])
ylabel(c2,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% LH LFP-HbT XCorr [SSP-SAP] NREM
subplot(3,2,3);
imagesc(data.SSP_SAP.NREM.adjLH.meanLFP_lags,data.SSP_SAP.NREM.adjLH.meanLFP_F,data.SSP_SAP.NREM.adjLH.meanHbTvLFPxcVals)
title({'SSP-SAP UnRx LH NREM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c5 = colorbar;
caxis([-0.4,0.4])
ylabel(c5,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% RH LFP-HbT XCorr [SSP-SAP] NREM
subplot(3,2,4);
imagesc(data.SSP_SAP.NREM.adjRH.meanLFP_lags,data.SSP_SAP.NREM.adjRH.meanLFP_F,data.SSP_SAP.NREM.adjRH.meanHbTvLFPxcVals)
title({'SSP-SAP treated RH NREM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c4 = colorbar;
caxis([-0.4,0.4])
ylabel(c4,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% LH LFP-HbT XCorr [Blank-SAP] NREM
subplot(3,2,5);
imagesc(data.Blank_SAP.NREM.adjLH.meanLFP_lags,data.Blank_SAP.NREM.adjLH.meanLFP_F,data.Blank_SAP.NREM.adjLH.meanHbTvLFPxcVals)
title({'Blank-SAP UnRx LH NREM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c5 = colorbar;
caxis([-0.4,0.4])
ylabel(c5,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% RH LFP-HbT XCorr [Blank-SAP] NREM
subplot(3,2,6);
imagesc(data.Blank_SAP.NREM.adjRH.meanLFP_lags,data.Blank_SAP.NREM.adjRH.meanLFP_F,data.Blank_SAP.NREM.adjRH.meanHbTvLFPxcVals)
title({'Blank-SAP treated RH NREM','LFP-[HbT] XCorr'})
xticks([-sleepLag*freq,-sleepLag*freq/2,0,sleepLag*freq/2,sleepLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-sleepLag*freq,sleepLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c6 = colorbar;
caxis([-0.4,0.4])
ylabel(c6,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'XCorr_NREM']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'XCorr_NREM'])
end
%% REM XCorr
summaryFigure3 = figure;
sgtitle('REM LFP-HbT Cross-correlation')
%% LH LFP-HbT XCorr [C57BL6J] rest
subplot(3,2,1);
imagesc(data.C57BL6J.REM.adjLH.meanLFP_lags,data.C57BL6J.REM.adjLH.meanLFP_F,data.C57BL6J.REM.adjLH.meanHbTvLFPxcVals)
title({'C57BL6J LH REM','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c1 = colorbar;
caxis([-0.1,0.1])
ylabel(c1,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% RH LFP-HbT XCorr [C57BL6J] REM
subplot(3,2,2);
imagesc(data.C57BL6J.REM.adjRH.meanLFP_lags,data.C57BL6J.REM.adjRH.meanLFP_F,data.C57BL6J.REM.adjRH.meanHbTvLFPxcVals)
title({'C57BL6J RH REM','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c2 = colorbar;
caxis([-0.1,0.1])
ylabel(c2,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% LH LFP-HbT XCorr [SSP-SAP] REM
subplot(3,2,3);
imagesc(data.SSP_SAP.REM.adjLH.meanLFP_lags,data.SSP_SAP.REM.adjLH.meanLFP_F,data.SSP_SAP.REM.adjLH.meanHbTvLFPxcVals)
title({'SSP-SAP UnRx LH REM','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c5 = colorbar;
caxis([-0.1,0.1])
ylabel(c5,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% RH LFP-HbT XCorr [SSP-SAP] REM
subplot(3,2,4);
imagesc(data.SSP_SAP.REM.adjRH.meanLFP_lags,data.SSP_SAP.REM.adjRH.meanLFP_F,data.SSP_SAP.REM.adjRH.meanHbTvLFPxcVals)
title({'SSP-SAP treated RH REM','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c4 = colorbar;
caxis([-0.1,0.1])
ylabel(c4,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% LH LFP-HbT XCorr [Blank-SAP] REM
subplot(3,2,5);
imagesc(data.Blank_SAP.REM.adjLH.meanLFP_lags,data.Blank_SAP.REM.adjLH.meanLFP_F,data.Blank_SAP.REM.adjLH.meanHbTvLFPxcVals)
title({'Blank-SAP UnRx LH REM','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c5 = colorbar;
caxis([-0.1,0.1])
ylabel(c5,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% RH LFP-HbT XCorr [Blank-SAP] REM
subplot(3,2,6);
imagesc(data.Blank_SAP.REM.adjRH.meanLFP_lags,data.Blank_SAP.REM.adjRH.meanLFP_F,data.Blank_SAP.REM.adjRH.meanHbTvLFPxcVals)
title({'Blank-SAP treated RH REM','LFP-[HbT] XCorr'})
xticks([-restLag*freq,-restLag*freq/2,0,restLag*freq/2,restLag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-restLag*freq,restLag*freq])
xlabel('Lags (s)')
ylabel('Freq (Hz)')
ylim([1,100])
c6 = colorbar;
caxis([-0.1,0.1])
ylabel(c6,{'Corr. coefficient';'LFP vs. \Delta[HbT] (\muM)'},'rotation',-90,'VerticalAlignment','bottom')
axis xy
axis square
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure3,[dirpath 'XCorr_REM']);
    set(summaryFigure3,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'XCorr_REM'])
end

end
