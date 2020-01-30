%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Calculate the average correlation coefficient of the different behavioral states
%________________________________________________________________________________________________________________________
%
%   Inputs: none
%
%   Outputs: Generates summary figures saved to C: drive Documents folder
%
%   Last Revised: Oct 1st, 2019
%________________________________________________________________________________________________________________________
clear
clc

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111'};
driveLetters = {'M','M','M','M','M','M','M','M','M'};
behavFields = {'Whisk','Rest','NREM','REM'};
colorbrewer_setA_colorA = [(31/256) (120/256) (180/256)];
colorbrewer_setA_colorB = [(51/256) (160/256) (44/256)];
colorbrewer_setA_colorC = [(255/256) (140/256) (0/256)];
colorbrewer_setA_colorD = [(255/256) (0/256) (115/256)];

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\Turner_Manuscript_Summer2020\' animalID '\Bilateral Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        data.(behavField).CBV_HbT.meanLH(a,1) = mean(AnalysisResults.MeanCBV.(behavField).CBV_HbT.adjLH);
        data.(behavField).CBV_HbT.meanRH(a,1) = mean(AnalysisResults.MeanCBV.(behavField).CBV_HbT.adjRH);
        data.(behavField).CBV_HbT.allLH{a,1} = AnalysisResults.MeanCBV.(behavField).CBV_HbT.adjLH;
        data.(behavField).CBV_HbT.allRH{a,1} = AnalysisResults.MeanCBV.(behavField).CBV_HbT.adjRH;
        data.(behavField).HR(a,1) = mean(AnalysisResults.MeanHR.(behavField));
    end
end

% take the mean and standard deviation of each set of signals
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    data.(behavField).CBV_HbT.Comb = cat(1,data.(behavField).CBV_HbT.meanLH,data.(behavField).CBV_HbT.meanRH);
    data.(behavField).CBV_HbT.catAllLH = [];
    data.(behavField).CBV_HbT.catAllRH = [];
    for h = 1:length(data.(behavField).CBV_HbT.allLH)
        data.(behavField).CBV_HbT.catAllLH = cat(1,data.(behavField).CBV_HbT.catAllLH,data.(behavField).CBV_HbT.allLH{h,1});
        data.(behavField).CBV_HbT.catAllRH = cat(1,data.(behavField).CBV_HbT.catAllRH,data.(behavField).CBV_HbT.allRH{h,1});
    end
    data.(behavField).CBV_HbT.allComb = cat(1,data.(behavField).CBV_HbT.catAllLH,data.(behavField).CBV_HbT.catAllRH);
    data.(behavField).CBV_HbT.meanCBV = mean(data.(behavField).CBV_HbT.Comb);
    data.(behavField).CBV_HbT.stdCBV = std(data.(behavField).CBV_HbT.Comb,0,1);
    data.(behavField).meanHR = mean(data.(behavField).HR);
    data.(behavField).stdHR = std(data.(behavField).HR,0,1);
end

%% summary figure(s)
summaryFigure = figure;
sgtitle('Mean Hemodynamics and Heart Rate')
xIndsA = ones(1,length(animalIDs)*2);
%% CBV HbT
subplot(1,3,1);
s1 = scatter(xIndsA*1,data.Whisk.CBV_HbT.Comb,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.CBV_HbT.meanCBV,data.Whisk.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
s2 = scatter(xIndsA*2,data.Rest.CBV_HbT.Comb,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Rest.CBV_HbT.meanCBV,data.Rest.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
s3 = scatter(xIndsA*3,data.NREM.CBV_HbT.Comb,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.NREM.CBV_HbT.meanCBV,data.NREM.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
s4 = scatter(xIndsA*4,data.REM.CBV_HbT.Comb,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.REM.CBV_HbT.meanCBV,data.REM.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
title('Mean \DeltaHbT (\muM)')
ylabel('\DeltaHbT (\muM)')
legend([s1,s2,s3,s4],'Whisking','Awake Rest','NREM','REM','Location','NorthWest')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')

subplot(1,3,2);
edges = -25:15:130;
[curve1] = SmoothHistogramBins_IOS(data.Whisk.CBV_HbT.allComb,edges);
[curve2] = SmoothHistogramBins_IOS(data.Rest.CBV_HbT.allComb,edges);
[curve3] = SmoothHistogramBins_IOS(data.NREM.CBV_HbT.allComb,edges);
[curve4] = SmoothHistogramBins_IOS(data.REM.CBV_HbT.allComb,edges);
before = findall(gca);
fnplt(curve1);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorD)
hold on
before = findall(gca);
fnplt(curve2);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorA)
before = findall(gca);
fnplt(curve3);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorB)
before = findall(gca);
fnplt(curve4);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorC)
title('\DeltaHbT (\muM) Distribution')
xlabel('\DeltaHbT (\muM)')
ylabel('Probability')
axis square
set(gca,'box','off')
axis tight

%% Heart rate
subplot(1,3,3)
xIndsB = ones(1,length(animalIDs));
scatter(xIndsB*1,data.Whisk.HR,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD,'jitter','on','jitterAmount',0.25);
hold on
e5 = errorbar(1,data.Whisk.meanHR,data.Whisk.stdHR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
scatter(xIndsB*2,data.Rest.HR,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA,'jitter','on','jitterAmount',0.25);
e6 = errorbar(2,data.Rest.meanHR,data.Rest.stdHR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
scatter(xIndsB*3,data.NREM.HR,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB,'jitter','on','jitterAmount',0.25);
e7 = errorbar(3,data.NREM.meanHR,data.NREM.stdHR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
scatter(xIndsB*4,data.REM.HR,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC,'jitter','on','jitterAmount',0.25);
e8 = errorbar(4,data.REM.meanHR,data.REM.stdHR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
title('Mean Heart Rate')
ylabel('Frequency (Hz)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields) + 1])
set(gca,'box','off')

% save figure(s)
dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\';
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Hemodynamics and Heart Rate']);
