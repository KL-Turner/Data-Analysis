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
driveLetters = {'E','E','E','F','F','F','D','D','D'};
behavFields = {'Whisk','Rest','NREM','REM'};
colorbrewer_setA_colorA = [(31/256) (120/256) (180/256)];
colorbrewer_setA_colorB = [(51/256) (160/256) (44/256)];
colorbrewer_setA_colorC = [(255/256) (140/256) (0/256)];
colorbrewer_setA_colorD = [(255/256) (0/256) (115/256)];

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        data.(behavField).meanLH(a,1) = mean(AnalysisResults.MeanCBV.(behavField).CBV_HbT.adjLH);
        data.(behavField).meanRH(a,1) = mean(AnalysisResults.MeanCBV.(behavField).CBV_HbT.adjRH);
        data.(behavField).allLH{a,1} = AnalysisResults.MeanCBV.(behavField).CBV_HbT.adjLH;
        data.(behavField).allRH{a,1} = AnalysisResults.MeanCBV.(behavField).CBV_HbT.adjRH;
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
scatter(xIndsA*1,data.Whisk.CBV_HbT.Comb,'MarkerEdgeColor',colorbrewer_setA_colorF,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.CBV_HbT.meanCBV,data.Whisk.CBV_HbT.stdCBV,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
e1.Color = 'black';
scatter(xIndsA*2,data.Rest.CBV_HbT.Comb,'MarkerEdgeColor',colorbrewer_setA_colorA,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Rest.CBV_HbT.meanCBV,data.Rest.CBV_HbT.stdCBV,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
e2.Color = 'black';
scatter(xIndsA*3,data.NREM.CBV_HbT.Comb,'MarkerEdgeColor',colorbrewer_setA_colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.NREM.CBV_HbT.meanCBV,data.NREM.CBV_HbT.stdCBV,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
e3.Color = 'black';
scatter(xIndsA*4,data.REM.CBV_HbT.Comb,'MarkerEdgeColor',colorbrewer_setA_colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.REM.CBV_HbT.meanCBV,data.REM.CBV_HbT.stdCBV,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
e4.Color = 'black';
title('Mean \DeltaHbT (\muM)')
title('\DeltaHbT (\muM)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])

subplot(2,2,4);
edges = -20:2.5:125;
h5 = histogram(data.Whisk.CBV_HbT.allComb,edges,'Normalization','probability');
hold on
p5 = plot(conv(h5.BinEdges,[0.5 0.5],'valid'),sgolayfilt((h5.BinCounts/sum(h5.BinCounts)),2,3));
h6 = histogram(data.Rest.CBV_HbT.(baselineType).allComb,edges,'Normalization','probability');
p6 = plot(conv(h6.BinEdges,[0.5 0.5],'valid'),sgolayfilt((h6.BinCounts/sum(h6.BinCounts)),2,3));
h7 = histogram(data.NREM.CBV_HbT.allComb,edges,'Normalization','probability');
p7 = plot(conv(h7.BinEdges,[0.5 0.5],'valid'),sgolayfilt((h7.BinCounts/sum(h7.BinCounts)),2,17));
h8 = histogram(data.REM.CBV_HbT.allComb,edges,'Normalization','probability');
p8 = plot(conv(h8.BinEdges,[0.5 0.5],'valid'),sgolayfilt((h8.BinCounts/sum(h8.BinCounts)),2,17));
axis square
ylim([0 0.5])

%% Heart rate
xIndsB = ones(1,length(animalIDs));
scatter(xIndsB*1,data.Whisk.(baselineType).HR,'MarkerEdgeColor',colorbrewer_setA_colorF,'jitter','on','jitterAmount',0.25);
hold on
e5 = errorbar(1,data.Whisk.(baselineType).meanHR,data.Whisk.(baselineType).stdHR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
e5.Color = 'black';
scatter(xIndsB*2,data.Rest.(baselineType).HR,'MarkerEdgeColor',colorbrewer_setA_colorA,'jitter','on','jitterAmount',0.25);
e6 = errorbar(2,data.Rest.(baselineType).meanHR,data.Rest.(baselineType).stdHR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
e6.Color = 'black';
scatter(xIndsB*3,data.NREM.HR,'MarkerEdgeColor',colorbrewer_setA_colorB,'jitter','on','jitterAmount',0.25);
e7 = errorbar(3,data.NREM.meanHR,data.NREM.stdHR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
e7.Color = 'black';
scatter(xIndsB*4,data.REM.HR,'MarkerEdgeColor',colorbrewer_setA_colorC,'jitter','on','jitterAmount',0.25);
e8 = errorbar(4,data.REM.meanHR,data.REM.stdHR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
e8.Color = 'black';
title('Mean Heart Rate')
ylabel('Frequency (Hz)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
legend([e1 e2 e3 e4 e5 e6],'Whisk','Rest','NREM','REM','Unstim','All')
axis square
xlim([0 length(behavFields)+1])

% save figure(s)
dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\';
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Hemodynamics and Heart Rate']);

