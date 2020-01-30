%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised: Oct 1st, 2019
%________________________________________________________________________________________________________________________
clear
clc

animalIDs = {'T101','T102','T103','T105','T108','T109','T110','T111'};
driveLetters = {'M','M','M','M','M','M','M','M','M'};
behavFields = {'Rest','NREM','REM','Whisk'};
corrCoeff_dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
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
        for c = 1:length(corrCoeff_dataTypes)
            powerSpec_dataType = corrCoeff_dataTypes{1,c};
            data.(behavField).(powerSpec_dataType).R{a,1} = AnalysisResults.CorrCoeff.(behavField).(powerSpec_dataType).R;
            data.(behavField).(powerSpec_dataType).meanRs(a,1) = AnalysisResults.CorrCoeff.(behavField).(powerSpec_dataType).meanR;
        end
    end
end

% concatenate the data and take mean/STD
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(corrCoeff_dataTypes)
        powerSpec_dataType = corrCoeff_dataTypes{1,f};
        data.(behavField).(powerSpec_dataType).catR = [];
        for z = 1:length(data.(behavField).(powerSpec_dataType).R) 
            data.(behavField).(powerSpec_dataType).catR = cat(1,data.(behavField).(powerSpec_dataType).catR,data.(behavField).(powerSpec_dataType).R{z,1});
        end
        data.(behavField).(powerSpec_dataType).meanR = mean(data.(behavField).(powerSpec_dataType).meanRs,1);
        data.(behavField).(powerSpec_dataType).stdR = std(data.(behavField).(powerSpec_dataType).meanRs,0,1);
    end
end

%% summary figure(s)
summaryFigure = figure;
sgtitle('Pearson''s Correlation Coefficients')
xInds = ones(1,length(animalIDs));
% CBV HbT
p1 = subplot(4,3,1);
s1 = scatter(xInds*1,data.Whisk.CBV_HbT.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.CBV_HbT.meanR,data.Whisk.CBV_HbT.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
s2 = scatter(xInds*2,data.Rest.CBV_HbT.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Rest.CBV_HbT.meanR,data.Rest.CBV_HbT.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
s3 = scatter(xInds*3,data.NREM.CBV_HbT.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NREM.CBV_HbT.meanR,data.NREM.CBV_HbT.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
s4 = scatter(xInds*4,data.REM.CBV_HbT.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.REM.CBV_HbT.meanR,data.REM.CBV_HbT.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
title('\DeltaHbT (\muM)')
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
legend([s1,s2,s3,s4],'Whisking','Awake Rest','NREM','REM','Location','SouthEast')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')

% Delta-band power
p2 = subplot(4,3,2);
scatter(xInds*1,data.Whisk.deltaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.deltaBandPower.meanR,data.Whisk.deltaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
scatter(xInds*2,data.Rest.deltaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Rest.deltaBandPower.meanR,data.Rest.deltaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
scatter(xInds*3,data.NREM.deltaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NREM.deltaBandPower.meanR,data.NREM.deltaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
scatter(xInds*4,data.REM.deltaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.REM.deltaBandPower.meanR,data.REM.deltaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
title('Delta-band [1-4 Hz]')
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')

% Theta-band power
p3 = subplot(4,3,3);
scatter(xInds*1,data.Whisk.thetaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.thetaBandPower.meanR,data.Whisk.thetaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
scatter(xInds*2,data.Rest.thetaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Rest.thetaBandPower.meanR,data.Rest.thetaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
scatter(xInds*3,data.NREM.thetaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NREM.thetaBandPower.meanR,data.NREM.thetaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
scatter(xInds*4,data.REM.thetaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.REM.thetaBandPower.meanR,data.REM.thetaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
title('Theta-band [4-10 Hz]')
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')

% Alpha-band power
p4 = subplot(4,3,7);
scatter(xInds*1,data.Whisk.alphaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.alphaBandPower.meanR,data.Whisk.alphaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
scatter(xInds*2,data.Rest.alphaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Rest.alphaBandPower.meanR,data.Rest.alphaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
scatter(xInds*3,data.NREM.alphaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NREM.alphaBandPower.meanR,data.NREM.alphaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
scatter(xInds*4,data.REM.alphaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.REM.alphaBandPower.meanR,data.REM.alphaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
title('Alpha-band [10-13 Hz]')
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')

% Beta-band power
p5 = subplot(4,3,8);
scatter(xInds*1,data.Whisk.betaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.betaBandPower.meanR,data.Whisk.betaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
scatter(xInds*2,data.Rest.betaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Rest.betaBandPower.meanR,data.Rest.betaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
scatter(xInds*3,data.NREM.betaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NREM.betaBandPower.meanR,data.NREM.betaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
scatter(xInds*4,data.REM.betaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.REM.betaBandPower.meanR,data.REM.betaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
title('Beta-band [13-30 Hz]')
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')

% Gamma-band power
p6 = subplot(4,3,9);
scatter(xInds*1,data.Whisk.gammaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.gammaBandPower.meanR,data.Whisk.gammaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
scatter(xInds*2,data.Rest.gammaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Rest.gammaBandPower.meanR,data.Rest.gammaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
scatter(xInds*3,data.NREM.gammaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NREM.gammaBandPower.meanR,data.NREM.gammaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
scatter(xInds*4,data.REM.gammaBandPower.meanRs,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.REM.gammaBandPower.meanR,data.REM.gammaBandPower.stdR,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
title('Gamma-band [30-100 Hz]')
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')

linkaxes([p1,p2,p3,p4,p5,p6],'xy')

%% summary figure(s)
edges = -0.1:0.025:1;
% CBV HbT
q1 = subplot(4,3,4);
[curve1] = SmoothHistogramBins_IOS(data.Whisk.CBV_HbT.catR,edges);
[curve2] = SmoothHistogramBins_IOS(data.Rest.CBV_HbT.catR,edges);
[curve3] = SmoothHistogramBins_IOS(data.NREM.CBV_HbT.catR,edges);
[curve4] = SmoothHistogramBins_IOS(data.REM.CBV_HbT.catR,edges);
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
title('\DeltaHbT (\muM)')
xlabel({'Corr. Coefficient';'Left hem vs. Right hem'})
ylabel('Probability')
axis square
set(gca,'box','off')

% Delta-band power
q2 = subplot(4,3,5);
[curve5] = SmoothHistogramBins_IOS(data.Whisk.deltaBandPower.catR,edges);
[curve6] = SmoothHistogramBins_IOS(data.Rest.deltaBandPower.catR,edges);
[curve7] = SmoothHistogramBins_IOS(data.NREM.deltaBandPower.catR,edges);
[curve8] = SmoothHistogramBins_IOS(data.REM.deltaBandPower.catR,edges);
before = findall(gca);
fnplt(curve5);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorD)
hold on
before = findall(gca);
fnplt(curve6);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorA)
before = findall(gca);
fnplt(curve7);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorB)
before = findall(gca);
fnplt(curve8);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorC)
title('Delta-band [1-4 Hz]')
xlabel({'Corr. Coefficient';'Left hem vs. Right hem'})
ylabel('Probability')
axis square
set(gca,'box','off')

% Theta-band power
q3 = subplot(4,3,6);
[curve9] = SmoothHistogramBins_IOS(data.Whisk.thetaBandPower.catR,edges);
[curve10] = SmoothHistogramBins_IOS(data.Rest.thetaBandPower.catR,edges);
[curve11] = SmoothHistogramBins_IOS(data.NREM.thetaBandPower.catR,edges);
[curve12] = SmoothHistogramBins_IOS(data.REM.thetaBandPower.catR,edges);
before = findall(gca);
fnplt(curve9);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorD)
hold on
before = findall(gca);
fnplt(curve10);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorA)
before = findall(gca);
fnplt(curve11);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorB)
before = findall(gca);
fnplt(curve12);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorC)
title('Theta-band [4-10 Hz]')
xlabel({'Corr. Coefficient';'Left hem vs. Right hem'})
ylabel('Probability')
axis square
set(gca,'box','off')

% Alpha-band power
q4 = subplot(4,3,10);
[curve13] = SmoothHistogramBins_IOS(data.Whisk.alphaBandPower.catR,edges);
[curve14] = SmoothHistogramBins_IOS(data.Rest.alphaBandPower.catR,edges);
[curve15] = SmoothHistogramBins_IOS(data.NREM.alphaBandPower.catR,edges);
[curve16] = SmoothHistogramBins_IOS(data.REM.alphaBandPower.catR,edges);
before = findall(gca);
fnplt(curve13);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorD)
hold on
before = findall(gca);
fnplt(curve14);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorA)
before = findall(gca);
fnplt(curve15);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorB)
before = findall(gca);
fnplt(curve16);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorC)
title('Alpha-band [10-13 Hz]')
xlabel({'Corr. Coefficient';'Left hem vs. Right hem'})
ylabel('Probability')
axis square
set(gca,'box','off')

% Beta-band power
q5 = subplot(4,3,11);
[curve17] = SmoothHistogramBins_IOS(data.Whisk.betaBandPower.catR,edges);
[curve18] = SmoothHistogramBins_IOS(data.Rest.betaBandPower.catR,edges);
[curve19] = SmoothHistogramBins_IOS(data.NREM.betaBandPower.catR,edges);
[curve20] = SmoothHistogramBins_IOS(data.REM.betaBandPower.catR,edges);
before = findall(gca);
fnplt(curve17);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorD)
hold on
before = findall(gca);
fnplt(curve18);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorA)
before = findall(gca);
fnplt(curve19);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorB)
before = findall(gca);
fnplt(curve20);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorC)
title('Beta-band [13-30 Hz]')
xlabel({'Corr. Coefficient';'Left hem vs. Right hem'})
ylabel('Probability')
axis square
set(gca,'box','off')

% Gamma-band power
q6 = subplot(4,3,12);
[curve21] = SmoothHistogramBins_IOS(data.Whisk.gammaBandPower.catR,edges);
[curve22] = SmoothHistogramBins_IOS(data.Rest.gammaBandPower.catR,edges);
[curve23] = SmoothHistogramBins_IOS(data.NREM.gammaBandPower.catR,edges);
[curve24] = SmoothHistogramBins_IOS(data.REM.gammaBandPower.catR,edges);
before = findall(gca);
fnplt(curve21);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorD)
hold on
before = findall(gca);
fnplt(curve22);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorA)
before = findall(gca);
fnplt(curve23);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorB)
before = findall(gca);
fnplt(curve24);
added = setdiff(findall(gca),before);
set(added,'Color',colorbrewer_setA_colorC)
title('Gamma-band [30-100 Hz]')
xlabel({'Corr. Coefficient';'Left hem vs. Right hem'})
ylabel('Probability')
axis square
set(gca,'box','off')

linkaxes([q1,q2,q3,q4,q5,q6],'xy')


%% save figure(s)
dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\';
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Peason''s Correlation Coefficients']);

