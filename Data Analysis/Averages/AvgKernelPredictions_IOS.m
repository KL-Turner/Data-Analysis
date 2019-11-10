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
behavFields = {'Whisk','Rest','NREM','REM','Unstim','All'};
baselineTypes = {'manualSelection','setDuration','entireDuration'};
corrCoeff_dataTypes = {'CBV','CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','muaPower'};
colorbrewer_setA_colorA = [0.520000 0.520000 0.510000];
colorbrewer_setA_colorB = [(31/256) (120/256) (180/256)];
colorbrewer_setA_colorC = [(255/256) (0/256) (115/256)];
colorbrewer_setA_colorD = [(51/256) (160/256) (44/256)];
colorbrewer_setA_colorE = [(255/256) (140/256) (0/256)];
colorbrewer_setA_colorF = [0.750000 0.000000 1.000000];

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        for c = 1:length(corrCoeff_dataTypes)
            corrCoeff_dataType = corrCoeff_dataTypes{1,c};
            if strcmp(behavField,'Rest') == true || strcmp(behavField,'Whisk') == true
                for d = 1:length(baselineTypes)
                    baselineType = baselineTypes{1,d};
                    data.(behavField).(corrCoeff_dataType).(baselineType).R(a,1) = AnalysisResults.CorrCoeff.(behavField).(corrCoeff_dataType).(baselineType).meanR;
                    data.(behavField).(corrCoeff_dataType).(baselineType).stdR(a,1) = AnalysisResults.CorrCoeff.(behavField).(corrCoeff_dataType).(baselineType).stdR;
                    data.(behavField).(corrCoeff_dataType).(baselineType).allCombR{a,1} = AnalysisResults.CorrCoeff.(behavField).(corrCoeff_dataType).(baselineType).R;
                end
            else
                data.(behavField).(corrCoeff_dataType).R(a,1) = AnalysisResults.CorrCoeff.(behavField).(corrCoeff_dataType).meanR;
                data.(behavField).(corrCoeff_dataType).stdR(a,1) = AnalysisResults.CorrCoeff.(behavField).(corrCoeff_dataType).stdR;
                data.(behavField).(corrCoeff_dataType).allCombR{a,1} = AnalysisResults.CorrCoeff.(behavField).(corrCoeff_dataType).R;
            end
        end
    end
end

% take the mean and standard deviation of each set of signals
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(corrCoeff_dataTypes)
        corrCoeff_dataType = corrCoeff_dataTypes{1,f};
        if strcmp(behavField,'Rest') || strcmp(behavField,'Whisk') == true
            for g = 1:length(baselineTypes)
                baselineType = baselineTypes{1,g};
                catR = [];
                data.(behavField).(corrCoeff_dataType).(baselineType).meanR = mean(data.(behavField).(corrCoeff_dataType).(baselineType).R,1);
                data.(behavField).(corrCoeff_dataType).(baselineType).stdMeanR = std(data.(behavField).(corrCoeff_dataType).(baselineType).R,0,1);
                data.(behavField).(corrCoeff_dataType).(baselineType).meanStdR = mean(data.(behavField).(corrCoeff_dataType).(baselineType).stdR,1);
                for z = 1:length(data.(behavField).(corrCoeff_dataType).(baselineType).allCombR)
                    catR = cat(1,catR,data.(behavField).(corrCoeff_dataType).(baselineType).allCombR{z,1});
                end
                data.(behavField).(corrCoeff_dataType).(baselineType).allComb = catR;
            end
        else
            catR = [];
            data.(behavField).(corrCoeff_dataType).meanR = mean(data.(behavField).(corrCoeff_dataType).R,1);
            data.(behavField).(corrCoeff_dataType).stdMeanR = std(data.(behavField).(corrCoeff_dataType).R,0,1);
            data.(behavField).(corrCoeff_dataType).meanStdR = mean(data.(behavField).(corrCoeff_dataType).stdR,1);
            for z = 1:length(data.(behavField).(corrCoeff_dataType).allCombR)
                catR = cat(1,catR,data.(behavField).(corrCoeff_dataType).allCombR{z,1});
            end
            data.(behavField).(corrCoeff_dataType).allComb = catR;
        end
    end
end

%% summary figure(s)
for h = 1:length(baselineTypes)
    baselineType = baselineTypes{1,h};
    xInds = ones(1,length(animalIDs));
    summaryFigure = figure;
    sgtitle({['L/R Pearson''s correlation coefficients for cortical data - ' baselineType],' '})
    %% CBV
    subplot(2,4,1);
    scatter(xInds*1,data.Whisk.CBV.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorF,'jitter','on', 'jitterAmount',0.25);
    hold on
    e1 = errorbar(1,data.Whisk.CBV.(baselineType).meanR,data.Whisk.CBV.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e1.Color = 'black';
    e2 = errorbar(1,data.Whisk.CBV.(baselineType).meanR,data.Whisk.CBV.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e2.Color = colorbrewer_setA_colorF;
    
    scatter(xInds*2,data.Rest.CBV.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
    e3 = errorbar(2,data.Rest.CBV.(baselineType).meanR,data.Rest.CBV.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e3.Color = 'black';
    e4 = errorbar(2,data.Rest.CBV.(baselineType).meanR,data.Rest.CBV.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e4.Color = colorbrewer_setA_colorA;
    
    scatter(xInds*3,data.NREM.CBV.R,'MarkerEdgeColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
    e5 = errorbar(3,data.NREM.CBV.meanR,data.NREM.CBV.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e5.Color = 'black';
    e6 = errorbar(3,data.NREM.CBV.meanR,data.NREM.CBV.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e6.Color = colorbrewer_setA_colorB;
    
    scatter(xInds*4,data.REM.CBV.R,'MarkerEdgeColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
    e7 = errorbar(4,data.REM.CBV.meanR,data.REM.CBV.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e7.Color = 'black';
    e8 = errorbar(4,data.REM.CBV.meanR,data.REM.CBV.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e8.Color = colorbrewer_setA_colorC;
    
    scatter(xInds*5,data.Unstim.CBV.R,'MarkerEdgeColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
    e9 = errorbar(5,data.Unstim.CBV.meanR,data.Unstim.CBV.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e9.Color = 'black';
    e10 = errorbar(5,data.Unstim.CBV.meanR,data.Unstim.CBV.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e10.Color = colorbrewer_setA_colorD;
    
    scatter(xInds*6,data.All.CBV.R,'MarkerEdgeColor',colorbrewer_setA_colorE,'jitter','on', 'jitterAmount',0.25);
    e11 = errorbar(6,data.All.CBV.meanR,data.All.CBV.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e11.Color = 'black';
    e12 = errorbar(6,data.All.CBV.meanR,data.All.CBV.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e12.Color = colorbrewer_setA_colorE;
    
    title('CBV Reflectance')
    ylabel('Correlation')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    legend([e1,e3,e5,e7,e9,e11],'Whisk','Rest','NREM','REM','Unstim','All')
    axis square
    ylim([0 1])
    xlim([0 length(behavFields)+1])
    
    %% CBV HbT
    subplot(2,4,2);
    scatter(xInds*1,data.Whisk.CBV_HbT.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorF,'jitter','on', 'jitterAmount',0.25);
    hold on
    e1 = errorbar(1,data.Whisk.CBV_HbT.(baselineType).meanR,data.Whisk.CBV_HbT.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e1.Color = 'black';
    e2 = errorbar(1,data.Whisk.CBV_HbT.(baselineType).meanR,data.Whisk.CBV_HbT.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e2.Color = colorbrewer_setA_colorF;
    
    scatter(xInds*2,data.Rest.CBV_HbT.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
    e3 = errorbar(2,data.Rest.CBV_HbT.(baselineType).meanR,data.Rest.CBV_HbT.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e3.Color = 'black';
    e4 = errorbar(2,data.Rest.CBV_HbT.(baselineType).meanR,data.Rest.CBV_HbT.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e4.Color = colorbrewer_setA_colorA;
    
    scatter(xInds*3,data.NREM.CBV_HbT.R,'MarkerEdgeColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
    e5 = errorbar(3,data.NREM.CBV_HbT.meanR,data.NREM.CBV_HbT.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e5.Color = 'black';
    e6 = errorbar(3,data.NREM.CBV_HbT.meanR,data.NREM.CBV_HbT.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e6.Color = colorbrewer_setA_colorB;
    
    scatter(xInds*4,data.REM.CBV_HbT.R,'MarkerEdgeColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
    e7 = errorbar(4,data.REM.CBV_HbT.meanR,data.REM.CBV_HbT.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e7.Color = 'black';
    e8 = errorbar(4,data.REM.CBV_HbT.meanR,data.REM.CBV_HbT.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e8.Color = colorbrewer_setA_colorC;
    
    scatter(xInds*5,data.Unstim.CBV_HbT.R,'MarkerEdgeColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
    e9 = errorbar(5,data.Unstim.CBV_HbT.meanR,data.Unstim.CBV_HbT.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e9.Color = 'black';
    e10 = errorbar(5,data.Unstim.CBV_HbT.meanR,data.Unstim.CBV_HbT.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e10.Color = colorbrewer_setA_colorD;
    
    scatter(xInds*6,data.All.CBV_HbT.R,'MarkerEdgeColor',colorbrewer_setA_colorE,'jitter','on', 'jitterAmount',0.25);
    e11 = errorbar(6,data.All.CBV_HbT.meanR,data.All.CBV_HbT.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e11.Color = 'black';
    e12 = errorbar(6,data.All.CBV_HbT.meanR,data.All.CBV_HbT.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e12.Color = colorbrewer_setA_colorE;
    
    title('CBV HbT')
    ylabel('Correlation')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis square
    ylim([0 1])
    xlim([0 length(behavFields)+1])
    
    %% Delta-band power
    subplot(2,4,3);
    scatter(xInds*1,data.Whisk.deltaBandPower.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorF,'jitter','on', 'jitterAmount',0.25);
    hold on
    e1 = errorbar(1,data.Whisk.deltaBandPower.(baselineType).meanR,data.Whisk.deltaBandPower.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e1.Color = 'black';
    e2 = errorbar(1,data.Whisk.deltaBandPower.(baselineType).meanR,data.Whisk.deltaBandPower.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e2.Color = colorbrewer_setA_colorF;
    
    scatter(xInds*2,data.Rest.deltaBandPower.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
    e3 = errorbar(2,data.Rest.deltaBandPower.(baselineType).meanR,data.Rest.deltaBandPower.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e3.Color = 'black';
    e4 = errorbar(2,data.Rest.deltaBandPower.(baselineType).meanR,data.Rest.deltaBandPower.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e4.Color = colorbrewer_setA_colorA;
    
    scatter(xInds*3,data.NREM.deltaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
    e5 = errorbar(3,data.NREM.deltaBandPower.meanR,data.NREM.deltaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e5.Color = 'black';
    e6 = errorbar(3,data.NREM.deltaBandPower.meanR,data.NREM.deltaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e6.Color = colorbrewer_setA_colorB;
    
    scatter(xInds*4,data.REM.deltaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
    e7 = errorbar(4,data.REM.deltaBandPower.meanR,data.REM.deltaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e7.Color = 'black';
    e8 = errorbar(4,data.REM.deltaBandPower.meanR,data.REM.deltaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e8.Color = colorbrewer_setA_colorC;
    
    scatter(xInds*5,data.Unstim.deltaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
    e9 = errorbar(5,data.Unstim.deltaBandPower.meanR,data.Unstim.deltaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e9.Color = 'black';
    e10 = errorbar(5,data.Unstim.deltaBandPower.meanR,data.Unstim.deltaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e10.Color = colorbrewer_setA_colorD;
    
    scatter(xInds*6,data.All.deltaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorE,'jitter','on', 'jitterAmount',0.25);
    e11 = errorbar(6,data.All.deltaBandPower.meanR,data.All.deltaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e11.Color = 'black';
    e12 = errorbar(6,data.All.deltaBandPower.meanR,data.All.deltaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e12.Color = colorbrewer_setA_colorE;
    
    title('Delta-band power')
    ylabel('Correlation')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis square
    ylim([0 1])
    xlim([0 length(behavFields)+1])
    
    %% Theta-band power
    subplot(2,4,4);
    scatter(xInds*1,data.Whisk.thetaBandPower.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorF,'jitter','on', 'jitterAmount',0.25);
    hold on
    e1 = errorbar(1,data.Whisk.thetaBandPower.(baselineType).meanR,data.Whisk.thetaBandPower.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e1.Color = 'black';
    e2 = errorbar(1,data.Whisk.thetaBandPower.(baselineType).meanR,data.Whisk.thetaBandPower.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e2.Color = colorbrewer_setA_colorF;
    
    scatter(xInds*2,data.Rest.thetaBandPower.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
    e3 = errorbar(2,data.Rest.thetaBandPower.(baselineType).meanR,data.Rest.thetaBandPower.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e3.Color = 'black';
    e4 = errorbar(2,data.Rest.thetaBandPower.(baselineType).meanR,data.Rest.thetaBandPower.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e4.Color = colorbrewer_setA_colorA;
    
    scatter(xInds*3,data.NREM.thetaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
    e5 = errorbar(3,data.NREM.thetaBandPower.meanR,data.NREM.thetaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e5.Color = 'black';
    e6 = errorbar(3,data.NREM.thetaBandPower.meanR,data.NREM.thetaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e6.Color = colorbrewer_setA_colorB;
    
    scatter(xInds*4,data.REM.thetaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
    e7 = errorbar(4,data.REM.thetaBandPower.meanR,data.REM.thetaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e7.Color = 'black';
    e8 = errorbar(4,data.REM.thetaBandPower.meanR,data.REM.thetaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e8.Color = colorbrewer_setA_colorC;
    
    scatter(xInds*5,data.Unstim.thetaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
    e9 = errorbar(5,data.Unstim.thetaBandPower.meanR,data.Unstim.thetaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e9.Color = 'black';
    e10 = errorbar(5,data.Unstim.thetaBandPower.meanR,data.Unstim.thetaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e10.Color = colorbrewer_setA_colorD;
    
    scatter(xInds*6,data.All.thetaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorE,'jitter','on', 'jitterAmount',0.25);
    e11 = errorbar(6,data.All.thetaBandPower.meanR,data.All.thetaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e11.Color = 'black';
    e12 = errorbar(6,data.All.thetaBandPower.meanR,data.All.thetaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e12.Color = colorbrewer_setA_colorE;
    
    title('Theta-band power')
    ylabel('Correlation')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis square
    ylim([0 1])
    xlim([0 length(behavFields)+1])
    
    %% Alpha-band power
    subplot(2,4,5);
    scatter(xInds*1,data.Whisk.alphaBandPower.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorF,'jitter','on', 'jitterAmount',0.25);
    hold on
    e1 = errorbar(1,data.Whisk.alphaBandPower.(baselineType).meanR,data.Whisk.alphaBandPower.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e1.Color = 'black';
    e2 = errorbar(1,data.Whisk.alphaBandPower.(baselineType).meanR,data.Whisk.alphaBandPower.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e2.Color = colorbrewer_setA_colorF;
    
    scatter(xInds*2,data.Rest.alphaBandPower.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
    e3 = errorbar(2,data.Rest.alphaBandPower.(baselineType).meanR,data.Rest.alphaBandPower.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e3.Color = 'black';
    e4 = errorbar(2,data.Rest.alphaBandPower.(baselineType).meanR,data.Rest.alphaBandPower.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e4.Color = colorbrewer_setA_colorA;
    
    scatter(xInds*3,data.NREM.alphaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
    e5 = errorbar(3,data.NREM.alphaBandPower.meanR,data.NREM.alphaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e5.Color = 'black';
    e6 = errorbar(3,data.NREM.alphaBandPower.meanR,data.NREM.alphaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e6.Color = colorbrewer_setA_colorB;
    
    scatter(xInds*4,data.REM.alphaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
    e7 = errorbar(4,data.REM.alphaBandPower.meanR,data.REM.alphaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e7.Color = 'black';
    e8 = errorbar(4,data.REM.alphaBandPower.meanR,data.REM.alphaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e8.Color = colorbrewer_setA_colorC;
    
    scatter(xInds*5,data.Unstim.alphaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
    e9 = errorbar(5,data.Unstim.alphaBandPower.meanR,data.Unstim.alphaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e9.Color = 'black';
    e10 = errorbar(5,data.Unstim.alphaBandPower.meanR,data.Unstim.alphaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e10.Color = colorbrewer_setA_colorD;
    
    scatter(xInds*6,data.All.alphaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorE,'jitter','on', 'jitterAmount',0.25);
    e11 = errorbar(6,data.All.alphaBandPower.meanR,data.All.alphaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e11.Color = 'black';
    e12 = errorbar(6,data.All.alphaBandPower.meanR,data.All.alphaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e12.Color = colorbrewer_setA_colorE;
    
    title('Alpha-band power')
    ylabel('Correlation')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis square
    ylim([0 1])
    xlim([0 length(behavFields)+1])
    
    %% Beta-band power
    subplot(2,4,6);
    scatter(xInds*1,data.Whisk.betaBandPower.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorF,'jitter','on', 'jitterAmount',0.25);
    hold on
    e1 = errorbar(1,data.Whisk.betaBandPower.(baselineType).meanR,data.Whisk.betaBandPower.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e1.Color = 'black';
    e2 = errorbar(1,data.Whisk.betaBandPower.(baselineType).meanR,data.Whisk.betaBandPower.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e2.Color = colorbrewer_setA_colorF;
    
    scatter(xInds*2,data.Rest.betaBandPower.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
    e3 = errorbar(2,data.Rest.betaBandPower.(baselineType).meanR,data.Rest.betaBandPower.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e3.Color = 'black';
    e4 = errorbar(2,data.Rest.betaBandPower.(baselineType).meanR,data.Rest.betaBandPower.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e4.Color = colorbrewer_setA_colorA;
    
    scatter(xInds*3,data.NREM.betaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
    e5 = errorbar(3,data.NREM.betaBandPower.meanR,data.NREM.betaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e5.Color = 'black';
    e6 = errorbar(3,data.NREM.betaBandPower.meanR,data.NREM.betaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e6.Color = colorbrewer_setA_colorB;
    
    scatter(xInds*4,data.REM.betaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
    e7 = errorbar(4,data.REM.betaBandPower.meanR,data.REM.betaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e7.Color = 'black';
    e8 = errorbar(4,data.REM.betaBandPower.meanR,data.REM.betaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e8.Color = colorbrewer_setA_colorC;
    
    scatter(xInds*5,data.Unstim.betaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
    e9 = errorbar(5,data.Unstim.betaBandPower.meanR,data.Unstim.betaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e9.Color = 'black';
    e10 = errorbar(5,data.Unstim.betaBandPower.meanR,data.Unstim.betaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e10.Color = colorbrewer_setA_colorD;
    
    scatter(xInds*6,data.All.betaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorE,'jitter','on', 'jitterAmount',0.25);
    e11 = errorbar(6,data.All.betaBandPower.meanR,data.All.betaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e11.Color = 'black';
    e12 = errorbar(6,data.All.betaBandPower.meanR,data.All.betaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e12.Color = colorbrewer_setA_colorE;
    
    title('Beta-band power')
    ylabel('Correlation')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis square
    ylim([0 1])
    xlim([0 length(behavFields)+1])
    
    %% Gamma-band power
    subplot(2,4,7);
    scatter(xInds*1,data.Whisk.gammaBandPower.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorF,'jitter','on', 'jitterAmount',0.25);
    hold on
    e1 = errorbar(1,data.Whisk.gammaBandPower.(baselineType).meanR,data.Whisk.gammaBandPower.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e1.Color = 'black';
    e2 = errorbar(1,data.Whisk.gammaBandPower.(baselineType).meanR,data.Whisk.gammaBandPower.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e2.Color = colorbrewer_setA_colorF;
    
    scatter(xInds*2,data.Rest.gammaBandPower.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
    e3 = errorbar(2,data.Rest.gammaBandPower.(baselineType).meanR,data.Rest.gammaBandPower.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e3.Color = 'black';
    e4 = errorbar(2,data.Rest.gammaBandPower.(baselineType).meanR,data.Rest.gammaBandPower.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e4.Color = colorbrewer_setA_colorA;
    
    scatter(xInds*3,data.NREM.gammaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
    e5 = errorbar(3,data.NREM.gammaBandPower.meanR,data.NREM.gammaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e5.Color = 'black';
    e6 = errorbar(3,data.NREM.gammaBandPower.meanR,data.NREM.gammaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e6.Color = colorbrewer_setA_colorB;
    
    scatter(xInds*4,data.REM.gammaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
    e7 = errorbar(4,data.REM.gammaBandPower.meanR,data.REM.gammaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e7.Color = 'black';
    e8 = errorbar(4,data.REM.gammaBandPower.meanR,data.REM.gammaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e8.Color = colorbrewer_setA_colorC;
    
    scatter(xInds*5,data.Unstim.gammaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
    e9 = errorbar(5,data.Unstim.gammaBandPower.meanR,data.Unstim.gammaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e9.Color = 'black';
    e10 = errorbar(5,data.Unstim.gammaBandPower.meanR,data.Unstim.gammaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e10.Color = colorbrewer_setA_colorD;
    
    scatter(xInds*6,data.All.gammaBandPower.R,'MarkerEdgeColor',colorbrewer_setA_colorE,'jitter','on', 'jitterAmount',0.25);
    e11 = errorbar(6,data.All.gammaBandPower.meanR,data.All.gammaBandPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e11.Color = 'black';
    e12 = errorbar(6,data.All.gammaBandPower.meanR,data.All.gammaBandPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e12.Color = colorbrewer_setA_colorE;
    
    title('Gamma-band power')
    ylabel('Correlation')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis square
    ylim([0 1])
    xlim([0 length(behavFields)+1])
    
    %% MUA power
    subplot(2,4,8);
    scatter(xInds*1,data.Whisk.muaPower.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorF,'jitter','on', 'jitterAmount',0.25);
    hold on
    e1 = errorbar(1,data.Whisk.muaPower.(baselineType).meanR,data.Whisk.muaPower.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e1.Color = 'black';
    e2 = errorbar(1,data.Whisk.muaPower.(baselineType).meanR,data.Whisk.muaPower.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e2.Color = colorbrewer_setA_colorF;
    
    scatter(xInds*2,data.Rest.muaPower.(baselineType).R,'MarkerEdgeColor',colorbrewer_setA_colorA,'jitter','on', 'jitterAmount',0.25);
    e3 = errorbar(2,data.Rest.muaPower.(baselineType).meanR,data.Rest.muaPower.(baselineType).stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e3.Color = 'black';
    e4 = errorbar(2,data.Rest.muaPower.(baselineType).meanR,data.Rest.muaPower.(baselineType).meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e4.Color = colorbrewer_setA_colorA;
    
    scatter(xInds*3,data.NREM.muaPower.R,'MarkerEdgeColor',colorbrewer_setA_colorB,'jitter','on', 'jitterAmount',0.25);
    e5 = errorbar(3,data.NREM.muaPower.meanR,data.NREM.muaPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e5.Color = 'black';
    e6 = errorbar(3,data.NREM.muaPower.meanR,data.NREM.muaPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e6.Color = colorbrewer_setA_colorB;
    
    scatter(xInds*4,data.REM.muaPower.R,'MarkerEdgeColor',colorbrewer_setA_colorC,'jitter','on', 'jitterAmount',0.25);
    e7 = errorbar(4,data.REM.muaPower.meanR,data.REM.muaPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e7.Color = 'black';
    e8 = errorbar(4,data.REM.muaPower.meanR,data.REM.muaPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e8.Color = colorbrewer_setA_colorC;
    
    scatter(xInds*5,data.Unstim.muaPower.R,'MarkerEdgeColor',colorbrewer_setA_colorD,'jitter','on', 'jitterAmount',0.25);
    e9 = errorbar(5,data.Unstim.muaPower.meanR,data.Unstim.muaPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e9.Color = 'black';
    e10 = errorbar(5,data.Unstim.muaPower.meanR,data.Unstim.muaPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e10.Color = colorbrewer_setA_colorD;
    
    scatter(xInds*6,data.All.muaPower.R,'MarkerEdgeColor',colorbrewer_setA_colorE,'jitter','on', 'jitterAmount',0.25);
    e11 = errorbar(6,data.All.muaPower.meanR,data.All.muaPower.stdMeanR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e11.Color = 'black';
    e12 = errorbar(6,data.All.muaPower.meanR,data.All.muaPower.meanStdR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e12.Color = colorbrewer_setA_colorE;
    
    title('MUA power')
    ylabel('Correlation')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis square
    ylim([0 1])
    xlim([0 length(behavFields)+1])
    
    % save figure(s)
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\HRF Response Function Predictions\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure, [dirpath baselineType '_AverageKernelPredictions']);
end
