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

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110'};
driveLetters = {'E','E','E','F','F','F','D','D'};
behavFields = {'Whisk','Rest','NREM','REM','Unstim','All'};
baselineTypes = {'manualSelection','setDuration','entireDuration'};
colorbrewer_setA_colorA = [0.520000 0.520000 0.510000];
colorbrewer_setA_colorB = [(31/256) (120/256) (180/256)];
colorbrewer_setA_colorC = [(255/256) (0/256) (115/256)];
colorbrewer_setA_colorD = [(51/256) (160/256) (44/256)];
colorbrewer_setA_colorE = [(255/256) (140/256) (0/256)];
colorbrewer_setA_colorF = [(255/256) (140/256) (0/256)];

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        if strcmp(behavField,'Whisk') == true || strcmp(behavField,'Rest') == true
            for d = 1:length(baselineTypes)
                baselineType = baselineTypes{1,d};
                data.(behavField).(baselineType).HR(a,1) = mean(AnalysisResults.MeanHR.(behavField).(baselineType));
            end
        else
            data.(behavField).HR(a,1) = mean(AnalysisResults.MeanHR.(behavField));
        end
    end
end

% take the mean and standard deviation of each set of signals
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    if strcmp(behavField,'Whisk') == true || strcmp(behavField,'Rest') == true
        for g = 1:length(baselineTypes)
            baselineType = baselineTypes{1,g};
            data.(behavField).(baselineType).meanHR = mean(data.(behavField).(baselineType).HR);
            data.(behavField).(baselineType).stdHR = std(data.(behavField).(baselineType).HR,0,1);
        end
    else
        data.(behavField).meanHR = mean(data.(behavField).HR);
        data.(behavField).stdHR = std(data.(behavField).HR,0,1);
    end
end

%% summary figure(s)
for h = 1:length(baselineTypes)
    baselineType = baselineTypes{1,h};
    xInds = ones(1,length(animalIDs));
    summaryFigure = figure;
    %% CBV
    scatter(xInds*1,data.Whisk.(baselineType).HR,'MarkerEdgeColor',colorbrewer_setA_colorF,'jitter','on','jitterAmount',0.25);
    hold on
    e1 = errorbar(1,data.Whisk.(baselineType).meanHR,data.Whisk.(baselineType).stdHR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorF);
    e1.Color = 'black';
    
    scatter(xInds*2,data.Rest.(baselineType).HR,'MarkerEdgeColor',colorbrewer_setA_colorA,'jitter','on','jitterAmount',0.25);
    hold on
    e2 = errorbar(2,data.Rest.(baselineType).meanHR,data.Rest.(baselineType).stdHR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorA);
    e2.Color = 'black';
    
    scatter(xInds*3,data.NREM.HR,'MarkerEdgeColor',colorbrewer_setA_colorB,'jitter','on','jitterAmount',0.25);
    hold on
    e3 = errorbar(3,data.NREM.meanHR,data.NREM.stdHR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorB);
    e3.Color = 'black';
    
    scatter(xInds*4,data.REM.HR,'MarkerEdgeColor',colorbrewer_setA_colorC,'jitter','on','jitterAmount',0.25);
    hold on
    e4 = errorbar(4,data.REM.meanHR,data.REM.stdHR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorC);
    e4.Color = 'black';
    
    scatter(xInds*5,data.Unstim.HR,'MarkerEdgeColor',colorbrewer_setA_colorD,'jitter','on','jitterAmount',0.25);
    hold on
    e5 = errorbar(5,data.Unstim.meanHR,data.Unstim.stdHR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD);
    e5.Color = 'black';
    
    scatter(xInds*6,data.All.HR,'MarkerEdgeColor',colorbrewer_setA_colorE,'jitter','on','jitterAmount',0.25);
    hold on
    e6 = errorbar(6,data.All.meanHR,data.All.stdHR,'o','MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE);
    e6.Color = 'black';
    
    title(['Mean Heart Rate - ' baselineType])
    ylabel('Freq (Hz)')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    legend([e1 e2 e3 e4 e5 e6],'Whisk','Rest','NREM','REM','Unstim','All')
    axis square
    xlim([0 length(behavFields)+1])
    
    % save figure(s)
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Heart Rate\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure, [dirpath baselineType '_AverageHeartRate']);
end
