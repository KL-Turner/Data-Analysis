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
behavFields = {'Contra','Whisk','Rest'};
neuralBands = {'gammaBandPower','muaPower'};
colorbrewer_setA_colorA = [(31/256) (120/256) (180/256)];
colorbrewer_setA_colorB = [(51/256) (160/256) (44/256)];
colorbrewer_setA_colorC = [(255/256) (140/256) (0/256)];
colorbrewer_setA_colorD = [(255/256) (0/256) (115/256)];
colorbrewer_setA_colorE = [(192/256) (0/256) (256/256)];

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for c = 1:length(neuralBands)
        neuralBand = neuralBands{1,c};
        for b = 1:length(behavFields)
            behavior = behavFields{1,b};
            % Gamma function kernels
            data.(neuralBand).(behavior).adjLH.gammaTimeVec{a,1} = AnalysisResults.HRFs.(neuralBand).adjLH.(behavior).gammaTimeVec;
            data.(neuralBand).(behavior).adjRH.gammaTimeVec{a,1}  = AnalysisResults.HRFs.(neuralBand).adjLH.(behavior).gammaTimeVec;
            data.(neuralBand).(behavior).adjLH.gammaFunc{a,1} = AnalysisResults.HRFs.(neuralBand).adjRH.(behavior).gammaFunc;
            data.(neuralBand).(behavior).adjRH.gammaFunc{a,1} = AnalysisResults.HRFs.(neuralBand).adjRH.(behavior).gammaFunc;
            % Impulse response function kernels
            data.(neuralBand).(behavior).adjLH.IRtimeVec{a,1} = AnalysisResults.HRFs.(neuralBand).adjLH.(behavior).IRtimeVec;
            data.(neuralBand).(behavior).adjRH.IRtimeVec{a,1}  = AnalysisResults.HRFs.(neuralBand).adjLH.(behavior).IRtimeVec;
            data.(neuralBand).(behavior).adjLH.IR{a,1} = AnalysisResults.HRFs.(neuralBand).adjRH.(behavior).IR;
            data.(neuralBand).(behavior).adjRH.IR{a,1} = AnalysisResults.HRFs.(neuralBand).adjRH.(behavior).IR;
            % Behavior-derived R2 predictions
            data.(neuralBand).(behavior).adjLH.MedR2(a,1) = AnalysisResults.HRFs.(neuralBand).adjLH.(behavior).Med_IndR2;
            data.(neuralBand).(behavior).adjRH.MedR2(a,1) = AnalysisResults.HRFs.(neuralBand).adjRH.(behavior).Med_IndR2;
            data.(neuralBand).(behavior).adjLH.AveR2(a,1) = AnalysisResults.HRFs.(neuralBand).adjLH.(behavior).AveR2;
            data.(neuralBand).(behavior).adjRH.AveR2(a,1) = AnalysisResults.HRFs.(neuralBand).adjRH.(behavior).AveR2;
            % Sleep R2 predictions for each behavior's kernel
            data.(neuralBand).NREM.(behavior).adjLH.MedR2(a,1) = AnalysisResults.HRFs.(neuralBand).adjLH.NREM.(behavior).Med_IndR2;
            data.(neuralBand).NREM.(behavior).adjRH.MedR2(a,1) = AnalysisResults.HRFs.(neuralBand).adjRH.NREM.(behavior).Med_IndR2;
            data.(neuralBand).REM.(behavior).adjLH.MedR2(a,1) = AnalysisResults.HRFs.(neuralBand).adjLH.REM.(behavior).Med_IndR2;
            data.(neuralBand).REM.(behavior).adjRH.MedR2(a,1) = AnalysisResults.HRFs.(neuralBand).adjRH.REM.(behavior).Med_IndR2;
        end
    end
end

% take the mean and standard deviation of each set of signals
for c = 1:length(neuralBands)
    neuralBand = neuralBands{1,c};
    for b = 1:length(behavFields)
        behavior = behavFields{1,b};
        % average IR functions
        data.(neuralBand).(behavior).Comb.IR = cat(1,cell2mat(data.(neuralBand).(behavior).adjLH.IR),cell2mat(data.(neuralBand).(behavior).adjRH.IR));
        data.(neuralBand).(behavior).Comb.IRtimeVec = cat(1,cell2mat(data.(neuralBand).(behavior).adjLH.IRtimeVec),cell2mat(data.(neuralBand).(behavior).adjRH.IRtimeVec));
        data.(neuralBand).(behavior).meanIR = mean(data.(neuralBand).(behavior).Comb.IR,1);
        data.(neuralBand).(behavior).stdIR = std(data.(neuralBand).(behavior).Comb.IR,0,1);
        data.(neuralBand).(behavior).meanIRtimeVec = mean(data.(neuralBand).(behavior).Comb.IRtimeVec);
        data.(neuralBand).(behavior).stdIRtimeVec = std(data.(neuralBand).(behavior).Comb.IRtimeVec,0,1);
        % Behavior-derived R2 predictions
        data.(neuralBand).(behavior).Comb.MedR2 = cat(1,data.(neuralBand).(behavior).adjLH.MedR2,data.(neuralBand).(behavior).adjRH.MedR2);
        data.(neuralBand).(behavior).Comb.AveR2 = cat(1,data.(neuralBand).(behavior).adjLH.AveR2,data.(neuralBand).(behavior).adjRH.AveR2);
        data.(neuralBand).(behavior).meanMedR2 = mean(data.(neuralBand).(behavior).Comb.MedR2);
        data.(neuralBand).(behavior).stdMedR2 = std(data.(neuralBand).(behavior).Comb.MedR2,0,1);
        data.(neuralBand).(behavior).meanAveR2 = mean(data.(neuralBand).(behavior).Comb.AveR2);
        data.(neuralBand).(behavior).stdAveR2 = std(data.(neuralBand).(behavior).Comb.AveR2,0,1);
        % NREM
        data.(neuralBand).NREM.(behavior).Comb.MedR2 = cat(1,data.(neuralBand).NREM.(behavior).adjLH.MedR2,data.(neuralBand).NREM.(behavior).adjRH.MedR2);
        data.(neuralBand).NREM.(behavior).meanMedR2 = mean(data.(neuralBand).NREM.(behavior).Comb.MedR2);
        data.(neuralBand).NREM.(behavior).stdMedR2 = std(data.(neuralBand).NREM.(behavior).Comb.MedR2,0,1);
        % REM
        data.(neuralBand).REM.(behavior).Comb.MedR2 = cat(1,data.(neuralBand).REM.(behavior).adjLH.MedR2,data.(neuralBand).REM.(behavior).adjRH.MedR2);
        data.(neuralBand).REM.(behavior).meanMedR2 = mean(data.(neuralBand).REM.(behavior).Comb.MedR2);
        data.(neuralBand).REM.(behavior).stdMedR2 = std(data.(neuralBand).REM.(behavior).Comb.MedR2,0,1);
    end
end

%% summary figure(s)
summaryFigure = figure;
sgtitle('HRF Kernels and Median R^2 Predictions')
xIndsA = ones(1,length(animalIDs)*2);
%%
subplot(2,2,1)
p1 = plot(data.gammaBandPower.Contra.meanIRtimeVec,data.gammaBandPower.Contra.meanIR,'color',colorbrewer_setA_colorE,'LineWidth',2);
hold on
p2 = plot(data.gammaBandPower.Whisk.meanIRtimeVec,data.gammaBandPower.Whisk.meanIR,'color',colorbrewer_setA_colorD,'LineWidth',2);
p3 = plot(data.gammaBandPower.Rest.meanIRtimeVec,data.gammaBandPower.Rest.meanIR,'color',colorbrewer_setA_colorA,'LineWidth',2);
title('Mean IR Function')
xlabel('HRF Time (s)')
ylabel({'Gamma-band [30-100 Hz] derived';'HRF amplitude (A.U.)'})
legend([p1,p2,p3],'Sensory-evoked','Volitional whisk','Awake rest','Location','NorthWest')
axis square
set(gca,'box','off')

subplot(2,2,2)
p1 = plot(data.muaPower.Contra.meanIRtimeVec,data.muaPower.Contra.meanIR,'color',colorbrewer_setA_colorE,'LineWidth',2);
hold on
p2 = plot(data.muaPower.Whisk.meanIRtimeVec,data.muaPower.Whisk.meanIR,'color',colorbrewer_setA_colorD,'LineWidth',2);
p3 = plot(data.muaPower.Rest.meanIRtimeVec,data.muaPower.Rest.meanIR,'color',colorbrewer_setA_colorA,'LineWidth',2);
title('Mean IR Function')
xlabel('HRF Time (s)')
ylabel({'MUA [0.3-3 KHz] derived';'HRF amplitude (A.U.)'})
legend([p1,p2,p3],'Sensory-evoked','Volitional whisk','Awake rest','Location','NorthWest')
axis square
set(gca,'box','off')
% gamma derived
subplot(2,2,3);
s1 = scatter(xIndsA*1,data.gammaBandPower.Contra.Comb.AveR2,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.gammaBandPower.Contra.Comb.meanAveR2,data.gammaBandPower.Contra.Comb.stdAveR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
s2 = scatter(xIndsA*2,data.gammaBandPower.Contra.Comb.MedR2,'MarkerEdgeColor',colorbrewer_setA_colorE,'MarkerFaceColor',colorbrewer_setA_colorE,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.gammaBandPower.Contra.Comb.meanMedR2,data.gammaBandPower.Contra.Comb.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
s3 = scatter(xIndsA*3,data.gammaBandPower.Whisk.Comb.AveR2,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.gammaBandPower.Whisk.Comb.meanAveR2,data.gammaBandPower.Whisk.Comb.stdAveR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
s4 = scatter(xIndsA*4,data.gammaBandPower.Whisk.Comb.MedR2,'MarkerEdgeColor',colorbrewer_setA_colorD,'MarkerFaceColor',colorbrewer_setA_colorD,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.gammaBandPower.Whisk.Comb.meanMedR2,data.gammaBandPower.Whisk.Comb.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
s5 = scatter(xIndsA*5,data.gammaBandPower.Rest.Comb.MedR2,'MarkerEdgeColor',colorbrewer_setA_colorA,'MarkerFaceColor',colorbrewer_setA_colorA,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.gammaBandPower.Rest.Comb.meanMedR2,data.gammaBandPower.Rest.Comb.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
title({'Gamma-band [30-100 Hz] derived';'Gamma function kernel R^2 predictions'})
ylabel('R^2')
legend([s1,s2,s3,s4,s5],'Stimulus-evoked AvgData','Stimulus-evoked MedData','Volition whisk AvgData','Volitional whisk MedData','Awake Rest MedData','Location','NorthEast')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 6])
set(gca,'box','off')
% mua derived
subplot(2,2,4);
s1 = scatter(xIndsA*1,data.muaPower.Contra.Comb.AveR2,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorE,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.muaPower.Contra.Comb.meanAveR2,data.muaPower.Contra.Comb.stdAveR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
s2 = scatter(xIndsA*2,data.muaPower.Contra.Comb.MedR2,'MarkerEdgeColor',colorbrewer_setA_colorE,'MarkerFaceColor',colorbrewer_setA_colorE,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.muaPower.Contra.Comb.meanMedR2,data.muaPower.Contra.Comb.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
s3 = scatter(xIndsA*3,data.muaPower.Whisk.Comb.AveR2,'MarkerEdgeColor','k','MarkerFaceColor',colorbrewer_setA_colorD,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.muaPower.Whisk.Comb.meanAveR2,data.muaPower.Whisk.Comb.stdAveR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
s4 = scatter(xIndsA*4,data.muaPower.Whisk.Comb.MedR2,'MarkerEdgeColor',colorbrewer_setA_colorD,'MarkerFaceColor',colorbrewer_setA_colorD,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.muaPower.Whisk.Comb.meanMedR2,data.muaPower.Whisk.Comb.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
s5 = scatter(xIndsA*5,data.muaPower.Rest.Comb.MedR2,'MarkerEdgeColor',colorbrewer_setA_colorA,'MarkerFaceColor',colorbrewer_setA_colorA,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.muaPower.Rest.Comb.meanMedR2,data.muaPower.Rest.Comb.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
title({'MUA [0.3-3 kHz] derived';'Gamma function kernel R^2 predictions'})
ylabel('R^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 6])
set(gca,'box','off')

% save figure(s)
dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\';
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Kernel Predictions']);
