clear
clc

cd('G:\T48\Combined Imaging\');
load('T48_AnalysisAverages.mat');
sleepMean(:, 1) = AnalysisAverages.CorrCoeff.Sleep.mean;
sleepSTD(:, 1) = AnalysisAverages.CorrCoeff.Sleep.std;
restMean(:, 1) = AnalysisAverages.CorrCoeff.Rest.mean;
restSTD(:, 1) = AnalysisAverages.CorrCoeff.Rest.std;

cd('G:\T49\Combined Imaging\');
load('T49_AnalysisAverages.mat');
sleepMean(:, 2) = AnalysisAverages.CorrCoeff.Sleep.mean;
sleepSTD(:, 2) = AnalysisAverages.CorrCoeff.Sleep.std;
restMean(:, 2) = AnalysisAverages.CorrCoeff.Rest.mean;
restSTD(:, 2) = AnalysisAverages.CorrCoeff.Rest.std;

cd('E:\T52\Combined Imaging\');
load('52_AnalysisAverages.mat');
sleepMean(:, 3) = AnalysisAverages.CorrCoeff.Sleep.mean;
sleepSTD(:, 3) = AnalysisAverages.CorrCoeff.Sleep.std;
restMean(:, 3) = AnalysisAverages.CorrCoeff.Rest.mean;
restSTD(:, 3) = AnalysisAverages.CorrCoeff.Rest.std;

cd('E:\T61\Combined Imaging\');
load('T61_AnalysisAverages.mat');
sleepMean(:, 4) = AnalysisAverages.CorrCoeff.Sleep.mean;
sleepSTD(:, 4) = AnalysisAverages.CorrCoeff.Sleep.std;
restMean(:, 4) = AnalysisAverages.CorrCoeff.Rest.mean;
restSTD(:, 4) = AnalysisAverages.CorrCoeff.Rest.std;

cd('E:\T62\Combined Imaging\');
load('T62_AnalysisAverages.mat');
sleepMean(:, 5) = AnalysisAverages.CorrCoeff.Sleep.mean;
sleepSTD(:, 5) = AnalysisAverages.CorrCoeff.Sleep.std;
restMean(:, 5) = AnalysisAverages.CorrCoeff.Rest.mean;
restSTD(:, 5) = AnalysisAverages.CorrCoeff.Rest.std;

cd('E:\T64\Combined Imaging\');
load('T64_AnalysisAverages.mat');
sleepMean(:, 6) = AnalysisAverages.CorrCoeff.Sleep.mean;
sleepSTD(:, 6) = AnalysisAverages.CorrCoeff.Sleep.std;
restMean(:, 6) = AnalysisAverages.CorrCoeff.Rest.mean;
restSTD(:, 6) = AnalysisAverages.CorrCoeff.Rest.std;

figure;
subplot(1, 2, 1)
errorbar(1:length(restMean), restMean, restSTD, 'ok')
title('Correlation Coefficient during Rest')
ylabel('Correlation Coefficient')
xlabel('Animal number)')
xlim([0 7])
ylim([0 1])

subplot(1, 2, 2)
errorbar(1:length(sleepMean), sleepMean, sleepSTD, 'sk')
title('Correlation Coefficient during Rest')
ylabel('Correlation Coefficient')
xlabel('Animal number)')
xlim([0 7])
ylim([0 1])
