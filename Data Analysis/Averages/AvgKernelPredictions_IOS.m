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
neuralBands = {'gammaBandPower','muaPower'};
hemDataTypes = {'adjLH','adjRH'};
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
    for z = 1:length(neuralBands)
        neuralBand = neuralBands{1,z};
        data.(neuralBand).adjLH.kernelT{a,1} = AnalysisResults.HRFs.(neuralBand).adjLH.kernelT;
        data.(neuralBand).adjRH.kernelT{a,1} = AnalysisResults.HRFs.(neuralBand).adjRH.kernelT;
        data.(neuralBand).adjLH.bestKernel{a,1} = AnalysisResults.HRFs.(neuralBand).adjLH.bestKernel;
        data.(neuralBand).adjRH.bestKernel{a,1} = AnalysisResults.HRFs.(neuralBand).adjRH.bestKernel;
    end
end

%% summary figure
summaryFigure = figure;
sgtitle('Hemodynamic Response Functions')
% Gamma-band derived kernels
subplot(1,2,1);
for a = 1:length(data.gammaBandPower.adjLH.bestKernel)
    plot(data.gammaBandPower.adjLH.kernelT{a,1},data.gammaBandPower.adjLH.bestKernel{a,1})
    hold on
    plot(data.gammaBandPower.adjRH.kernelT{a,1},data.gammaBandPower.adjRH.bestKernel{a,1})
end
title({'Gamma-band power [30-100 Hz]';'derived kernels'})
xlabel('Time (s)')
ylabel('A.U.')
axis square
set(gca,'box','off')
axis tight

% MUA derived kernels
subplot(1,2,2);
for a = 1:length(data.muaPower.adjLH.bestKernel)
    plot(data.muaPower.adjLH.kernelT{a,1},data.muaPower.adjLH.bestKernel{a,1})
    hold on
    plot(data.muaPower.adjRH.kernelT{a,1},data.muaPower.adjRH.bestKernel{a,1})
end
title({'MUA power [0.3-3 kHz]';'derived kernels'})
xlabel('Time (s)')
ylabel('A.U.')
axis square
set(gca,'box','off')
axis tight







