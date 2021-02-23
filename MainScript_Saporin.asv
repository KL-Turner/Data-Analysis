function [] = MainScript_Saporin()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
% Purpose: Generates KLT's main and supplemental figs for Turner et al. eLife2020
%
% Scripts used to pre-process the original data are located in the folder "Pre-Processing Scripts".
% Functions that are used in both the analysis and pre-processing are located in the analysis folder.
%________________________________________________________________________________________________________________________

clear; clc; close all;
%% make sure the code repository and data are present in the current directory
currentFolder = pwd;
addpath(genpath(currentFolder));
fileparts = strsplit(currentFolder,filesep);
if ismac
    rootFolder = fullfile(filesep,fileparts{1:end});
    delim = '/';
else
    rootFolder = fullfile(fileparts{1:end});
    delim = '\';
end
% add root folder to Matlab's working directory
addpath(genpath(rootFolder))
%% run the data analysis. The progress bars will show the analysis progress
rerunAnalysis = 'n';
saveFigs = 'y';
if exist('AnalysisResults.mat','file') ~= 2 || strcmp(rerunAnalysis,'y') == true
    multiWaitbar('Analyzing behavioral hemodynamics',0,'Color','R'); pause(0.25);
    multiWaitbar('Analyzing coherence',0,'Color','Y'); pause(0.25);
    multiWaitbar('Analyzing neural-hemo coherence',0,'Color','R'); pause(0.25);
    multiWaitbar('Analyzing power spectra',0,'Color','R'); pause(0.25);
    multiWaitbar('Analyzing Pearson''s correlation coefficients',0,'Color','Y'); pause(0.25);
    multiWaitbar('Analyzing cross correlation',0,'Color','R'); pause(0.25);
    multiWaitbar('Analyzing evoked responses',0,'Color','Y'); pause(0.25);
    % run analysis and output a structure containing all the analyzed data
    [AnalysisResults] = AnalyzeData(rootFolder);
    multiWaitbar('CloseAll');
else
    disp('Loading analysis results and generating figures...'); disp(' ')
    load('AnalysisResults.mat')
end
%% generate figures 
[AnalysisResults] = WhiskEvoked_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = StimEvoked_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Coherence_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = PowerSpec_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = PearsonsCorr_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = XCorr_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = MeanHbT_Saporin(rootFolder,saveFigs,delim,AnalysisResults); %#ok<*NASGU> 
%% fin.
disp('MainScript Analysis - Complete'); disp(' ')
end

function [AnalysisResults] = AnalyzeData(rootFolder)
% IOS animal IDs
animalIDs = {'T135','T141','T142','T144','T151','T155','T156','T157','T159'};
saveFigs = 'y';
if exist('AnalysisResults.mat','file') == 2
    load('AnalysisResults.mat')
else
    AnalysisResults = [];
end
%% Analyze the hemodynamic signal [HbT] during different arousal states (IOS)
runFromStart = 'n';
for aa = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,aa})) == false || isfield(AnalysisResults.(animalIDs{1,aa}),'MeanCBV') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeMeanCBV(animalIDs{1,aa},rootFolder,AnalysisResults);
    end
    multiWaitbar('Analyzing behavioral hemodynamics','Value',aa/length(animalIDs));
end
%% Analyze the spectral coherence between bilateral hemodynamic [HbT] and neural signals (IOS)
runFromStart = 'n';
for bb = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'Coherence') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeCoherence(animalIDs{1,bb},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar('Analyzing coherence','Value',bb/length(animalIDs));
end
%% Analyze the spectral coherence between neural-hemodynamic [HbT] signals (IOS)
runFromStart = 'n';
for cc = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,cc})) == false || isfield(AnalysisResults.(animalIDs{1,cc}),'NeuralHemoCoherence') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeNeuralHemoCoherence(animalIDs{1,cc},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar('Analyzing neural-hemo coherence','Value',cc/length(animalIDs));
end
%% Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
runFromStart = 'n';
for dd = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,dd})) == false || isfield(AnalysisResults.(animalIDs{1,dd}),'PowerSpectra') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzePowerSpectrum(animalIDs{1,dd},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar('Analyzing power spectra','Value',dd/length(animalIDs));
end
%% Analyze Pearson's correlation coefficient between bilateral hemodynamic [HbT] and neural signals (IOS)
runFromStart = 'n';
for ee = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,ee})) == false || isfield(AnalysisResults.(animalIDs{1,ee}),'CorrCoeff') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeCorrCoeffs(animalIDs{1,ee},rootFolder,AnalysisResults);
    end
    multiWaitbar('Analyzing Pearson''s correlation coefficients','Value',ee/length(animalIDs));
end
%% Analyze the cross-correlation between neural activity and hemodynamics [HbT] (IOS)
runFromStart = 'n';
for ff = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,ff})) == false || isfield(AnalysisResults.(animalIDs{1,ff}),'XCorr') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeXCorr(animalIDs{1,ff},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar('Analyzing cross correlation','Value',ff/length(animalIDs));
end
%% Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
runFromStart = 'n';
for gg = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,gg})) == false || isfield(AnalysisResults.(animalIDs{1,gg}),'EvokedAvgs') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeEvokedResponses(animalIDs{1,gg},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar('Analyzing evoked responses','Value',gg/length(animalIDs));
end
%% fin.
disp('Loading analysis results and generating figures...'); disp(' ')
end
