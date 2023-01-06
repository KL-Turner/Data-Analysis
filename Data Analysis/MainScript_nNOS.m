function [] = MainScript_nNOS()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panels for Turner et al. nNOS Manuscript
%
% Functions used to pre-process the original data are located in the folder "Pre-Processing Functions"
% Functions used to analyze data for figures are located in the folder "Data Analysis Functions"
% Functions optained from 3rd party are located in the folder "Shared Functions"
%________________________________________________________________________________________________________________________

clear; clc; close all;
% verify code repository and data are in the current directory/added path
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
% analysis subfunctions
runAnalysis = false;
if runAnalysis == true
    AnalyzeBehavioralHbT_Handler(rootFolder,delim,false)
    AnalyzeBilateralCoherence_Handler(rootFolder,delim,false)
    AnalyzeNeuralHemoCoherence_Handler(rootFolder,delim,false)
    AnalyzePowerSpectrum_Handler(rootFolder,delim,false)
    AnalyzePearsonCorrelation_Handler(rootFolder,delim,false)
    AnalyzePowerSpectrum_LFP_Handler(rootFolder,delim,false)
    AnalyzeCrossCorrelation_Handler(rootFolder,delim,false)
    AnalyzeEvokedResponses_Handler(rootFolder,delim,false)
    AnalyzeArousalTransitions_Handler(rootFolder,delim,false)
    multiWaitbar('close all')
end
%% generate figures
disp('Loading analysis results and generating figures...'); disp(' ')
saveFigs = true;
DiaphoraseCellCounts_Figures(rootFolder,saveFigs,delim);

% BilateralCoherence_Figures(rootFolder,saveFigs,delim);
% DiaphoraseCellCounts_Bilateral_IOS(rootFolder,saveFigs,delim);
% WhiskEvoked_Bilateral_IOS(rootFolder,saveFigs,delim);
% WhiskingBehavior_Bilateral_IOS(rootFolder,saveFigs,delim);
% StimEvoked_Bilateral_IOS(rootFolder,saveFigs,delim);
% StimEvoked_PulseTrain_IOS(rootFolder,saveFigs,delim);
% StimEvoked_PulseTrain_2PLSM(rootFolder,saveFigs,delim);
% BaselineShift_2PLSM(rootFolder,saveFigs,delim);
% CrossCorrelation_Bilateral_IOS(rootFolder,saveFigs,delim);
% NeuralHemoCoherence_Bilateral_IOS(rootFolder,saveFigs,delim);
% PowerSpectrum_Bilateral_IOS(rootFolder,saveFigs,delim);
% ArousalTransitions_Bilateral_IOS(rootFolder,saveFigs,delim);
% PowerSpectrumLFP_Bilateral_IOS(rootFolder,saveFigs,delim);
% PearsonsCorr_Bilateral_IOS(rootFolder,saveFigs,delim);
% ArousalStateHemodynamics_Bilateral_IOS(rootFolder,saveFigs,delim);
% SleepAmounts_Bilateral_IOS(rootFolder,saveFigs,delim)

end
