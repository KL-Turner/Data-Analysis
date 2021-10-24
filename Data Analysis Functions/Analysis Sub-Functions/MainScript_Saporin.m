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
%% verify code repository and data are in the current directory/added path
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
%% analysis subfunctions
runAnalysis = true;
if runAnalysis == true
%     AnalyzeBehavioralHbT_Handler(rootFolder,delim,false)
%     AnalyzeBilateralCoherence_Handler(rootFolder,delim,false)
%     AnalyzeNeuralHemoCoherence_Handler(rootFolder,delim,false)
%     AnalyzePowerSpectrum_Handler(rootFolder,delim,false)
%     AnalyzePearsonCorrelation_Handler(rootFolder,delim,false)
%     AnalyzeCrossCorrelation_Handler(rootFolder,delim,false)
%     AnalyzeEvokedResponsesA_Handler(rootFolder,delim,false)
%     AnalyzeEvokedResponsesB_Handler(rootFolder,delim,false)
%     AnalyzeVesselEvokedResponses_Handler(rootFolder,delim,false)
%     AnalyzeArousalTransitions_Handler(rootFolder,delim,false)
%     AnalyzeArousalStateProbability_Handler(rootFolder,delim,false)
%     AnalyzeWhiskingBehavior_Handler(rootFolder,delim,false)
%     AnalyzePowerSpectrum_LFP_Handler(rootFolder,delim,false)
%     AnalyzeVesselBaselineShift_Handler(rootFolder,delim,true)
end
%% generate figures
disp('Loading analysis results and generating figures...'); disp(' ')
saveFigs = true;
% DiaphoraseCellCounts_Bilateral_IOS(rootFolder,saveFigs,delim);
% WhiskEvoked_Bilateral_IOS(rootFolder,saveFigs,delim);
% WhiskingBehavior_Bilateral_IOS(rootFolder,saveFigs,delim);
% StimEvoked_Bilateral_IOS(rootFolder,saveFigs,delim);
% StimEvoked_PulseTrain_IOS(rootFolder,saveFigs,delim);
% StimEvoked_PulseTrain_2PLSM(rootFolder,saveFigs,delim);
BaselineShift_2PLSM(rootFolder,saveFigs,delim);
% CrossCorrelation_Bilateral_IOS(rootFolder,saveFigs,delim);
% BilateralCoherence_Bilateral_IOS(rootFolder,saveFigs,delim);
% NeuralHemoCoherence_Bilateral_IOS(rootFolder,saveFigs,delim);
% PowerSpectrum_Bilateral_IOS(rootFolder,saveFigs,delim);
% ArousalTransitions_Bilateral_IOS(rootFolder,saveFigs,delim);
% PowerSpectrumLFP_Bilateral_IOS(rootFolder,saveFigs,delim);
% PearsonsCorr_Bilateral_IOS(rootFolder,saveFigs,delim);
% ArousalStateHemodynamics_Bilateral_IOS(rootFolder,saveFigs,delim);
% SleepAmounts_Bilateral_IOS(rootFolder,saveFigs,delim)

end
