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
    AnalyzeBehavioralHbT_Handler(rootFolder,delim,'n')
    AnalyzeBilateralCoherence_Handler(rootFolder,delim,'n')
    AnalyzeNeuralHemoCoherence_Handler(rootFolder,delim,'n')
    AnalyzePowerSpectrum_Handler(rootFolder,delim,'n')
    AnalyzePearsonCorrelation_Handler(rootFolder,delim,'n')
    AnalyzeCrossCorrelation_Handler(rootFolder,delim,'n')
    AnalyzeEvokedResponses_Handler(rootFolder,delim,'n')
end
%% generate figures
disp('Loading analysis results and generating figures...'); disp(' ')
saveFigs = true;
% DiaphoraseCellCounts_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% WhiskEvoked_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% StimEvoked_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% StimEvoked2_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% StimEvoked3_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% Coherence_Saporin(rootFolder,saveFigs,delim);
% NeuralHemoCoherence_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% WhiskHemoCoherence_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% PowerSpec_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% PowerSpec2_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% PearsonsCorr_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% XCorr_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% MeanHbT_Saporin(rootFolder,saveFigs,delim,AnalysisResults); %#ok<*NASGU>

end
