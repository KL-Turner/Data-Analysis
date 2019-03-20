%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
clc;
clear;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')

[animalID, ~, ~, EventData, RestData, RestingBaselines, SpectrogramData, ~, ComparisonData] = LoadDataStructs_2P;

mergedDirectory = dir('*_MergedData.mat');
mergedDataFiles = {mergedDirectory.name}';
mergedDataFiles = char(mergedDataFiles);

disp('Block [0] structs loaded.'); disp(' ')

%% BLOCK PURPOSE: [1] Whisking evoked averages
disp('Analyzing Block [1] Analyzing the whisking-evoked responses.'); disp(' ')
[ComparisonData] = AnalyzeEvokedResponses_2P(animalID,mergedDataFiles, RestingBaselines, EventData, SpectrogramData, ComparisonData);

%% BLOCK PURPOSE: [2] Power spectrum analysis
disp('Analyzing Block [2] Analyzing the vessel and whisking power spectra.'); disp(' ')
[ComparisonData] = AnalyzePowerSpectrum_2P(animalID, mergedDataFiles, ComparisonData);

%% BLOCK PURPOSE: [3] Coherence between vessel diameter and whisk acceleration
disp('Analyzing Block [3] Analyzing the coherence between vessel diameter and whisking acceleration.'); disp(' ')
[ComparisonData] = AnalyzeCoherence_2P(animalID, mergedDataFiles, ComparisonData);

%% BLOCK PURPOSE: [4] Cross correlation between whisk acceleration and vessel diameter
disp('Analyzing Block [4] Analyzing the cross correlation between vessel diameter and whisker acceleration.'); disp(' ')
[ComparisonData] = AnalyzeXCorr_2P(animalID, mergedDataFiles, ComparisonData);

disp('Two Photon Data Analysis - Complete.'); disp(' ')
