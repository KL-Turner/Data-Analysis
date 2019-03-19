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
[ComparisonData] = AnalyzeEvokedResponses_2P(animalID, RestingBaselines, EventData, SpectrogramData, ComparisonData);

%% BLOCK PURPOSE: [2] Power spectrum analysis
[ComparisonData] = AnalyzePowerSpectrum_2P(animalID, mergedDataFiles, ComparisonData);

%% BLOCK PURPOSE: [3] Coherence between vessel diameter and whisk acceleration
[ComparisonData] = AnalyzeCoherence_2P(animalID, mergedDataFiles, ComparisonData);
