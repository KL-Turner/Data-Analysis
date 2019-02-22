%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: 1) Analyze the spectrogram for each file and normalizing by the resting baseline.
%            2) Creating single trial figures for each 5 minute session.
%            3) Run Sleep scoring analysis.
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

[animalID, ~, ~, ~, RestData, RestingBaselines, SpectrogramData, ~, ~] = LoadDataStructs_2P;

mergedDirectory = dir('*_MergedData.mat');
mergedDataFiles = {mergedDirectory.name}';
mergedDataFiles = char(mergedDataFiles);

disp('Block [0] structs loaded.'); disp(' ')

%% BLOCK PURPOSE: [1] Create spectrograms for each file
disp('Analyzing Block [1] Analyzing the spectrogram for each file and normalizing by the resting baseline.'); disp(' ')
[SpectrogramData] = CreateTrialSpectrograms_2P(animalID, mergedDataFiles, SpectrogramData);

% Find spectrogram baselines for each day
[RestingBaselines] = CalculateSpectrogramBaselines_2P(animalID, RestingBaselines, SpectrogramData);

% Normalize spectrogram by baseline
[SpectrogramData] = NormalizeSpectrograms_2P(animalID, RestingBaselines, SpectrogramData);

%% BLOCK PURPOSE: [2] Single Trial Checks
disp('Analyzing Block [2] Creating single trial figures for each 5 minute session.'); disp(' ')
CreateSingleTrialFigs_2P(mergedDataFiles, RestingBaselines, SpectrogramData)


% %% BLOCK PURPOSE: [3] Sleep scoring
% disp('Analyzing Block [3] Running Sleep scoring analysis.'); disp(' ')
% electrodeInput = input('Which electrode(s) would you like to use to sleep score? (L, R, B): ', 's'); disp(' ')
% [SleepData] = SleepScore(animal, hem, rawDataFiles, procDataFiles, electrodeInput, RestingBaselines, SpectrogramData);

disp('Stage five analysis - complete'); disp(' ')
