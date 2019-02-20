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

[animal, hem, ~, ~, ~, ~, RestingBaselines, SpectrogramData, ~, ~] = LoadDataStructs();

dataTypes = {'LH', 'RH'};
windowCamFiles = ls('*_WindowCam.bin');
procDataFiles = ls('*_ProcData.mat');
rawDataFiles = ls('*_RawData.mat');

disp('Block [0] structs loaded.'); disp(' ')

%% BLOCK PURPOSE: [1] Create spectrograms for each file
disp('Analyzing Block [1] Analyzing the spectrogram for each file.'); disp(' ')
for neuralDT = 1:length(dataTypes)
    dataType = dataTypes{neuralDT};
    [SpectrogramData] = CreateTrialSpectrograms(animal, dataType, rawDataFiles, SpectrogramData);
end

% Find spectrogram baselines for each day
for neuralDT = 1:length(dataTypes)
    dataType = dataTypes{neuralDT};
    [RestingBaselines] = CalculateSpectrogramBaselines(animal, dataType, RestingBaselines, SpectrogramData);
end

% Normalize spectrogram by baseline
for neuralDT = 1:length(dataTypes)
    dataType = dataTypes{neuralDT};
    [SpectrogramData] = NormalizeSpectrograms(animal, dataType, RestingBaselines, SpectrogramData);
end

%% BLOCK PURPOSE: [2] Single Trial Checks
disp('Analyzing Block [2] Creating single trial figures for each 5 minute session.'); disp(' ')
for pDF = 1:size(procDataFiles, 1)
    procDataFile = procDataFiles(pDF, :);
    disp(['Creating single trial figure for file number ' num2str(pDF) ' of ' num2str(size(procDataFiles, 1)) '...']); disp(' ')
    CreateSingleTrialFigs(procDataFile, RestingBaselines, SpectrogramData)
    close all
end

%% BLOCK PURPOSE: [3] Sleep scoring
disp('Analyzing Block [3] Running Sleep scoring analysis...'); disp(' ')
electrodeInput = input('Which electrode(s) would you like to use to sleep score? (L, R, B): ', 's'); disp(' ')
[SleepData] = SleepScore(animal, hem, rawDataFiles, procDataFiles, electrodeInput, RestingBaselines, SpectrogramData);

disp('Stage five analysis - complete'); disp(' ')
