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

[animal, hem, ~, ~, EventData, RestData, RestingBaselines, SleepData, ComparisonData] = LoadDataStructs();

infusionStatement = input('Is this an infusion trial? (y/n): ', 's'); disp(' ')
params.Infusion = infusionStatement;

targetMinutes = input('What is the target minute mark?: '); disp(' ')
params.targetMinutes = targetMinutes;

params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;

dataTypes = {'LH', 'RH'};

windowCamFileStruct = dir('*_WindowCam.bin');
windowCamFiles = {windowCamFileStruct.name}';
windowCamFiles = char(windowCamFiles);

procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFiles = char(procDataFiles);

rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFiles = char(rawDataFiles);

disp('Block [0] structs loaded.'); disp(' ')

%% BLOCK PURPOSE: [1] Stimulus and whisking evoked averages
disp('Analyzing Block [1] Analyzing the stimulus and whisking-evoked responses.'); disp(' ')
for dT = 1:length(dataTypes)
    dataType = dataTypes{dT};
    [ComparisonData] = AnalyzeEvokedResponses(animal, dataType, params, RestData, EventData, SpectrogramData, ComparisonData);
end 

%% BLOCK PURPOSE: [2] Cross correlation
disp('Analyzing Block [2] Analzying the cross-correlation between CBV and LFP.'); disp(' ')
for dT = 1:length(dataTypes)
    CBVdataType = ([dataTypes{dT} '_Electrode']);
    neuralDataType = dataTypes{dT};
    [ComparisonData] = AnalyzeXCorr(animal, CBVdataType, neuralDataType, params, RestData, RestingBaselines, SpectrogramData, SleepData, ComparisonData);
end

%% BLOCK PURPOSE: [3] Coherence
disp('Analyzing Block [3] Analyzing the coherence between L/R CBV and Gamma signals.'); disp(' ')
[ComparisonData] = AnalyzeCoherence(animal, params, procDataFiles, RestingBaselines, RestData, SleepData, ComparisonData);

%% BLOCK PURPOSE: [4] Power Spectrum
disp('Analyzing Block [4] Analyzing the power spectrums for CBV and Gamma signals.'); disp(' ')
[ComparisonData] = AnalyzePowerSpectrum(animal, params, procDataFiles, RestingBaselines, RestData, SleepData, ComparisonData);

%% BLOCK PURPOSE: [4] Analyze mean CBV and STD/variance
for dT = 1:length(dataTypes)
    dataType = dataTypes{dT};
    [ComparisonData] = AnalyzeCBV_STD(animal, dataType, params, procDataFiles, RestingBaselines, RestData, SleepData, ComparisonData);
end

%% BLOCK PURPOSE: [5] Correlation Coefficients between behaviors
[ComparisonData] = AnalyzeCorrCoeffs(animal, params, procDataFiles, RestingBaselines, RestData, SleepData, ComparisonData);

%% BLOCK PURPOSE: [6] Hemodynamic response functions
for dT = 1:length(dataTypes)
    dataType = dataTypes{dT};
    disp(['Generating ' num2str(dataType) ' HRF...']); disp(' ')
    [ComparisonData] = AnalyzeAwakeHRF(targetMinutes, dataType, procDataFiles, RestingBaselines, ComparisonData);
end

[ComparisonData] = GenerateHRFTable(procDataFiles, RestingBaselines, SleepData, ComparisonData);

disp('Data Analysis - Complete.'); disp(' ')

