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
close all
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')

% Load or create the AnalysisResults.mat structure into the Workspace
resultsDataFileStruct = dir('*_AnalysisResults.mat');
resultsDataFile = {resultsDataFileStruct.name}';
resultsDataFileID = char(resultsDataFile);
if exist(resultsDataFileID)
    load(resultsDataFileID)
else
    AnalysisResults = [];
end

%% BLOCK PURPOSE: [1] Stimulus and whisking evoked averages
disp('Analyzing Block [1] Analyzing the whisking-evoked and stimulus-evoked CBV/neural response.'); disp(' ')
evoked_dataTypes = {'LH', 'RH'};
params.targetMinutes = 30;
for dT = 1:length(evoked_dataTypes)
    evoked_dataType = evoked_dataTypes{dT};
    [AnalysisResults] = AnalyzeEvokedResponses_IOS(evoked_dataType, params, AnalysisResults);
end 

%% BLOCK PURPOSE: [2] Cross correlation
disp('Analyzing Block [2] Analzying the cross-correlation between CBV and LFP.'); disp(' ')
xcorr_CBVdataTypes = {'LH', 'RH', 'LH_Electrode', 'RH_Electrode'};
xcorr_neuralDataTypes = {'cortical_LH', 'cortical_RH', 'cortical_LH', 'cortical_RH'};
baselineType = 'manualSelection';
params.targetMinutes = 30;
params.minTime.Rest = 10;
params.minTime.NREM = 55;
params.minTime.REM = 55;

for a = 1:length(xcorr_CBVdataTypes)
    xcorr_CBVdataType = xcorr_CBVdataTypes{a};
    xcorr_neuralDataType = xcorr_neuralDataTypes{a};
    [AnalysisResults] = AnalyzeXCorr_IOS(xcorr_CBVdataType, xcorr_neuralDataType, baselineType, params, AnalysisResults);
end

%% BLOCK PURPOSE: [3] Coherence
disp('Analyzing Block [3] Analyzing the coherence between L/R CBV and neural-band signals.'); disp(' ')
coherr_dataTypes = {'CBV', 'deltaBandPower', 'thetaBandPower', 'alphaBandPower', 'betaBandPower', 'gammaBandPower'};
params.targetMinutes = 30;
params.minTime.Rest = 10;
params.minTime.NREM = 55;
params.minTime.REM = 55;

for b = 1:length(coherr_dataTypes)
    coherr_dataType = coherr_dataTypes{b};
[AnalysisResults] = AnalyzeCoherence_IOS(coherr_dataType, params, AnalysisResults);
end

%% BLOCK PURPOSE: [4] Power Spectra
disp('Analyzing Block [3] Analyzing the power spectra of CBV and neural-band signals.'); disp(' ')
powerspec_dataTypes = {'CBV', 'deltaBandPower', 'thetaBandPower', 'alphaBandPower', 'betaBandPower', 'gammaBandPower'};
params.targetMinutes = 30;
params.minTime.Rest = 10;
params.minTime.NREM = 55;
params.minTime.REM = 55;

for b = 1:length(powerspec_dataTypes)
    powerspec_dataType = powerspec_dataTypes{b};
    [AnalysisResults] = AnalyzePowerSpectrum_IOS(powerspec_dataType, params, AnalysisResults);
end

%% BLOCK PURPOSE: [4] Analyze mean CBV and STD/variance
% for dT = 1:length(dataTypes)
%     dataType = dataTypes{dT};
%     [ComparisonData] = AnalyzeCBV_STD(animal, dataType, params, procDataFiles, RestingBaselines, RestData, SleepData, ComparisonData);
% end

% %% BLOCK PURPOSE: [5] Correlation Coefficients between behaviors
% [AnalysisResults] = AnalyzeCorrCoeffs(animal, params, procDataFiles, RestingBaselines, RestData, SleepData, AnalysisResults);

%% BLOCK PURPOSE: [6] Hemodynamic response functions
% for dT = 1:length(dataTypes)
%     dataType = dataTypes{dT};
%     disp(['Generating ' num2str(dataType) ' HRF...']); disp(' ')
%     [ComparisonData] = AnalyzeAwakeHRF(targetMinutes, dataType, procDataFiles, RestingBaselines, ComparisonData);
% end

% [ComparisonData] = GenerateHRFTable(procDataFiles, RestingBaselines, SleepData, ComparisonData);

disp('Data Analysis - Complete.'); disp(' ')

