function [] = DataAnalysis_IOS()
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
% clc;
% clear;
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
disp('Analyzing Block [1] Analyzing the whisking-evoked and stimulus-evoked Hemodynamic/neural response.'); disp(' ')
evoked_dataTypes = {'LH','RH'};
params.targetMinutes = 30;
for dT = 1:length(evoked_dataTypes)
    evoked_dataType = evoked_dataTypes{dT};
    [AnalysisResults] = AnalyzeEvokedResponses_IOS(evoked_dataType,params,AnalysisResults);
end 
close all

%% BLOCK PURPOSE: [2] Cross correlation
disp('Analyzing Block [2] Analzying the cross-correlation between hemodynamics and neural data.'); disp(' ')
xcorr_CBVdataTypes = {'LH','RH','LH_Electrode','RH_Electrode'};
xcorr_neuralDataTypes = {'cortical_LH','cortical_RH','cortical_LH','cortical_RH'};
baselineType = 'manualSelection';
params.targetMinutes = 30;   % minutes
params.minTime.Rest = 10;   % seconds
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 30;   % seconds
[AnalysisResults] = AnalyzeXCorr_IOS(xcorr_CBVdataTypes,xcorr_neuralDataTypes,baselineType,params,AnalysisResults);
close all

%% BLOCK PURPOSE: [3] Coherence
disp('Analyzing Block [3] Analyzing the coherence between L/R CBV and neural-band signals.'); disp(' ')
coherr_dataTypes = {'CBV','CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','muaPower'};
baselineType = 'manualSelection';
params.targetMinutes = 30;   % minutes
params.minTime.Rest = 10;   % seconds
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 30;   % seconds
[AnalysisResults] = AnalyzeCoherence_IOS(coherr_dataTypes,baselineType,params, AnalysisResults);
close all

%% BLOCK PURPOSE: [4] Power Spectra
disp('Analyzing Block [4] Analyzing the power spectra of CBV and neural-band signals.'); disp(' ')
powerspec_dataTypes =  {'CBV','CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','muaPower'};
baselineType = 'manualSelection';
params.targetMinutes = 30;   % minutes
params.minTime.Rest = 10;   % seconds
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 30;   % seconds
[AnalysisResults] = AnalyzePowerSpectrum_IOS(powerspec_dataTypes,baselineType,params,AnalysisResults);
close all

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

