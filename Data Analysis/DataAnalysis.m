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


% %% BLOCK PURPOSE: [1] Stimulus and whisking evoked averages
% % disp('Analyzing Block [1] Analyzing the stimulus and whisking-evoked responses.'); disp(' ')
% % for dT = 1:length(dataTypes)
% %     dataType = dataTypes{dT};
% %     [ComparisonData] = AnalyzeEvokedResponses(animal, dataType, params, RestData, EventData, SpectrogramData, ComparisonData);
% % end 
% 
%% BLOCK PURPOSE: [2] Cross correlation
% disp('Analyzing Block [2] Analzying the cross-correlation between CBV and LFP.'); disp(' ')
xcorr_CBVdataTypes = {'LH', 'RH', 'LH_Electrode', 'RH_Electrode'};
xcorr_neuralDataTypes = {'LH', 'RH', 'LH', 'RH'};
params.targetMinutes = 30;

params.minTime.Rest = 10;
params.minTime.NREM = 60;
params.minTime.REM = 60;

for dT = 1:length(xcorr_CBVdataTypes)
    xcorr_CBVdataType = xcorr_CBVdataTypes{dT};
    xcorr_neuralDataType = xcorr_CBVdataTypes{dT};
    [AnalysisResults] = AnalyzeXCorr(xcorr_CBVdataType, xcorr_neuralDataType, params, AnalysisResults);
end

%% BLOCK PURPOSE: [3] Coherence
% disp('Analyzing Block [3] Analyzing the coherence between L/R CBV and Gamma signals.'); disp(' ')
% [ComparisonData] = AnalyzeCoherence(animal, params, procDataFiles, RestingBaselines, RestData, SleepData, ComparisonData);

%% BLOCK PURPOSE: [4] Power Spectrum
% disp('Analyzing Block [4] Analyzing the power spectrums for CBV and Gamma signals.'); disp(' ')
% [ComparisonData] = AnalyzePowerSpectrum(animal, params, procDataFiles, RestingBaselines, RestData, SleepData, ComparisonData);

%% BLOCK PURPOSE: [4] Analyze mean CBV and STD/variance
% for dT = 1:length(dataTypes)
%     dataType = dataTypes{dT};
%     [ComparisonData] = AnalyzeCBV_STD(animal, dataType, params, procDataFiles, RestingBaselines, RestData, SleepData, ComparisonData);
% end

%% BLOCK PURPOSE: [5] Correlation Coefficients between behaviors
[AnalysisResults] = AnalyzeCorrCoeffs(animal, params, procDataFiles, RestingBaselines, RestData, SleepData, AnalysisResults);

%% BLOCK PURPOSE: [6] Hemodynamic response functions
% for dT = 1:length(dataTypes)
%     dataType = dataTypes{dT};
%     disp(['Generating ' num2str(dataType) ' HRF...']); disp(' ')
%     [ComparisonData] = AnalyzeAwakeHRF(targetMinutes, dataType, procDataFiles, RestingBaselines, ComparisonData);
% end

% [ComparisonData] = GenerateHRFTable(procDataFiles, RestingBaselines, SleepData, ComparisonData);

disp('Data Analysis - Complete.'); disp(' ')

