function [] = DataAnalysis_IOS()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Step through various different data analysis and save the results in a summary structure.
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
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

% %% BLOCK PURPOSE: [1] Stimulus and whisking evoked averages
% disp('Analyzing Block [1] Analyzing the whisking-evoked and stimulus-evoked hemodynamic and neural responses.'); disp(' ')
% evoked_dataTypes = {'adjLH','adjRH'};
% [AnalysisResults] = AnalyzeEvokedResponses_IOS(evoked_dataTypes,AnalysisResults);
% 
% %% BLOCK PURPOSE: [2] Cross correlation
% disp('Analyzing Block [2] Analzying the cross-correlation between hemodynamics and neural data.'); disp(' ')
% xcorr_dataTypes = {'adjLH','adjRH'};
% params.minTime.Rest = 10;   % seconds
% params.minTime.NREM = 30;   % seconds
% params.minTime.REM = 30;   % seconds
% [AnalysisResults] = AnalyzeXCorr_IOS(xcorr_dataTypes,params,AnalysisResults);
%  
% %% BLOCK PURPOSE: [3] Coherence
% disp('Analyzing Block [3] Analyzing the coherence between L/R hemodynamic and neural data.'); disp(' ')
% coherr_dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
% params.minTime.Rest = 10;   % seconds
% params.minTime.NREM = 30;   % seconds
% params.minTime.REM = 30;   % seconds
% [AnalysisResults] = AnalyzeCoherence_IOS(coherr_dataTypes,params, AnalysisResults);
% 
% %% BLOCK PURPOSE: [4] Power Spectra
% disp('Analyzing Block [4] Analyzing the power spectra of hemodynamic and neural data.'); disp(' ')
% powerspec_dataTypes =  {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
% params.minTime.Rest = 10;   % seconds
% params.minTime.NREM = 30;   % seconds
% params.minTime.REM = 30;   % seconds
% [AnalysisResults] = AnalyzePowerSpectrum_IOS(powerspec_dataTypes,params,AnalysisResults);
% 
% %% BLOCK PURPOSE: [5] Pearson's correlation coefficient
% disp('Analyzing Block [5] Analyzing the Pearson''s correlation coefficient between bilateral hemodynamic and neural data.'); disp(' ')
% corrCoeff_dataTypes =  {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
% params.minTime.Rest = 10;   % seconds
% [AnalysisResults] = AnalyzeCorrCoeffs_IOS(corrCoeff_dataTypes,params,AnalysisResults);

% %% BLOCK PURPOSE: [6] Mean CBV values
% disp('Analyzing Block [6] Analyzing the mean CBV during different behaviors.'); disp(' ')
% params.minTime.Rest = 10;   % seconds
% [AnalysisResults] = AnalyzeMeanCBV_IOS(params,AnalysisResults);
% 
% %% BLOCK PURPOSE: [7] Mean heart rate values
% disp('Analyzing Block [7] Analyzing the mean heart rate during different behaviors.'); disp(' ')
% params.minTime.Rest = 10;   % seconds
% [AnalysisResults] = AnalyzeMeanHeartRate_IOS(params,AnalysisResults);

%% BLOCK PURPOSE: [8] Hemodynamic response functions
disp('Analyzing Block [8] Analyzing the hemodynamic response function and predictability of awake data.'); disp(' ')
HRF_hemDataTypes = {'adjLH','adjRH'};
HRF_neuralBands = {'gammaBandPower','muaPower'};
behaviors = {'Contra','Whisk','Rest'};
for a = 1:length(HRF_hemDataTypes)
    hemDataType = HRF_hemDataTypes{1,a};
    for b = 1:length(HRF_neuralBands)
        neuralBand = HRF_neuralBands{1,b};
        for c = 1:length(behaviors)
            behavior = behaviors{1,c};
            [AnalysisResults] = CalculateHRFDeconvolution_IOS(neuralBand,hemDataType,behavior,AnalysisResults);
            [AnalysisResults] = EvaluateCBVPredictionAccuracy_IOS(neuralBand,hemDataType,behavior,AnalysisResults);
        end
    end
end

disp('Data Analysis - Complete.'); disp(' ')

