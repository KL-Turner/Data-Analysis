function [] = MainScript_nNOS()
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panels for Turner et al. nNOS Manuscript
%
% Functions used to pre-process the original data are located in the folder "Pre-Processing Functions"
% Functions used to analyze data for figures are located in the folder "Data Analysis Functions"
% Functions optained from 3rd party are located in the folder "Shared Functions"
%----------------------------------------------------------------------------------------------------------
% zap;
cd('F:\NO Project\')
% save figures?
saveFigs = true;
% verify code repository and data are in the current directory/added path
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
% analysis subfunctions
addpath(genpath(rootFolder))
runAnalysis = false;
if runAnalysis == true
    %% evoked responses
    AnalyzeEvokedResponses_Ephys_Handler(rootFolder,delim,false)
    AnalyzeEvokedResponses_GCaMP_Handler(rootFolder,delim,false)
    AnalyzeEvokedResponses_2P_Handler(rootFolder,delim,false)
    AnalyzeEvokedResponses_Pulse_Handler(rootFolder,delim,false)

    %% HbT, HbO, HbR, GCaMP
    AnalyzeIntrinsicSignals_Ephys_Handler(rootFolder,delim,false)
    AnalyzeIntrinsicSignals_GCaMP_Handler(rootFolder,delim,false)
    AnalyzeHemoGFPRelationship_EGFP_Handler(rootFolder,delim,false)
    AnalyzeArousalDerivative_Ephys_Handler(rootFolder,delim,false)

    %% arteriole diameter and baseline shift
    AnalyzeArterioleBaseline_2P_Handler(rootFolder,delim,false)
    AnalyzeArterioleDiameter_2P_Handler(rootFolder,delim,false)

    %% power spectral density
    AnalyzePowerSpectrum_Ephys_Handler(rootFolder,delim,false)
    AnalyzePowerSpectrum_GCaMP_Handler(rootFolder,delim,false)
    AnalyzePowerSpectrum_LFP_Handler(rootFolder,delim,false)
    AnalyzePowerSpectrum_2P_Handler(rootFolder,delim,false)

    %% coherence between hemipheres
    AnalyzeBilateralCoherence_Ephys_Handler(rootFolder,delim,false)
    AnalyzeBilateralCoherence_GCaMP_Handler(rootFolder,delim,false)

    %% coherence between neural-hemo
    AnalyzeNeuralHemoCoherence_Ephys_Handler(rootFolder,delim,false)
    AnalyzeNeuralHemoCoherence_GCaMP_Handler(rootFolder,delim,false)

    %% Pearson's correlation between hemispheres
    AnalyzePearsonCorrelation_Ephys_Handler(rootFolder,delim,false)
    AnalyzePearsonCorrelation_GCaMP_Handler(rootFolder,delim,false)

    %% cross correlation between neural-hemo
    AnalyzeCrossCorrelation_Ephys_Handler(rootFolder,delim,false)
    AnalyzeCrossCorrelation_GCaMP_Handler(rootFolder,delim,false)
    AnalyzeCrossCorrelation_EGFP_Handler(rootFolder,delim,false)

    %% pupil analysis
    AnalyzePupilArea_Ephys_Handler(rootFolder,delim,false)
    AnalyzePupilEvokedResponses_Ephys_Handler(rootFolder,delim,false)
    AnalyzePupilCoherence_Ephys_Handler(rootFolder,delim,false)
    AnalyzePupilCrossCorrelation_Ephys_Handler(rootFolder,delim,false)
    AnalyzePupilInterBlinkInterval_Ephys_Handler(rootFolder,delim,false)

    %% state probability
    AnalyzeArousalStateProbability_Ephys_Handler(rootFolder,delim,false)
    AnalyzeArousalStateProbability_GCaMP_Handler(rootFolder,delim,false)

    %% sleep model accuracy
    AnalyzeModelAccuracy_Ephys_Handler(rootFolder,delim,false)
    AnalyzeModelAccuracy_GCaMP_Handler(rootFolder,delim,false)

    %% quanitification of whisking
    AnalyzeWhiskingBehavior_Ephys_Handler(rootFolder,delim,false)
    AnalyzeWhiskingBehavior_GCaMP_Handler(rootFolder,delim,false)

    %% arousal state transitions
    AnalyzeArousalTransitions_Ephys_Handler(rootFolder,delim,false)
    AnalyzeArousalTransitions_GCaMP_Handler(rootFolder,delim,false)
end
%{
Ephys exp are bilateral windows/stereotrodes SIBF @ 568 nm, sleep included
GCaMP exp are full bilateral windows SIBF + FC @ 480, 530, 630 nm, sleep included
2P, Pulse, Running spectroscopy are all animals from same data set, RH window only, awake only
Open field and fMRI are awake only
%}
disp('Loading analysis results and generating figures...'); disp(' ')
% main figure panels
% Fig1_nNOS(rootFolder,saveFigs,delim)
% Fig2_nNOS(rootFolder,saveFigs,delim)
% Fig3_nNOS(rootFolder,saveFigs,delim)
% Fig4_nNOS(rootFolder,saveFigs,delim)
% Fig5_nNOS(rootFolder,saveFigs,delim)
% Fig6_nNOS(rootFolder,saveFigs,delim)
% Fig7_nNOS(rootFolder,saveFigs,delim)
% Fig8_nNOS(rootFolder,saveFigs,delim)
% supplemental figure panels
FigS1_nNOS(rootFolder,saveFigs,delim)
% FigS2_nNOS()
% FigS3_nNOS()
% FigS4_nNOS()