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
zap;
cd('D:\NO Project\')
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
% add root folder to Matlab's working directory
addpath(genpath(rootFolder))
% analysis subfunctions
runAnalysis = true;
if runAnalysis == true
% % % % % % % %     % state probability
% % % % % % % %     AnalyzeArousalStateProbability_Ephys_Handler(rootFolder,delim,false)
% % % % % % % %     % quanitification of behavior
% % % % % % % %     AnalyzeWhiskingBehavior_Ephys_Handler(rootFolder,delim,false)
% % % % % % % %     % arousal state transitions
% % % % % % % %     AnalyzeArousalTransitions_Ephys_Handler(rootFolder,delim,false)
% % % % % % % %     % vessel baseline shif
% % % % % % % %     AnalyzeVesselBaselineShift_2P_Handler(rootFolder,delim,false)
% % % % % % % %     % stimulus and whisking evoked responses
% % % % % % % %     AnalyzeEvokedResponses_Ephys_Handler(rootFolder,delim,false)
% % % % % % % %     AnalyzeEvokedResponses_GCaMP_Handler(rootFolder,delim,false)
% % % % % % % %     AnalyzeEvokedResponses_2P_Handler(rootFolder,delim,false)
% % % % % % % %     AnalyzeEvokedResponses_Pulse_Handler(rootFolder,delim,false)
% % % % % % % %     % IOS signal analysis (HbT, GCaMP, Deoxy)
% % % % % % % %     AnalyzeIntrinsicSignals_Ephys_Handler(rootFolder,delim,false)
% % % % % % % %     AnalyzeIntrinsicSignals_GCaMP_Handler(rootFolder,delim,false)
% % % % % % % %     % power spectral density of IOS and second spectra for neural signals
% % % % % % % %     AnalyzePowerSpectrum_Ephys_Handler(rootFolder,delim,false)
% % % % % % % %     AnalyzePowerSpectrum_GCaMP_Handler(rootFolder,delim,false)
% % % % % % % %     AnalyzePowerSpectrum_LFP_Handler(rootFolder,delim,false)
% % % % % % % %     % coherence between hemipheres
% % % % % % % %     AnalyzeBilateralCoherence_Ephys_Handler(rootFolder,delim,false)
% % % % % % % %     AnalyzeBilateralCoherence_GCaMP_Handler(rootFolder,delim,false)
% % % % % % % %     % Pearson's correlation between hemispheres
% % % % % % % %     AnalyzePearsonCorrelation_Ephys_Handler(rootFolder,delim,false)
% % % % % % % %     AnalyzePearsonCorrelation_GCaMP_Handler(rootFolder,delim,false)
    % neural-hemo coherence within hemispheres
    AnalyzeNeuralHemoCoherence_Ephys_Handler(rootFolder,delim,false)
    AnalyzeNeuralHemoCoherence_GCaMP_Handler(rootFolder,delim,false)
    % cross correlation between neural and hemodynamic signals
    AnalyzeCrossCorrelation_Ephys_Handler(rootFolder,delim,false)
    AnalyzeCrossCorrelation_GCaMP_Handler(rootFolder,delim,false)
    % pupil analysis
    AnalyzePupilArea_Ephys_Handler(rootFolder,delim,false)
    AnalyzePupilEvokedResponses_Ephys_Handler(rootFolder,delim,false)
    AnalyzePupilCoherence_EPhys_Handler(rootFolder,delim,false)
    AnalyzePupilCrossCorrelation_EPhys_Handler(rootFolder,delim,false)
    % close wait bars
    multiWaitbar('close all')
end
%{
Ephys exp are bilateral windows/stereotrodes SIBF @ 568 nm, sleep included
GCaMP exp are full bilateral windows SIBF + FC @ 480, 530, 630 nm, sleep included
2P, Pulse, Running spectroscopy are all animals from same data set, RH window only, awake only
Open field, fMRI, and NPY DREADDs are awake only
%}
disp('Loading analysis results and generating figures...'); disp(' ')
% histology, IHC figures
DiaphoraseCellCounts_Figures(rootFolder,saveFigs,delim); % Dakota, Denver, Kyle, Nikki
% behavior analysis
OpenFieldBehavior_Figures(rootFolder,saveFigs,delim) % Shakhawat
ArousalStateProbability_Ephys_Figures(rootFolder,saveFigs,delim)
WhiskingBehavior_Ephys_Figures(rootFolder,saveFigs,delim);
ArousalTransitions_Bilateral_IOS(rootFolder,saveFigs,delim);
% IOS signals (HbT, GCaMP, HbO/R) and 2p baselines
IntrinsicSignals_Ephys_Figures(rootFolder,saveFigs,delim)
IntrinsicSignals_GCaMP_Figures(rootFolder,saveFigs,delim)
BaselineShift_2P_Figures(rootFolder,saveFigs,delim);
% stimulus evoked (HbT, GCaMP, HbO/R, gamma/MUA, 2p diameter, running)
StimEvoked_Pulse_Figures(rootFolder,saveFigs,delim);
StimEvoked_GCaMP_Figures(rootFolder,saveFigs,delim);
StimEvoked_2P_Figures(rootFolder,saveFigs,delim);
StimEvoked_Ephys_Figures(rootFolder,saveFigs,delim);
RunningSpectroscopy_Figures(rootFolder,saveFigs,delim) % Qingguang
% whisking evoked (HbT, GCaMP, HbO/R, gamma/MUA, 2p diameter)
WhiskEvoked_Pulse_Figures(rootFolder,saveFigs,delim);
WhiskEvoked_GCaMP_Figures(rootFolder,saveFigs,delim);
WhiskEvoked_2P_Figures(rootFolder,saveFigs,delim);
WhiskEvoked_Ephys_Figures(rootFolder,saveFigs,delim);
% power spectra (LFP, HbT, GCaMP, HbO/R, gamma, 2p diameter)
PowerSpectrum_LFP_Figures(rootFolder,saveFigs,delim);
PowerSpectrum_Ephys_Figures(rootFolder,saveFigs);
PowerSpectrum_GCaMP_Figures(rootFolder,saveFigs);
PowerSpectrum_2P_Figures(rootFolder,saveFigs);
% coherence (HbT, GCaMP, HbO/R, gamma)
BilateralCoherence_Ephys_Figures(rootFolder,saveFigs,delim)
BilateralCoherence_GCaMP_Figures(rootFolder,saveFigs,delim)
% Pearson's correlations (HbT, GCaMP, HbO/R, gamma)
PearsonsCorrelation_Ephys_Figures(rotFolder,saveFigs,delim)
PearsonsCorrelation_Ephys_Figures(rotFolder,saveFigs,delim)
% neural-hemo coherence (HbT, GCaMP, HbO/R, gamma)
NeuralHemoCoherence_Ephys_Figures(rootFolder,saveFigs,delim)
NeuralHemoCoherence_GCaMP_Figures(rootFolder,saveFigs,delim)
% cross correlation (HbT, GCaMP, HbO/R, gamma)
CrossCorrelation_Ephys_Figures(rootFolder,saveFigs,delim);
CrossCorrelation_GCaMP_Figures(rootFolder,saveFigs,delim);
% pupil area, stimulus evoked response, coherence/XC w/ HbT/gamma
AnalyzeBehavioralArea_Ephys_Handler(rootFolder,delim,false)
AnalyzeEvokedResponses_Ephys_Handler(rootFolder,delim,false)
AnalyzeCoherence_EPhys_Handler(rootFolder,delim,false)
AnalyzeCrossCorrelation_EPhys_Handler(rootFolder,delim,false)
%% fMRI experiments - Said, Qingqing, Nanyin
% BOLD stimulus evoked responses
% resting state functional connectivity mapping, barrels vs. FC
% bold & rsFC changes over time due to adaptation
%% NPY G(i) vs. G(q) experiments - Mike, Dakota IHC NPY w/ DREADD (mCherry), antiNK1R, antiNPY
% stimulus evoked 5 min 1/10 sec, 5 min 5 sec
% whisking evoked
% resting baseline diameter