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
% analysis subfunctions
addpath(genpath(rootFolder))
runAnalysis = false;
if runAnalysis == true
    %     %% evoked responses
    %     AnalyzeEvokedResponses_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzeEvokedResponses_GCaMP_Handler(rootFolder,delim,false)
    %     AnalyzeEvokedResponses_2P_Handler(rootFolder,delim,false)
    %     AnalyzeEvokedResponses_Pulse_Handler(rootFolder,delim,false)
    %
    %     %% HbT, HbO, HbR, GCaMP
    %     AnalyzeIntrinsicSignals_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzeIntrinsicSignals_GCaMP_Handler(rootFolder,delim,false)
    %     AnalyzeHemoGFPRelationship_EGFP_Handler(rootFolder,delim,false)
    %
    %     %% arteriole diameter and baseline shift
    %     AnalyzeArterioleBaseline_2P_Handler(rootFolder,delim,false)
    %     AnalyzeArterioleDiameter_2P_Handler(rootFolder,delim,false)
    %
    %     %% power spectral density
    %     AnalyzePowerSpectrum_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzePowerSpectrum_GCaMP_Handler(rootFolder,delim,false)
    %     AnalyzePowerSpectrum_LFP_Handler(rootFolder,delim,false)
    %     AnalyzePowerSpectrum_2P_Handler(rootFolder,delim,false)
    %
    %     %% coherence between hemipheres
    %     AnalyzeBilateralCoherence_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzeBilateralCoherence_GCaMP_Handler(rootFolder,delim,false)
    %
    %     %% coherence between neural-hemo
    AnalyzeNeuralHemoCoherence_Ephys_Handler(rootFolder,delim,true)
    %     AnalyzeNeuralHemoCoherence_GCaMP_Handler(rootFolder,delim,false)
    %
    %     %% Pearson's correlation between hemispheres
    %     AnalyzePearsonCorrelation_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzePearsonCorrelation_GCaMP_Handler(rootFolder,delim,false)
    %
    %     %% cross correlation between neural-hemo
    %     % AnalyzeCrossCorrelation_Ephys_Handler(rootFolder,delim,true)
    %     AnalyzeCrossCorrelation_GCaMP_Handler(rootFolder,delim,false)
    %
    %     %% pupil analysis
    %     AnalyzePupilArea_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzePupilEvokedResponses_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzePupilCoherence_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzePupilCrossCorrelation_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzePupilInterBlinkInterval_Ephys_Handler(rootFolder,delim,false)
    %
    %     %% state probability
    %     AnalyzeArousalStateProbability_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzeArousalStateProbability_GCaMP_Handler(rootFolder,delim,false)
    %
    %     %% sleep model accuracy
    %     AnalyzeModelAccuracy_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzeModelAccuracy_GCaMP_Handler(rootFolder,delim,false)
    %
    %     %% quanitification of whisking
    %     AnalyzeWhiskingBehavior_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzeWhiskingBehavior_GCaMP_Handler(rootFolder,delim,false)
    %
    %     %% arousal state transitions
    %     AnalyzeArousalTransitions_Ephys_Handler(rootFolder,delim,false)
    %     AnalyzeArousalTransitions_GCaMP_Handler(rootFolder,delim,false)

end
%{
Ephys exp are bilateral windows/stereotrodes SIBF @ 568 nm, sleep included
GCaMP exp are full bilateral windows SIBF + FC @ 480, 530, 630 nm, sleep included
2P, Pulse, Running spectroscopy are all animals from same data set, RH window only, awake only
Open field and fMRI are awake only
%}
disp('Loading analysis results and generating figures...'); disp(' ')
%% Single trial examples - DONE
SingleTrialExample_Awake_GCaMP_Figures(rootFolder,saveFigs,delim);
SingleTrialExample_Asleep_GCaMP_Figures(rootFolder,saveFigs,delim);

%% histology, behavior, controls - DONE
% DiaphoraseCellCounts_Figures(rootFolder,saveFigs,delim);
% OpenFieldBehavior_Figures(rootFolder,saveFigs,delim)
% HemoGFPRelationship_EGFP_Figures(rootFolder,saveFigs,delim)

%% arousal-state probability - DONE
% ArousalStateProb_Ephys_Figures(rootFolder,saveFigs,delim)
% ArousalStateProb_GCaMP_Figures(rootFolder,saveFigs,delim)

%% whisking behavior - DONE
% WhiskingBehavior_Ephys_Figures(rootFolder,saveFigs,delim);
% WhiskingBehavior_GCaMP_Figures(rootFolder,saveFigs,delim);

%% pupil size, stimulus evoked response, coherence/XC w/ HbT/gamma - DONE
% PupilArea_Ephys_Figures(rootFolder,saveFigs,delim)
% PupilStimEvoked_Ephys_Figures(rootFolder,saveFigs,delim)
% PupilCoherence_Ephys_Figures(rootFolder,saveFigs,delim)
% PupilCrossCorrelation_Ephys_Figures(rootFolder,saveFigs,delim)
% PupilInterBlinkInterval_Ephys_Figures(rootFolder,saveFigs,delim)

%% arousal-state transitions - DONE
% ArousalTransitions_Ephys_Figures(rootFolder,saveFigs,delim);
% ArousalTransitions_GCaMP_SI_Figures(rootFolder,saveFigs,delim);
% ArousalTransitions_GCaMP_FC_Figures(rootFolder,saveFigs,delim);

%% IOS signals (HbT, GCaMP, HbO/R) - DONE
% IntrinsicSignals_Ephys_Figures(rootFolder,saveFigs,delim)
% IntrinsicSignals_GCaMP_SI_Figures(rootFolder,saveFigs,delim)
% IntrinsicSignals_GCaMP_FC_Figures(rootFolder,saveFigs,delim)

%% 2p and baseline shift - DONE
% ArterioleDiameter_2P_Figures(rootFolder,saveFigs,delim);
% BaselineShift_2P_Figures(rootFolder,saveFigs,delim);

%% stimulus evoked (HbT, GCaMP, HbO/R, gamma/MUA, 2p diameter, running) - DONE
% StimEvoked_Ephys_Figures(rootFolder,saveFigs,delim);
% StimEvoked_GCaMP_SI_Figures(rootFolder,saveFigs,delim);
% StimEvoked_GCaMP_FC_Figures(rootFolder,saveFigs,delim);
% StimEvoked_2P_Figures(rootFolder,saveFigs,delim);
% StimEvoked_Pulse_Figures(rootFolder,saveFigs,delim);
% RunningSpectroscopy_Figures(rootFolder,saveFigs,delim)

%% whisking evoked (HbT, GCaMP, HbO/R, gamma/MUA, 2p diameter) - DONE
% VolitionalWhisk_Ephys_Figures(rootFolder,saveFigs,delim);
% VolitionalWhisk_GCaMP_SI_Figures(rootFolder,saveFigs,delim);
% VolitionalWhisk_GCaMP_FC_Figures(rootFolder,saveFigs,delim);
% VolitionalWhisk_2P_Figures(rootFolder,saveFigs,delim);
% VolitionalWhisk_Pulse_Figures(rootFolder,saveFigs,delim);

%% power spectra (LFP, HbT, GCaMP, HbO/R, gamma, 2p diameter) - DONE
% PowerSpectrum_LFP_Figures(rootFolder,saveFigs,delim);
% PowerSpectrum_Ephys_Figures(rootFolder,saveFigs,delim);
% PowerSpectrum_GCaMP_SI_Figures(rootFolder,saveFigs,delim);
% PowerSpectrum_GCaMP_FC_Figures(rootFolder,saveFigs,delim);
% PowerSpectrum_2P_Figures(rootFolder,saveFigs,delim);

%% bilateral coherence (HbT, GCaMP, HbO/R, gamma) - DONE
% BilateralCoherence_Ephys_Figures(rootFolder,saveFigs,delim)
% BilateralCoherence_GCaMP_SI_Figures(rootFolder,saveFigs,delim)
% BilateralCoherence_GCaMP_FC_Figures(rootFolder,saveFigs,delim)

%% Pearson's correlations (HbT, GCaMP, HbO/R, gamma) - DONE
% PearsonsCorrelation_Ephys_Figures(rootFolder,saveFigs,delim)
% PearsonsCorrelation_GCaMP_SI_Figures(rootFolder,saveFigs,delim)
% PearsonsCorrelation_GCaMP_FC_Figures(rootFolder,saveFigs,delim)

%% neural-hemo coherence (HbT, GCaMP, HbO/R, gamma) - DONE
% NeuralHemoCoherence_Ephys_Figures(rootFolder,saveFigs,delim)
% NeuralHemoCoherence_GCaMP_SI_Figures(rootFolder,saveFigs,delim)
% NeuralHemoCoherence_GCaMP_FC_Figures(rootFolder,saveFigs,delim)

%% cross correlation (HbT, GCaMP, HbO/R, gamma) - DONE
% CrossCorrelation_Ephys_Figures(rootFolder,saveFigs,delim);
% CrossCorrelation_GCaMP_SI_Figures(rootFolder,saveFigs,delim);
% CrossCorrelation_GCaMP_FC_Figures(rootFolder,saveFigs,delim);

%% BOLD
% Zhang data TBD