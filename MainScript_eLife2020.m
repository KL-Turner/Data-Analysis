function [] = MainScript_eLife2020()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
% Purpose: Generates KLT's main and supplemental figs for Turner et al. eLife2020
%
% Scripts used to pre-process the original data are located in the folder "Pre-Processing Scripts".
% Functions that are used in both the analysis and pre-processing are located in the analysis folder.
%________________________________________________________________________________________________________________________

clear; clc; close all;
%% make sure the code repository and data are present in the current directory
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
%% run the data analysis. The progress bars will show the analysis progress
rerunAnalysis = 'n';
saveFigs = 'n';
if exist('AnalysisResults.mat','file') ~= 2 || strcmp(rerunAnalysis,'y') == true
    multiWaitbar_eLife2020('Analyzing sleep probability',0,'Color','B'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing behavioral distributions',0,'Color','W'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing behavioral heart rate' ,0,'Color','B'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing behavioral transitions',0,'Color','W'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing vessel behavioral transitions',0,'Color','B'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing behavioral hemodynamics',0,'Color','W'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing behavioral vessel diameter',0,'Color','B'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing laser doppler flow',0,'Color','W'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing coherence',0,'Color','B'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing neural-hemo coherence',0,'Color','W'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing power spectra',0,'Color','B'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing vessel power spectra',0,'Color','W'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing Pearson''s correlation coefficients',0,'Color','B'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing cross correlation',0,'Color','W'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing model cross validation distribution',0,'Color','B'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing evoked responses',0,'Color','W'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing vessel evoked responses',0,'Color','B'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing CBV-Gamma relationship',0,'Color','W'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing HbT-Sleep probability',0,'Color','B'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing TwoP-Sleep probability',0,'Color','W'); pause(0.25);
    multiWaitbar_eLife2020('Analyzing arteriole durations',0,'Color','B'); pause(0.25);
    % run analysis and output a structure containing all the analyzed data
    [AnalysisResults] = AnalyzeData_eLife2020(rootFolder);
    multiWaitbar_eLife2020('CloseAll');
else
    disp('Loading analysis results and generating figures...'); disp(' ')
    load('AnalysisResults.mat')
end
%% supplemental figure panels
[AnalysisResults] = Fig8_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig7_S3_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig7_S2_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig7_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig6_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig5_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig4_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig3_S5_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig3_S4_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig3_S3_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig3_S2_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig3_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig2_S2_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig2_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig1_S9_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig1_S8_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig1_S7_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig1_S6_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig1_S5_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig1_S4_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig1_S3_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig1_S2_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig1_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
%% supplemental tables
[AnalysisResults] = TableS12_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = TableS11_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = TableS10_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = TableS9_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = TableS8_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = TableS7_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = TableS6_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = TableS5_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
% TableS4 - text only, no figure
% TableS3 - text only, no figure
[AnalysisResults] = TableS2_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = TableS1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
%% main figure panels
[AnalysisResults] = Fig8_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig7_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig6_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig5_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig4_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig3_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig2_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Fig1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
%% tables
[AnalysisResults] = Table5_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Table4_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Table3_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Table2_eLife2020(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Table1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults); %#ok<NASGU>
%% fin.
disp('MainScript Analysis - Complete'); disp(' ')
end

function [AnalysisResults] = AnalyzeData_eLife2020(rootFolder)
% IOS animal IDs
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
% 2PLSM animal IDs
TwoP_animalIDs = {'T115','T116','T117','T118','T125','T126'};
saveFigs = 'y';
if exist('AnalysisResults.mat','file') == 2
    load('AnalysisResults.mat')
else
    AnalysisResults = [];
end
%% Block [1] Analyze the arousal-state probability of trial duration and resting events (IOS)
runFromStart = 'n';
for aa = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,aa})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,aa}),'SleepProbability') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeAwakeProbability_eLife2020(IOS_animalIDs{1,aa},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing sleep probability','Value',aa/length(IOS_animalIDs));
end
%% Block [2] Analyze the arousal-state distribution of different behavioral measurements (IOS)
runFromStart = 'n';
for bb = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,bb})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,bb}),'BehaviorDistributions') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeBehavioralDistributions_eLife2020(IOS_animalIDs{1,bb},rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing behavioral distributions','Value',bb/length(IOS_animalIDs));
end
%% Block [3] Analyze the heart rate during different arousal-states (IOS)
runFromStart = 'n';
for cc = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,cc})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,cc}),'MeanHR') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeMeanHeartRate_eLife2020(IOS_animalIDs{1,cc},rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing behavioral heart rate','Value',cc/length(IOS_animalIDs));
end
%% Block [4] Analyze the transitions between different arousal-states (IOS)
runFromStart = 'n';
for dd = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,dd})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,dd}),'Transitions') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeTransitionalAverages_eLife2020(IOS_animalIDs{1,dd},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing behavioral transitions','Value',dd/length(IOS_animalIDs));
end
%% Block [5] Analyze the transitions between different arousal-states (2PLSM)
runFromStart = 'n';
for ee = 1:length(TwoP_animalIDs)
    if isfield(AnalysisResults,(TwoP_animalIDs{1,ee})) == false || isfield(AnalysisResults.(TwoP_animalIDs{1,ee}),'Transitions') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeVesselTransitionalAverages_eLife2020(TwoP_animalIDs{1,ee},rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing vessel behavioral transitions','Value',ee/length(TwoP_animalIDs));
end
%% Block [6] Analyze the hemodynamic signal [HbT] during different arousal states (IOS)
runFromStart = 'n';
for ff = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,ff})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,ff}),'MeanCBV') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeMeanCBV_eLife2020(IOS_animalIDs{1,ff},rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing behavioral hemodynamics','Value',ff/length(IOS_animalIDs));
end
%% Block [7] Analyze the arteriole diameter D/D during different arousal states (2PLSM)
runFromStart = 'n';
for gg = 1:length(TwoP_animalIDs)
    if isfield(AnalysisResults,(TwoP_animalIDs{1,gg})) == false || isfield(AnalysisResults.(TwoP_animalIDs{1,gg}),'MeanVesselDiameter') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeMeanVesselDiameter_eLife2020(TwoP_animalIDs{1,gg},rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing behavioral vessel diameter','Value',gg/length(TwoP_animalIDs));
end
%% Block [8] Analyze the laser Doppler flowmetry during different arousal states (IOS)
runFromStart = 'n';
for hh = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,hh})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,hh}),'LDFlow') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeLaserDoppler_eLife2020(IOS_animalIDs{1,hh},rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing laser doppler flow','Value',hh/length(IOS_animalIDs));
end
%% Block [9] Analyze the spectral coherence between bilateral hemodynamic [HbT] and neural signals (IOS)
runFromStart = 'n';
for jj = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,jj})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,jj}),'Coherence') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeCoherence_eLife2020(IOS_animalIDs{1,jj},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing coherence','Value',jj/length(IOS_animalIDs));
end
%% Block [10] Analyze the spectral coherence between neural-hemodynamic [HbT] signals (IOS)
runFromStart = 'n';
for jj = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,jj})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,jj}),'NeuralHemoCoherence') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeNeuralHemoCoherence_eLife2020(IOS_animalIDs{1,jj},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing neural-hemo coherence','Value',jj/length(IOS_animalIDs));
end
%% Block [11] Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
runFromStart = 'n';
for kk = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,kk})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,kk}),'PowerSpectra') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzePowerSpectrum_eLife2020(IOS_animalIDs{1,kk},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing power spectra','Value',kk/length(IOS_animalIDs));
end
%% Block [12] Analyze the spectral power of arteriole diameter D/D (2PLSM)
runFromStart = 'n';
for ll = 1:length(TwoP_animalIDs)
    if isfield(AnalysisResults,(TwoP_animalIDs{1,ll})) == false || isfield(AnalysisResults.(TwoP_animalIDs{1,ll}),'PowerSpectra') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeVesselPowerSpectrum_eLife2020(TwoP_animalIDs{1,ll},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing vessel power spectra','Value',ll/length(TwoP_animalIDs));
end
%% Block [13] Analyze Pearson's correlation coefficient between bilateral hemodynamic [HbT] and neural signals (IOS)
runFromStart = 'n';
for mm = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,mm})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,mm}),'CorrCoeff') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeCorrCoeffs_eLife2020(IOS_animalIDs{1,mm},rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing Pearson''s correlation coefficients','Value',mm/length(IOS_animalIDs));
end
%% Block [14] Analyze the cross-correlation between neural activity and hemodynamics [HbT] (IOS)
runFromStart = 'n';
for nn = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,nn})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,nn}),'XCorr') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeXCorr_eLife2020(IOS_animalIDs{1,nn},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing cross correlation','Value',nn/length(IOS_animalIDs));
end
%% Block [15] Analyze the out-of-bag error (model accuracy) of each random forest classification model (IOS)
runFromStart = 'n';
for oo = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,oo})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,oo}),'ModelAccuracy') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeModelAccuracy_eLife2020(IOS_animalIDs{1,oo},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing model cross validation distribution','Value',oo/length(IOS_animalIDs));
end
%% Block [16] Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
runFromStart = 'n';
for pp = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,pp})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,pp}),'EvokedAvgs') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeEvokedResponses_eLife2020(IOS_animalIDs{1,pp},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing evoked responses','Value',pp/length(IOS_animalIDs));
end
%% Block [17] Analyze the whisking-evoked arteriole D/D responses (2PLSM)
runFromStart = 'n';
for qq = 1:length(TwoP_animalIDs)
    if isfield(AnalysisResults,(TwoP_animalIDs{1,qq})) == false || isfield(AnalysisResults.(TwoP_animalIDs{1,qq}),'EvokedAvgs') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeVesselEvokedResponses_eLife2020(TwoP_animalIDs{1,qq},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing vessel evoked responses','Value',qq/length(TwoP_animalIDs));
end
%% Block [18] Analyze the relationship between gamma-band power and hemodynamics [HbT] (IOS)
runFromStart = 'n';
for qq = 1:length(IOS_animalIDs)
    if isfield(AnalysisResults,(IOS_animalIDs{1,qq})) == false || isfield(AnalysisResults.(IOS_animalIDs{1,qq}),'HbTvsGamma') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeCBVGammaRelationship_eLife2020(IOS_animalIDs{1,qq},rootFolder,AnalysisResults);
    end
    multiWaitbar_eLife2020('Analyzing CBV-Gamma relationship','Value',qq/length(IOS_animalIDs));
end
%% Block [19] Analyze the probability of arousal-state classification based on hemodynamic [HbT] changes (IOS)
runFromStart = 'n';
if isfield(AnalysisResults,'HbTSleepProbability') == false || strcmp(runFromStart,'y') == true
    [AnalysisResults] = AnalyzeHbTSleepProbability_eLife2020(IOS_animalIDs,rootFolder,AnalysisResults);
end
multiWaitbar_eLife2020('Analyzing HbT-Sleep probability','Value',1/length(1));
%% Block [20] Analyze the probability of arousal-state classification based on arteriole D/D changes (2PLSM)
runFromStart = 'n';
if isfield(AnalysisResults,'TwoPSleepProbability') == false || strcmp(runFromStart,'y') == true
    [AnalysisResults] = AnalyzeTwoPSleepProbability_eLife2020(TwoP_animalIDs,rootFolder,AnalysisResults);
end
multiWaitbar_eLife2020('Analyzing TwoP-Sleep probability','Value',1/length(1));
%% Block [21] Analyze the time of each arousal-state data per artery (2PLSM)
runFromStart = 'n';
if isfield(AnalysisResults,'ArterioleDurations') == false || strcmp(runFromStart,'y') == true
    [AnalysisResults] = AnalyzeArterioleDurations_eLife2020(TwoP_animalIDs,rootFolder,AnalysisResults);
end
multiWaitbar_eLife2020('Analyzing arteriole durations','Value',1/length(1));
%% fin.
disp('Loading analysis results and generating figures...'); disp(' ')

end
