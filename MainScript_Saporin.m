function [] = MainScript_Saporin()
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
rerunAnalysis = 'y';
saveFigs = 'y';
if exist('AnalysisResults.mat','file') ~= 2 || strcmp(rerunAnalysis,'y') == true
    multiWaitbar('Analyzing behavioral hemodynamics',0,'Color','P'); pause(0.25);
    multiWaitbar('Analyzing coherence',0,'Color','B'); pause(0.25);
    multiWaitbar('Analyzing neural-hemo coherence',0,'Color','G'); pause(0.25);
    multiWaitbar('Analyzing power spectra',0,'Color','P'); pause(0.25);
    multiWaitbar('Analyzing Pearson''s correlation coefficients',0,'Color','B'); pause(0.25);
    multiWaitbar('Analyzing cross correlation',0,'Color','G'); pause(0.25);
    multiWaitbar('Analyzing evoked responses',0,'Color','P'); pause(0.25);
    multiWaitbar('Analyzing gamma-HbT relationship',0,'Color','B'); pause(0.25);
    multiWaitbar('Analyzing power spectra2',0,'Color','G'); pause(0.25);
    multiWaitbar('Analyzing whisk-hemo coherence',0,'Color','P'); pause(0.25);
    [AnalysisResults] = AnalyzeData(rootFolder);
    multiWaitbar('CloseAll');
else
    disp('Loading analysis results and generating figures...'); disp(' ')
    load('AnalysisResults.mat','-mat')
end
%% generate figures
% [AnalysisResults] = DiaphoraseCellCounts_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = WhiskEvoked_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = StimEvoked_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = Coherence_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = NeuralHemoCoherence_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = WhiskHemoCoherence_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = PowerSpec_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = PowerSpec2_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = PearsonsCorr_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = XCorr_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = NeuralHemoLinearity_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
[AnalysisResults] = MeanHbT_Saporin(rootFolder,saveFigs,delim,AnalysisResults); %#ok<*NASGU>
%% fin.
disp('MainScript Analysis - Complete'); disp(' ')
end

function [AnalysisResults] = AnalyzeData(rootFolder)
% IOS animal IDs
expGroups = {'C57BL6J','SSP-SAP','Blank-SAP'};
saveFigs = 'y';
if exist('AnalysisResults.mat','file') == 2
    load('AnalysisResults.mat','-mat')
else
    AnalysisResults = [];
end
% determine waitbar length
waitBarLength = 0;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name}, '.'));
    animalIDs = {folderList.name};
    waitBarLength = waitBarLength + length(animalIDs);
end
%% Analyze the hemodynamic signal [HbT] during different arousal states (IOS)
runFromStart = 'n';
cc = 1;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name}, '.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'MeanCBV') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzeMeanCBV(animalIDs{1,bb},expGroups{1,aa},rootFolder,AnalysisResults);
        end
        multiWaitbar('Analyzing behavioral hemodynamics','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
%% Analyze the spectral coherence between bilateral hemodynamic [HbT] and neural signals (IOS)
runFromStart = 'n';
cc = 1;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'Coherence') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzeCoherence(animalIDs{1,bb},expGroups{1,aa},rootFolder,AnalysisResults);
        end
        multiWaitbar('Analyzing coherence','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
%% Analyze the spectral coherence between neural-hemodynamic [HbT] signals (IOS)
runFromStart = 'n';
cc = 1;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name}, '.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'NeuralHemoCoherence') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzeNeuralHemoCoherence(animalIDs{1,bb},expGroups{1,aa},rootFolder,AnalysisResults);
        end
        multiWaitbar('Analyzing neural-hemo coherence','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
%% Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
runFromStart = 'n';
cc = 1;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name}, '.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'PowerSpectra') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzePowerSpectrum(animalIDs{1,bb},expGroups{1,aa},rootFolder,AnalysisResults);
        end
        multiWaitbar('Analyzing power spectra','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
%% Analyze Pearson's correlation coefficient between bilateral hemodynamic [HbT] and neural signals (IOS)
runFromStart = 'n';
cc = 1;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name}, '.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'CorrCoeff') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzeCorrCoeffs(animalIDs{1,bb},expGroups{1,aa},rootFolder,AnalysisResults);
        end
        multiWaitbar('Analyzing Pearson''s correlation coefficients','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
%% Analyze the cross-correlation between neural activity and hemodynamics [HbT] (IOS)
runFromStart = 'n';
cc = 1;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name}, '.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'XCorr') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzeXCorr(animalIDs{1,bb},expGroups{1,aa},saveFigs,rootFolder,AnalysisResults);
        end
        multiWaitbar('Analyzing cross correlation','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
%% Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
runFromStart = 'n';
cc = 1;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name}, '.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'EvokedAvgs') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzeEvokedResponses(animalIDs{1,bb},expGroups{1,aa},saveFigs,rootFolder,AnalysisResults);
        end
        multiWaitbar('Analyzing evoked responses','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
%% Analyze the relationship between gamma-band power and hemodynamics [HbT] (IOS)
runFromStart = 'n';
cc = 1;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name}, '.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'HbTvsGamma') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzeCBVGammaRelationship(animalIDs{1,bb},expGroups{1,aa},rootFolder,AnalysisResults);
        end
        multiWaitbar('Analyzing gamma-HbT relationship','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
%% Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
runFromStart = 'n';
cc = 1;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name}, '.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'PowerSpectra2') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzePowerSpectrum2(animalIDs{1,bb},expGroups{1,aa},rootFolder,AnalysisResults);
        end
        multiWaitbar('Analyzing power spectra2','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
%% Analyze the spectral coherence between whisk-hemodynamic [HbT] signals (IOS)
runFromStart = 'n';
cc = 1;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name}, '.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'WhiskHemoCoherence') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzeWhiskHemoCoherence(animalIDs{1,bb},expGroups{1,aa},rootFolder,AnalysisResults);
        end
        multiWaitbar('Analyzing whisk-hemo coherence','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
%% fin.
disp('Loading analysis results and generating figures...'); disp(' ')
end