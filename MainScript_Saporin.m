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
%% verify code repository and data are in the current directory/added path
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
%% analysis subfunctions
runAnalysis = true;
if runAnalysis == true
    AnalyzeBehavioralHbT_Handler(rootFolder,delim,'n')
    AnalyzeBilateralCoherence_Handler(rootFolder,delim,'n')
    AnalyzeNeuralHemoCoherence_Handler(rootFolder,delim,'n')
    AnalyzePowerSpectrum_Handler(rootFolder,delim,'n')
end
%% generate figures
disp('Loading analysis results and generating figures...'); disp(' ')
saveFigs = true;
% DiaphoraseCellCounts_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% WhiskEvoked_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% StimEvoked_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% StimEvoked2_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% StimEvoked3_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% Coherence_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% NeuralHemoCoherence_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% WhiskHemoCoherence_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% PowerSpec_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% PowerSpec2_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% PearsonsCorr_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% XCorr_Saporin(rootFolder,saveFigs,delim,AnalysisResults);
% MeanHbT_Saporin(rootFolder,saveFigs,delim,AnalysisResults); %#ok<*NASGU>

end

% 
% %% Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
% setName = 'IOS Set A';
% expGroups = {'C57BL6J','SSP-SAP','Blank-SAP'};
% % determine waitbar length
% waitBarLength = 0;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]);
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     waitBarLength = waitBarLength + length(animalIDs);
% end
% % run analysis for each animal in the group
% runFromStart = 'n';
% cc = 1;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]); 
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     for bb = 1:length(animalIDs)
%         if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'PowerSpectra') == false || strcmp(runFromStart,'y') == true
%             [AnalysisResults] = AnalyzePowerSpectrum(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,AnalysisResults);
%         end
%         multiWaitbar('Analyzing power spectra','Value',cc/waitBarLength);
%         cc = cc + 1;
%     end
% end
% %% Analyze Pearson's correlation coefficient between bilateral hemodynamic [HbT] and neural signals (IOS)
% setName = 'IOS Set A';
% expGroups = {'C57BL6J','SSP-SAP','Blank-SAP'};
% % determine waitbar length
% waitBarLength = 0;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]);
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     waitBarLength = waitBarLength + length(animalIDs);
% end
% % run analysis for each animal in the group
% runFromStart = 'n';
% cc = 1;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]); 
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     for bb = 1:length(animalIDs)
%         if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'CorrCoeff') == false || strcmp(runFromStart,'y') == true
%             [AnalysisResults] = AnalyzeCorrCoeffs(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,AnalysisResults);
%         end
%         multiWaitbar('Analyzing Pearson''s correlation coefficients','Value',cc/waitBarLength);
%         cc = cc + 1;
%     end
% end
% %% Analyze the cross-correlation between neural activity and hemodynamics [HbT] (IOS)
% setName = 'IOS Set A';
% expGroups = {'C57BL6J','SSP-SAP','Blank-SAP'};
% % determine waitbar length
% waitBarLength = 0;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]);
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     waitBarLength = waitBarLength + length(animalIDs);
% end
% % run analysis for each animal in the group
% runFromStart = 'n';
% cc = 1;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]); 
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     for bb = 1:length(animalIDs)
%         if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'XCorr') == false || strcmp(runFromStart,'y') == true
%             [AnalysisResults] = AnalyzeXCorr(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,AnalysisResults);
%         end
%         multiWaitbar('Analyzing cross correlation','Value',cc/waitBarLength);
%         cc = cc + 1;
%     end
% end
% %% Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
% setName = 'IOS Set A';
% expGroups = {'C57BL6J','SSP-SAP','Blank-SAP'};
% % determine waitbar length
% waitBarLength = 0;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]);
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     waitBarLength = waitBarLength + length(animalIDs);
% end
% % run analysis for each animal in the group
% runFromStart = 'n';
% cc = 1;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]); 
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     for bb = 1:length(animalIDs)
%         if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'EvokedAvgs') == false || strcmp(runFromStart,'y') == true
%             [AnalysisResults] = AnalyzeEvokedResponses(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,AnalysisResults);
%         end
%         multiWaitbar('Analyzing evoked responses (IOS Set A)','Value',cc/waitBarLength);
%         cc = cc + 1;
%     end
% end
% %% Analyze the relationship between gamma-band power and hemodynamics [HbT] (IOS)
% setName = 'IOS Set A';
% expGroups = {'C57BL6J','SSP-SAP','Blank-SAP'};
% % determine waitbar length
% waitBarLength = 0;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]);
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     waitBarLength = waitBarLength + length(animalIDs);
% end
% % run analysis for each animal in the group
% runFromStart = 'n';
% cc = 1;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]); 
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     for bb = 1:length(animalIDs)
%         if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'HbTvsGamma') == false || strcmp(runFromStart,'y') == true
%             [AnalysisResults] = AnalyzeCBVGammaRelationship(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,AnalysisResults);
%         end
%         multiWaitbar('Analyzing gamma-HbT relationship','Value',cc/waitBarLength);
%         cc = cc + 1;
%     end
% end
% 
% %% Analyze the spectral coherence between whisk-hemodynamic [HbT] signals (IOS)
% setName = 'IOS Set A';
% expGroups = {'C57BL6J','SSP-SAP','Blank-SAP'};
% % determine waitbar length
% waitBarLength = 0;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]);
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     waitBarLength = waitBarLength + length(animalIDs);
% end
% % run analysis for each animal in the group
% runFromStart = 'n';
% cc = 1;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]); 
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     for bb = 1:length(animalIDs)
%         if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'WhiskHemoCoherence') == false || strcmp(runFromStart,'y') == true
%             [AnalysisResults] = AnalyzeWhiskHemoCoherence(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,AnalysisResults);
%         end
%         multiWaitbar('Analyzing whisk-hemo coherence','Value',cc/waitBarLength);
%         cc = cc + 1;
%     end
% end
% %% Analyze the arteriole diameter D/D during different arousal states (2PLSM)
% expGroups = {'SSP-SAP','Blank-SAP'};
% setName = '2PLSM Set B';
% % determine waitbar length
% waitBarLength = 0;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]);
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     waitBarLength = waitBarLength + length(animalIDs);
% end
% % run analysis for each animal in the group
% runFromStart = 'n';
% cc = 1;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]); 
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     for bb = 1:length(animalIDs)
%         if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'MeanVesselDiameter') == false || strcmp(runFromStart,'y') == true
%             [AnalysisResults] = AnalyzeMeanVesselDiameter(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,AnalysisResults);
%         end
%         multiWaitbar('Analyzing behavioral arteriole diameter','Value',cc/waitBarLength);
%         cc = cc + 1;
%     end
% end
% %% Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
% expGroups = {'SSP-SAP','Blank-SAP'};
% setName = '2PLSM Set B';
% % determine waitbar length
% waitBarLength = 0;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]);
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     waitBarLength = waitBarLength + length(animalIDs);
% end
% % run analysis for each animal in the group
% runFromStart = 'n';
% cc = 1;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]); 
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     for bb = 1:length(animalIDs)
%         if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'EvokedAvgs_2P') == false || strcmp(runFromStart,'y') == true
%             [AnalysisResults] = AnalyzeVesselEvokedResponses(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,AnalysisResults);
%         end
%         multiWaitbar('Analyzing evoked responses (2PLSM Set B)','Value',cc/waitBarLength);
%         cc = cc + 1;
%     end
% end
% %% Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
% expGroups = {'SSP-SAP','Blank-SAP'};
% setName = 'IOS Set B';
% % determine waitbar length
% waitBarLength = 0;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]);
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     waitBarLength = waitBarLength + length(animalIDs);
% end
% % run analysis for each animal in the group
% runFromStart = 'n';
% cc = 1;
% for aa = 1:length(expGroups)
%     folderList = dir([expGroups{1,aa} delim setName]); 
%     folderList = folderList(~startsWith({folderList.name}, '.'));
%     animalIDs = {folderList.name};
%     for bb = 1:length(animalIDs)
%         if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'EvokedAvgs') == false || strcmp(runFromStart,'y') == true
%             [AnalysisResults] = AnalyzeEvokedResponses(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,AnalysisResults);
%         end
%         multiWaitbar('Analyzing evoked responses (IOS Set B)','Value',cc/waitBarLength);
%         cc = cc + 1;
%     end
% end
% %% fin.
% disp('Loading analysis results and generating figures...'); disp(' ')
% 
% end
