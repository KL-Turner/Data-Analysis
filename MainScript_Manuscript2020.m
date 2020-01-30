function [] = MainScript_Manuscript2020()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Generates KLT's main and supplemental figs for the 2020 sleep paper. 
%
%            Scripts used to pre-process the original data are located in the folder "Pre-processing-scripts".
%            Functions that are used in both the analysis and pre-processing are located in the analysis folder.
%________________________________________________________________________________________________________________________

clear; clc; close all
%% Make sure the current directory is 'TurnerFigs-Neuron2020' and that the code repository is present.
currentFolder = pwd;
addpath(genpath(currentFolder));
fileparts = strsplit(currentFolder,filesep);
if ismac
    rootfolder = fullfile(filesep,fileparts{1:end});
else
    rootfolder = fullfile(fileparts{1:end});
end
% Add root folder to Matlab's working directory.
addpath(genpath(rootfolder))

%% Run the data analysis. The progress bars will show the analysis progress.
dataSummary = dir('ComparisonData.mat');
% If the analysis structure has already been created, load it and skip the analysis.
if ~isempty(dataSummary)
    load(dataSummary.name);
    disp('Loading analysis results and generating figures...'); disp(' ')
else
    multiWaitbar_Neuron2020('Analyzing cross correlation',0,'Color',[0.720000 0.530000 0.040000]); pause(0.25);
    multiWaitbar_Neuron2020('Analyzing coherence',0,'Color',[0.720000 0.530000 0.040000]); pause(0.25);
    % Run analysis and output a structure containing all the analyzed data.
    [ComparisonData] = AnalyzeData_Neuron2020;
    multiWaitbar_Neuron2020('CloseAll');
end

%% Informational figures with function dependencies for the various analysis and the time per vessel.
% To view individual summary figures, change the value of line 72 to false. You will then be prompted to manually select
% any number of figures (CTL-A for all) inside any of the five folders. You can only do one animal at a time.
functionNames = {'MainScript_Neuron2020','StageOneProcessing_Neuron2020','StageTwoProcessing_Neuron2020','StageThreeProcessing_Neuron2020'};
functionList = {};
for a = 1:length(functionNames)
    [functionList] = GetFuncDependencies_Neuron2020(a,functionNames{1,a},functionNames,functionList);
end
DetermineVesselStatistics_Neuron2020(ComparisonData);

%% Individual figures can be re-run after the analysis has completed.
Fig7_Angle_Neuron2020(ComparisonData)
disp('MainScript Analysis - Complete'); disp(' ')
end

function [ComparisonData] = AnalyzeData_Neuron2020()
animalIDs = {'T72','T73','T74','T75','T76','T80','T81','T82','T83'};   % list of animal IDs
ComparisonData = [];   % pre-allocate the results structure as empty

%% BLOCK PURPOSE: [1] Analyze the cross-correlation between abs(whisker acceleration) and vessel diameter.
for b = 1:length(animalIDs)
    [ComparisonData] = AnalyzeXCorr_Neuron2020(animalIDs{1,b},ComparisonData);
    multiWaitbar_Neuron2020('Analyzing cross correlation','Value',b/length(animalIDs));
end

% BLOCK PURPOSE: [2] Analyze the spectral coherence between abs(whisker acceleration) and vessel diameter.
for c = 1:length(animalIDs)
    [ComparisonData] = AnalyzeCoherence_Neuron2020(animalIDs{1,c},ComparisonData);
    multiWaitbar_Neuron2020('Analyzing coherence','Value',c/length(animalIDs));
end

answer = questdlg('Would you like to save the analysis results structure?','','yes','no','yes');
if strcmp(answer,'yes')
    save('ComparisonData.mat','ComparisonData')
end
disp('Loading analysis results and generating figures...'); disp(' ')

end
