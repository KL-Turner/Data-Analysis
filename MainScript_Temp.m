function [] = MainScript_APOE()
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
if exist('AnalysisResults_APOE.mat','file') ~= 2 || strcmp(rerunAnalysis,'y') == true
    multiWaitbar('Analyzing whisking-evoked responses',0,'Color','P'); pause(0.25);
    [AnalysisResults] = AnalyzeData(rootFolder);
    multiWaitbar('CloseAll');
else
    disp('Loading analysis results and generating figures...'); disp(' ')
    load('AnalysisResults_APOE.mat','-mat')
end
%% generate figures
[AnalysisResults] = WhiskEvoked_APOE4(rootFolder,saveFigs,delim,AnalysisResults); %#ok<*NASGU>
%% fin.
disp('MainScript Analysis - Complete'); disp(' ')
end

function [AnalysisResults_APOE] = AnalyzeData(rootFolder)
% IOS animal IDs
expGroups = {'Dural','Capillary'};
% saveFigs = 'y';
if exist('AnalysisResults_Wenke.mat','file') == 2
    load('AnalysisResults_Wenke.mat','-mat')
else
    AnalysisResults_APOE = [];
end
% determine waitbar length
waitBarLength = 0;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    waitBarLength = waitBarLength + length(animalIDs);
end
%% Analyze the whisking-evoked arteriole diameter
runFromStart = 'n';
cc = 1;
for aa = 1:length(expGroups)
    folderList = dir(expGroups{1,aa});
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(AnalysisResults_APOE,(animalIDs{1,bb})) == false || isfield(AnalysisResults_APOE.(animalIDs{1,bb}),'EvokedAvgs') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults_APOE] = AnalyzeVesselEvokedResponses(animalIDs{1,bb},expGroups{1,aa},rootFolder,AnalysisResults_APOE);
        end
        multiWaitbar('Analyzing whisking-evoked responses','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
%% fin.
disp('Loading analysis results and generating figures...'); disp(' ')
end