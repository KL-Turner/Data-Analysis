function [] = MainScript_SSS()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panels for Garborg et al. Manuscript in preparation
%
% Functions used to pre-process the original data are located in the folder "Pre-Processing Functions"
% Functions used to analyze data for figures are located in the folder "Data Analysis Functions"
% Functions optained from 3rd party are located in the folder "Shared Functions"
%________________________________________________________________________________________________________________________

clear; clc;
% Verify code repository and data are in the current directory/added path
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
% Add root folder to Matlab's working directory
addpath(genpath(rootFolder))
multiWaitbar('CloseAll');
% Analysis subfunctions
runAnalysis = true;
if runAnalysis == true
    dataLocation = [rootFolder delim 'Analysis Structures'];
    cd(dataLocation)
    AnalyzeEvokedResponses_Handler_GarborgTBD(rootFolder,delim,false)
    AnalyzeSagittalSinusGFP_Handler_GarborgTBD(rootFolder,delim,false)
    multiWaitbar('CloseAll');
    cd(rootFolder)
end
% Main figures
disp('Loading analysis results and generating figures...'); disp(' ')
saveFigs = false;
AwakeEvokedResponses_GarborgTBD(rootFolder,saveFigs,delim)
% Figure Panel 8 is schematic diagram
% if exist('VideoS1.mp4','file') ~= 2
 
end
