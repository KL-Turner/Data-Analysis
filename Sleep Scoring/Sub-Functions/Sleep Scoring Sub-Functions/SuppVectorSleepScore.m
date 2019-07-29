%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised: July 26th, 2019
%________________________________________________________________________________________________________________________

%% Character list of all ProcData files
clear; clc; close all
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

% Character name for the RestingBaselines structure
baselinesFileStruct = dir('*_RestingBaselines.mat');
baselinesFile = {baselinesFileStruct.name}';
baselinesFileID = char(baselinesFile);
load(baselinesFileID)

%%
AddSleepParameters_SVM(procDataFileIDs, RestingBaselines)

%%
CreateTrainingDataSet_SVM(procDataFileIDs, RestingBaselines)













