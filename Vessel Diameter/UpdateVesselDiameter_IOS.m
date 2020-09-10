%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

clear; clc; close all
% character list of all RawData files
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat'); 
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
%% for window cam files
% draw ROIs along vein
DrawVesselROIs_IOS(procDataFileIDs)
% calculate FWHM for each file
CalcVesselDiameterFWHM_IOS(procDataFileIDs)
%% for processed data files
UpdateVesselProcDataFiles_IOS(procDataFileIDs)
