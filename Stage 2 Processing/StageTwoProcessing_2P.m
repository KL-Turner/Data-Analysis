%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: 1) Pull vessel notes from Excel sheet.
%            2) Analyze vessel diameter.
%            3) Analyze neural bands, force sensors, and whisker motion.
%            4) Correct LabVIEW time offset.
%            5) Combine LabVIEWData and MScan Data files to ceate MergedData. 
%________________________________________________________________________________________________________________________
%
%   Inputs: Scan the workspace for filenames *_LabVIEWData.mat' and '*_MScanData.mat'
%
%   Outputs: Updated fields in LabVIEWData.mat, new/updated fields of MScanData.mat and MergedData.mat.       
%
%   Last Revised: February 21st, 2019    
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
clc;
clear;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')

[~, ~, ~, ~, ~, ~, ~, ~] = LoadDataStructs_2P;

msExcel_File = uigetfile('*.xlsx');
Analyze2PDataNotes(msExcel_File);

disp('Block [0] structs loaded.'); disp(' ')

%% BLOCK PURPOSE: [1] Use ms Excel sheet to create MScanData.mat files with vessel information
disp('Analyzing Block [1] Pulling vessel notes from Excel sheet.'); disp(' ')
Analyze2PDataNotes(msExcel_File);

%% BLOCK PURPOSE: [2] Analyze vessel diameter and add it to MScanData.mat
disp('Analyzing Block [2] Analyzing vessel diameter.'); disp(' ')
mscanDirectory = dir('*_MScanData.mat');
mscanDataFiles = {mscanDirectory.name}';
mscanDataFiles = char(mscanDataFiles);
Analyze2PDiameter(mscanDirectory);

%% BLOCK PURPOSE: [3] Process neural, whiskers, and force sensor data
disp('Analyzing Block [3] Analyzing neural bands, force sensors, and whiskers.'); disp(' ')
labviewDirectory = dir('*_LabVIEWData.mat');
labviewDataFiles = {labviewDirectory.name}';
labviewDataFiles = char(labviewDataFiles);
Process2PDataFiles(labviewDataFiles, mscanDataFiles)

%% BLOCK PURPOSE: [4] Analyze vessel diameter and add it to MScanData.mat
disp('Analyzing Block [4] Correcting LabVIEW time offset.'); disp(' ')
CorrectLabVIEWOffset(labviewDataFiles, mscanDataFiles)

%% BLOCK PURPOSE: [5] Analyze vessel diameter and add it to MScanData.mat
disp('Analyzing Block [5] Combing LabVIEWData and MScan Data files to ceate MergedData.'); disp(' ')
CombineLabVIEWMScanFiles(labviewDataFiles, mscanDataFiles)

disp('Two Photon Stage Two Processing - Complete.'); disp(' ')
