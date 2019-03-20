%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpse:  1) Additions to the MergedData structure including flags and scores.
%            2) A RestData.mat structure with periods of rest.
%            3) A EventData.mat structure with event-related information.
%            4) Find the resting baseline for vessel diameter and neural data.
%________________________________________________________________________________________________________________________
%
%   Inputs: MergedData files, followed by newly created RestData and EventData structs for normalization 
%           by the unique day and vessel's resting baseline.
%
%   Outputs: 1) Additions to the MergedData structure including behavioral flags and scores.
%            2) A RestData.mat structure with periods of rest.
%            3) A EventData.mat structure with event-related information.
%            4) A Baselines.mat structure with resting baselines.
%
%   Last Revised: February 21st, 2019
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
clc;
clear;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')

[animalID, ~, ~, EventData, RestData, ~, ~, ~, ~] = LoadDataStructs_2P; %#ok<ASGLU>

mergedDirectory = dir('*_MergedData.mat');
mergedDataFiles = {mergedDirectory.name}';
mergedDataFiles = char(mergedDataFiles);
dataTypes = {'Vessel_Diameter', 'DeltaBand_Power', 'ThetaBand_Power', 'AlphaBand_Power', 'BetaBand_Power', 'GammaBand_Power', 'MUA_Power'};

disp('Block [0] structs loaded.'); disp(' ')

%% BLOCK PURPOSE: [1] Categorize data 
disp('Analyzing Block [1] Categorizing behavioral data, adding flags to MergedData structures.'); disp(' ')
for fileNumber = 1:size(mergedDataFiles, 1)
    fileName = mergedDataFiles(fileNumber, :);
    disp(['Analyzing file ' num2str(fileNumber) ' of ' num2str(size(mergedDataFiles, 1)) '...']); disp(' ')
    CategorizeData_2P(fileName)
end

%% BLOCK PURPOSE: [2] Create RestData data structure
disp('Analyzing Block [2] Creating RestData struct for vessels and neural data.'); disp(' ')
[RestData] = ExtractRestingData_2P(mergedDataFiles, dataTypes);
    
%% BLOCK PURPOSE: [3] Create EventData data structure
disp('Analyzing Block [3] Creating EventData struct for vessels and neural data.'); disp(' ')
[EventData] = ExtractEventTriggeredData_2P(mergedDataFiles, dataTypes);

%% BLOCK PURPOSE: [4] Create Baselines data structure
disp('Analyzing Block [4] Finding the resting baseline for vessel diameter and neural data.'); disp(' ')
targetMinutes = 15;
[RestingBaselines] = CalculateRestingBaselines_2P(animalID, targetMinutes, RestData);

%% BLOCK PURPOSE: [5]: Determine vessel statistics
disp('Analyzing Block [5] Finding the different vessel statistics for each animal.'); disp(' ')



disp('Two Photon Stage Three Processing - Complete.'); disp(' ')
