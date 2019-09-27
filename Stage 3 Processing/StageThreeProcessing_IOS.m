%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: 1) Categorize data using previously processed ProcData data structures, add 'flags'  
%            2) Create RestData structure that contains periods of rest.
%            3) Create EventData structure that contains periods after stimuli and whisks.
%            4) Uses periods when animal is not being stimulated or moving to establish a 
%               baseline for a given session of imaging.
%            5) Normalizes the different data structures.
%________________________________________________________________________________________________________________________
%
%   Inputs: 1) Select all _ProcData files from all days. Follow the command window prompts.
%           2) Select one single _RawData file for the animal information. The ProcData files are
%              already in the list and will be used to run the function.
%           3) No inputs. ProcData files already loaded.
%           4) No inputs. RestData.mat is already loaded.
%           5) No inputs. RestData.mat and EventData.mat are already loaded.
%
%   Outputs: 1) Additions to the ProcData structure including flags and scores.
%            2) A RestData.mat structure with periods of rest.
%            3) A EventData.mat structure with event-related information.
%            4) Baselines.mat containing the baselines for individual resting periods.
%            5) Creates NormData in the rest/event structures.
%
%   Last Revised: October 5th, 2018
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
clc;
clear;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')

% Character list of all RawData files
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);

% Character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
[animalID, ~, ~] = GetFileInfo_IOS(procDataFileIDs(1,:));

targetMinutes = 30;
baselineType = 'manualSelection';
dataTypes = {'CBV', 'cortical_LH', 'cortical_RH', 'hippocampus', 'EMG'};
updatedDataTypes = {'CBV', 'CBV_HbT', 'cortical_LH', 'cortical_RH', 'hippocampus', 'EMG'};
neuralDataTypes = {'cortical_LH', 'cortical_RH', 'hippocampus'};

%% BLOCK PURPOSE: [1] Categorize data 
disp('Analyzing Block [1] Categorizing data.'); disp(' ')
for a = 1:size(procDataFileIDs, 1)
    procDataFile = procDataFileIDs(a, :);
    disp(['Analyzing file ' num2str(a) ' of ' num2str(size(procDataFileIDs, 1)) '...']); disp(' ')
    CategorizeData_IOS(procDataFile)
end

%% BLOCK PURPOSE: [2] Create RestData data structure
disp('Analyzing Block [2] Create RestData struct for CBV and neural data.'); disp(' ')
[RestData] = ExtractRestingData_IOS(procDataFileIDs, dataTypes);

%% BLOCK PURPOSE: [4] Create Baselines data structure
disp('Analyzing Block [4] Create Baselines struct for CBV and neural data.'); disp(' ')
trialDuration_sec = 900;
[RestingBaselines] = CalculateRestingBaselines_IOS(animalID, targetMinutes, trialDuration_sec, RestData);

%% BLOCK PURPOSE: [5] Manually select files for custom baseline calculation
disp('Analyzing Block [5] Manually select files for custom baseline calculation.'); disp(' ')
[RestingBaselines] = CalculateManualRestingBaselines_IOS(animalID, procDataFileIDs, RestData, RestingBaselines);

%% BLOCK PURPOSE [6] Add delta HbT field to each processed data file
disp('Analyzing Block [6] Adding delta HbT to each ProcData file.'); disp(' ')
UpdateTotalHemoglobin_IOS(procDataFileIDs, RestingBaselines, baselineType)

%% BLOCK PURPOSE: [7] Re-create the RestData structure now that HbT is available
disp('Analyzing Block [7] Creating RestData struct for CBV and neural data.'); disp(' ')
[RestData] = ExtractRestingData_IOS(procDataFileIDs, updatedDataTypes);

%% BLOCK PURPOSE: [8] Create the EventData structure for CBV and neural data
disp('Analyzing Block [8] Create EventData struct for CBV and neural data.'); disp(' ')
[EventData] = ExtractEventTriggeredData_IOS(procDataFileIDs, updatedDataTypes);

%% BLOCK PURPOSE: [9] Normalize RestData and EventData structures by the resting baseline
disp('Analyzing Block [9] Normalizing RestData and EventData structures by the resting baseline.'); disp(' ')
[RestData] = NormBehavioralDataStruct_IOS(RestData, RestingBaselines, baselineType);
save([animalID '_RestData.mat'], 'RestData')

[EventData] = NormBehavioralDataStruct_IOS(EventData, RestingBaselines, baselineType);
save([animalID '_EventData.mat'], 'EventData')

%% BLOCK PURPOSE: [10] Analyze the spectrogram for each session.
disp('Analyzing Block [10] Analyzing the spectrogram for each file and normalizing by the resting baseline.'); disp(' ')
CreateTrialSpectrograms_IOS(rawDataFileIDs, neuralDataTypes);

% Find spectrogram baselines for each day
specDirectory = dir('*_SpecData.mat');
specDataFiles = {specDirectory.name}';
specDataFileIDs = char(specDataFiles);
[RestingBaselines] = CalculateSpectrogramBaselines_IOS(animalID, neuralDataTypes, trialDuration_sec, specDataFileIDs, RestingBaselines, baselineType);

% Normalize spectrogram by baseline
NormalizeSpectrograms_IOS(specDataFileIDs, neuralDataTypes, RestingBaselines);

%% BLOCK PURPOSE [11] Generate single trial figures
% disp('Analyzing Block [11] Gennerating single trial summary figures'); disp(' ')
% saveFigs = 'y';
% GenerateSingleFigures_IOS(procDataFileIDs, RestingBaselines, baselineType, saveFigs)

disp('Stage Three Processing - Complete.'); disp(' ')
