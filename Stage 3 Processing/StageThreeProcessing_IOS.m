%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 1) Categorize behavioral (rest,whisk,stim) data using previously processed data structures, add 'flags'  
%            2) Create a temporary RestData structure that contains periods of rest - use this for initial figures
%            3) Analyze neural data and create different spectrograms for each file's electrodes
%            4) Uses periods when animal is not being stimulated or moving to establish an initial baseline
%            5) Manually select awake files for a slightly different baseline not based on hard time vals
%            6) Use the best baseline to convert reflectance changes to total hemoglobin
%            7) Re-create the RestData structure now that we can deltaHbT
%            8) Create an EventData structure looking at the different data types after whisking or stimulation
%            9) Apply the resting baseline to each data type to create a percentage change 
%            10) Use the time indeces of the resting baseline file to apply a percentage change to the spectrograms
%            11) Use the time indeces of the resting baseline file to create a reflectance pixel-based baseline
%            12) Generate a summary figure for all of the analyzed and processed data
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
[animalID,~,~] = GetFileInfo_IOS(procDataFileIDs(1,:));

targetMinutes = 30;
timeOverride = 'y';
curDir = cd;
dirBreaks = strfind(curDir,'\');
curFolder = curDir(dirBreaks(end) + 1:end);
if strcmp(curFolder,'Bilateral Imaging') == true
    imagingType = 'bilateral';
elseif strcmp(curFolder,'Single Hemisphere') == true
    imagingType = 'single';
end
dataTypes = {'CBV','cortical_LH','cortical_RH','hippocampus','EMG','flow'};
updatedDataTypes = {'CBV','CBV_HbT','cortical_LH','cortical_RH','hippocampus','EMG','flow'};
neuralDataTypes = {'cortical_LH','cortical_RH','hippocampus'};

%% BLOCK PURPOSE: [1] Categorize data 
disp('Analyzing Block [1] Categorizing data.'); disp(' ')
for a = 1:size(procDataFileIDs, 1)
    procDataFile = procDataFileIDs(a, :);
    disp(['Analyzing file ' num2str(a) ' of ' num2str(size(procDataFileIDs, 1)) '...']); disp(' ')
    CategorizeData_IOS(procDataFile)
end

%% BLOCK PURPOSE: [2] Create RestData data structure
disp('Analyzing Block [2] Create RestData struct for CBV and neural data.'); disp(' ')
[RestData] = ExtractRestingData_IOS(procDataFileIDs,dataTypes,imagingType);

%% BLOCK PURPOSE: [3] Analyze the spectrogram for each session.
disp('Analyzing Block [3] Analyzing the spectrogram for each file and normalizing by the resting baseline.'); disp(' ')
CreateTrialSpectrograms_IOS(rawDataFileIDs,neuralDataTypes);

%% BLOCK PURPOSE: [4] Create Baselines data structure
disp('Analyzing Block [4] Create Baselines struct for CBV and neural data.'); disp(' ')
baselineType = 'setDuration';
trialDuration_sec = 900;
[RestingBaselines] = CalculateRestingBaselines_IOS(animalID,targetMinutes,trialDuration_sec,RestData);

% Find spectrogram baselines for each day
specDirectory = dir('*_SpecData.mat');
specDataFiles = {specDirectory.name}';
specDataFileIDs = char(specDataFiles);
[RestingBaselines] = CalculateSpectrogramBaselines_IOS(animalID,neuralDataTypes,trialDuration_sec,specDataFileIDs,RestingBaselines,baselineType);

% Normalize spectrogram by baseline
NormalizeSpectrograms_IOS(specDataFileIDs,neuralDataTypes,RestingBaselines);

%% BLOCK PURPOSE: [5] Manually select files for custom baseline calculation
disp('Analyzing Block [5] Manually select files for custom baseline calculation.'); disp(' ')
if strcmp(timeOverride,'n') == true
    [RestingBaselines] = CalculateManualRestingBaselines_IOS(animalID,procDataFileIDs,RestData,RestingBaselines,imagingType);
else
    [RestingBaselines] = CalculateManualRestingBaselinesWithTimeIndeces_IOS(imagingType);
end

%% BLOCK PURPOSE [6] Add delta HbT field to each processed data file
disp('Analyzing Block [6] Adding delta HbT to each ProcData file.'); disp(' ')
updatedBaselineType = 'manualSelection';
UpdateTotalHemoglobin_IOS(procDataFileIDs, RestingBaselines,updatedBaselineType,imagingType)

%% BLOCK PURPOSE: [7] Re-create the RestData structure now that HbT is available
disp('Analyzing Block [7] Creating RestData struct for CBV and neural data.'); disp(' ')
[RestData] = ExtractRestingData_IOS(procDataFileIDs,updatedDataTypes,imagingType);

%% BLOCK PURPOSE: [8] Create the EventData structure for CBV and neural data
disp('Analyzing Block [8] Create EventData struct for CBV and neural data.'); disp(' ')
[EventData] = ExtractEventTriggeredData_IOS(procDataFileIDs,updatedDataTypes,imagingType);

%% BLOCK PURPOSE: [9] Normalize RestData and EventData structures by the resting baseline
disp('Analyzing Block [9] Normalizing RestData and EventData structures by the resting baseline.'); disp(' ')
[RestData] = NormBehavioralDataStruct_IOS(RestData,RestingBaselines,updatedBaselineType);
save([animalID '_RestData.mat'],'RestData','-v7.3')

[EventData] = NormBehavioralDataStruct_IOS(EventData,RestingBaselines,updatedBaselineType);
save([animalID '_EventData.mat'],'EventData','-v7.3')

%% BLOCK PURPOSE: [10] Analyze the spectrogram baseline for each session.
disp('Analyzing Block [10] Analyzing the spectrogram for each file and normalizing by the resting baseline.'); disp(' ')
% Find spectrogram baselines for each day
specDirectory = dir('*_SpecData.mat');
specDataFiles = {specDirectory.name}';
specDataFileIDs = char(specDataFiles);
[RestingBaselines] = CalculateSpectrogramBaselines_IOS(animalID,neuralDataTypes,trialDuration_sec,specDataFileIDs,RestingBaselines,updatedBaselineType);

% Normalize spectrogram by baseline
NormalizeSpectrograms_IOS(specDataFileIDs,neuralDataTypes,RestingBaselines);

% Create a structure with all spectrograms for convenient analysis further downstream
CreateAllSpecDataStruct_IOS(animalID,neuralDataTypes)

%% BLOCK PURPOSE: [11] Generate pixel baseline from WindowCam.mat files
% disp('Analyzing Block [11] Generating pixel-based resting baselines for reflectance data'); disp(' ')
% if strcmp(imagingType,'single') == true
%     [RestingBaselines] = CalculatePixelBaselines_IOS(procDataFileIDs,RestingBaselines,baselineType);
% end

%% BLOCK PURPOSE [12] Generate single trial figures
disp('Analyzing Block [12] Generating single trial summary figures'); disp(' ')
% saveFigs = 'y';
% GenerateSingleFigures_IOS(procDataFileIDs,RestingBaselines,updatedBaselineType,saveFigs,imagingType)

disp('Stage Three Processing - Complete.'); disp(' ')
