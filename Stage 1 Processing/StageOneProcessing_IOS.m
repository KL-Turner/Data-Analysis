function StageOneProcessing_IOS(fileNames, trackWhiskers)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Data acquired during trials must be in a form that MATLAB can work with easily. This code converts the
%            various forms of data listed below into MATLAB structures that can be easily manipulated.
%
%            .bin - Cameras
%            .tdms - Digital and Analog Data
%            .tdms_index - Index for the LabVIEW data in the .tdms file
%________________________________________________________________________________________________________________________
%
%   Inputs: fileNames - [cell array] list of filames with the extension _WhiskerCam.bin.
%           TrackWhiskers - [binary] tells code whether to track the whiskers or not.
%
%   Outputs: A processed RawData file for each filename that is saved to the current directory.
%
%   Last Revised: October 3rd, 2018
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Note: This function can be run independently of a main script.
% Clear the workspace variables and command window.
clc;
clear;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')

% If there are no inputs to the function, it asks the user to load all files with a _WhiskerCam.bin extension.
if nargin == 0
    fileNames = uigetfile('*_WhiskerCam.bin', 'MultiSelect', 'on');   % CTL-A to select all files
end

% Default setting - if you automatically play the function, it will track the whiskers.
% To debug the code without taking the time to track whiskers, set trackWhiskers = 0.
if nargin < 2
    trackWhiskers = 0;
end

% Control for single file instead of a list.
if iscell(fileNames) == 0
    fileName = fileNames;
    fileNames = 1;
end

%% BLOCK PURPOSE: [1] Preparing to create RawData files.
disp('Analyzing Block [1] Preparing to create RawData file.'); disp(' ')
% Load in each file one at a time, looping through the list.
for fileNumber = 1:length(fileNames)
    disp(['Analyzing file ' num2str(fileNumber) ' of ' num2str(length(fileNames)) '...']); disp(' ')
    
    % Adapt to list or single file. The purpose of this is control the way uigetfile handles an instance of a
    % single file input (character string) vs. multiple files, which it puts in cells.
    if iscell(fileNames) == 1
        indFile = fileNames{fileNumber};
    else
        indFile = fileName;
    end
    
    % Pull out the file ID for the file - this is the numerical string after the animal name/hem
    [~, ~, ~, fileID] = GetFileInfo_IOS(indFile);
    
    % Determine if a RawData file has already been created for this file. If it has, skip it.
    fileExist = ls(['*' fileID '_RawData.mat']);
    if not(isempty(fileExist))
        disp('File already exists. Continuing...'); disp(' ')
    end
    
    %% BLOCK PURPOSE: [2] Import .tdms data (All channels);
    disp('Analyzing Block [2] Importing .tdms data from all channels.'); disp(' ')
    trialData = ReadInTDMSWhiskerTrials_IOS([fileID '.tdms']);
    
    dataRow = strcmp(trialData.Data.Names, 'Neural_LH');   % Left hem neural data
    Neural_LH = trialData.Data.Vals(dataRow,:) / str2double(trialData.amplifierGain);
    
    dataRow = strcmp(trialData.Data.Names, 'Neural_RH');   % Right hem neural data
    Neural_RH = trialData.Data.Vals(dataRow,:) / str2double(trialData.amplifierGain);
    
    dataRow = strcmp(trialData.Data.Names, 'EMG');   % EMG data
    EMG = trialData.Data.Vals(dataRow,:) / str2double(trialData.amplifierGain);
    
%     dataRow = strcmp(trialData.Data.Names, 'Respiration');   % EMG data
%     Respiration = trialData.Data.Vals(dataRow,:) / str2double(trialData.amplifierGain);
    
    dataRow = strcmp(trialData.Data.Names, 'Force_Sensor');   % Force sensor data
    Force_Sensor = trialData.Data.Vals(dataRow,:);
    
    dataRow = strcmp(trialData.Data.Names, 'Solenoid_LeftPad');   % Left whisker pad puffs
    Left_Pad_Solenoid = gt(trialData.Data.Vals(dataRow,:), 0.5)*1;   % Amplitude is 1
    
    dataRow = strcmp(trialData.Data.Names, 'Solenoid_RightPad');   % Right whisker pad puffs
    Right_Pad_Solenoid = gt(trialData.Data.Vals(dataRow,:), 0.5)*2;   % Amplitude is 2
    
%     dataRow = strcmp(trialData.Data.Names, 'Solenoid_Tail');   % Puffs to the tail
%     Tail_Solenoid = gt(trialData.Data.Vals(dataRow,:), 0.5)*3;   % Amplitude is 3
    
    dataRow = strcmp(trialData.Data.Names, 'Solenoid_Auditory');   % Auditory puffs (in distance)
    Auditory_Solenoid = gt(trialData.Data.Vals(dataRow,:), 0.5)*4;   % Amplitude is 4
    
    % Add together puffing arrays (no two puffs occur simultaneously)
    Solenoids = Left_Pad_Solenoid + Right_Pad_Solenoid + Auditory_Solenoid;
    
    %% BLOCK PURPOSE: [3] Start Whisker tracker
    disp('Analyzing Block [3] Starting whisker tracker.'); disp(' ')
    if trackWhiskers   % Logical statement (if trackWhiskers == 1)
        [WhiskerAngle] = WhiskerTrackerParallel_IOS(fileID);
        inds = isnan(WhiskerAngle) == 1;
        WhiskerAngle(inds) = [];
    else
        WhiskerAngle = [];
    end
    
    %% BLOCK PURPOSE: [4] Evaluate Data and Save
    disp('Analyzing Block [4] Evaluating data to save to RawData file.'); disp(' ')
    % Notes - all variables are descriptive
    RawData.Notes.experimenter = trialData.experimenter;
    RawData.Notes.animalID = trialData.animalID;
    RawData.Notes.imagedHemisphere = trialData.imagedHemisphere;
    RawData.Notes.solenoidPressure_PSI = str2double(trialData.solenoidPressure_PSI);
    RawData.Notes.isofluraneTime_Military = str2double(trialData.isofluraneTime_Military);
    RawData.Notes.sessionID = trialData.sessionID;
    RawData.Notes.amplifierGain = str2double(trialData.amplifierGain);
    RawData.Notes.CBVCamSamplingRate = str2double(trialData.CBVCamSamplingRate);
    RawData.Notes.whiskerCamSamplingRate = str2double(trialData.whiskerCamSamplingRate);
    RawData.Notes.webCamSamplingRate = str2double(trialData.webCamSamplingRate);
    RawData.Notes.analogSamplingRate = str2double(trialData.analogSamplingRate);
    RawData.Notes.trialDuration_Seconds = str2double(trialData.trialDuration_Seconds);
    RawData.Notes.CBVCamPixelHeight = str2double(trialData.CBVCamPixelHeight);
    RawData.Notes.CBVCamPixelWidth = str2double(trialData.CBVCamPixelWidth);
    RawData.Notes.CBVCamBitDepth = str2double(trialData.CBVCamBitDepth);
    RawData.Notes.whiskerCamPixelHeight = str2double(trialData.whiskerCamPixelHeight);
    RawData.Notes.whiskerCamPixelWidth = str2double(trialData.whiskerCamPixelWidth);
    RawData.Notes.CBVCamExposureTime_Microseconds = str2double(trialData.CBVCamExposureTime_Microseconds);
    RawData.Notes.CBVCamBinning = trialData.CBVCamBinning;
    RawData.Notes.numberDroppedWhiskerCamFrames = str2double(trialData.numberDroppedWhiskerCamFrames);
    RawData.Notes.droppedWhiskerCamFrameIndex = trialData.droppedWhiskerCamFrameIndex;
    
    % Data
    RawData.Data.Neural_LH = Neural_LH;
    RawData.Data.Neural_RH = Neural_RH;
    RawData.Data.EMG = EMG;
    RawData.Data.Respiration = Respiration;
    RawData.Data.Force_Sensor = Force_Sensor;
    RawData.Data.Solenoids = Solenoids;
    RawData.Data.WhiskerAngle = WhiskerAngle;
    
    disp(['File Created. Saving RawData File ' num2str(fileNumber) '...']); disp(' ')
    save([trialData.animalID '_' trialData.imagedHemisphere '_' fileID '_RawData'], 'RawData')
end

disp('Stage One Processing - Complete.'); disp(' ')

end
