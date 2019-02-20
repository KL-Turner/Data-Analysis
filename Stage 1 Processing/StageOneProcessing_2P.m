function StageOneProcessing_2P(fileNames, trackWhiskers)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
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
%   Inputs: fileNames - [cell array] list of filames with the extension _dalsa.bin.
%           TrackWhiskers - [binary] tells code whether to track the whiskers or not.
%
%   Outputs: A processed RawData file for each filename that is saved to the current directory.
%
%   Last Revised: January 18th, 2019
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
    trackWhiskers = 1;
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
    [~, ~, ~, fileID] = GetFileInfo(indFile);
    
    % Determine if a RawData file has already been created for this file. If it has, skip it.
    fileExist = ls(['*' fileID '_RawData.mat']);
    if not(isempty(fileExist))
        disp('File already exists. Continuing...'); disp(' ')
    end
    
    %% BLOCK PURPOSE: [2] Import .tdms data (All channels);
    disp('Analyzing Block [2] Importing .tdms data from all channels.'); disp(' ')
    trialData = ReadInTDMSWhiskerTrials_2P([fileID '.tdms']);
    
    dataRow = strcmp(trialData.Data.Names, 'Force_Sensor');   % Force sensor data
    Force_Sensor = trialData.Data.Vals(dataRow,:);
    
    %% BLOCK PURPOSE: [3] Start Whisker tracker
    disp('Analyzing Block [3] Starting whisker tracker.'); disp(' ')
    if trackWhiskers   % Logical statement (if trackWhiskers == 1)
        [WhiskerAngle] = WhiskerTrackerParallel(fileID);
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
    RawData.Notes.isofluraneTime_Military = str2double(trialData.isofluraneTime_Military);
    RawData.Notes.sessionID = trialData.sessionID;
    RawData.Notes.amplifierGain = str2double(trialData.amplifierGain);
    RawData.Notes.whiskerCamSamplingRate = str2double(trialData.whiskerCamSamplingRate);
    RawData.Notes.analogSamplingRate = str2double(trialData.analogSamplingRate);
    RawData.Notes.trialDuration_Seconds = str2double(trialData.trialDuration_Seconds);
    RawData.Notes.whiskerCamPixelHeight = str2double(trialData.whiskerCamPixelHeight);
    RawData.Notes.whiskerCamPixelWidth = str2double(trialData.whiskerCamPixelWidth);
    RawData.Notes.numberDroppedWhiskerCamFrames = str2double(trialData.numberDroppedWhiskerCamFrames);
    RawData.Notes.droppedWhiskerCamFrameIndex = trialData.droppedWhiskerCamFrameIndex;
    
    % Data
    RawData.Data.Force_Sensor = Force_Sensor;
    RawData.Data.WhiskerAngle = WhiskerAngle;
    
    disp(['File Created. Saving RawData File ' num2str(fileNumber) '...']); disp(' ')
    save([trialData.animalID '_' trialData.imagedHemisphere '_' fileID '_RawData'], 'RawData')
end

disp('Stage One Processing - Complete.'); disp(' ')

end
