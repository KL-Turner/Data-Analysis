function StageOneProcessing_SlowOscReview2019(fileNames, trackWhiskers)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Data acquired during trials must be in a form that Matlab can work with easily. This code converts the
%            various forms of data listed below into MATLAB structures that can be easily manipulated.
%
%            .bin - Cameras
%            .tdms - Digital and Analog Data
%            .tdms_index - Index for the LabVIEW data in the .tdms file
%________________________________________________________________________________________________________________________
%
%   Inputs: fileNames - [cell array] list of filames with the extension '_WhiskerCam.bin'.
%           TrackWhiskers - [binary] tells code whether to track the whiskers or not.
%
%   Outputs: A processed LabVIEWData file for each filename that is saved to the current directory.
%
%   Last Revised: March 21st, 2019
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Note: This function can be run independently of a main script
% Clear the workspace variables and command window
clc;
clear;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')

% If there are no inputs to the function, it asks the user to load all files with a '_WhiskerCam.bin' extension
if nargin == 0
    fileNames = uigetfile('*_WhiskerCam.bin', 'MultiSelect', 'on');   % CTL-A to select all files
end

% Default setting - if you automatically play the function, it will track the whiskers
% To debug the code without taking the time to track whiskers, set trackWhiskers = 0
if nargin < 2
    trackWhiskers = 1;
end

%% BLOCK PURPOSE: [1] Preparing to create LabVIEWData files.
disp('Analyzing Block [1] Preparing to create LabVIEWData file(s).'); disp(' ')
% Load in each file one at a time, looping through the list
for a = 1:length(fileNames)
    disp(['Analyzing file ' num2str(a) ' of ' num2str(length(fileNames)) '...']); disp(' ')
    % Adapt to list or single file. The purpose of this is control the way uigetfile handles an instance of a
    % single file input (character string) vs. multiple files, which it puts in cells
    if iscell(fileNames) == 1
        indFile = fileNames{a};
    else
        indFile = fileName;
    end
    
    % Pull out the file ID for the file - this is the numerical string after the animal name/hemisphere
    [~, ~, ~, fileID] = GetFileInfo_SlowOscReview2019(indFile);
    
    % Determine if a LabVIEWData file has already been created for this file. If it has, skip it
    fileExist = ls(['*' fileID '_LabVIEWData.mat']);
    if isempty(fileExist)
        
        
        %% BLOCK PURPOSE: [2] Import .tdms data (All channels).
        disp('Analyzing Block [2] Importing .tdms data from all channels.'); disp(' ')
        trialData = ReadInTDMSWhiskerTrials_SlowOscReview2019([fileID '.tdms']);
        
        dataRow = strcmp(trialData.data.names, 'Force_Sensor');   % Force sensor data
        forceSensor = trialData.data.vals(dataRow,:);
        
        %% BLOCK PURPOSE: [3] Start Whisker tracker.
        disp('Analyzing Block [3] Starting whisker tracking.'); disp(' ')
        if trackWhiskers
            [whiskerAngle] = WhiskerTrackerParallel_SlowOscReview2019(fileID);
            inds = isnan(whiskerAngle) == 1;
            whiskerAngle(inds) = [];
        else
            whiskerAngle = [];
        end
        
        %% BLOCK PURPOSE: [4] Save the notes and data.
        disp('Analyzing Block [4] Evaluating data to save to LabVIEWData file.'); disp(' ')
        % notes - all variables are descriptive
        LabVIEWData.notes.experimenter = trialData.experimenter;
        LabVIEWData.notes.animalID = trialData.animalID;
        LabVIEWData.notes.imagedHemisphere = trialData.imagedHemisphere;
        LabVIEWData.notes.isofluraneTime_Military = str2double(trialData.isofluraneTime_Military);
        LabVIEWData.notes.sessionID = trialData.sessionID;
        LabVIEWData.notes.amplifierGain = str2double(trialData.amplifierGain);
        LabVIEWData.notes.whiskerCamSamplingRate_Hz = str2double(trialData.whiskerCamSamplingRate_Hz);
        LabVIEWData.notes.analogSamplingRate_Hz = str2double(trialData.analogSamplingRate_Hz);
        LabVIEWData.notes.trialDuration_Seconds = str2double(trialData.trialDuration_Seconds);
        LabVIEWData.notes.whiskerCamPixelHeight = str2double(trialData.whiskerCamPixelHeight);
        LabVIEWData.notes.whiskerCamPixelWidth = str2double(trialData.whiskerCamPixelWidth);
        LabVIEWData.notes.numberDroppedWhiskerCamFrames = str2double(trialData.numberDroppedWhiskerCamFrames);
        LabVIEWData.notes.droppedWhiskerCamFrameIndex = trialData.droppedWhiskerCamFrameIndex;
        
        % Data
        LabVIEWData.data.forceSensor = forceSensor;
        LabVIEWData.data.whiskerAngle = whiskerAngle;
        
        % Checklist for analysis steps - debugging purposes
        LabVIEWData.notes.checklist.processData = false;
        LabVIEWData.notes.checklist.offsetCorrect = false;
        
        disp(['File Created. Saving LabVIEWData File ' num2str(a) '...']); disp(' ')
        save([trialData.animalID '_' trialData.imagedHemisphere '_' fileID '_LabVIEWData'], 'LabVIEWData')
    else
        disp('File already exists. Continuing...'); disp(' ')
    end
end

disp('Two Photon Stage One Processing - Complete.'); disp(' ')

end
