function StageOneProcessing_IOS(fileNames, trackWhiskers)
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
%   Outputs: A processed LabVIEWRawData file for each filename that is saved to the current directory.
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

%% BLOCK PURPOSE: [1] Preparing to create LabVIEWRawData files.
disp('Analyzing Block [1] Preparing to create LabVIEWRawData file(s).'); disp(' ')
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
    [~, ~, ~, fileID] = GetFileInfo_IOS(indFile);
    
    % Determine if a LabVIEWRawData file has already been created for this file. If it has, skip it
    fileExist = ls(['*' fileID '_LabVIEWRawData.mat']);
    if isempty(fileExist)
        
        %% BLOCK PURPOSE: [2] Import .tdms data (All channels).
        disp('Analyzing Block [2] Importing .tdms data from all channels.'); disp(' ')
        trialData = ReadInTDMSWhiskerTrials_IOS([fileID '.tdms']);
        
        dataRow = strcmp(trialData.data.names, 'Neural_LH');  
        cortical_LH = trialData.data.vals(dataRow,:) / str2double(trialData.amplifierGain);
        
        dataRow = strcmp(trialData.data.names, 'Force_Sensor'); 
        forceSensor = trialData.data.vals(dataRow,:);
        
        dataRow = strcmp(trialData.data.names, 'LPadSol'); 
        LPadSol = gt(trialData.data.vals(dataRow,:), 0.5)*1;   % Amplitude is 1 
        
        dataRow = strcmp(trialData.data.names, 'RPadSol'); 
        RPadSol = gt(trialData.data.vals(dataRow,:), 0.5)*2;   % Amplitude is 1
        
        dataRow = strcmp(trialData.data.names, 'AudSol'); 
        AudSol = gt(trialData.data.vals(dataRow,:), 0.5)*3;   % Amplitude is 1
        
        dataRow = strcmp(trialData.data.names, 'Respiration'); 
        hippocampus = trialData.data.vals(dataRow,:) / str2double(trialData.amplifierGain);
        
        dataRow = strcmp(trialData.data.names, 'EMG');
        EMG = trialData.data.vals(dataRow,:) / str2double(trialData.amplifierGain);
        
        dataRow = strcmp(trialData.data.names, 'Neural_RH'); 
        cortical_RH = trialData.data.vals(dataRow,:) / str2double(trialData.amplifierGain);
               
        solenoids = LPadSol + RPadSol + AudSol;

        %% BLOCK PURPOSE: [3] Start Whisker tracker.
        disp('Analyzing Block [3] Starting whisker tracking.'); disp(' ')
        if trackWhiskers == true
            [whiskerAngle] = WhiskerTrackerParallel_IOS(fileID);
            inds = isnan(whiskerAngle) == 1;
            whiskerAngle(inds) = [];
        else
            whiskerAngle = [];
        end
        
        %% BLOCK PURPOSE: [4] Save the notes and data.
        disp('Analyzing Block [4] Evaluating data to save to LabVIEWRawData file.'); disp(' ')
        % notes - all variables are descriptive
        LabVIEWRawData.notes.experimenter = trialData.experimenter;
        LabVIEWRawData.notes.animalID = trialData.animalID;
        LabVIEWRawData.notes.hemisphere = trialData.hemisphere;
        LabVIEWRawData.notes.solenoidPSI = str2double(trialData.solenoidPSI);
        LabVIEWRawData.notes.isofluraneTime = str2double(trialData.isofluraneTime);
        LabVIEWRawData.notes.sessionID = trialData.sessionID;
        LabVIEWRawData.notes.amplifierGain = str2double(trialData.amplifierGain);
        LabVIEWRawData.notes.CBVCamSamplingRate = str2double(trialData.CBVCamSamplingRate);
        LabVIEWRawData.notes.whiskerCamSamplingRate = str2double(trialData.whiskerCamSamplingRate);
        LabVIEWRawData.notes.webCamSamplingRate = str2double(trialData.webCamSamplingRate);
        LabVIEWRawData.notes.pupilCamSamplingRate = str2double(trialData.pupilCamSamplingRate);
        LabVIEWRawData.notes.analogSamplingRate = str2double(trialData.analogSamplingRate);
        LabVIEWRawData.notes.trialDuration_sec = str2double(trialData.trialDuration_Sec);
        LabVIEWRawData.notes.CBVCamPixelWidth = str2double(trialData.CBVCamPixelWidth);
        LabVIEWRawData.notes.CBVCamPixelHeight = str2double(trialData.CBVCamPixelHeight);
        LabVIEWRawData.notes.CBVCamBitDepth = str2double(trialData.CBVCamBitDepth);
        LabVIEWRawData.notes.pupilCamPixelWidth = str2double(trialData.pupilCamPixelWidth);
        LabVIEWRawData.notes.pupilCamPixelHeight = str2double(trialData.pupilCamPixelHeight);
        LabVIEWRawData.notes.whiskerCamPixelHeight = str2double(trialData.whiskerCamPixelHeight);
        LabVIEWRawData.notes.whiskerCamPixelWidth = str2double(trialData.whiskerCamPixelWidth);
        LabVIEWRawData.notes.CBVCamExposureTime_microsec = str2double(trialData.CBVCamExposureTime_microsec);
        LabVIEWRawData.notes.CBVCamBinning = str2double(trialData.CBVCamBinning);
        LabVIEWRawData.notes.numberDroppedPupilCamFrames = str2double(trialData.numberDroppedPupilCamFrames);
        LabVIEWRawData.notes.droppedPupilCamFrameIndex = trialData.droppedPupilCamFrameIndex;
        LabVIEWRawData.notes.numberDroppedWhiskCamFrames = str2double(trialData.numberDroppedWhiskCamFrames);
        LabVIEWRawData.notes.droppedWhiskCamFrameIndex = trialData.droppedWhiskCamFrameIndex;
        
        % Data
        LabVIEWRawData.data.cortical_LH = cortical_LH;
        LabVIEWRawData.data.cortical_RH = cortical_RH;
        LabVIEWRawData.data.hippocampus = hippocampus;
        LabVIEWRawData.data.forceSensor = forceSensor;
        LabVIEWRawData.data.EMG = EMG;
        LabVIEWRawData.data.whiskerAngle = whiskerAngle;
        LabVIEWRawData.data.solenoids = solenoids;
     
        disp(['File Created. Saving LabVIEWRawData File ' num2str(a) '...']); disp(' ')
        save([trialData.animalID '_' trialData.hemisphere '_' fileID '_LabVIEWRawData'], 'LabVIEWRawData')
    else
        disp('File already exists. Continuing...'); disp(' ')
    end
end

disp('IOS Stage One Processing - Complete.'); disp(' ')

end
