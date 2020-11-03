function StageOneProcessing_IOS_eLife2020(fileNames,trackWhiskers)
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

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Note: This function can be run independently of a main script
% Clear the workspace variables and command window
clc;
clear;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')
% If there are no inputs to the function, it asks the user to load all files with a '_WhiskerCam.bin' extension
if nargin == 0
    fileNames = uigetfile('*_WhiskerCam.bin','MultiSelect','on');   % CTL-A to select all files
end
% Default setting - if you automatically play the function, it will track the whiskers
% To debug the code without taking the time to track whiskers, set trackWhiskers = 0
if nargin < 2
    trackWhiskers = 1;
end
% Prompt user if the laser doppler was aquired for this day of imaging
ldInput = input('Was laser doppler acquired during this trial? (y/n): ','s'); disp(' ')
p2Input = input('Is this IOS imaging for 2P data? (y/n): ','s'); disp(' ')

%% BLOCK PURPOSE: [1] Preparing to create RawData files.
disp('Analyzing Block [1] Preparing to create RawData file(s).'); disp(' ')
% Load in each file one at a time, looping through the list
for a = 1:length(fileNames)
    disp(['Analyzing WhiskerCam file (' num2str(a) ' of ' num2str(length(fileNames)) ')']); disp(' ')
    % Adapt to list or single file. The purpose of this is control the way uigetfile handles an instance of a
    % single file input (character string) vs. multiple files, which it puts in cells
    if iscell(fileNames) == true
        indFile = fileNames{a};
    else
        indFile = fileNames;
    end
    % Pull out the file ID for the file - this is the numerical string after the animal name/hemisphere
    [~,~,fileID] = GetFileInfo_IOS_eLife2020(indFile);
    % Determine if a RawData file has already been created for this file. If it has, skip it
    fileExist = ls(['*' fileID '_RawData.mat']);
    if isempty(fileExist)
        %% BLOCK PURPOSE: [2] Import .tdms data (All channels).
        disp('Analyzing Block [2] Importing .tdms data from all channels.'); disp(' ')
        trialData = ReadInTDMSWhiskerTrials_IOS_eLife2020([fileID '.tdms']);
        % Left, Right, and hippocampal electrodes
        dataRow = strcmp(trialData.data.names,'Cortical_LH');  
        cortical_LH = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
        dataRow = strcmp(trialData.data.names,'Cortical_RH');
        cortical_RH = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
        dataRow = strcmp(trialData.data.names,'Hippocampus');
        hippocampus = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
        % Left, Right, Auditory solenoids. Combine the arrays together.
        dataRow = strcmp(trialData.data.names,'LPadSol'); 
        LPadSol = gt(trialData.data.vals(dataRow,:),0.5)*1;   % ID amplitude is 1 
        dataRow = strcmp(trialData.data.names,'RPadSol'); 
        RPadSol = gt(trialData.data.vals(dataRow,:),0.5)*2;   % ID amplitude is 2
        dataRow = strcmp(trialData.data.names,'AudSol'); 
        AudSol = gt(trialData.data.vals(dataRow,:),0.5)*3;   % ID amplitude is 3
        solenoids = LPadSol + RPadSol + AudSol;
        % Force sensor and EMG
        dataRow = strcmp(trialData.data.names,'Force_Sensor'); 
        forceSensor = trialData.data.vals(dataRow,:);
        dataRow = strcmp(trialData.data.names,'EMG');
        EMG = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
        % Laser doppler
        if strcmp(ldInput,'y') == true
            trialData2 = ReadInTDMSWhiskerTrials_LD_IOS_eLife2020([fileID '_LD.tdms']);
            % LD backscatter
            dataRow = strcmp(trialData2.data.names,'LD_BackScatter');
            backScatter = trialData2.data.vals(dataRow,:);          
            % LD Flow
            dataRow = strcmp(trialData2.data.names,'LD_Flow');
            flow = trialData2.data.vals(dataRow,:);
        end
        
        %% BLOCK PURPOSE: [3] Start Whisker tracker.
        disp('Analyzing Block [3] Starting whisker tracking.'); disp(' ')
        if trackWhiskers == true
            [whiskerAngle] = WhiskerTrackerParallel_IOS_eLife2020(fileID);
            inds = isnan(whiskerAngle) == 1;
            whiskerAngle(inds) = [];
        else
            whiskerAngle = [];
        end
        
        %% BLOCK PURPOSE: [4] Save the notes and data.
        disp('Analyzing Block [4] Evaluating data to save to RawData file.'); disp(' ')
        % notes - all variables are descriptive
        RawData.notes.experimenter = trialData.experimenter;
        RawData.notes.animalID = trialData.animalID;
        RawData.notes.hemisphere = trialData.hemisphere;
        RawData.notes.solenoidPSI = str2double(trialData.solenoidPSI);
        RawData.notes.isofluraneTime = str2double(trialData.isofluraneTime);
        RawData.notes.sessionID = trialData.sessionID;
        RawData.notes.amplifierGain = str2double(trialData.amplifierGain);
        RawData.notes.CBVCamSamplingRate = str2double(trialData.CBVCamSamplingRate);
        RawData.notes.whiskCamSamplingRate = str2double(trialData.whiskCamSamplingRate);
        RawData.notes.webCamSamplingRate = str2double(trialData.webCamSamplingRate);
        RawData.notes.pupilCamSamplingRate = str2double(trialData.pupilCamSamplingRate);
        RawData.notes.analogSamplingRate = str2double(trialData.analogSamplingRate);
        RawData.notes.trialDuration_sec = str2double(trialData.trialDuration_sec);
        RawData.notes.CBVCamPixelWidth = str2double(trialData.CBVCamPixelWidth);
        RawData.notes.CBVCamPixelHeight = str2double(trialData.CBVCamPixelHeight);
        RawData.notes.CBVCamBitDepth = str2double(trialData.CBVCamBitDepth);
        RawData.notes.pupilCamPixelWidth = str2double(trialData.pupilCamPixelWidth);
        RawData.notes.pupilCamPixelHeight = str2double(trialData.pupilCamPixelHeight);
        RawData.notes.whiskCamPixelHeight = str2double(trialData.whiskCamPixelHeight);
        RawData.notes.whiskCamPixelWidth = str2double(trialData.whiskCamPixelWidth);
        RawData.notes.CBVCamExposureTime_microsec = str2double(trialData.CBVCamExposureTime_microsec);
        RawData.notes.CBVCamBinning = trialData.CBVCamBinning;
        RawData.notes.droppedPupilCamFrameIndex = trialData.droppedPupilCamFrameIndex;
        RawData.notes.droppedWhiskCamFrameIndex = trialData.droppedWhiskCamFrameIndex;
        % Data
        RawData.data.cortical_LH = cortical_LH;
        RawData.data.cortical_RH = cortical_RH;
        RawData.data.hippocampus = hippocampus;
        RawData.data.forceSensor = forceSensor;
        RawData.data.EMG = EMG;
        RawData.data.whiskerAngle = whiskerAngle;
        RawData.data.solenoids = solenoids;
        if strcmp(ldInput,'y') == true
            RawData.data.backScatter = backScatter;
            RawData.data.flow = flow;
        end
        disp(['File Created. Saving RawData File ' num2str(a) '...']); disp(' ')
        save([trialData.animalID '_' fileID '_RawData'],'RawData')
    else
        disp('File already exists. Continuing...'); disp(' ')
    end
end
% add ROI pixels to RawData file if desired - for 2P data only
if strcmp(p2Input,'y') == true
    for a = 1:length(fileNames)
        disp(['Analyzing WindowCam file (' num2str(a) ' of ' num2str(length(fileNames)) ')']); disp(' ')
        % Adapt to list or single file. The purpose of this is control the way uigetfile handles an instance of a
        % single file input (character string) vs. multiple files, which it puts in cells
        if iscell(fileNames) == true
            indFile = fileNames{a};
        else
            indFile = fileNames;
        end
        % Pull out the file ID for the file - this is the numerical string after the animal name/hemisphere
        [~,~,fileID] = GetFileInfo_IOS_eLife2020(indFile);
        ExtractImageMatrixFor2PData_IOS_eLife2020(fileID);
    end
end
disp('IOS Stage One Processing - Complete.'); disp(' ')

end
