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
% Clear the workspace variables and command window
clc;
clear;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')
% Asks the user to load all files with a '_WhiskerCam.bin' extension
fileNames = uigetfile('*_WhiskerCam.bin','MultiSelect','on');   % CTL-A to select all files
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
    [~,~,fileID] = GetFileInfo_IOS(indFile);
    % Determine if a RawData file has already been created for this file. If it has, skip it
    fileExist = ls(['*' fileID '_RawData.mat']);
    if isempty(fileExist)
        %% BLOCK PURPOSE: [2] Import .tdms data (All channels).
        disp('Analyzing Block [2] Importing .tdms data from all channels.'); disp(' ')
        trialData = ReadInTDMSWhiskerTrials_IOS([fileID '.tdms']);
        trialData2 = ReadInTDMSWhiskerTrials_LD_IOS([fileID '_LD.tdms']);
        % Left, Right, and hippocampal electrodes
        dataRow = strcmp(trialData.data.names,'Cortical_LH');  
        cortical_LH = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
        dataRow = strcmp(trialData.data.names,'Cortical_RH');
        cortical_RH = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
        dataRow = strcmp(trialData.data.names,'Hippocampus');
        hippocampus = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
        % Left, Right, Auditory solenoids. Combine the arrays together.
        dataRow = strcmp(trialData.data.names,'LPadSol'); 
        LPadSol = gt(trialData.data.vals(dataRow,:),0.5)*1;    % ID amplitude is 1 
        dataRow = strcmp(trialData.data.names,'RPadSol'); 
        RPadSol = gt(trialData.data.vals(dataRow,:),0.5)*2;    % ID amplitude is 2
        dataRow = strcmp(trialData.data.names,'AudSol'); 
        AudSol = gt(trialData.data.vals(dataRow,:),0.5)*3;     % ID amplitude is 3
        dataRow = strcmp(trialData2.data.names,'Opto_LED');
        OptoLED = gt(trialData2.data.vals(dataRow,:),0.5)*4;   % ID amplitude is 4
        analogLengthA = length(LPadSol);
        analogLengthB = length(OptoLED);
        if analogLengthA > analogLengthB
            stimulations = LPadSol(1:analogLengthB) + RPadSol(1:analogLengthB) + AudSol(1:analogLengthB) + OptoLED;
        else
            stimulations = LPadSol + RPadSol + AudSol + OptoLED(1:analogLengthA);
        end
        % Force sensor and EMG
        dataRow = strcmp(trialData.data.names,'Force_Sensor'); 
        forceSensor = trialData.data.vals(dataRow,:);
        dataRow = strcmp(trialData.data.names,'EMG');
        EMG = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
        % Laser doppler
        dataRow = strcmp(trialData2.data.names,'LD_BackScatter');
        backScatter = trialData2.data.vals(dataRow,:);
        dataRow = strcmp(trialData2.data.names,'LD_Flow');
        flow = trialData2.data.vals(dataRow,:);        
        %% BLOCK PURPOSE: [3] Start Whisker tracker.
        disp('Analyzing Block [3] Starting whisker tracking.'); disp(' ')
        [whiskerAngle] = WhiskerTrackerParallel_IOS(fileID);
        inds = isnan(whiskerAngle) == 1;
        whiskerAngle(inds) = [];
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
        RawData.notes.LEDpower_mW = trialData.LEDpower_mW;
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
        RawData.notes.solenoidDutyCycle = str2double(trialData.Sol_DutyCycle);
        RawData.notes.solenoidFreq = str2double(trialData.Sol_Freq);
        RawData.notes.solenoidDuration_sec = str2double(trialData.Sol_Duration_sec);
        RawData.notes.LEDdutyCycle = str2double(trialData.LED_DutyCycle);
        RawData.notes.LEDfreq = str2double(trialData.LED_Freq);
        RawData.notes.LEDduration_sec = str2double(trialData.LED_Duration_sec);
        RawData.notes.interstim_sec = str2double(trialData.Interstim_sec);
        RawData.notes.stimOffset_sec = str2double(trialData.Stim_Offset_sec);
        % Data
        RawData.data.cortical_LH = cortical_LH;
        RawData.data.cortical_RH = cortical_RH;
        RawData.data.hippocampus = hippocampus;
        RawData.data.forceSensor = forceSensor;
        RawData.data.EMG = EMG;
        RawData.data.whiskerAngle = whiskerAngle;
        RawData.data.stimulations = stimulations;
        RawData.data.backScatter = backScatter;
        RawData.data.flow = flow;
        disp(['File Created. Saving RawData File ' num2str(a) '...']); disp(' ')
        save([trialData.animalID '_' fileID '_RawData'],'RawData')
    else
        disp('File already exists. Continuing...'); disp(' ')
    end
end
disp('IOS Stage One Processing - Complete.'); disp(' ')
