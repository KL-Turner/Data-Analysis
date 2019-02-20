%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: (1) Verify the whisker stimulus pattern by plotting each whisker stimulus for each day
%            (2) Create a plot of the windows for each day along with the associated ROI that was drawn
%            (3) Check for movement artifacts in the first X seconds of the first file of each day
%            (4) Check the quality of the electrodes
%________________________________________________________________________________________________________________________
%
%   Inputs:            
%
%   Outputs:        
%
%   Last Revised: October 5th, 2018
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
clc;
clear;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')

[animal, ~, ~, ~, ~, ~, ~, ~, ~, ~] = LoadDataStructs();

windowCamFiles = ls('*_WindowCam.bin');
[~, ~, fileDates, ~] = GetFileInfo(windowCamFiles);
[uniqueDays, ~, DayID] = GetUniqueDays(fileDates);
firstFileOfDay = cell(1, length(uniqueDays));
for uD = 1:length(uniqueDays)
    FileInd = DayID == uD;
    dayFilenames = windowCamFiles(FileInd,:);
    firstFileOfDay(uD) = {dayFilenames(1,:)};
end
procDataFiles = ls('*_ProcData.mat');

disp('Block [0] structs loaded.'); disp(' ')

%% BLOCK PURPOSE: [1] Compare windows from the various days.
% Obtain the first file name for each unique day of imaging
disp('Analyzing Block [1] Compare windows from different days.'); disp(' ')
for fFOD = 1:length(firstFileOfDay)
    firstFile = firstFileOfDay{fFOD};
    [~, ~, fileDate, ~] = GetFileInfo(firstFile);
    strDay = ConvertDate(fileDate); 
    imageWidth = 256;
    imageHeight = 256;
    GetSingleCBVFrame(firstFile, imageWidth, imageHeight, animal, strDay);
    UniqueDayROIs(firstFile, imageWidth, imageHeight, animal, strDay);
end

%% BLOCK PURPOSE: [2] Check for Movement Artifacts.
% Obtain the first file name for each unique day of imaging
disp('Analyzing Block [2] Check for movement artifacts from each day.'); disp(' ')
for fFOD = 1:length(firstFileOfDay)
    firstFile = firstFileOfDay{fFOD};
    [~, ~, fileDate, ~] = GetFileInfo(firstFile);
    strDay = ConvertDate(fileDate); 
    imageHeight = 256;
    imageWidth = 256;
    frameInds = 1:900;
    [imageStack] = GetCBVFrameSubset(firstFile, imageHeight, imageWidth, frameInds);
    CheckMovementArtifacts(imageStack);
    prompt = msgbox('Click when finished checking for movement artifacts');
    waitfor(prompt);
end

% implay_SingleTrial

%% BLOCK PURPOSE: [3] Check the left and right electrode quality
disp('Analyzing Block [3] Check the electrode quality.'); disp(' ')
for fFOD = 1:length(firstFileOfDay)
    firstFile = firstFileOfDay{fFOD}(1:15);
    rawDataFile = ([animal '_both_' firstFile '_rawdata.mat']);
    procDataFile = ([animal '_both_' firstFile '_ProcData.mat']);
    CheckElectrodeQuality(rawDataFile, procDataFile);
end

close all;
disp('Stage Four Processing - Complete.'); disp(' ') 
