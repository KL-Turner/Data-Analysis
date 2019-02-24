function [] = CreateBilateralROI(animal, hem, rawDataFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs: 
%
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

%% Draw ROIs
ROINames = {'LH', 'RH', 'LH_Electrode', 'RH_Electrode'};
DrawROIs(animal, hem, ROINames);

%% ROI Trace
disp('Loading in files for the Left Hemisphere...'); disp(' ')
ROIname1 = 'LH';
AddROIIntensityToRawdataFile(ROIname1, rawDataFiles)   % Add CBV data to RawData file after drawing the ROI

disp('Loading in files for the Right Hemisphere...'); disp(' ')
ROIname2 = 'RH';
AddROIIntensityToRawdataFile(ROIname2, rawDataFiles)   % Add CBV data to RawData file after drawing the ROI

disp('Loading in files for the Left Hemisphere...'); disp(' ')
ROIname3 = 'LH_Electrode';
AddROIIntensityToRawdataFile(ROIname3, rawDataFiles)   % Add CBV data to RawData file after drawing the ROI

disp('Loading in files for the Right Hemisphere...'); disp(' ')
ROIname4 = 'RH_Electrode';
AddROIIntensityToRawdataFile(ROIname4, rawDataFiles)   % Add CBV data to RawData file after drawing the ROI

close all;
disp('All ROIs added to RawData files - Complete'); disp(' ')
disp('To create a new ROI, delete the current ROIs.mat file'); disp(' ')

end
