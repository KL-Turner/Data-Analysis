%_______________________________________________________________________________________________
% Written by Kevin L. Turner, Jr. 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%_______________________________________________________________________________________________
%
%   Purpose: 
%_______________________________________________________________________________________________
%
%   Inputs: 
%
%   Outputs: 
%________________________________________________________________________________________________

%% Clear workspace, command window, close all figures
clear;
clc;

%% Load in a single 5 minute dalsa file
filename = uigetfile('*_WindowCam.bin', 'Pick a single .dalsa CBV file for the movie');
image_height = 256;     % Default image value                                                                                                                
image_width = 256;      % Default image value
Fs = 30;                % Sampling rate: note this should be 60 for GCaMP experiments

Start_Time = input('Input the start time in seconds (0 - 300): '); disp(' ')  % Integers only
End_Time = input('Input the end time in seconds: '); disp(' ')                % Integers only

Frame_Start_Time = floor(Start_Time)*Fs;      % Convert to frame # and round down to nearest integer
Frame_End_Time = floor(End_Time)*Fs;         
FrameInds = Frame_Start_Time:Frame_End_Time;  % Range of implay

CBV_Frames = GetCBVFrameSubset_IOS(filename, image_height, image_width, FrameInds);  % Obtain subset of desired frames
Baseline_Frame = mean(CBV_Frames, 3);   % Average in the 3rd dimension (time) to create a baseline

%% Create a figure to show the baseline frame in color and grey-scale
% figure; 
% imagesc(Baseline_Frame);    % Show scaled baseline frame
% title('Color-scaled baseline frame (max 4096)');
% xlabel('width (pixels)');
% ylabel('height (pixels)');
% colorbar;
% axis square;
% 
% figure;
% imagesc(Baseline_Frame);
% colormap('gray');
% title('Grey-scale baseline frame');
% xlabel('width (pixels)');
% ylabel('height (pixels)');
% axis square;

%% Create implay movie for the desired timeframe
N_cbv_frames = zeros(size(CBV_Frames));     % Pre-allocate space for normalized frames

for N = 1:size(CBV_Frames, 3)
    disp(['Adding frame numer ' num2str(N) ' of ' num2str(size(CBV_Frames, 3)) '.']); disp(' ')
    N_cbv_frames(:,:, N) = CBV_Frames(:,:, N)./ Baseline_Frame;     % Normalize by baseline
    imagesc(CBV_Frames(:,:, N))
    imagesc(N_cbv_frames(:,:, N) )
end

min = 0.95;
max = 1.05;
ImplayAutoColorMap_IOS(N_cbv_frames, min, max);

Prompt = msgbox('Warning: May need to edit implay colormap range for best results');
