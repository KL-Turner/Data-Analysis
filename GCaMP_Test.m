clear; clc; close all;
% user inputs for file information
rawDataFileID = uigetfile('*_RawData.mat','MultiSelect','off');
load(rawDataFileID)
[~,~,fileID] = GetFileInfo_IOS(rawDataFileID);
windowCamFileID = [fileID '_WindowCam.bin'];
Fs = RawData.notes.CBVCamSamplingRate;
imageHeight = RawData.notes.CBVCamPixelHeight;                                                                                                            
imageWidth = RawData.notes.CBVCamPixelWidth;
trialDuration = RawData.notes.trialDuration_sec;
pixelsPerFrame = imageWidth*imageHeight;
% open the file, get file size, back to the begining
fid = fopen(windowCamFileID);
fseek(fid,0,'eof');
fileSize = ftell(fid);
fseek(fid,0,'bof');
% identify the number of frames to read
nFramesToRead = floor(fileSize/(pixelsPerFrame*2));
skippedPixels = pixelsPerFrame*2; 
% loop over each frame
imageStack = NaN*ones(imageHeight,imageWidth,nFramesToRead);
for n = 1:nFramesToRead
    fseek(fid,(n - 1)*skippedPixels,'bof');
    z = fread(fid,pixelsPerFrame,'*int16','b');
    img = reshape(z,imageHeight,imageWidth);
    imageStack(:,:,n) = rot90(img',2);
end
fclose('all');
%% Create implay movie for the desired timeframe
handle = implay(imageStack,Fs);
handle.Visual.ColorMap.UserRange = 1; 
handle.Visual.ColorMap.UserRangeMin = 0; 
handle.Visual.ColorMap.UserRangeMax = 2000; % max is 4096
%% Draw ROI over L/R hemispheres
