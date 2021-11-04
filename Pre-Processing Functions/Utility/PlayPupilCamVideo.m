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
%   Last Revised:
%________________________________________________________________________________________________________________________

zap;
% User inputs for file information
procDataFileID = uigetfile('*_ProcData.mat','MultiSelect','off');
load(procDataFileID)
[~,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
pupilCamFileID = [fileID '_PupilCam.bin'];
trialDuration = ProcData.notes.trialDuration_sec;
fid = fopen(pupilCamFileID); % reads the binary file in to the work space
fseek(fid,0,'eof'); % find the end of the video frame
fileSize = ftell(fid); % calculate file size
fseek(fid,0,'bof'); % find the begining of video frames
imageHeight = ProcData.notes.pupilCamPixelHeight; % how many pixels tall is the frame
imageWidth = ProcData.notes.pupilCamPixelWidth; % how many pixels wide is the frame
samplingRate = ProcData.notes.pupilCamSamplingRate;
pixelsPerFrame = imageWidth*imageHeight;
skippedPixels = pixelsPerFrame;
nFramesToRead = floor(fileSize/(pixelsPerFrame));
imageStack = zeros(imageHeight,imageWidth,nFramesToRead);
for dd = 1:nFramesToRead
    fseek(fid,(dd - 1)*skippedPixels,'bof');
    z = fread(fid,pixelsPerFrame,'*uint8','b');
    img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
    imageStack(:,:,dd) = flip(imrotate(img,-90),2);
end
fclose('all');
sliceViewer(imageStack)
