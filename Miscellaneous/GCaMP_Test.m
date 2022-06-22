clear; clc;
% user inputs for file information
windowCamFileID = uigetfile('*_WindowCam.bin','MultiSelect','off');
Fs = 30;
imageHeight = 128;                                                                                                            
imageWidth = 128;
trialDuration = 120;
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
% fclose('all');
%% Create implay movie for the desired timeframe
% handle = implay(imageStack,Fs);
% handle.Visual.ColorMap.UserRange = 1; 
% handle.Visual.ColorMap.UserRangeMin = 0; 
% handle.Visual.ColorMap.UserRangeMax = 2000; % max is 4096
%% Draw ROI over L/R hemispheres
% create/load pre-existing ROI file with the coordinates
% ROIFileDir = dir('*_ROIs.mat');
% if isempty(ROIFileDir) == true
%     ROIs = [];
% else
%     ROIFileName = {ROIFileDir.name}';
%     ROIFileID = char(ROIFileName);
%     load(ROIFileID);
% end
%% Draw ROIs
% create figure of the image frame
roiFig = figure;
imagesc(imageStack(:,:,3))
colormap(gray)
axis image
xlabel('Caudal')
caxis([0,2000])
% draw ROI over the cement
disp('Please select your region of interest for the left hemisphere'); disp(' ')
[~,LH_xi,LH_yi] = roipoly;
ROIs.LH_xi = LH_xi;
ROIs.LH_yi = LH_yi;
disp('Please select your region of interest for the right hemisphere'); disp(' ')
[~,RH_xi,RH_yi] = roipoly;
ROIs.RH_xi = RH_xi;
ROIs.RH_yi = RH_yi;
close(roiFig)
%% extract CBV/GCaMP frames
cbvFrames = imageStack(:,:,3:2:end - 1);
gcampFrames = imageStack(:,:,4:2:end);
numFrames = size(cbvFrames,3);
% figure; 
% imagesc(cbvFrames(:,:,3)); 
% caxis([0,2000]); 
% colormap gray; 
% axis image
% LH - apply image mask to each frame of the stack
LH_cbvRefl = zeros(1,numFrames);
LH_gcampRefl = zeros(1,numFrames);
for n = 1:numFrames
    LH_mask = roipoly(cbvFrames(:,:,1),ROIs.LH_xi,ROIs.LH_yi);
    LH_cbvMask = LH_mask.*double(cbvFrames(:,:,n));
    LH_gcampMask = LH_mask.*double(gcampFrames(:,:,n));
    LH_cbvRefl(n) = mean(nonzeros(LH_cbvMask));
    LH_gcampRefl(n) = mean(nonzeros(LH_gcampMask));
end
LH_cbvRefl_norm = (LH_cbvRefl - mean(LH_cbvRefl))./mean(LH_cbvRefl);
LH_gcampRefl_norm = (LH_gcampRefl - mean(LH_gcampRefl))./mean(LH_gcampRefl);
% RH - apply image mask to each frame of the stack
RH_cbvRefl = zeros(1,numFrames);
RH_gcampRefl = zeros(1,numFrames);
for n = 1:numFrames
    RH_mask = roipoly(cbvFrames(:,:,1),ROIs.RH_xi,ROIs.RH_yi);
    RH_cbvMask = RH_mask.*double(cbvFrames(:,:,n));
    RH_gcampMask = RH_mask.*double(gcampFrames(:,:,n));
    RH_cbvRefl(n) = mean(nonzeros(RH_cbvMask));
    RH_gcampRefl(n) = mean(nonzeros(RH_gcampMask));
end
RH_cbvRefl_norm = (RH_cbvRefl - mean(RH_cbvRefl))./mean(RH_cbvRefl);
RH_gcampRefl_norm = (RH_gcampRefl - mean(RH_gcampRefl))./mean(RH_gcampRefl);
%% generate figure
figure; 
plot(detrend(LH_cbvRefl_norm,'constant'))
hold on
plot(detrend(RH_cbvRefl_norm,'constant'))
plot(detrend(LH_gcampRefl_norm,'constant'))
plot(detrend(RH_gcampRefl_norm,'constant'))
legend('LH CBV','RH CBV','LH GCaMP7s','RH GCaMP7s')
ylabel('Percent (%)')
%% display correlation between cbv and gcamp
LH_R = corrcoef(LH_cbvRefl,LH_gcampRefl);
RH_R = corrcoef(LH_cbvRefl,LH_gcampRefl);
disp(['LH R: ' num2str(LH_R(2))]);
disp(['RH R: ' num2str(RH_R(2))]); disp(' ')
