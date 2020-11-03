%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate supplemental videos for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_eLife2020
%________________________________________________________________________________________________________________________

clear; clc; close all;
%% information and data for first example
% dataLocation = [rootFolder '\Summary Figures and Structures\Supplemental Movies\'];
% cd(dataLocation)
fs1 = 30;   % cbv/pupil camera
fs2 = 150;  % whisker camera
colorRfcAwake = [(0/256),(64/256),(64/256)];
colorRfcNREM = [(0/256),(174/256),(239/256)];
colorRfcREM = [(190/256),(30/256),(45/256)];
%% Take frames from the CBV camera file for baseline/awake/rem examples
exampleRawDataFileID_A = 'T123_200301_14_48_14_RawData.mat';
load(exampleRawDataFileID_A,'-mat')
exampleCBVcamFileID_A = '200301_14_48_14_WindowCam.bin';
baseStartTime = 580;
baseEndTime = 600;
baseFrameIndex = baseStartTime*fs1:baseEndTime*fs1;
awakeStartTime = 400;
awakeEndTime = 600;
awakeFrameIndex = awakeStartTime*fs1:awakeEndTime*fs1;
awakeFrameWhiskIndex = awakeStartTime*fs2:awakeEndTime*fs2;
nremStartTime = 40;
nremEndTime = 100;
nremFrameIndex = nremStartTime*fs1:nremEndTime*fs1;
nremFrameWhiskIndex = nremStartTime*fs2:nremEndTime*fs2;
remStartTime = 100;
remEndTime = 300;
remFrameIndex = remStartTime*fs1:remEndTime*fs1;
remFrameWhiskIndex = remStartTime*fs2:remEndTime*fs2;
% Obtain subset of desired frames
cbvImageHeight = RawData.notes.CBVCamPixelHeight;
cbvImageWidth = RawData.notes.CBVCamPixelWidth;
baseCBVframes = GetCBVFrameSubset_IOS_eLife2020(exampleCBVcamFileID_A,cbvImageHeight,cbvImageWidth,baseFrameIndex);
awakeCBVframes = GetCBVFrameSubset_IOS_eLife2020(exampleCBVcamFileID_A,cbvImageHeight,cbvImageWidth,awakeFrameIndex);
nremCBVframes = GetCBVFrameSubset_IOS_eLife2020(exampleCBVcamFileID_A,cbvImageHeight,cbvImageWidth,nremFrameIndex);
remCBVframes = GetCBVFrameSubset_IOS_eLife2020(exampleCBVcamFileID_A,cbvImageHeight,cbvImageWidth,remFrameIndex);
%% Use baseline frames to establish ROIs and resting baseline
% draw a rectangular ROI around the window to remove outside pixels
boxFig = figure;
imagesc(baseCBVframes(:,:,1))
colormap gray
caxis([0,2^12])
axis image
boxROI = drawrectangle;
boxPosition = round(boxROI.Vertices);
boxX = boxPosition(:,1);
boxY = boxPosition(:,2);
boxMask = poly2mask(boxX,boxY,cbvImageWidth,cbvImageHeight);
close(boxFig)
% take values from within a square ROI within the window
boxWidth = abs(boxPosition(1,1) - boxPosition(3,1));
boxHeight = abs(boxPosition(1,2) - boxPosition(2,2));
% baseline frames
for a = 1:size(baseCBVframes,3)
    baseCBVframe = baseCBVframes(:,:,a);
    baseBoxVals = baseCBVframe(boxMask);
    baseBoxFrames(:,:,a) = reshape(baseBoxVals,boxHeight,boxWidth); %#ok<*SAGROW>
end
% awake frames
for a = 1:size(awakeCBVframes,3)
    awakeCBVframe = awakeCBVframes(:,:,a);
    awakeBoxVals = awakeCBVframe(boxMask);
    awakeBoxFrames(:,:,a) = reshape(awakeBoxVals,boxHeight,boxWidth);
end
% nrem frames
for a = 1:size(nremCBVframes,3)
    nremCBVframe = nremCBVframes(:,:,a);
    nremBoxVals = nremCBVframe(boxMask);
    nremBoxFrames(:,:,a) = reshape(nremBoxVals,boxHeight,boxWidth);
end
% rem frames
for a = 1:size(remCBVframes,3)
    remCBVframe = remCBVframes(:,:,a);
    remBoxVals = remCBVframe(boxMask);
    remBoxFrames(:,:,a) = reshape(remBoxVals,boxHeight,boxWidth);
end
% set values from outside the window to NaN
windowFig = figure;
imagesc(baseBoxFrames(:,:,1))
colormap gray
caxis([0,2^12])
axis image
windowMask = roipoly;
close(windowFig)
% baseline frames
for w = 1:size(baseBoxFrames,3)
    baseWindowFrame = baseBoxFrames(:,:,w);
    baseWindowFrame(~windowMask) = NaN;
    baseWindowFrames(:,:,w) = baseWindowFrame;
end
% awake frames
for w = 1:size(awakeBoxFrames,3)
    awakeWindowFrame = awakeBoxFrames(:,:,w);
    awakeWindowFrame(~windowMask) = NaN;
    awakeWindowFrames(:,:,w) = awakeWindowFrame;
end
% nrem frames
for w = 1:size(nremBoxFrames,3)
    nremWindowFrame = nremBoxFrames(:,:,w);
    nremWindowFrame(~windowMask) = NaN;
    nremWindowFrames(:,:,w) = nremWindowFrame;
end
% rem frames
for w = 1:size(remBoxFrames,3)
    remWindowFrame = remBoxFrames(:,:,w);
    remWindowFrame(~windowMask) = NaN;
    remWindowFrames(:,:,w) = remWindowFrame;
end
%% Normalize each data set by the baseline frame
baselineFrame = mean(baseWindowFrames,3);
% awake data
for a = 1:size(awakeWindowFrames,3)
    awakeCBVImageStack(:,:,a) = ((awakeWindowFrames(:,:,a)));% - baselineFrame));%./(baselineFrame)).*100;
end
% nrem data
for a = 1:size(nremWindowFrames,3)
    nremCBVImageStack(:,:,a) = ((nremWindowFrames(:,:,a)));% - baselineFrame));%./(baselineFrame)).*100;
end
% rem data
for a = 1:size(remWindowFrames,3)
    remCBVImageStack(:,:,a) = ((remWindowFrames(:,:,a)));% - baselineFrame));%./(baselineFrame)).*100;
end
%% figure
figure;
sgtitle('IOS camera')
subplot(2,2,1)
imagesc(baseBoxFrames(:,:,1))
title('base image')
colormap(gca,'gray')
c0 = colorbar;
ylabel(c0,'Pixel intensity','rotation',-90,'VerticalAlignment','bottom')
caxis([0,2^12])
axis image
subplot(2,2,2)
imagesc(awakeCBVImageStack(:,:,1))
title('Awake reflectance')
colormap(gca,'gray')
c0 = colorbar;
ylabel(c0,'\DeltaPixel intensity','rotation',-90,'VerticalAlignment','bottom')
caxis([0,2^12])
axis image
subplot(2,2,3)
imagesc(nremCBVImageStack(:,:,662))
title('NREM reflectance')
colormap(gca,'gray')
c0 = colorbar;
ylabel(c0,'\DeltaPixel intensity','rotation',-90,'VerticalAlignment','bottom')
caxis([0,2^12])
axis image
subplot(2,2,4)
imagesc(remCBVImageStack(:,:,1175))
title('REM reflectance')
colormap(gca,'gray')
c0 = colorbar;
ylabel(c0,'\DeltaPixel intensity','rotation',-90,'VerticalAlignment','bottom')
caxis([0,2^12])
axis image
        
% handle = implay(remCBVImageStack,30);
% handle.Visual.ColorMap.UserRange = 1; 
% handle.Visual.ColorMap.UserRangeMin = -10; 
% handle.Visual.ColorMap.UserRangeMax = 10;
