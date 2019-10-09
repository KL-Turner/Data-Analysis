%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpse:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised:
%________________________________________________________________________________________________________________________

clear
clc

% User inputs for file information
rawDataFileID = uigetfile('*_RawData.mat','MultiSelect','off');
[animalID,~,fileID] = GetFileInfo_IOS(rawDataFileID);
windowCamFileID = [animalID '_' fileID '_WindowCam.bin'];
pupilCamFileID = [animalID '_' fileID '_PupilCam.bin'];
whiskCamFileID = [animalID '_' fileID '_WhiskerCam.bin'];
procDataFileID = [animalID '_' fileID '_ProcData.mat'];
specDataFileID = [animalID '_' fileID '_SpecData.mat'];
baselineFileID = [animalID '_RestingBaselines.mat'];

% load all relevent file structures
load(rawDataFileID)
load(procDataFileID)
load(specDataFileID)
load(baselineFileID)

% establish each video file's information
trialDuration = RawData.notes.trialDuration_sec;
disp(['Trial is ' num2str(trialDuration) ' seconds long.']); disp(' ')
startTime = input('Input the desired start time (sec): '); disp(' ')
endTime = input('Input the desired end time (sec): '); disp(' ')

if startTime >= trialDuration || startTime < 0
    disp(['A start time of  ' num2str(startTime) ' is not a valid input']); disp(' ')
    return
elseif endTime > trialDuration || endTime <= startTime || endTime <= 0
    disp(['An end time of  ' num2str(startTime) ' is not a valid input']); disp(' ')
    return
end

%% CBV reflectance movie
imageHeight = RawData.notes.CBVCamPixelHeight;                                                                                                            
imageWidth = RawData.notes.CBVCamPixelWidth;
Fs = RawData.notes.CBVCamSamplingRate;

frameStart = floor(startTime)*Fs;
frameEnd = floor(endTime)*Fs;         
frameInds = frameStart:frameEnd;

 % Obtain subset of desired frames - normalize by an artificial baseline
frames = GetCBVFrameSubset_IOS(windowCamFileID,imageHeight,imageWidth,frameInds);
baselineFrame = mean(frames,3);
normFrames = zeros(size(frames));
for a = 1:size(frames,3)
    disp(['Creating CBV image stack: (' num2str(a) '/' num2str(size(frames,3)) ')']); disp(' ')
    normFrames(:,:,a) = frames(:,:,a)./baselineFrame;   % Normalize by baseline
end

%% Whisker 
imageHeight = RawData.notes.whiskCamPixelHeight;                                                                                                            
imageWidth = RawData.notes.whiskCamPixelWidth;
Fs = RawData.notes.whiskCamSamplingRate;

frameStart = floor(startTime)*Fs;
frameEnd = floor(endTime)*Fs;         
frameInds = frameStart:frameEnd;

pixelsPerFrame = imageWidth*imageHeight;
skippedPixels = pixelsPerFrame*2; % Multiply by two because there are 16 bits (2 bytes) per pixel
fid = fopen(whiskCamFileID);
fseek(fid,0,'eof');
fileSize = ftell(fid);
fseek(fid,0,'bof');
nFramesToRead = length(frameInds);
imageStack = zeros(imageHeight,imageWidth,nFramesToRead);
for a = 1:nFramesToRead
    disp(['Creating image stack: (' num2str(a) '/' num2str(nFramesToRead) ')']); disp(' ')
    fseek(fid,frameInds(a)*skippedPixels,'bof');
    z = fread(fid,pixelsPerFrame,'*uint8','b');
    img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
    imageStack(:,:,a) = flip(imrotate(img,-90),2);
end
fclose('all');

handle = implay(imageStack,Fs);
handle.Visual.ColorMap.UserRange = 1; 
handle.Visual.ColorMap.UserRangeMin = min(img(:)); 
handle.Visual.ColorMap.UserRangeMax = max(img(:));

%% Pupil
imageHeight = RawData.notes.pupilCamPixelHeight;                                                                                                            
imageWidth = RawData.notes.pupilCamPixelWidth;
Fs = RawData.notes.pupilCamSamplingRate;

frameStart = floor(startTime)*Fs;
frameEnd = floor(endTime)*Fs;         
frameInds = frameStart:frameEnd;

pixelsPerFrame = imageWidth*imageHeight;
skippedPixels = pixelsPerFrame*2;   % Multiply by two because there are 16 bits (2 bytes) per pixel
fid = fopen(pupilCamFileID);
fseek(fid,0,'eof');
fileSize = ftell(fid);
fseek(fid,0,'bof');
nFramesToRead = length(frameInds);
imageStack = zeros(imageWidth,imageHeight,nFramesToRead);
for a = 1:nFramesToRead
    disp(['Creating image stack: (' num2str(a) '/' num2str(nFramesToRead) ')']); disp(' ')
    fseek(fid, frameInds(a)*skippedPixels,'bof');
    z = fread(fid, pixelsPerFrame,'*uint8','b');
    img = reshape(z(1:pixelsPerFrame),imageHeight,imageWidth);
    imageStack(:,:,a) = flip(imrotate(img,-90),2);
end
fclose('all');

handle = implay(imageStack, Fs);
handle.Visual.ColorMap.UserRange = 1; 
handle.Visual.ColorMap.UserRangeMin = min(img(:)); 
handle.Visual.ColorMap.UserRangeMax = max(img(:));




















































%% movie file comparing rgb with original data
outputVideo = VideoWriter('PresentationMotionVideo_RGBvsOriginalDepth.avi');
fps = 15;   % default fps from video acquisition
speedUp = 2;   % speed up by factor of
outputVideo.FrameRate = fps*speedUp;
open(outputVideo);
fig = figure('Position', get(0, 'Screensize'));
for a = 1:size(rgbStack,1)
    subplot(1,2,1)
    imshow(rgbStack{a,1})
    subplot(1,2,2)
    imagesc(originalStack{a,1});
    colormap jet
    caxis([0 .52])
    axis image
    axis off
    currentFrame = getframe(fig);
    writeVideo(outputVideo, currentFrame);
end
close(outputVideo)
close(fig)

%% movie file comparing rgb with processed data
outputVideo = VideoWriter('PresentationMotionVideo_RGBvsProcDepth.avi');
fps = 15;   % default fps from video acquisition
speedUp = 2;   % speed up by factor of
outputVideo.FrameRate = fps*speedUp;
open(outputVideo);
fig = figure('Position', get(0, 'Screensize'));
for a = 1:size(rgbStack,1)
    subplot(1,2,1)
    imshow(rgbStack{a,1})
    subplot(1,2,2)
    imagesc(procStack(:,:,a));
    colormap jet
    caxis(SuppData.caxis)
    axis image
    axis off
    currentFrame = getframe(fig);
    writeVideo(outputVideo, currentFrame);
end
close(outputVideo)
close(fig)

%% movie file showing motion and height tracking
outputVideo = VideoWriter('PresentationMotionVideo_Results.avi');
fps = 15;   % default fps from video acquisition
speedUp = 2;   % speed up by factor of
outputVideo.FrameRate = fps*speedUp;
open(outputVideo);
fig = figure('Position', get(0, 'Screensize'));
avg20Height = NaN(1,length(binStack));
max_caxis = SuppData.caxis;
maxVal = max_caxis(2);
x = [];
y = [];
distanceTraveled = 0;
distancePath = NaN(1,length(binStack));
binWidth_inches = 14;
distancePerPixel = (binWidth_inches/SuppData.binWidth)*2.54;   % in to cm
for a = 1:size(binStack,3)
    %% Motion
    imageA = binStack(:,:,a);
    [yA,xA] = ndgrid(1:size(imageA,1), 1:size(imageA,2));
    centroidA = mean([xA(logical(imageA)), yA(logical(imageA))]);
    x = horzcat(x,centroidA(1));
    y = horzcat(y,centroidA(2));
    
    
    if a > 1
        imageB = binStack(:,:,a-1);
        [yB,xB] = ndgrid(1:size(imageB,1), 1:size(imageB,2));
        centroidB = mean([xB(logical(imageB)), yB(logical(imageB))]);
        
        centroidCoord = [centroidB; centroidA];
        d = pdist(centroidCoord, 'euclidean');
        if isnan(d) == true
            d = 0;
        end
            distanceTraveled = distanceTraveled+d;
    end
    distancePath(1,a) = distanceTraveled;
           
    %% Rearing
    depthImg = procStack(:,:,a);
    maxInds = depthImg == maxVal;
    depthImg(maxInds) = NaN;
    validPix = imcomplement(isnan(depthImg));
    pixelVec = depthImg(validPix);
    ascendPixelVals = sort(pixelVec(:),'ascend');
    twentyPercentile = ascendPixelVals(1:ceil(length(ascendPixelVals)*0.2));
    avg20Height(1,a) = mean(twentyPercentile);
   
    %% figure
    subplot(2,2,[1,3])
    imagesc(procStack(:,:,a));
    colormap jet
    caxis(SuppData.caxis)
    axis image
    axis off
    hold on
    scatter(centroidA(1), centroidA(2), 'MarkerEdgeColor', 'white', 'MarkerFaceColor', 'white')
    plot(x,y,'Color','w','LineWidth',2)
    
    subplot(2,2,2)
    plot((1:size(binStack,3))/fps, 100*avg20Height, 'k')
    set(gca, 'YDir','reverse')
    title('Distance from camera')
    ylabel('Height (cm)')
    xlabel('Time (sec)')
    xlim([0 70])
    ylim([40 48])
    
    subplot(2,2,4)
    plot((1:size(binStack,3))/fps, distancePath.*distancePerPixel, 'k')
    title('Distance traveled')
    ylabel('Distance (cm)')
    xlabel('Time (sec)')
    xlim([0 70])

    currentFrame = getframe(fig);
    writeVideo(outputVideo, currentFrame);
end
close(outputVideo)
close(fig)
