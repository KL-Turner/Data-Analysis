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
[animalID,fileDate,fileID] = GetFileInfo_IOS(rawDataFileID);
strDay = ConvertDate_IOS(fileDate);
windowCamFileID = [fileID '_WindowCam.bin'];
pupilCamFileID = [fileID '_PupilCam.bin'];
whiskCamFileID = [fileID '_WhiskerCam.bin'];
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
cbvImageHeight = RawData.notes.CBVCamPixelHeight;                                                                                                            
cbvImageWidth = RawData.notes.CBVCamPixelWidth;
cbvFs = RawData.notes.CBVCamSamplingRate;

cbvFrameStart = floor(startTime)*cbvFs;
cbvFrameEnd = floor(endTime)*cbvFs;         
cbvFrameInds = cbvFrameStart:cbvFrameEnd;

% Obtain subset of desired frames - normalize by an artificial baseline
cbvFrames = GetCBVFrameSubset_IOS(windowCamFileID,cbvImageHeight,cbvImageWidth,cbvFrameInds);
% draw a rectangular ROI around the window to remove outside pixels
boxFig = figure;
imagesc(cbvFrames(:,:,1))
colormap gray
caxis([0 2^12])
axis image
boxROI = drawrectangle;
boxPosition = round(boxROI.Vertices);
boxX = boxPosition(:,1);
boxY = boxPosition(:,2);
boxMask = poly2mask(boxX,boxY,cbvImageWidth,cbvImageHeight);
close(boxFig)

boxWidth = abs(boxPosition(1,1) - boxPosition(3,1));
boxHeight = abs(boxPosition(1,2) - boxPosition(2,2));

for a = 1:size(cbvFrames,3)
    cbvFrame = cbvFrames(:,:,a);
    boxVals = cbvFrame(boxMask);
    boxCBVFrames(:,:,a) = reshape(boxVals,boxHeight,boxWidth);
end
% 
% windowFig = figure;
% imagesc(boxCBVFrames(:,:,1))
% colormap gray
% caxis([0 2^12])
% axis image
% windowMask = roipoly;
% close(windowFig)
% 
% for w = 1:size(boxCBVFrames,3)
%     windowFrame = boxCBVFrames(:,:,w);
%     windowFrame(~windowMask) = NaN;
%     windowFrames(:,:,w) = windowFrame;
% end

baselineFrame = mean(boxCBVFrames,3);
for a = 1:size(boxCBVFrames,3)
    normCBVFrames(:,:,a) = ((boxCBVFrames(:,:,a) - baselineFrame)./(baselineFrame)).*100;
end

% cbvHandle = implay(normCBVFrames,cbvFs);
% cbvHandle.Visual.ColorMap.UserRange = 1; 
% cbvHandle.Visual.ColorMap.UserRangeMin = -10; 
% cbvHandle.Visual.ColorMap.UserRangeMax = 10;
electrodeInput = input('Input the cortical hemisphere (LH/RH): ','s'); disp(' ')
electrodeHem = ['cortical_' electrodeInput];

%% Whisker 
whiskImageHeight = RawData.notes.whiskCamPixelHeight;                                                                                                            
whiskImageWidth = RawData.notes.whiskCamPixelWidth;
whiskFs = RawData.notes.whiskCamSamplingRate;

whiskFrameStart = floor(startTime)*whiskFs;
whiskFrameEnd = floor(endTime)*whiskFs;         
whiskFrameInds = whiskFrameStart:whiskFrameEnd;

whiskerPixelsPerFrame = whiskImageWidth*whiskImageHeight;
whiskerSkippedPixels = whiskerPixelsPerFrame*2; % Multiply by two because there are 16 bits (2 bytes) per pixel
whiskFid = fopen(whiskCamFileID);
fseek(whiskFid,0,'eof');
whiskNFramesToRead = length(whiskFrameInds);
whiskImageStack = zeros(whiskImageHeight,whiskImageWidth,whiskNFramesToRead);
for a = 1:whiskNFramesToRead
    fseek(whiskFid,whiskFrameInds(a)*whiskerSkippedPixels,'bof');
    whiskZ = fread(whiskFid,whiskerPixelsPerFrame,'*uint8','b');
    whiskImg = reshape(whiskZ(1:whiskerPixelsPerFrame),whiskImageWidth,whiskImageHeight);
    whiskImageStack(:,:,a) = flip(imrotate(whiskImg,-90),2);
end
fclose('all');

c = 1;
for b = 1:size(whiskImageStack,3)
    if rem(b,5) == 1
        dsWhiskImageStack(:,:,c) = whiskImageStack(:,:,b);
        c = c + 1;
    end
end

% whiskHandle = implay(dsWhiskImageStack,30);
% whiskHandle.Visual.ColorMap.UserRange = 1; 
% whiskHandle.Visual.ColorMap.UserRangeMin = min(whiskImg(:)); 
% whiskHandle.Visual.ColorMap.UserRangeMax = max(whiskImg(:));

%% Pupil
pupilImageHeight = RawData.notes.pupilCamPixelHeight;                                                                                                            
pupilImageWidth = RawData.notes.pupilCamPixelWidth;
pupilFs = RawData.notes.pupilCamSamplingRate;

pupilFrameStart = floor(startTime)*pupilFs;
pupilFrameEnd = floor(endTime)*pupilFs;         
pupilFrameInds = pupilFrameStart:pupilFrameEnd;

pupilPixelsPerFrame = pupilImageWidth*pupilImageHeight;
pupilSkippedPixels = pupilPixelsPerFrame*2;   % Multiply by two because there are 16 bits (2 bytes) per pixel
pupilFid = fopen(pupilCamFileID);
fseek(pupilFid,0,'eof');
pupilFileSize = ftell(pupilFid);
fseek(pupilFid,0,'bof');
pupilNFramesToRead = length(pupilFrameInds);
pupilImageStack = zeros(pupilImageWidth,pupilImageHeight,pupilNFramesToRead);
for a = 1:pupilNFramesToRead
    fseek(pupilFid,pupilFrameInds(a)*pupilSkippedPixels,'bof');
    pupilZ = fread(pupilFid,pupilPixelsPerFrame,'*uint8','b');
    pupilImg = reshape(pupilZ(1:pupilPixelsPerFrame),pupilImageHeight,pupilImageWidth);
    pupilImageStack(:,:,a) = flip(imrotate(pupilImg,-90),2);
end
fclose('all');

% pupilHandle = implay(pupilImageStack,pupilFs);
% pupilHandle.Visual.ColorMap.UserRange = 1; 
% pupilHandle.Visual.ColorMap.UserRangeMin = min(pupilImg(:)); 
% pupilHandle.Visual.ColorMap.UserRangeMax = max(pupilImg(:));

%% Reflectance data
[D,C] = butter(3,1/(ProcData.notes.CBVCamSamplingRate/2),'low');
CBV = ProcData.data.CBV.Barrels;
normCBV = (CBV - RestingBaselines.setDuration.CBV.Barrels.(strDay))./(RestingBaselines.setDuration.CBV.Barrels.(strDay));
filtCBV = filtfilt(D,C,normCBV)*100;

%% HbT data

%% Neural data
T = SpecData.(electrodeHem).fiveSec.T;
F = SpecData.(electrodeHem).fiveSec.F;
normS = SpecData.(electrodeHem).fiveSec.normS;

%% movie file comparing rgb with original data
outputVideo = VideoWriter([animalID '_' fileID '_PresentationVideo.avi']);
fps = 30;   % default fps from video acquisition
speedUp = 1;   % speed up by factor of
outputVideo.FrameRate = fps*speedUp;
open(outputVideo);
cbvVector = NaN(1,length(filtCBV));
fig = figure('Position',get(0,'Screensize'));
T1 = (0:1/6:(endTime-startTime));
T2 = 1;
neuralMatrix = NaN(length(F),length(T1));
for a = 1:size(pupilImageStack,3)
    sgtitle({[animalID ' ' strrep(fileID,'_',' ') ' Presentation Example'],' '})
    % window movie
    s1 = subplot(3,3,1);
    imagesc(normCBVFrames(:,:,a))
    title('Pixel reflectance')
    colormap(gca,'gray')
    colorbar
    caxis([-7.5 7.5])
    axis image
    axis off
    % pupil movie
    s2 = subplot(3,3,2);
    imagesc(pupilImageStack(:,:,a))
    title('Pupil camera')
    colormap(gca,'gray')
    caxis([min(pupilImg(:)) max(pupilImg(:))])
    axis image
    axis off
    % whisker movie
    s3 = subplot(3,3,3);
    imagesc(dsWhiskImageStack(:,:,a))
    title('Whisker camera')
    colormap(gca,'gray')
    caxis([min(whiskImg(:)) max(whiskImg(:))])
    axis image
    axis off
    % reflectance or HbT
    s4 = subplot(3,3,4:6);
    pos = get(gcf, 'Position'); % gives x left, y bottom, width, height
    x = pos(1);
    y = pos(2);
    w = pos(3);
    h = pos(4);
    cbvVector(1,a) = filtCBV(1,a + startTime*cbvFs);
    plot((1:length(cbvVector))/cbvFs,cbvVector,'k')
    title('Mean pixel (CBV) reflectance')
    xlabel('Time (sec)')
    ylabel('\DeltaR/R')
    xlim([0 (size(pupilImageStack,3))/pupilFs])
    % neural data
    s5 = subplot(3,3,7:9);
    if rem(a,5) == 1
        neuralMatrix(:,T2) = normS(:,a + startTime*6);
        T2 = T2 + 1;
    end
    semilog_imagesc_IOS(T1,F,neuralMatrix,'y')
    title('Cortical spectrogram')
    xlabel('Time (sec)')
    ylabel('Freq (Hz)')
    xlim([0 (size(pupilImageStack,3))/pupilFs])
    colormap(gca,'parula')
    cbar = colorbar;
    caxis([-1 2])
    axis xy
    set(gca,'TickLength',[0, 0])
    set(gca,'box','off')
    
    s4Pos = get(s4,'position');
    s5Pos = get(s5,'position');
    s5Pos(3:4) = s4Pos(3:4);
    set(s5,'position',s5Pos);
    
    currentFrame = getframe(fig);    writeVideo(outputVideo,currentFrame);
end
close(outputVideo)
close(fig)
