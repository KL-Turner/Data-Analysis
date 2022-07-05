%% approximate blink frames for example
% 8760-8770
% 14215-14216
% 22392-22400
% 23732-23734
exampleRawDataFileID = 'T123_200304_14_32_00_RawData.mat';
load(exampleRawDataFileID,'-mat')
pupilCamFileID = '200304_14_32_00_PupilCam.bin';
% Relevant information from RawData file's notes
trialDuration = RawData.notes.trialDuration_sec;
imageHeight = RawData.notes.pupilCamPixelHeight;                                                                                                            
imageWidth = RawData.notes.pupilCamPixelWidth;
pupilFs = RawData.notes.pupilCamSamplingRate;
% Input time indeces for video file
startTime = 0;
endTime = 900;
% Index binary file to desired frames
frameStart = floor(startTime)*pupilFs;
frameEnd = floor(endTime)*pupilFs;
frameInds = frameStart:frameEnd;
pixelsPerFrame = imageWidth*imageHeight;
skippedPixels = pixelsPerFrame;   % Multiply by two because there are 16 bits (2 bytes) per pixel
fid = fopen(pupilCamFileID);
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
% Play movie
handle = implay(imageStack,pupilFs);
handle.Visual.ColorMap.UserRange = 1; 
handle.Visual.ColorMap.UserRangeMin = min(img(:)); 
handle.Visual.ColorMap.UserRangeMax = max(img(:));
% Constants
theangles = (1:1:180);
RadonThresh = 0.5;
PupilThresh = 0.35;
BlinkThresh =  75;
% Empty Structures
Pupil_Area(1:size(imageStack,3)) = NaN; % Area of pupil
Pupil_Major(1:size(imageStack,3)) = NaN; % Length of major axis of pupil
Pupil_Minor(1:size(imageStack,3)) = NaN; % Length of minor axis of pupil
Pupil_Centroid(1:size(imageStack,3),2) = NaN; % Center of pupil
PupilBoundary(1:size(imageStack,1),1:size(imageStack,2),1:size(imageStack,3)) = NaN;
PupilPix=cell(1,size(imageStack,3));
% Select ROI
WorkingImg = imcomplement(uint8(imageStack(:,:,1))); % grab frame from image stack
fprintf('Draw roi around eye\n')
[BW] = roipoly(WorkingImg);
fprintf('Running\n')
tic
% Framewise pupil area measurement
parfor framenum = 1:size(imageStack,3)
    WorkingImg = imcomplement(uint8(imageStack(:,:,framenum))); % grab frame from image stack
    FiltImg = medfilt2(WorkingImg,[5,5]); % median filter image
    ThreshImg  = uint8(double(FiltImg).*BW); % Only look at pixel values in ROI
    HoldImg = double(ThreshImg);
    HoldImg(HoldImg == 0) = NaN;
    AvgPix = mean(HoldImg(:),'omitnan');
    MnsubImg = HoldImg - AvgPix; % Mean subtrack ROI pixels
    MnsubImg(isnan(MnsubImg)) = 0;
    MnsubImg(MnsubImg < 0) = 0; % set all negative pixel values to 0;
    RadPupil = radon(MnsubImg); % transform movie frame in to radon space
    minPupil = min(RadPupil,[],1);
    minMat = repmat(minPupil,size(RadPupil,1),1);
    MaxMat = repmat(max((RadPupil - minMat),[],1),size(RadPupil,1),1);
    NormPupil = (RadPupil - minMat)./MaxMat; % Normalize each projection angle to its min and max values. Each value should now be between [0 1]
    ThreshPupil = NormPupil;
    ThreshPupil(NormPupil >= RadonThresh) = 1;
    ThreshPupil(NormPupil < RadonThresh) = 0; % Binarize radon projection
    RadonPupil = iradon(double(ThreshPupil > RadonThresh*max(ThreshPupil(:))),theangles,'linear','Hamming',size(WorkingImg,2)); % transform back to image space
    [Pupil_Pix,Pupil_Boundary] = bwboundaries(RadonPupil > PupilThresh*max(RadonPupil(:)),'noholes'); % find area corresponding to pupil on binary image
    numPixels = cellfun(@length,Pupil_Pix);
    [~,idx] = max(numPixels);
    FillPupil = Pupil_Boundary;
    FillPupil(FillPupil ~= idx) = 0;
    FillPupil(FillPupil == idx) = 1;
    FillPupil = imfill(FillPupil,'holes'); % fill any subthreshold pixels inside the pupil boundary
    if framenum == 1
        CheckPupil = labeloverlay(uint8(imageStack(:,:,framenum)),FillPupil);
        figure;
        imshow(CheckPupil);
        title('Detected pupil vs video frame');
    end
    area_filled = regionprops(FillPupil,'FilledArea','Image','FilledImage','Centroid','MajorAxisLength','MinorAxisLength');
    Pupil_Area(framenum) = area_filled.FilledArea;
    Pupil_Major(framenum) = area_filled.MajorAxisLength;
    Pupil_Minor(framenum) = area_filled.MinorAxisLength;
    Pupil_Centroid(framenum,:) = area_filled.Centroid;
    PupilBoundary(:,:,framenum) = FillPupil;
    PupilPix{framenum} = Pupil_Pix{idx};
    Hold = labeloverlay(uint8(imageStack(:,:,framenum)),FillPupil);
    Overlay(:,:,:,framenum) = Hold;  
end
PupilTracker.Pupil_Area = Pupil_Area;
PupilTracker.Pupil_Major = Pupil_Major;
PupilTracker.Pupil_Minor = Pupil_Minor;
PupilTracker.Pupil_Centroid = Pupil_Centroid;
PupilTracker.Pupil_Boundary = PupilBoundary;
PupilTracker.Pupil_Pix = PupilPix;
PupilTracker.Overlay = Overlay;
PupilTracker.Eye_ROI = BW;
toc
Blinks = find(abs(diff(PupilTracker.Pupil_Area)) >= BlinkThresh) + 1;
PlotPupilArea = PupilTracker.Pupil_Area;
PlotPupilArea(Blinks) = NaN;
BlinkTimes(1:size(PupilTracker.Pupil_Area,2)) = NaN;
BlinkTimes(Blinks) = 1;
BlinkTimes = BlinkTimes*(1.1*max(PupilTracker.Pupil_Area(:)));
PupilTime = (1:length(PupilTracker.Pupil_Area))/30;
figure;
plot(PupilTime,PlotPupilArea);
hold on; 
scatter(PupilTime,BlinkTimes);
xlabel('Time (sec)');
ylabel('Pupil area (pixels)');
title('Pupil area changes');
legend('Pupil area','Eyes closed');
