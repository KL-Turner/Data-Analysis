function [data] = FuncRunPupilTracker(procDataFileID)
%________________________________________________________________________________________________________________________
% Written by Kyle W. Gheres & Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Track changes in pupil area and detect periods of blinking
%________________________________________________________________________________________________________________________

load(procDataFileID)
% draw an ROI around the eye
[~,~,fileID] = GetFileInfo_IOS(procDataFileID);
pupilCamFileID = [fileID '_PupilCam.bin'];
fid = fopen(pupilCamFileID); % reads the binary file in to the work space
fseek(fid,0,'eof'); % find the end of the video frame
imageHeight = ProcData.notes.pupilCamPixelHeight; %#ok<*NODEF>
imageWidth = ProcData.notes.pupilCamPixelWidth;
pixelsPerFrame = imageWidth*imageHeight;
skippedPixels = pixelsPerFrame;
roiImage = zeros(imageHeight,imageWidth,1);
fseek(fid,1*skippedPixels,'bof'); % read .bin File to roiImage
z = fread(fid,pixelsPerFrame,'*uint8','b');
img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
roiImage(:,:,1) = flip(imrotate(img,-90),2);
roiImage = uint8(roiImage); % convert double floating point data to unsignned 8bit integers
data.workingImg = imcomplement(roiImage); % grab frame from image stack
disp('Draw ROI around eye'); disp(' ')
eyeFigure = figure;
title('Draw ROI around eye')
[eyeROI,data.x12,data.y12] = roipoly(data.workingImg);
close(eyeFigure)
% model the distribution of pixel intensities as a gaussian to estimate/isolate the population of pupil pixels
threshSet = 4.5; % StD beyond mean intensity to binarize image for pupil tracking
medFiltParams = [5,5]; % [x,y] dimensions for 2d median filter of images
pupilHistEdges = 1:1:256; % camera data is unsigned 8bit integers. Ignore 0 values
filtImg = medfilt2(data.workingImg,medFiltParams); % median filter image
threshImg = uint8(double(filtImg).*eyeROI); % only look at pixel values in ROI
[phat,~] = mle(reshape(threshImg(threshImg ~= 0),1,numel(threshImg(threshImg ~= 0))),'distribution','Normal');
intensityThresh = phat(1) + (threshSet*phat(2)); % set threshold as 4.5 sigma above population mean estimated from MLE
testFig = figure;
data.pupilHist = histogram(threshImg((threshImg ~= 0)),'BinEdges',pupilHistEdges,'Normalization','Probability');
data.theFit = pdf('normal',data.pupilHist.BinEdges,phat(1),phat(2)); % generate distribution from mle fit of data
data.normFit = data.theFit./sum(data.theFit); % normalize fit so sum of gaussian ==1
data.pupilHistEdges = pupilHistEdges;
data.threshImg = threshImg;
hold on;
plot(data.pupilHist.BinEdges,data.normFit,'r','LineWidth',2);
xline(intensityThresh,'--m','LineWidth',1);
title('Histogram of image pixel intensities')
xlabel('Pixel intensities');
ylabel('Bin Counts');
legend({'Normalized Bin Counts','MLE fit of data','Starting 4.5 StD ROI threshold'},'Location','northwest');
xlim([0,256]);
axis square
% figure for verifying pupil threshold
testImg = threshImg;
testImg(threshImg >= intensityThresh) = 1;
testImg(threshImg < intensityThresh) = 0;
testThresh = labeloverlay(roiImage(:,:,1),testImg);
threshFig = figure;
imagesc(testThresh);
colormap gray
axis off
axis square
% check threshold
threshOK = false;
while threshOK == false
    disp(['Intensity threshold: ' num2str(intensityThresh)]); disp (' ')
    threshCheck = input('Is pupil threshold value ok? (y/n): ','s'); disp(' ')
    if strcmp(threshCheck,'y') == true
        threshOK = true;
    else
        intensityThresh = input('Manually set pupil intensity threshold: '); disp(' ')
        testImg(threshImg >= intensityThresh) = 1;
        testImg(threshImg < intensityThresh) = 0;
        testThresh = labeloverlay(roiImage(:,:,1),testImg);
        imagesc(testThresh);
        colormap gray
        axis image
        axis off
        title('Pixels above threshold');
    end
end
data.intensityThresh = intensityThresh;
close(testFig)
close(threshFig)
% run pupil/blink tracking on all data files
theAngles = 1:1:180; % projection angles measured during radon transform of pupil
radonThresh = 0.05; % arbitrary threshold used to clean up radon transform above values == 1 below == 0
pupilThresh = 0.25; % arbitrary threshold used to clean up inverse radon transform above values == 1 below == 0
blinkThresh = 0.35; % arbitrary threshold used to binarize data for blink detection above values == 1 below == 0
medFiltParams = [5,5]; % [x,y] dimensions for 2d median filter of images
fid = fopen(pupilCamFileID); % reads the binary file in to the work space
fseek(fid,0,'eof'); % find the end of the video frame
fileSize = ftell(fid); % calculate file size
fseek(fid,0,'bof'); % find the begining of video frames
imageHeight = ProcData.notes.pupilCamPixelHeight; % how many pixels tall is the frame
imageWidth = ProcData.notes.pupilCamPixelWidth; % how many pixels wide is the frame
pixelsPerFrame = imageWidth*imageHeight;
skippedPixels = pixelsPerFrame;
nFramesToRead = floor(fileSize/(pixelsPerFrame));
imageStack = zeros(200,200,nFramesToRead);
% read .bin file to imageStack
for dd = 1:nFramesToRead
    fseek(fid,(dd - 1)*skippedPixels,'bof');
    z = fread(fid,pixelsPerFrame,'*uint8','b');
    img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
    imageStack(:,:,dd) = flip(imrotate(img,-90),2);
end
% convert double floating point data to unsignned 8bit integers
imageStack = uint8(imageStack);
% grab frame from image stack
data.workingImg = imcomplement(imageStack(:,:,2));
% pre-allocate empty structures
pupilArea(1:size(imageStack,3)) = NaN; % area of pupil
pupilMajor(1:size(imageStack,3)) = NaN; % length of major axis of pupil
pupilMinor(1:size(imageStack,3)) = NaN; % length of minor axis of pupil
pupilCentroid(1:size(imageStack,3),2) = NaN; % center of pupil
pupilBoundary(1:size(imageStack,1),1:size(imageStack,2),1:size(imageStack,3)) = NaN;
procStart = tic;
disp(['Running pupil tracker for: ' pupilCamFileID]); disp(' ')
imageFrames = gpuArray(imageStack);
roiInt(1:size(imageFrames,3)) = NaN;
roiInt = gpuArray(roiInt);
correctedFlag=false;
last_frame_ok=1;
for frameNum = 1:size(imageStack,3)
    filtImg = medfilt2(imcomplement(imageFrames(:,:,frameNum)),medFiltParams);
    % only look at pixel values in ROI
    threshImg = uint8(double(filtImg).*eyeROI);
    roiInt_temp = sum(threshImg,1);
    roiInt(frameNum) = sum(roiInt_temp,2);
    isoPupil = threshImg;
    isoPupil(isoPupil < intensityThresh) = 0;
    isoPupil(isoPupil >= intensityThresh) = 1;
    isoPupil = medfilt2(isoPupil,medFiltParams);
    radPupil = radon(isoPupil);
    minPupil = min(radPupil,[],1);
    minMat = repmat(minPupil,size(radPupil,1),1);
    maxMat = repmat(max((radPupil - minMat),[],1),size(radPupil,1),1);
    % normalize each projection angle to its min and max values. Each value should now be between [0 1]
    normPupil = (radPupil - minMat)./maxMat;
    threshPupil = normPupil;
    % binarize radon projection
    threshPupil(normPupil >= radonThresh) = 1;
    threshPupil(normPupil < radonThresh) = 0;
    % transform back to image space
    radonPupil = gather(iradon(double(threshPupil),theAngles,'linear','Hamming',size(data.workingImg,2)));
    if frameNum == 1
        data.saveRadonImg = radonPupil;
    end
    % find area corresponding to pupil on binary image
    [~,pupilBoundaries] = bwboundaries(radonPupil > pupilThresh*max(radonPupil(:)),8,'noholes');
    fillPupil = pupilBoundaries;
    % fill any subthreshold pixels inside the pupil boundary
    fillPupil = imfill(fillPupil,8,'holes');
    areaFilled = regionprops(fillPupil,'FilledArea','Image','FilledImage','Centroid','MajorAxisLength','MinorAxisLength');
    if frameNum > 1
        if abs((roiInt(frameNum) - roiInt(frameNum - 1))/roiInt(frameNum)) >= blinkThresh % Exclude fitting frame during blinks
            areaFilled = [];
            fillPupil(:) = 0;
        end
    end
    if ~isempty(areaFilled) %Is an pupil identified
        if size(areaFilled,1) > 1 %Is the pupil fragmented in to multiple ROI
            clear theArea areaLogical
            for num = 1:size(areaFilled,1)
                theArea(num) = areaFilled(num).FilledArea; %#ok<*SAGROW>
            end
            maxArea = max(theArea);%Find the ROI with the largest area
            areaLogical = theArea == maxArea;
            areaFilled = areaFilled(areaLogical);
            areaFilled = areaFilled(1);
            % check for aberrant pupil diameter changes
            if frameNum > 1
                fracChange = (maxArea - pupilArea(frameNum-1))/pupilArea(frameNum-1); %frame-wise fractional change
                volFlag = fracChange <- 0.1; % does the change exceed a 10% reduction in pupil size
                if ~isnan(pupilArea(frameNum - 1)) %does the current frame follow a blink
                    if volFlag==true
                        if correctedFlag==false
                            last_frame_ok=frameNum-1;
                            correctedFlag=true;
                        end
                        % correct aberrant diameters by altering radon threshold
                        pupilSweep = pupilThresh - (0:0.01:pupilThresh);
                        for sweepNum = 2:size(pupilSweep,2)
                            if  volFlag == true
                                sweepArea=[];
                                [~,sweepBoundaries] = bwboundaries(medfilt2(radonPupil,[7 7]) > pupilSweep(sweepNum)*max(radonPupil(:)),8,'noholes');
                                fillSweep = sweepBoundaries;
                                fillSweep =imfill(fillSweep,8,'holes');
                                areaSweep = regionprops(fillSweep,'FilledArea','Image','FilledImage','Centroid','MajorAxisLength','MinorAxisLength');
                                for num = 1:size(areaSweep,1)
                                    sweepArea(num) = areaSweep(num).FilledArea;
                                end
                                maxSweep = max(sweepArea);
                                sweepLogical = sweepArea == maxSweep;
                                if sum(sweepLogical) > 1
                                    x = true;
                                    for aa = 1:length(sweepLogical)
                                        if x == true && sweepLogical(1,aa) == 1
                                            sweepLogical(1,aa) = 1;
                                            x = false;
                                        else
                                            sweepLogical(1,aa) = 0;
                                        end
                                    end
                                end
                                fracChange = (maxSweep - pupilArea(last_frame_ok)/pupilArea(last_frame_ok));
                                volFlag = fracChange < -0.1;
                            end
                        end
                        % this can be used to insert NaN if the change is > 10%
                        if  abs(fracChange) < 0.1 %changed to only fill data withing a +/- 10% change in area
                            %    fillPupil(:) = 0;
                            %    areaFilled.FilledArea = NaN;
                            %    areaFilled.MajorAxisLength = NaN;
                            %    areaFilled.MinorAxisLength = NaN;
                            %    areaFilled.Centroid = [NaN,NaN];
                            %else
                            fillPupil = fillSweep;
                            areaFilled = areaSweep(sweepLogical);
                        end
                        if ~exist('correctedFrames','var')
                            frameInd = 1;
                            correctedFrames(frameInd) = frameNum;
                        else
                            frameInd = frameInd + 1;
                            correctedFrames(frameInd) = frameNum;
                        end
                    else
                        correctedFlag=false;
                        last_frame_ok=frameNum;
                    end
                end
            end
        else
            if frameNum > 1
                % check for aberrant pupil diameter changes
                fracChange = (areaFilled.FilledArea - pupilArea(frameNum-1))/pupilArea(frameNum-1);
                volFlag = fracChange < -0.1;
                if ~isnan(pupilArea(frameNum - 1))
                    if volFlag==true
                        if correctedFlag==false
                            last_frame_ok=frameNum-1;
                            correctedFlag=true;
                        end
                        % correct aberrant diameters by altering radon threshold
                        pupilSweep = pupilThresh - (0:0.01:pupilThresh);
                        for sweepNum = 2:size(pupilSweep,2)
                            if volFlag == true
                                sweepArea = [];
                                [~,sweepBoundaries] = bwboundaries(medfilt2(radonPupil,[7 7]) > pupilSweep(sweepNum)*max(radonPupil(:)),8,'noholes');
                                fillSweep = sweepBoundaries;
                                fillSweep = imfill(fillSweep,8,'holes');
                                areaSweep = regionprops(fillSweep,'FilledArea','Image','FilledImage','Centroid','MajorAxisLength','MinorAxisLength');
                                for num = 1:size(areaSweep,1)
                                    sweepArea(num) = areaSweep(num).FilledArea; %#ok<*AGROW>
                                end
                                maxSweep = max(sweepArea);
                                sweepLogical = sweepArea == maxSweep;
                                fracChange = (maxSweep - pupilArea(last_frame_ok))/pupilArea(last_frame_ok);
                                volFlag = fracChange < -0.1;
                            end
                        end
                        % this can be used to insert NaN if the change is > 10%
                        if abs(fracChange) < 0.1
                            %    fillPupil(:) = 0;
                            %    areaFilled.FilledArea = NaN;
                            %    areaFilled.MajorAxisLength = NaN;
                            %    areaFilled.MinorAxisLength = NaN;
                            %    areaFilled.Centroid = [NaN,NaN];
                            % else
                            fillPupil = fillSweep;
                            areaFilled = areaSweep(sweepLogical);
                        end
                        if ~exist('correctedFrames','var')
                            frameInd = 1;
                            correctedFrames(frameInd)=frameNum;
                        else
                            frameInd = frameInd + 1;
                            correctedFrames(frameInd) = frameNum;
                        end
                    else
                        correctedFlag=false;
                        last_frame_ok=frameNum;
                    end
                end
            end
        end
        pupilArea(frameNum) = areaFilled.FilledArea;
        pupilMajor(frameNum) = areaFilled.MajorAxisLength;
        pupilMinor(frameNum) = areaFilled.MinorAxisLength;
        pupilCentroid(frameNum,:) = areaFilled.Centroid;
        pupilBoundary(:,:,frameNum) = fillPupil;
        holdMat = labeloverlay(imageStack(:,:,frameNum),fillPupil,'Transparency',0.8);
        if size(holdMat,3) == 1
            overlay(:,:,:,frameNum) = repmat(holdMat,1,1,3);
        else
            overlay(:,:,:,frameNum) = holdMat;
        end
    else
        pupilArea(frameNum) = NaN;
        pupilMajor(frameNum) = NaN;
        pupilMinor(frameNum) = NaN;
        pupilCentroid(frameNum,:) = NaN;
        pupilBoundary(:,:,frameNum) = fillPupil;
        holdMat = labeloverlay(imageStack(:,:,frameNum),fillPupil);
        overlay(:,:,:,frameNum) = repmat(holdMat,1,1,3);
    end
end
data.pupilArea = pupilArea;
data.pupilMajor = pupilMajor;
data.pupilMinor = pupilMinor;
data.pupilCentroid = pupilCentroid;
data.eyeROI = eyeROI;
data.roiIntensity = gather(roiInt);
proceEnd = toc(procStart);
procMin = proceEnd/60;
minText = num2str(procMin);
procSec = round(str2double(minText(2:end))*60,0);
secText = num2str(procSec);
disp(['File processing time: ' minText(1) ' min ' secText ' seconds']); disp(' ')
blinks = find((abs(diff(data.roiIntensity))./data.roiIntensity(2:end)) >= blinkThresh) + 1;
data.blinkFrames = overlay(:,:,:,blinks);
data.blinkInds = blinks;
data.overlay = overlay(:,:,:,1);
% patch NaNs
data.nanPupilPatch = fillmissing(data.pupilArea,'movmedian',31);
% patch sudden spikes
diffArea = abs(diff(data.nanPupilPatch));
% threshold for interpolation
threshold = 500;
diffIndex = diffArea > threshold;
[linkedDiffIndex] = LinkBinaryEvents_IOS(gt(diffIndex,0),[30,0]);
edgeFound = false;
% identify edges for interpolation
xx = 1;
for aa = 1:length(linkedDiffIndex)
    if edgeFound == false
        if linkedDiffIndex(1,aa) == 1
            startEdge(xx,1) = aa;
            edgeFound = true;
        end
    elseif edgeFound == true
        if linkedDiffIndex(1,aa) == 0
            endEdge(xx,1) = aa;
            xx = xx + 1;
            edgeFound = false;
        end
    end
end
% fill from start:ending edges of rapid pupil fluctuations that weren't NaN
nanPupilArea = pupilArea;
for aa = 1:length(startEdge)
    startTime = startEdge(aa,1);
    endTime = endEdge(aa,1);
    nanPupilArea(startTime:endTime) = NaN;
    data.patchedPupilArea = fillmissing(nanPupilArea,'spline');
end

end


