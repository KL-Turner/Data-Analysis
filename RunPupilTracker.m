%________________________________________________________________________________________________________________________
% Written by Kyle W. Gheres & Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Track changes in pupil area and detect periods of blinking
%________________________________________________________________________________________________________________________

clear; clc;
procDataFileID = uigetfile('*_ProcData.mat');
load(procDataFileID)
% draw an ROI around the eye
[~,~,fileID] = GetFileInfo(procDataFileID);
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
workingImg = imcomplement(roiImage); % grab frame from image stack
disp('Draw ROI around eye'); disp(' ')
eyeFigure = figure;
title('Draw ROI around eye')
[eyeROI,x12,y12] = roipoly(workingImg);
close(eyeFigure)
% model the distribution of pixel intensities as a gaussian to estimate/isolate the population of pupil pixels
threshSet = 4.5; % StD beyond mean intensity to binarize image for pupil tracking
medFiltParams = [5,5]; % [x,y] dimensions for 2d median filter of images
pupilHistEdges = 1:1:256; % camera data is unsigned 8bit integers. Ignore 0 values
filtImg = medfilt2(workingImg,medFiltParams); % median filter image
threshImg = uint8(double(filtImg).*eyeROI); % only look at pixel values in ROI
[phat,~] = mle(reshape(threshImg(threshImg ~= 0),1,numel(threshImg(threshImg ~= 0))),'distribution','Normal');
intensityThresh = phat(1) + (threshSet*phat(2)); % set threshold as 4.5 sigma above population mean estimated from MLE
figure
pupilHist = histogram(threshImg((threshImg ~= 0)),'BinEdges',pupilHistEdges,'Normalization','Probability');
theFit = pdf('normal',pupilHist.BinEdges,phat(1),phat(2)); % generate distribution from mle fit of data
normFit = theFit./sum(theFit); % normalize fit so sum of gaussian ==1
hold on;
plot(pupilHist.BinEdges,normFit,'r','LineWidth',2);
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
figure;
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
samplingRate = ProcData.notes.pupilCamSamplingRate;
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
workingImg = imcomplement(imageStack(:,:,2));
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
    radonPupil = gather(iradon(double(threshPupil),theAngles,'linear','Hamming',size(workingImg,2)));
    if frameNum == 1
        saveRadonImg = radonPupil;
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
            data.overlay(:,:,:,frameNum) = repmat(holdMat,1,1,3);
        else
            data.overlay(:,:,:,frameNum) = holdMat;
        end
    else
        pupilArea(frameNum) = NaN;
        pupilMajor(frameNum) = NaN;
        pupilMinor(frameNum) = NaN;
        pupilCentroid(frameNum,:) = NaN;
        pupilBoundary(:,:,frameNum) = fillPupil;
        holdMat = labeloverlay(imageStack(:,:,frameNum),fillPupil);
        data.overlay(:,:,:,frameNum) = repmat(holdMat,1,1,3);
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
data.blinkFrames = data.overlay(:,:,:,blinks);
plotPupilArea = data.pupilArea;
plotPupilArea((blinks - 1:blinks + 1)) = NaN;
blinkTimes(1:size(data.pupilArea,2)) = NaN;
blinkTimes(blinks) = 1;
data.blinkInds = blinks;
% patch NaNs
data.nanPupilPatch = fillmissing(data.pupilArea,'movmedian',31);
% patch sudden spikes
diffArea = abs(diff(data.nanPupilPatch));
% threshold for interpolation
threshold = 500;
diffIndex = diffArea > threshold;
[linkedDiffIndex] = LinkBinaryEvents(gt(diffIndex,0),[30,0]);
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
% subplot for eye ROI
figure
subplot(2,2,1)
imagesc(workingImg)
hold on;
x1 = plot(x12,y12,'color','r','LineWidth',1');
title('ROI to measure changes in pupil area')
legend(x1,'eye ROI')
colormap gray
axis image
axis off
% subplot for ROI histrogram and threshold
subplot(2,2,2)
pupilHist = histogram(threshImg((threshImg ~= 0)),'BinEdges',pupilHistEdges,'Normalization','Probability');
theFit = pdf('normal',pupilHist.BinEdges,phat(1),phat(2)); % generate distribution from mle fit of data
normFit = theFit./sum(theFit); % normalize fit so sum of gaussian ==1
hold on;
plot(pupilHist.BinEdges,normFit,'r','LineWidth',2);
xline(intensityThresh,'--c','LineWidth',1);
title('Histogram of image pixel intensities')
xlabel('Pixel intensities');
ylabel('Bin Counts');
legend({'Normalized Bin Counts','MLE fit of data','Pupil intensity threshold'},'Location','northwest');
xlim([0,256]);
axis square
% subplot for radon transform
subplot(2,2,3)
imagesc(saveRadonImg)
title('Radon transform back to image space')
colormap gray
axis image
axis off
% subplot for measured pupil area
subplot(2,2,4)
imagesc(data.overlay(:,:,:,1));
title('Calculated pupil area')
colormap gray
axis image
% patch falsely rapid transitions
figure
sgtitle('Pupil area changes');
subplot(2,1,1)
hold on;
nanVals = isnan(data.pupilArea);
nanInds = find(nanVals == 1);
for aa = 1:length(nanInds)
    x1 = xline(nanInds(1,aa)/samplingRate,'g');
end
difInds = find(diffIndex == 1);
for bb = 1:length(difInds)
    x2 = xline(difInds(1,bb)/samplingRate,'b');
end
p1 = plot((1:length(data.pupilArea))/samplingRate,data.pupilArea,'k','LineWidth',1);
xlabel('Time (sec)');
ylabel('Pupil area (pixels)');
legend([p1,x1,x2],'Original pupil area','ROI loss','ROI inaccuracy');
set(gca,'box','off')
axis tight
subplot(2,1,2)
plot((1:length(data.patchedPupilArea))/samplingRate,medfilt1(data.patchedPupilArea,7),'k','LineWidth',1);
xlabel('Time (sec)');
ylabel('Pupil area (pixels)');
legend('Patched & medfilt pupil area');
set(gca,'box','off')
axis tight

function [animalID,fileDate,fileID] = GetFileInfo(fileName)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Identify important aspects of a file name and output each individually
%________________________________________________________________________________________________________________________

% Identify the extension
extInd = strfind(fileName(1,:),'.');
extension = fileName(1,extInd + 1:end);
% Identify the underscores
fileBreaks = strfind(fileName(1,:),'_');
switch extension
    case 'bin'
        animalID = [];
        fileDate = fileName(:,1:fileBreaks(1) - 1);
        fileID = fileName(:,1:fileBreaks(4) - 1);
    case 'mat'
        % Use the known format to parse
        animalID = fileName(:,1:fileBreaks(1) - 1);
        if numel(fileBreaks) > 3
            fileDate = fileName(:,fileBreaks(1) + 1:fileBreaks(2) - 1);
            fileID = fileName(:,fileBreaks(1) + 1:fileBreaks(5) - 1);
        else
            fileDate = [];
            fileID = [];
        end
end

end

function [linkedWF] = LinkBinaryEvents(binWF,dCrit)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
% Purpose: Link binary events that occur within a certain amount of time
%________________________________________________________________________________________________________________________

% Identify Edges, control for trial start/stop
dBinWF = diff(gt(binWF,0));
upInd = find(dBinWF == 1);
downInd = find(dBinWF == -1);
if binWF(end) > 0
    downInd = [downInd,length(binWF)];
end
if binWF(1) > 0
    upInd = [1,upInd];
end
% Link periods of bin_wf==0 together if less than dCrit(1). Calculate time between events
brkTimes = upInd(2:length(upInd)) - downInd(1:(length(downInd) - 1));
% Identify times less than user-defined period
sub_dCritDowns = find(lt(brkTimes,dCrit(1)));
% Link any identified breaks together
if isempty(sub_dCritDowns) == 0
    for d = 1:length(sub_dCritDowns)
        start = downInd(sub_dCritDowns(d));
        stop = upInd(sub_dCritDowns(d) + 1);
        binWF(start:stop) = 1;
    end
end
% Link periods of bin_wf==1 together if less than dCrit(2)
hitimes = downInd - upInd;
blips = find(lt(hitimes,dCrit(2)) == 1);
if isempty(blips) == 0
    for b = 1:length(blips)
        start = upInd(blips(b));
        stop = downInd(blips(b));
        binWF(start:stop) = 0;
    end
end
linkedWF = binWF;

end
