%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Track pupil diameter and blinking
%________________________________________________________________________________________________________________________

clear; clc; close all
% Character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% create/load pre-existing ROI file with the coordinates
ROIFileDir = dir('*_PupilData.mat');
if isempty(ROIFileDir) == true
    PupilData = [];
    PupilData.EyeROI = [];
else
    ROIFileName = {ROIFileDir.name}';
    ROIFileID = char(ROIFileName);
    load(ROIFileID);
end
% establish the number of unique days based on file IDs
[~,fileDates,~] = GetFileInfo_IOS(procDataFileIDs);
[uniqueDays,~,DayID] = GetUniqueDays_IOS(fileDates);
firstFileOfDay = cell(1,length(uniqueDays));
for aa = 1:length(uniqueDays)
    FileInd = DayID == aa;
    dayFilenames = procDataFileIDs(FileInd,:);
    firstFileOfDay(aa) = {dayFilenames(1,:)};
end
% Create the desired window ROI for each day if it doesn't yet exist
for bb = 1:length(firstFileOfDay)
    firstFile = firstFileOfDay{1,bb};
    load(firstFile)
    [animalID,fileDate,fileID] = GetFileInfo_IOS(firstFile);
    strDay = ConvertDate_IOS(fileDate);
    if ~isfield(PupilData.EyeROI,(strDay))
        pupilCamFileID = [fileID '_PupilCam.bin'];
        fid = fopen(pupilCamFileID); % reads the binary file in to the work space
        fseek(fid,0,'eof'); % find the end of the video frame
        pupilHistEdges = (1:1:256); % camera data is unsigned 8bit integers. Ignore 0 values
        imageHeight = ProcData.notes.pupilCamPixelHeight;
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
        disp('Draw roi around eye'); disp(' ')
        eyeFigure = figure;
        annotation('textbox',[0.4,0.9,0.1,0.1],'String','Draw ROI around eye','FitBoxToText','on','LineStyle','none','FontSize',16);
        [eyeROI] = roipoly(workingImg);
        close(eyeFigure)
        %% model the distribution of pixel intensities as a gaussian to estimate/isolate the population of pupil pixels
        threshSet = 4.5; % StD beyond mean intensity to binarize image for pupil tracking
        medFiltParams = [5,5]; % [x,y] dimensions for 2d median filter of images
        filtImg = medfilt2(workingImg,medFiltParams); % median filter image
        threshImg = uint8(double(filtImg).*eyeROI); % only look at pixel values in ROI
        [phat,pci] = mle(reshape(threshImg(threshImg ~= 0),1,numel(threshImg(threshImg ~= 0))),'distribution','Normal');
        %% figure for verifying pupil threshold
        pupilROIFig = figure;
        subplot(1,3,1)
        pupilHist = histogram(threshImg((threshImg ~= 0)),'BinEdges',pupilHistEdges);
        xlabel('Pixel intensities');
        ylabel('Bin Counts');
        title('Histogram of image pixel intensities')
        normCounts = pupilHist.BinCounts./sum(pupilHist.BinCounts); % normalizes bin count to total bin counts
        theFit = pdf('normal',pupilHist.BinEdges,phat(1),phat(2)); % generate distribution from mle fit of data
        normFit = theFit./sum(theFit); % normalize fit so sum of gaussian ==1
        intensityThresh = phat(1) + (threshSet*phat(2)); % set threshold as 4.5 sigma above population mean estimated from MLE
        testImg = threshImg;
        testImg(threshImg >= intensityThresh) = 1;
        testImg(threshImg < intensityThresh) = 0;
        testThresh = labeloverlay(roiImage(:,:,1),testImg);
        axis square
        subplot(1,3,2)
        plot(pupilHist.BinEdges(2:end),normCounts,'k','LineWidth',1);
        xlabel('Pixel intensities');
        ylabel('Normalized bin counts');
        title('Normalized histogram and MLE fit of histogram');
        hold on;
        plot(pupilHist.BinEdges,normFit,'r','LineWidth',2);
        xline(intensityThresh,'--c','LineWidth',1);
        legend({'Normalized Bin Counts','MLE fit of data','Pixel intensity threshold'},'Location','northwest');
        xlim([0,256]);
        axis square
        subplot(1,3,3)
        imshow(testThresh);
        title('Pixels above threshold');
        %% check threshold
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
                manualThreshFig = figure;
                imshow(testThresh);
                title('Pixels above threshold');
            end
        end
        close all
        PupilData.EyeROI.(strDay) = eyeROI;
        PupilData.Threshold.(strDay) = intensityThresh;
        save([animalID '_PupilData.mat'],'PupilData');
    end
end
%% run pupil/blink tracking on all data files
theAngles = 1:1:180; % projection angles measured during radon transform of pupil
pupilHistEdges = 1:1:256; % camera data is unsigned 8bit integers. Ignore 0 values
radonThresh = 0.05; % arbitrary threshold used to clean up radon transform above values ==1 below ==0
pupilThresh = 0.35; % arbitrary threshold used to clean up inverse radon transform above values ==1 below ==0
blinkThresh = 0.35; % arbitrary threshold used to binarize data for blink detection above values ==1 below ==0
medFiltParams = [5,5]; % [x,y] dimensions for 2d median filter of images
for cc = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(cc,:);
    load(procDataFileID)
    disp(['Checking pupil tracking status of file (' num2str(cc) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
    %     if isfield(ProcData.data,'Pupil') == false
    [animalID,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    pupilCamFileID = [fileID '_PupilCam.bin'];
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
    % read .bin file to imageStack
    for dd = 1:nFramesToRead
        fseek(fid,(dd - 1)*skippedPixels,'bof');
        z = fread(fid,pixelsPerFrame,'*uint8','b');
        img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
        imageStack(:,:,dd) = flip(imrotate(img,-90),2);
    end
    imageStack = uint8(imageStack); % convert double floating point data to unsignned 8bit integers
    workingImg = imcomplement(imageStack(:,:,2)); % grab frame from image stack
    % pre-allocate empty structures
    pupilArea(1:size(imageStack,3)) = NaN; % area of pupil
    pupilMajor(1:size(imageStack,3)) = NaN; % length of major axis of pupil
    pupilMinor(1:size(imageStack,3)) = NaN; % length of minor axis of pupil
    pupilCentroid(1:size(imageStack,3),2) = NaN; % center of pupil
    pupilBoundary(1:size(imageStack,1),1:size(imageStack,2),1:size(imageStack,3)) = NaN;
    % select ROI containing pupil
    eyeROI = PupilData.EyeROI.(strDay);
    % select pupil intensity threshold
    intensityThresh = PupilData.Threshold.(strDay);
    procStart = tic;
    disp(['Running pupil tracker for: ' pupilCamFileID]); disp(' ')
    imageFrames = gpuArray(imageStack);
    roiInt(1:size(imageFrames,3)) = NaN;
    roiInt = gpuArray(roiInt);
    for frameNum = 1:size(imageStack,3)
        filtImg = medfilt2(imcomplement(imageFrames(:,:,frameNum)),medFiltParams);
        threshImg = uint8(double(filtImg).*eyeROI); % only look at pixel values in ROI
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
        normPupil = (radPupil - minMat)./maxMat; % normalize each projection angle to its min and max values. Each value should now be between [0 1]
        threshPupil = normPupil;
        threshPupil(normPupil >= radonThresh) = 1;
        threshPupil(normPupil < radonThresh) = 0; % binarize radon projection
        radonPupil = gather(iradon(double(threshPupil),theAngles,'linear','Hamming',size(workingImg,2))); % transform back to image space
        [Pupil_Pix,pupilBoundaries] = bwboundaries(radonPupil > pupilThresh*max(radonPupil(:)),8,'noholes'); % find area corresponding to pupil on binary image
        fillPupil = pupilBoundaries;
        fillPupil = imfill(fillPupil,'holes'); % fill any subthreshold pixels inside the pupil boundary
        areaFilled = regionprops(fillPupil,'FilledArea','Image','FilledImage','Centroid','MajorAxisLength','MinorAxisLength');
        if size(areaFilled,1) > 1
            clear theArea areaLogical
            for num = 1:size(areaFilled,1)
                theArea(num) = areaFilled(num).FilledArea;
            end
            maxArea = max(theArea);
            areaLogical = theArea == maxArea;
            areaFilled = areaFilled(areaLogical);
        end
        if ~isempty(areaFilled)
            pupilArea(frameNum) = areaFilled.FilledArea;
            pupilMajor(frameNum) = areaFilled.MajorAxisLength;
            pupilMinor(frameNum) = areaFilled.MinorAxisLength;
            pupilCentroid(frameNum,:) = areaFilled.Centroid;
            pupilBoundary(:,:,frameNum) = fillPupil;
            holdMat = labeloverlay(imageStack(:,:,frameNum),fillPupil,'Transparency',0.8);
            overlay(:,:,:,frameNum) = holdMat;
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
    ProcData.data.Pupil.pupilArea = pupilArea;
    ProcData.data.Pupil.pupilMajor = pupilMajor;
    ProcData.data.Pupil.pupiMinor = pupilMinor;
    ProcData.data.Pupil.pupilCentroid = pupilCentroid;
    ProcData.data.Pupil.eyeROI = eyeROI;
    ProcData.data.Pupil.roiIntensity = gather(roiInt);
    proceEnd = toc(procStart);
    procMin = proceEnd/60;
    minText = num2str(procMin);
    procSec = round(str2double(minText(2:end))*60,0);
    secText = num2str(procSec);
    disp(['File processing time: ' minText(1) ' min ' secText ' seconds']); disp(' ')
    blinks = find((abs(diff(ProcData.data.Pupil.roiIntensity))./ProcData.data.Pupil.roiIntensity(2:end)) >= blinkThresh) + 1;
    ProcData.data.Pupil.blinkFrames = overlay(:,:,:,blinks);
    rowNum = ceil(size(ProcData.data.Pupil.blinkFrames,4)/4);
    plotPupilArea = ProcData.data.Pupil.pupilArea;
    plotPupilArea((blinks - 5:blinks + 30)) = NaN;
    blinkTimes(1:size(ProcData.data.Pupil.pupilArea,2)) = NaN;
    blinkTimes(blinks) = 1;
    ProcData.data.Pupil.blinkTimes = blinkTimes*(1.1*max(ProcData.data.Pupil.pupilArea(:)));
    ProcData.data.Pupil.blinkInds = blinks;
    % save data
    disp(['Saving pupil tracking results to: ' procDataFileID]); disp(' ')
    save(procDataFileID,'ProcData')
    %% visualize pupil diameter and blink times
    pupilFigure = figure;
    plot((1:length(ProcData.data.Pupil.pupilArea))/samplingRate,plotPupilArea,'k','LineWidth',1);
    hold on;
    scatter((1:length(ProcData.data.Pupil.pupilArea))/samplingRate,ProcData.data.Pupil.blinkTimes,50,'r','filled');
    xlabel('Time (sec)');
    ylabel('Pupil area (pixels)');
    title('Pupil area changes');
    legend('Pupil area','Eyes closed');
    set(gca,'box','off')
    % save the file to directory.
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Pupil Tracking/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(pupilFigure,[dirpath animalID '_' fileID '_PupilTracking']);
    close(pupilFigure)
    %     end
end
