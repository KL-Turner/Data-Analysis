function [ROIs] = CalculateROICorrelationMatrix(ROIs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Analyze the cross-correlation between the gamma-band power and each pixel to properly place a circular 1mm ROI
%________________________________________________________________________________________________________________________
% pull file information and camera frames
fileID = uigetfile('*_WindowCam.bin');
fileID2 = strrep(fileID,'_',' ');
animalID = 'T102';
rawDataFileID = [animalID '_' fileID(1:end - 13) 'RawData.mat'];
procDataFileID = [animalID '_' fileID(1:end - 13) 'ProcData.mat'];
load(rawDataFileID)
load(procDataFileID)
imageHeight = RawData.notes.CBVCamPixelHeight;
imageWidth = RawData.notes.CBVCamPixelWidth;
[frames] = ReadDalsaBinary_IOS(fileID,imageWidth,imageHeight);
% open figure for ROI drawing
windowFig = figure;
imagesc(frames{1})
title([animalID ' ' fileID2])
xlabel('Image size (pixels)')
ylabel('Image size (pixels)')
colormap gray
colorbar
axis image
caxis([0 2^RawData.notes.CBVCamBitDepth])
% draw ROI for the mask over the entire windows
hem = {'LH','RH'};
for a = 1:length(hem)
    isok = false;
    while isok == false
        disp(['Draw an ROI over the ' hem{1,a}]); disp(' ')
        [~,rect] = imcrop;
        hold on;
        ROIoutline = rectangle('Position',rect,'EdgeColor','r');
        checkMask = input('Is the ROI okay? (y/n): ','s'); disp(' ')
        if strcmp(checkMask,'y') == true
            isok = true;
            ROIs.([hem{1,a} '_windowMask']).rect = rect;
        end
        delete(ROIoutline);
    end
end
close(windowFig)
% extract the pixel values from the window ROIs
for b = 1:length(hem)
    imageMask = nan(size(frames{1}));
    rectMask = ROIs.([hem{1,b} '_windowMask']).rect;
    rectMask = round(rectMask);
    imageMask(rectMask(2):(rectMask(2) + rectMask(4)),rectMask(1):(rectMask(1) + rectMask(3))) = 1;
    for c = 1:length(frames) - 1
        frame = frames{1,c};
        frameHold = double(frame).*imageMask;
        imageStack.(hem{1,b})(:,c) = frameHold(~isnan(frameHold));
    end
end
% extract and process gamma band power
[B,A] = butter(3,2/(ProcData.notes.dsFs/2),'low');
LH_gammaBandPower = detrend(filtfilt(B,A,ProcData.data.cortical_LH.gammaBandPower - ProcData.data.cortical_LH.gammaBandPower(1)) + ProcData.data.cortical_LH.gammaBandPower(1),'constant');
RH_gammaBandPower = detrend(filtfilt(B,A,ProcData.data.cortical_RH.gammaBandPower - ProcData.data.cortical_RH.gammaBandPower(1)) + ProcData.data.cortical_RH.gammaBandPower(1),'constant');
% cross correlation
lagTime = 5;   % seconds
maxLag = lagTime*RawData.notes.CBVCamSamplingRate;
for d = 1:length(hem)
    hemisphere = hem{1,d};
    if strcmp(hemisphere,'LH') == true
        gammaBandArray = LH_gammaBandPower;
    elseif strcmp(hemisphere,'RH') == true
        gammaBandArray = RH_gammaBandPower;
    end
    % extract pixel values from each numel index in matrix image
    for e = 1:size(imageStack.(hemisphere),1)
        disp(['Analyzing matrix numel ' num2str(e) ' of ' num2str(size(imageStack.(hemisphere),1))]); disp(' ')
        pixelArray = imageStack.(hemisphere)(e,:);
        pixelArray = detrend(filtfilt(B,A,pixelArray - pixelArray(1)) + pixelArray(1),'constant');
        [xcorrVals,lags] = xcorr(pixelArray,gammaBandArray,maxLag,'coeff');
        validVals = xcorrVals(181:181 + 30);
        maxCorr = min(xcorrVals);
        if isnan(maxCorr) == true
            corrMatrix.(hemisphere)(1,e) = 0;
        else
            corrMatrix.(hemisphere)(1,e) = maxCorr;
        end
    end
end
% place circle along the most correlation region of each hemisphere
circRadius = 15;   % pixels to be 1 mm in diameter
for f = 1:length(hem)
    rectMask = ROIs.([hem{1,f} '_windowMask']).rect;
    rectMask = round(rectMask);
    imgWidth = rectMask(3) + 1;
    imgHeight = rectMask(4) + 1;
    corrImg = reshape(corrMatrix.(hem{1,f}),imgHeight,imgWidth);
    % generate image
    isok = false;
    while isok == false
        windowFig = figure;
        imagesc(corrImg)
        title([animalID ' ' hem{1,f} ' peak pixel correlations'])
        xlabel('Image size (pixels)')
        ylabel('Image size (pixels)')
        colormap parula
        colorbar
        axis image
        disp(['Move the ROI over the most correlated region for the ' hem{1,f}]); disp(' ')
        circ = drawcircle('Center',[0,0],'Radius',circRadius);
        checkCircle = input('Is the ROI okay? (y/n): ','s'); disp(' ')
        circPosition = round(circ.Center);
        if strcmp(checkCircle,'y') == true
            isok = true;
            rectBottomLeftCorner = [rectMask(1),rectMask(2) + rectMask(4)];
            circPositionEdit = [rectBottomLeftCorner(1) + circPosition(1),rectBottomLeftCorner(2) - circPosition(2)];
            ROIs.([hem{1,f} '_windowMask']).circPosition = circPositionEdit;
            ROIs.([hem{1,f} '_windowMask']).circRadius = circRadius;
        end
        delete(windowFig);
    end
end
% check final image
finalROIplacement = figure;
imagesc(frames{1})
hold on;
LH_ROI = drawcircle('Center',ROIs.LH_windowMask.circPosition,'Radius',ROIs.LH_windowMask.circRadius);
RH_ROI = drawcircle('Center',ROIs.RH_windowMask.circPosition,'Radius',ROIs.RH_windowMask.circRadius);
title([animalID ' final ROI placement'])
xlabel('Image size (pixels)')
ylabel('Image size (pixels)')
colormap gray
colorbar
axis image
caxis([0 2^RawData.notes.CBVCamBitDepth])
save(ROIs,'ROIs')

end
