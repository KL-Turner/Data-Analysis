function [imageMatrix] = ExtractImageMatrixFor2PData_IOS(fileID)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________

% pull file information and camera frames
fileID2 = strrep(fileID,'_',' ');
rawDataFileID = ls(['*' fileID '_RawData.mat']);
load(rawDataFileID)
windowCamFileID = [fileID '_WindowCam.bin'];
[~,fileDate,~] = GetFileInfo_IOS(rawDataFileID);
strDay = ConvertDate_IOS(fileDate);
animalID = RawData.notes.animalID;
[frames] = ReadDalsaBinary_IOS(animalID,windowCamFileID);
% ROI file
ROIFileDir = dir('*_ROIs.mat');
if isempty(ROIFileDir) == true
    ROIs = [];
else
    ROIFileName = {ROIFileDir.name}';
    ROIFileID = char(ROIFileName);
    load(ROIFileID);
end
% draw ROI if it doesn't exist
if isempty(ROIFileDir) == true
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
    % draw rectangular ROI around the image
    isok = false;
    while isok == false
        disp('Draw an ROI over the entire window'); disp(' ')
        [~,rect] = imcrop;
        hold on;
        ROIoutline = rectangle('Position',rect,'EdgeColor','r');
        checkMask = input('Is the ROI okay? (y/n): ','s'); disp(' ')
        if strcmp(checkMask,'y') == true
            isok = true;
            ROIs.(['IOS_' strDay]).rect = rect;
        end
        delete(ROIoutline);
    end
    close(windowFig)
    save([animalID '_ROIs.mat'],'ROIs')
end
% extract the pixel values from the window ROIs
imageMask = nan(size(frames{1}));
rectMask = ROIs.(['IOS_' strDay]).rect;
rectMask = round(rectMask);
imageWidth = rectMask(3) + 1;
imageHeight = rectMask(4) + 1;
imageMask(rectMask(2):(rectMask(2) + rectMask(4)),rectMask(1):(rectMask(1) + rectMask(3))) = 1;
imageMatrix = zeros(imageHeight,imageWidth,length(frames) - 1);
for c = 1:length(frames) - 1
    disp(['Reshaping image matrix frame (' num2str(c) '/' num2str(length(frames) - 1) ')']); disp(' ')
    frame = frames{1,c};
    frameHold = double(frame).*imageMask;
    imageArray = frameHold(~isnan(frameHold));
    RawData.data.imageMatrix(:,:,c) = reshape(imageArray,imageHeight,imageWidth);
end
disp('Saving RawData file...'); disp(' ')
save(rawDataFileID,'RawData','-v7.3')

end
