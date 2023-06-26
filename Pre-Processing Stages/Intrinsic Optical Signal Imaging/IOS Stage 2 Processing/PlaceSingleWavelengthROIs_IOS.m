function [ROIs] = PlaceSingleWavelengthROIs_IOS(animalID,fileID,ROIs,lensMag,imagingType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Place a circular 1 mm ROI over multi-wavelength IOS images
%________________________________________________________________________________________________________________________
strDay = ConvertDate_IOS(fileID);
fileDate = fileID(1:6);
% determine which ROIs to draw based on imaging type
if strcmpi(imagingType,'Single ROI (SI)') == true
    ROInames = {'barrels'};
elseif strcmpi(imagingType,'Single ROI (SSS)') == true
    ROInames = {'SSS','lSSS','rSSS'};
elseif strcmpi(imagingType,'Bilateral ROI (SI)') == true
    ROInames = {'LH','RH'};
elseif strcmpi(imagingType,'Bilateral ROI (SI,FC)') == true
    ROInames = {'LH','RH','fLH','fRH'};
end
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% character list of all WindowCam files
windowDataFileStruct = dir('*_PCO_Cam01.pcoraw');
windowDataFiles = {windowDataFileStruct.name}';
windowDataFileIDs = char(windowDataFiles);
bb = 1;
% double check file lists
for aa = 1:size(procDataFileIDs)
    procDataFileID = procDataFileIDs(aa,:);
    [~,procfileDate,~] = GetFileInfo_IOS(procDataFileID);
    windowDataFileID = windowDataFileIDs(aa,:);
    if strcmp(procfileDate,fileDate) == true
        procDataFileList(bb,:) = procDataFileID;
        windowDataFileList(bb,:) = windowDataFileID;
        bb = bb + 1;
    end
end
% go through each file and check the color o
load(procDataFileList(1,:));
nFramesToRead = 10;
% pre-allocate memory
frames = cell(1,nFramesToRead);
for n = 1:nFramesToRead
    frames{n} = imread(windowDataFileList(1,:),n);
end
roiFrame = frames{1};
% determine the proper size of the ROI based on camera/lens magnification
if strcmpi(lensMag,'0.75X') == true
    circRadius = 37; % pixels to be 1 mm in diameter
elseif strcmpi(lensMag,'1.0X') == true
    circRadius = 45;
elseif strcmpi(lensMag,'1.5X') == true
    circRadius = 60;
elseif strcmpi(lensMag,'2.0X') == true
    circRadius = 75;
elseif strcmpi(lensMag,'2.5X') == true
    circRadius = 90;
elseif strcmpi(lensMag,'3.0X') == true
    circRadius = 105;
end
% place circle along the most relevant region of each hemisphere
for ff = 1:length(ROInames)
    % generate image
    isok = false;
    while isok == false
        windowFig = figure;
        imagesc(roiFrame)
        title([animalID ' ' ROInames{1,ff} ' ROI'])
        xlabel('Image size (pixels)')
        ylabel('Image size (pixels)')
        colormap gray
        colorbar
        axis image
        disp(['Move the ROI over the desired region for ' ROInames{1,ff}]); disp(' ')
        drawnow
        circ = drawcircle('Center',[0,0],'Radius',circRadius,'Color','r');
        checkCircle = input('Is the ROI okay? (y/n): ','s'); disp(' ')
        circPosition = round(circ.Center);
        if strcmpi(checkCircle,'y') == true
            isok = true;
            ROIs.(strDay).(ROInames{1,ff}).circPosition = circPosition;
            ROIs.(strDay).(ROInames{1,ff}).circRadius = circRadius;
        end
        delete(windowFig);
    end
end
% check final image
fig = figure;
imagesc(roiFrame)
hold on;
for aa = 1:length(ROInames)
    drawcircle('Center',ROIs.([ROInames{1,aa} '_' strDay]).circPosition,'Radius',ROIs.([ROInames{1,aa} '_' strDay]).circRadius,'Color','r');
end
title([animalID ' final ROI placement'])
xlabel('Image size (pixels)')
ylabel('Image size (pixels)')
colormap gray
colorbar
axis image
clim([0,2^ProcData.notes.CBVCamBitDepth])
savefig(fig,[animalID '_' strDay '_ROIs.fig'])