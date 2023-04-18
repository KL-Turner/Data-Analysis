function [ROIs] = PlaceTriWavelengthROIs_IOS(animalID,fileID,ROIs,lensMag,imagingType)
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
for qq = 1:size(procDataFileList,1)
    disp(['Verifying first frame color from file (' num2str(qq) '/' num2str(size(procDataFileList,1)) ')']); disp(' ')
    load(procDataFileList(qq,:));
    nFramesToRead = 10;
    % pre-allocate memory
    frames = cell(1,nFramesToRead);
    for n = 1:nFramesToRead
        frames{n} = imread(windowDataFileList(qq,:),n);
    end
    gcampCheck = figure;
    frames = frames(1:end);
    for xx = 1:10
        subplot(2,5,xx)
        imagesc(frames{1,xx})
        axis image
        colormap gray
    end
    contCheck = false;
    while contCheck == false
        drawnow
        if isfield(ProcData.notes,'blueFrames') == true
            gcampFrames = ProcData.notes.blueFrames;
        else
            gcampFrames = input('Which index are blue LED frames (1,2,3): '); disp(' ')
        end
        if gcampFrames == 1
            if qq == 1
                roiFrame = frames{3};
            end
            ProcData.notes.blueFrames = 1;
            ProcData.notes.redFrames = 2;
            ProcData.notes.greenFrames = 3;
            save(procDataFileIDs(qq,:),'ProcData')
            contCheck = true;
        elseif gcampFrames == 2
            if qq == 1
                roiFrame = frames{1};
            end
            ProcData.notes.blueFrames = 2;
            ProcData.notes.redFrames = 3;
            ProcData.notes.greenFrames = 1;
            save(procDataFileIDs(qq,:),'ProcData')
            contCheck = true;
        elseif gcampFrames == 3
            if qq == 1
                roiFrame = frames{2};
            end
            ProcData.notes.blueFrames = 3;
            ProcData.notes.redFrames = 1;
            ProcData.notes.greenFrames = 2;
            save(procDataFileIDs(qq,:),'ProcData')
            contCheck = true;
        end
    end
    close(gcampCheck)
    fclose('all');
end
% determine the proper size of the ROI based on camera/lens magnification
if strcmpi(lensMag,'0.75X') == true
    circRadius = 37/2; % pixels to be 1 mm in diameter
elseif strcmpi(lensMag,'1.0X') == true
    circRadius = 45/2;
elseif strcmpi(lensMag,'1.5X') == true
    circRadius = 60/2;
elseif strcmpi(lensMag,'2.0X') == true
    circRadius = 75/2;
elseif strcmpi(lensMag,'2.5X') == true
    circRadius = 90/2;
elseif strcmpi(lensMag,'3.0X') == true
    circRadius = 105/2;
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
    drawcircle('Center',ROIs.(strDay).(ROInames{1,aa}).circPosition,'Radius',ROIs.(strDay).(ROInames{1,aa}).circRadius,'Color','r');
end
title([animalID ' final ROI placement'])
xlabel('Image size (pixels)')
ylabel('Image size (pixels)')
colormap gray
colorbar
axis image
caxis([0,2^ProcData.notes.CBVCamBitDepth])
savefig(fig,[animalID '_' strDay '_ROIs.fig'])