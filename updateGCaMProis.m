zap;
curDir = pwd;
% create/load pre-existing ROI file with the coordinates
cd ..
ROIFileDir = dir('*_UpdatedROIs.mat');
if isempty(ROIFileDir) == true
    UpdatedROIs = [];
else
    ROIFileName = {ROIFileDir.name}';
    ROIFileID = char(ROIFileName);
    load(ROIFileID);
end
cd(curDir)
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% load existing ROI structure if it exists
ROIFileDir = dir('*_ROIs.mat');
ROIFileName = {ROIFileDir.name}';
ROIFileID = char(ROIFileName);
load(ROIFileID);
qq = 1; stimFrames = [];
for zz = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(zz,:);
    load(procDataFileID)
    [animalID,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    windowCamFileID = [fileID '_WindowCam.bin'];
    % radius area
    circRadius = 15;
    circRadius = circRadius/2;
    % frames
    [frames] = ReadDalsaBinary_IOS(animalID,windowCamFileID);
    if ProcData.notes.greenFrames == 1
        cbvFrames = frames(1:3:end - 1);
        gcampFrames = frames(2:3:end);
        deoxyFrames = frames(3:3:end);
    elseif ProcData.notes.greenFrames == 2
        cbvFrames = frames(2:3:end);
        gcampFrames = frames(3:3:end);
        deoxyFrames = frames(1:3:end - 1);
    elseif ProcData.notes.greenFrames == 3
        cbvFrames = frames(3:3:end);
        gcampFrames = frames(1:3:end - 1);
        deoxyFrames = frames(2:3:end);
    end
    stimTimes = [ProcData.data.stimulations.LPadSol,ProcData.data.stimulations.RPadSol];
    % stimTimes = [100 200 300 400 500 600 700 800];
    if isempty(stimTimes) == false
        for aa = 1:length(stimTimes)
            stimTime = stimTimes(1,aa);
            frameTime = round(stimTime*10);
            baseline = [];
            for bb = 1:20
                baseline(:,:,bb) = cbvFrames{1,frameTime - bb};
            end
            baseFrame = mean(baseline,3);
            stimFrames(:,:,qq) = (cbvFrames{1,frameTime + 10} - baseFrame)./baseFrame;
            qq = qq + 1;
        end
    end
end
avgFrame = mean(stimFrames,3);
% place circle along the most relevant region of each hemisphere
ROInames = {'LH','RH'};
for ff = 1:length(ROInames)
    % generate image
    isok = false;
    while isok == false
        windowFig = figure;
        ROInames2 = {'LH','RH','frontalLH','frontalRH'};
        imagesc(avgFrame)
        hold on;
        for aa = 1:length(ROInames2)
            drawcircle('Center',ROIs.([ROInames2{1,aa} '_' strDay]).circPosition,'Radius',ROIs.([ROInames2{1,aa} '_' strDay]).circRadius,'Color','r');
        end
        title([animalID ' original ROI placement'])
        xlabel('Image size (pixels)')
        ylabel('Image size (pixels)')
        colormap jet
        colorbar
        axis image
        caxis([-.05,0.05])
        disp(['Move the ROI over the desired region for ' ROInames{1,ff}]); disp(' ')
        drawnow
        circ = drawcircle('Center',[0,0],'Radius',circRadius,'Color','r');
        checkCircle = input('Is the ROI okay? (y/n): ','s'); disp(' ')
        circPosition = round(circ.Center);
        if strcmpi(checkCircle,'y') == true
            isok = true;
            UpdatedROIs.(strDay).(ROInames{1,ff}).circPosition = circPosition;
            UpdatedROIs.(strDay).(ROInames{1,ff}).circRadius = circRadius;
        end
        delete(windowFig);
    end
end
UpdatedROIs.(strDay).fLH.circPosition = ROIs.(['frontalLH_' strDay]).circPosition;
UpdatedROIs.(strDay).fLH.circRadius = ROIs.(['frontalLH_' strDay]).circRadius;
UpdatedROIs.(strDay).fRH.circPosition = ROIs.(['frontalRH_' strDay]).circPosition;
UpdatedROIs.(strDay).fRH.circRadius = ROIs.(['frontalRH_' strDay]).circRadius;
% check final image
ROInames = {'LH','RH','fLH','fRH'};
fig = figure;
imagesc(avgFrame)
hold on;
for aa = 1:length(ROInames)
    drawcircle('Center',UpdatedROIs.(strDay).(ROInames{1,aa}).circPosition,'Radius',UpdatedROIs.(strDay).(ROInames{1,aa}).circRadius,'Color','b');
end
title([animalID ' updated ROI placement'])
xlabel('Image size (pixels)')
ylabel('Image size (pixels)')
colormap jet
colorbar
axis image
caxis([-.05,0.05])
savefig(fig,[animalID '_' strDay '_UpdatedROIs.fig'])
cd ..
save([animalID '_UpdatedROIs.mat'],'UpdatedROIs');
cd(curDir)