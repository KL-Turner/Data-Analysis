function [] = CalcVesselDiameterFWHM_IOS(procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________   

ROIFileDir = dir('*_ROIs.mat');
ROIFileName = {ROIFileDir.name}';
ROIFileID = char(ROIFileName);
load(ROIFileID);
for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    [animalID,~,fileID] = GetFileInfo_IOS(procDataFileID);
    windowCamFileID = [fileID '_WindowCam.bin'];
    load(procDataFileID)
    newFileID = [animalID '_' fileID '_VesselDiam.mat'];
    imageHeight = ProcData.notes.CBVCamPixelHeight;                                                                                                            
    imageWidth = ProcData.notes.CBVCamPixelWidth;
    samplingRate = ProcData.notes.CBVCamSamplingRate;
    trialDuration_sec = ProcData.notes.trialDuration_sec;
    frameInds = 1:samplingRate*trialDuration_sec;
    [imageStack] =  GetCBVFrameSubset_IOS(windowCamFileID,imageWidth,imageHeight,frameInds);
    for c = 1:length(ROIs.vesselIDs)
        vesselID = ROIs.vesselIDs{c,1};
        for d = 1:size(imageStack,3)
            currentImage = imageStack(:,:,d);
            dataLine = improfile(currentImage,ROIs.x_endpoints{c,1},ROIs.y_endpoints{c,1});
            % interpolate data_line
            dataLine = interp1(1:length(dataLine),dataLine,1:1/100:length(dataLine),'spline');
            % Find the half max value.
            halfMax = (min(dataLine) + max(dataLine))/2;
            % Find where the data first drops below half the max.
            index1 = find(dataLine <= halfMax,1,'first');
            % Find where the data last rises above half the max.
            index2 = find(dataLine <= halfMax,1,'last');
            VesselData.(vesselID)(d,1) = index2 - index1 + 1; % FWHM in indexes.
        end
    end
    disp(['Adding vessel FWHM analysis to ProcData file ' num2str(a) ' of ' num2str(size(procDataFileIDs,1))]); disp(' ')
    save(newFileID,'VesselData')
end

end