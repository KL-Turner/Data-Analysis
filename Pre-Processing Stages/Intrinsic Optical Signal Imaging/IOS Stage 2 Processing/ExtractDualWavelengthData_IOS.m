function [] = ExtractDualWavelengthData_IOS(ROIs,ROInames,rawDataFileIDs,procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Determine if each desired ROI is drawn, then go through each frame and extract the mean of valid pixels.
%________________________________________________________________________________________________________________________

for a = 1:size(rawDataFileIDs,1)
    rawDataFileID = rawDataFileIDs(a,:);
    procDataFileID = procDataFileIDs(a,:);
    disp(['Analyzing IOS ROIs from ProcData file (' num2str(a) '/' num2str(size(rawDataFileIDs,1)) ')']); disp(' ')
    [animalID,fileDate,fileID] = GetFileInfo_IOS(rawDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    load(rawDataFileID)
    load(procDataFileID)
    [frames] = ReadDalsaBinary_IOS(animalID,[fileID '_WindowCam.bin']);
    if ProcData.notes.greenFrames == 1
        cbvFrames = frames(1:2:end - 1);
        gcampFrames = frames(2:2:end);
    elseif ProcData.notes.greenFrames == 2
        cbvFrames = frames(2:2:end);
        gcampFrames = frames(1:2:end);
    end
    for b = 1:length(ROInames)
        ROIshortName = ROInames{1,b};
        ROIname = [ROInames{1,b} '_' strDay];
        disp(['Extracting ' ROIname ' ROI CBV data from ' rawDataFileID '...']); disp(' ')
        % draw circular ROIs based on XCorr for LH/RH/Barrels, then free-hand for cement ROIs
        maskFig = figure;
        imagesc(frames{1});
        axis image;
        colormap gray
        if any(strcmp(ROIshortName,{'LH','RH','frontalLH','frontalRH','Barrels'})) == true
            circROI = drawcircle('Center',ROIs.(ROIname).circPosition,'Radius',ROIs.(ROIname).circRadius);
            mask = createMask(circROI,frames{1});
            close(maskFig)
            cbvMeanIntensity = BinToIntensity_IOS(mask,cbvFrames);
            gcampMeanIntensity = BinToIntensity_IOS(mask,gcampFrames);
            RawData.data.CBV.(ROIname) = cbvMeanIntensity;
            RawData.data.GCaMP7s.(ROIname) = gcampMeanIntensity;
        elseif strcmp(ROIshortName,'Cement') == true
            mask = roipoly(frames{1},ROIs.(ROIname).xi,ROIs.(ROIname).yi);
            close(maskFig)
            cbvCementMeanIntensity = BinToIntensity_IOS(mask,cbvFrames);
            gcampCementMeanIntensity = BinToIntensity_IOS(mask,gcampFrames);
            RawData.data.CBV.(ROIname) = cbvCementMeanIntensity;
            RawData.data.GCaMP7s.(ROIname) = gcampCementMeanIntensity;
        end
    end
    save(rawDataFileID,'RawData')
end

end
