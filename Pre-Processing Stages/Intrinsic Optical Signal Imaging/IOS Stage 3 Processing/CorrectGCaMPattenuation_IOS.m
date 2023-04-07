function [] = CorrectGCaMPattenuation_IOS(procDataFileIDs,RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Converts reflectance values to changes in total hemoglobin using absorbance curves of hardware
%________________________________________________________________________________________________________________________

ROIFileDir = dir('*_UpdatedROIs.mat');
ROIFileName = {ROIFileDir.name}';
ROIFileID = char(ROIFileName);
load(ROIFileID);
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    disp(['Adding GCaMP correction to ProcData file (' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ')...']); disp(' ')
    load(procDataFileID)
    imagingWavelengths = ProcData.notes.imagingWavelengths;
    if any(strcmp(imagingWavelengths,{'Red, Green, & Blue','Lime, Green, & Blue','Green & Blue','Lime & Blue'})) == true
        [animalID,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(fileDate);
        windowCamFileID = [fileID '_WindowCam.bin'];
        [frames] = ReadDalsaBinary_IOS(animalID,windowCamFileID);
        if ProcData.notes.greenFrames == 1
            greenFrames = frames(1:3:end - 1);
            blueFrames = frames(2:3:end);
        elseif ProcData.notes.greenFrames == 2
            greenFrames = frames(2:3:end);
            blueFrames = frames(3:3:end);
        elseif ProcData.notes.greenFrames == 3
            greenFrames = frames(3:3:end);
            blueFrames = frames(1:3:end - 1);
        end
        greenImageStack = reshape(cell2mat(greenFrames),ProcData.notes.CBVCamPixelWidth,ProcData.notes.CBVCamPixelHeight,length(greenFrames));
        blueImageStack = reshape(cell2mat(blueFrames),ProcData.notes.CBVCamPixelWidth,ProcData.notes.CBVCamPixelHeight,length(blueFrames));
        scale = RestingBaselines.Pixel.green.(strDay)./RestingBaselines.Pixel.blue.(strDay);
        correctedStack = (blueImageStack./greenImageStack).*scale;
        % correct attentuation of somatosensory ROIs
        ROInames = fieldnames(UpdatedROIs.(strDay));
        for bb = 1:length(ROInames)
            roiName = ROInames{bb,1};
            maskFig = figure;
            imagesc(greenFrames{1});
            axis image;
            colormap gray
            circROI = drawcircle('Center',UpdatedROIs.(strDay).(roiName).circPosition,'Radius',UpdatedROIs.(strDay).(roiName).circRadius);
            mask = createMask(circROI,greenFrames{1});
            close(maskFig)
            for nn = 1:size(correctedStack,3)
                roiMask = mask.*double(correctedStack(:,:,nn));
                ProcData.data.GCaMP7s.(roiName)(1,nn) = mean(nonzeros(roiMask));
            end
        end
    end
    % save data
    save(procDataFileID,'ProcData')
end

end