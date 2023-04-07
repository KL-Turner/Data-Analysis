function [] = UpdateTotalHemoglobin_IOS(procDataFileIDs,RestingBaselines,baselineType,pixelWise)
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
    disp(['Adding changes in total hemoglobin to ProcData file (' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ')...']); disp(' ')
    load(procDataFileID)
    imagingWavelengths = ProcData.notes.imagingWavelengths;
    if any(strcmp(imagingWavelengths,{'Green','Green & Blue','Red, Green, & Blue'})) == true
        ledType = 'M530L3';
        bandfilterType = 'FB530-10';
        cutfilterType = 'MF525-39';
        [~,~,weightedcoeffHbT] = GetHbcoeffs_IOS(ledType,bandfilterType,cutfilterType);
    elseif any(strcmp(imagingWavelengths,{'Lime','Lime & Blue','Red, Lime, & Blue'})) == true
        ledType = 'M565L3';
        bandfilterType = 'FB570-10';
        cutfilterType = 'EO65160';
        [~,~,weightedcoeffHbT] = GetHbcoeffs_IOS(ledType,bandfilterType,cutfilterType);
    else
        weightedcoeffHbT = NaN;
    end
    conv2um = 1e6;
    [animalID,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    if strcmp(pixelWise,'y') == true
        windowCamFileID = [fileID '_WindowCam.bin'];
        [frames] = ReadDalsaBinary_IOS(animalID,windowCamFileID);
        if ProcData.notes.greenFrames == 1
            greenFrames = frames(1:3:end - 1);
        elseif ProcData.notes.greenFrames == 2
            greenFrames = frames(2:3:end);
        elseif ProcData.notes.greenFrames == 3
            greenFrames = frames(3:3:end);
        end
        greenImageStack = reshape(cell2mat(greenFrames),ProcData.notes.CBVCamPixelWidth,ProcData.notes.CBVCamPixelHeight,length(greenFrames));
        hbtImageStack = (log(greenImageStack./RestingBaselines.Pixel.green.(strDay))).*weightedcoeffHbT*conv2um;
        cbvFields = fieldnames(UpdatedROIs.(strDay));
        for bb = 1:length(cbvFields)
            cbvField = cbvFields{bb,1};
            maskFig = figure;
            imagesc(greenFrames{1});
            axis image;
            colormap gray
            circROI = drawcircle('Center',UpdatedROIs.(strDay).(cbvField).circPosition,'Radius',UpdatedROIs.(strDay).(cbvField).circRadius);
            mask = createMask(circROI,greenFrames{1});
            close(maskFig)
            for nn = 1:size(hbtImageStack,3)
                roiMask = mask.*double(hbtImageStack(:,:,nn));
                ProcData.data.HbT.(cbvField)(1,nn) = mean(nonzeros(roiMask));
            end
        end
    else
        cbvFields = fieldnames(ProcData.data.CBV);
        for bb = 1:length(cbvFields)
            cbvField = cbvFields{bb,1};
            ProcData.data.HbT.(cbvField) = (log(ProcData.data.CBV.(cbvField)/RestingBaselines.(baselineType).CBV.(cbvField).(strDay).mean))*weightedcoeffHbT*conv2um;
        end
    end
    save(procDataFileID,'ProcData')
end

end
