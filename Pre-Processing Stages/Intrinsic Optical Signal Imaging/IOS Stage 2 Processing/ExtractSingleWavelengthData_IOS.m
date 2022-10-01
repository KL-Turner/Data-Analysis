function [] = ExtractSingleWavelengthData_IOS(ROIs,ROInames,rawDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%
% Purpose: Determine if each desired ROI is drawn, then go through each frame and extract the mean of valid pixels.
%________________________________________________________________________________________________________________________

for a = 1:size(rawDataFileIDs,1)
    rawDataFile = rawDataFileIDs(a,:);
    disp(['Analyzing IOS ROIs from RawData file (' num2str(a) '/' num2str(size(rawDataFileIDs,1)) ')']); disp(' ')
    [animalID,fileDate,fileID] = GetFileInfo_IOS(rawDataFile);
    strDay = ConvertDate_IOS(fileDate);
    load(rawDataFile)
    [frames] = ReadDalsaBinary_IOS(animalID,[fileID '_WindowCam.bin']);
    for b = 1:length(ROInames)
        ROIshortName = ROInames{1,b};
        ROIname = [ROInames{1,b} '_' strDay];
        disp(['Extracting ' ROIname ' ROI CBV data from ' rawDataFile '...']); disp(' ')
        % draw circular ROIs based on XCorr for LH/RH/Barrels, then free-hand for cement ROIs
        maskFig = figure;
        imagesc(frames{1});
        axis image;
        colormap gray
        if any(strcmp(ROIshortName,{'LH','RH','frontalLH','frontalRH','Barrels'})) == true
            circROI = drawcircle('Center',ROIs.(ROIname).circPosition,'Radius',ROIs.(ROIname).circRadius);
            mask = createMask(circROI,frames{1});
            close(maskFig)
        else
            mask = roipoly(frames{1},ROIs.(ROIname).xi,ROIs.(ROIname).yi);
            close(maskFig)
        end
        meanIntensity = BinToIntensity_IOS(mask,frames);
        RawData.data.CBV.(ROIname) = meanIntensity;
    end
    save(rawDataFile,'RawData')
end

end
