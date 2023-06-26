function [] = ExtractTriWavelengthData_IOS(ROIs,ROInames,procDataFileIDs)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Analyzing IOS ROIs from ProcData file (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
    [~,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    load(procDataFileID)
    if isfield(ProcData.data,'reflectance') == false
        info = imfinfo([fileID '_PCO_Cam01.pcoraw']);
        numberOfPages = length(info);
        for k = 1:numberOfPages
            disp(num2str(k))
            imageStack(:,:,k) = imread([fileID '_PCO_Cam01.pcoraw'],k);
        end
        % separate image stack by wavelength
        if any(strcmp(imagingWavelengths,{'Red, Green, & Blue','Red, Lime, & Blue'})) == true
            reflFields = {'red','green','blue'};
            red = imageStack(:,:,1:end);
            green = imageStack(:,:,1:end);
            blue = imageStack(:,:,1:end);
        elseif any(strcmp(imagingWavelengths,{'Green & Blue','Lime & Blue'})) == true
            reflFields = {'green','blue'};
            green = imageStack(:,:,1:end);
            blue = imageStack(:,:,1:end);
        else
            reflFields = {'green'};
            green = imageStack(:,:,1:end);
        end
        for b = 1:length(ROInames)
            ROIshortName = ROInames{1,b};
            ROIname = ROInames{1,b};
            disp(['Extracting ' ROIname ' ROI IOS data from ' procDataFileID '...']); disp(' ')
            % draw circular ROIs based on XCorr for LH/RH/Barrels, then free-hand for cement ROIs
            maskFig = figure;
            imagesc(green(:,:,1));
            axis image;
            colormap gray
            for cc = 1:length(reflFields)
                reflField = reflFields{cc,1};
                if any(strcmp(ROIshortName,{'LH','RH','fLH','fRH','barrels'})) == true
                    circROI = drawcircle('Center',ROIs.(strDay).(ROIname).circPosition,'Radius',ROIs.(strDay).(ROIname).circRadius);
                    ROImask = createMask(circROI,reflField(:,:,1));
                    close(maskFig)
                    for aa = 1:size(reflField,3)
                        mask = ROImask.*double(reflField(:,:,aa));
                        refl(aa,1) = mean(nonzeros(mask));
                    end
                    ProcData.data.(reflField).(ROIname) = refl;
                else
                    ROImask = roipoly(double(reflField(:,:,1)),ROIs.(ROIname).xi,ROIs.(ROIname).yi);
                    close(maskFig)
                    for aa = 1:size(reflField,3)
                        mask = ROImask.*double(reflField(:,:,aa));
                        refl(aa,1) = mean(nonzeros(mask));
                    end
                    ProcData.data.(reflField).(ROIname) = refl;
                end
            end
        end
        save(procDataFileID,'ProcData')
    end
end