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
    if isfield(ProcData.data,'CBV') == false
        info = imfinfo([fileID '_PCO_Cam01.pcoraw']);
        numberOfPages = length(info);
        for k = 1:numberOfPages
            disp(num2str(k))
            imageStack(:,:,k) = imread([fileID '_PCO_Cam01.pcoraw'],k);
        end
        if ProcData.notes.greenFrames == 1
            cbvFrames = imageStack(:,:,1:3:end - 1);
        elseif ProcData.notes.greenFrames == 2
            cbvFrames = imageStack(:,:,2:3:end);
        elseif ProcData.notes.greenFrames == 3
            cbvFrames = imageStack(:,:,3:3:end);
        end
        for b = 1:length(ROInames)
            ROIshortName = ROInames{1,b};
            ROIname = ROInames{1,b};
            disp(['Extracting ' ROIname ' ROI CBV data from ' procDataFileID '...']); disp(' ')
            % draw circular ROIs based on XCorr for LH/RH/Barrels, then free-hand for cement ROIs
            maskFig = figure;
            imagesc(cbvFrames(:,:,1));
            axis image;
            colormap gray
            if any(strcmp(ROIshortName,{'LH','RH','fLH','fRH','barrels'})) == true
                circROI = drawcircle('Center',ROIs.(strDay).(ROIname).circPosition,'Radius',ROIs.(strDay).(ROIname).circRadius);
                ROImask = createMask(circROI,cbvFrames(:,:,1));
                close(maskFig)
                for aa = 1:size(cbvFrames,3)
                    mask = ROImask.*double(cbvFrames(:,:,aa));
                    refl(aa,1) = mean(nonzeros(mask));
                end
                ProcData.data.CBV.(ROIname) = refl;
            else
                ROImask = roipoly(double(cbvFrames(:,:,1)),ROIs.(ROIname).xi,ROIs.(ROIname).yi);
                close(maskFig)
                for aa = 1:size(cbvFrames,3)
                    mask = ROImask.*double(cbvFrames(:,:,aa));
                    refl(aa,1) = mean(nonzeros(mask));
                end
                ProcData.data.CBV.(ROIname) = refl;
            end
        end
        save(procDataFileID,'ProcData')
    end
end