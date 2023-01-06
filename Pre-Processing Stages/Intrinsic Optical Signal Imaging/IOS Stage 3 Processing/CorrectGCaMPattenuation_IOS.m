function [] = CorrectGCaMPattenuation_IOS(procDataFileIDs,RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Converts reflectance values to changes in total hemoglobin using absorbance curves of hardware
%________________________________________________________________________________________________________________________

for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    disp(['Adding GCaMP correction to ProcData file (' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ')...']); disp(' ')
    load(procDataFileID)
    imagingType = ProcData.notes.imagingType;
    imagingWavelengths = ProcData.notes.imagingWavelengths;
    if any(strcmp(imagingWavelengths,{'Red, Green, & Blue','Lime, Green, & Blue','Green & Blue','Lime & Blue'})) == true
        [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(fileDate);
        % correct attentuation of somatosensory ROIs
        % use imaging type to determine ROI names and typical lens magnification
        if strcmpi(imagingType,'Single ROI (SI)') == true
            ROInames = {'Barrels'};
        elseif strcmpi(imagingType,'Single ROI (SSS)') == true
            ROInames = {'SSS','lSSS','rSSS'};
        elseif strcmpi(imagingType,'Bilateral ROI (SI)') == true
            ROInames = {'LH','RH'};
        elseif strcmpi(imagingType,'Bilateral ROI (SI,FC)') == true
            ROInames = {'LH','RH','frontalLH','frontalRH'};
        end
        for bb = 1:length(ROInames)
            gcampField = ROInames{1,bb};
            scale = RestingBaselines.manualSelection.CBV.(gcampField).(strDay).mean/RestingBaselines.manualSelection.GCaMP7s.(gcampField).(strDay).mean;
            ProcData.data.GCaMP7s.(['cor' (gcampField)]) = (ProcData.data.GCaMP7s.(gcampField)./ProcData.data.CBV.(gcampField))*scale;
        end
    end
    % save data
    save(procDataFileID,'ProcData')
end

end