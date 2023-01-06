function [] = UpdateTotalHemoglobin_IOS(procDataFileIDs,RestingBaselines,baselineType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Converts reflectance values to changes in total hemoglobin using absorbance curves of hardware
%________________________________________________________________________________________________________________________

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
    [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    cbvFields = fieldnames(ProcData.data.CBV);
    for bb = 1:length(cbvFields)
        cbvField = cbvFields{bb,1};
        ProcData.data.CBV_HbT.(cbvField) = (log(ProcData.data.CBV.(cbvField)/RestingBaselines.(baselineType).CBV.(cbvField).(strDay).mean))*weightedcoeffHbT*conv2um;
    end
    save(procDataFileID,'ProcData')
end

end
