function [] = UpdateTotalHemoglobin_FP(procDataFileIDs,RestingBaselines,baselineType,imagingType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Converts reflectance values to changes in total hemoglobin using absorbance curves of hardware
%________________________________________________________________________________________________________________________

ledType = 'M565L3';
bandfilterType = 'FB570-10';
cutfilterType = 'EO65160';
conv2um = 1e6;
[~,~,weightedcoeffHbT] = getHbcoeffs_FP(ledType,bandfilterType,cutfilterType);
for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Adding changes in total hemoglobin to ProcData file (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')...']); disp(' ')
    load(procDataFileID)
    [~,fileDate,~] = GetFileInfo_FP(procDataFileID);
    strDay = ConvertDate_FP(fileDate);
    if strcmp(imagingType,'bilateral') == true
            cbvFields = {'LH','adjLH','RH','adjRH'};
    elseif strcmp(imagingType,'single') == true
            cbvFields = {'Barrels','adjBarrels'};
    end
    for b = 1:length(cbvFields)
        cbvField = cbvFields{1,b};
        ProcData.data.CBV_HbT.(cbvField) = (log(ProcData.data.CBV.(cbvField)/RestingBaselines.(baselineType).CBV.(cbvField).(strDay)))*weightedcoeffHbT*conv2um;
    end
    save(procDataFileID,'ProcData')
end

end
