function [] = UpdateTotalHemoglobin_IOS_eLife2020(procDataFileIDs,RestingBaselines,baselineType,imagingType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Converts reflectance values to changes in total hemoglobin using absorbance curves of hardware
%________________________________________________________________________________________________________________________

ledType = 'M530L3';
bandfilterType = 'FB530-10';
cutfilterType = 'EO46540';
conv2um = 1e6;
[~,~,weightedcoeffHbT] = getHbcoeffs_IOS_eLife2020(ledType,bandfilterType,cutfilterType);
for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Adding changes in total hemoglobin to ProcData file (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')...']); disp(' ')
    load(procDataFileID)
    [~,fileDate,~] = GetFileInfo_IOS_eLife2020(procDataFileID);
    strDay = ConvertDate_IOS_eLife2020(fileDate);
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
