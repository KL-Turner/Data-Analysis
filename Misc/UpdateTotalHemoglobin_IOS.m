function [] = UpdateTotalHemoglobin_IOS(procDataFileIDs, RestingBaselines, baselineType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised: July 27th, 2019
%________________________________________________________________________________________________________________________

ledType = 'M530L3';
bandfilterType = 'FB530-10';
cutfilterType = 'EO46540';
conv2um = 1e6;
[~,~,weightedcoeffHbT] = getHbcoeffs_IOS(ledType, bandfilterType, cutfilterType);

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Adding changes in total hemoglobin to ProcData file (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')...']); disp(' ')
    load(procDataFileID)
    [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    cbvFields = fieldnames(ProcData.data.CBV);
    for b = 1:length(cbvFields)
        cbvField = cbvFields{b,1};
        if strcmp(cbvField, 'Cement') == false
            ProcData.data.CBV_HbT.(cbvField) = (log(ProcData.data.CBV.(cbvField)/RestingBaselines.(baselineType).CBV.(cbvField).(strDay)))*weightedcoeffHbT*conv2um;
        end
    end
    save(procDataFileID, 'ProcData')
end


end
