function [] = CorrectGCaMPattenuation_IOS(procDataFileIDs,RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Converts reflectance values to changes in total hemoglobin using absorbance curves of hardware
%________________________________________________________________________________________________________________________

for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID)
    [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    LH_scale = RestingBaselines.manualSelection.CBV.adjLH.(strDay)/RestingBaselines.manualSelection.GCaMP7s.adjLH.(strDay);
    RH_scale = RestingBaselines.manualSelection.CBV.adjRH.(strDay)/RestingBaselines.manualSelection.GCaMP7s.adjRH.(strDay);
    ProcData.data.GCaMP7s.corLH = (ProcData.data.GCaMP7s.adjLH./ProcData.data.CBV.adjLH)*LH_scale;
    ProcData.data.GCaMP7s.corRH = (ProcData.data.GCaMP7s.adjRH./ProcData.data.CBV.adjRH)*RH_scale;
    save(procDataFileID,'ProcData')
end

end