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
    load(procDataFileID)
    [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    % correct attentuation of somatosensory ROIs
    LH_scale = RestingBaselines.manualSelection.CBV.LH.(strDay).mean/RestingBaselines.manualSelection.GCaMP7s.LH.(strDay).mean;
    RH_scale = RestingBaselines.manualSelection.CBV.RH.(strDay).mean/RestingBaselines.manualSelection.GCaMP7s.RH.(strDay).mean;
    ProcData.data.GCaMP7s.corLH = (ProcData.data.GCaMP7s.LH./ProcData.data.CBV.LH)*LH_scale;
    ProcData.data.GCaMP7s.corRH = (ProcData.data.GCaMP7s.RH./ProcData.data.CBV.RH)*RH_scale;
    % correct attentuation of frontal ROIs
    frontalLH_scale = RestingBaselines.manualSelection.CBV.frontalLH.(strDay).mean/RestingBaselines.manualSelection.GCaMP7s.frontalLH.(strDay).mean;
    frontalRH_scale = RestingBaselines.manualSelection.CBV.frontalRH.(strDay).mean/RestingBaselines.manualSelection.GCaMP7s.frontalRH.(strDay).mean;
    ProcData.data.GCaMP7s.corFrontalLH = (ProcData.data.GCaMP7s.frontalLH./ProcData.data.CBV.frontalLH)*frontalLH_scale;
    ProcData.data.GCaMP7s.corFrontalRH = (ProcData.data.GCaMP7s.frontalRH./ProcData.data.CBV.frontalRH)*frontalRH_scale;
    % save data
    save(procDataFileID,'ProcData')
end

end