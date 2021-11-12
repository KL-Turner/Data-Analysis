function [] = zScorePupilData_IOS(procDataFileID,RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Manually check pupil diameter
%________________________________________________________________________________________________________________________

load(procDataFileID);
if isfield(ProcData.data.Pupil,'zDiameter') == false %#ok<NODEF>
    [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
    if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
        strDay = ConvertDate_IOS(fileDate);
        areaMean = RestingBaselines.manualSelection.Pupil.mmArea.(strDay).mean;
        areaStd = RestingBaselines.manualSelection.Pupil.mmArea.(strDay).std;
        diameterMean = RestingBaselines.manualSelection.Pupil.mmDiameter.(strDay).mean;
        diameterStd = RestingBaselines.manualSelection.Pupil.mmDiameter.(strDay).std;
        pupilArea = ProcData.data.Pupil.mmArea;
        diameter = ProcData.data.Pupil.mmDiameter;
        ProcData.data.Pupil.zArea = (pupilArea - areaMean)./areaStd;
        ProcData.data.Pupil.zDiameter = (diameter - diameterMean)./diameterStd;
        save(procDataFileID,'ProcData')
    end
end

end