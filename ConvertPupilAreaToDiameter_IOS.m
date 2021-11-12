function [] = ConvertPupilAreaToDiameter_IOS(procDataFileID)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Manually check pupil diameter
%________________________________________________________________________________________________________________________

load(procDataFileID);
% if isfield(ProcData.data.Pupil,'pixelDiameter') == false %#ok<NODEF>
    pupilArea = ProcData.data.Pupil.pupilArea; 
    diameter = sqrt(pupilArea./pi)*2;
    ProcData.data.Pupil.diameter = diameter;
%     ProcData.data.Pupil.pupilMinor = ProcData.data.Pupil.pupiMinor;
%     ProcData.data.Pupil = rmfield(ProcData.data.Pupil,'pupiMinor');
    ProcData.data.Pupil.mmPerPixel = 0.018;
    ProcData.data.Pupil.mmDiameter = diameter.*ProcData.data.Pupil.mmPerPixel;
    ProcData.data.Pupil.mmArea = pupilArea.*(ProcData.data.Pupil.mmPerPixel^2);
    save(procDataFileID,'ProcData')
% end

end