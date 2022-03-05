function [] = CheckPupilDiameter_IOS(procDataFileID)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Manually check pupil diameter
%________________________________________________________________________________________________________________________

load(procDataFileID)
ProcData.data.Pupil.diameterCheckComplete = [];
if strcmpi(ProcData.data.Pupil.frameCheck,'y') == true %#ok<*NODEF>
    if isfield(ProcData.data.Pupil,'diameterCheckComplete') == false
        % load files and extract video information
        [animalID,~,fileID] = GetFileInfo_IOS(procDataFileID);
        [z,p,k] = butter(4,5/(ProcData.notes.pupilCamSamplingRate/2),'low');
        [sos,g] = zp2sos(z,p,k);
        pupilArea = ProcData.data.Pupil.pupilArea;
        try
            filtPupilArea = filtfilt(sos,g,pupilArea);
        catch
            filtPupilArea = pupilArea;
        end
        % request user input for this file
        check = false;
        while check == false
            % figure showing 10 frames
            diameterCheck = figure;
            plot((1:length(filtPupilArea))./ProcData.notes.pupilCamSamplingRate,filtPupilArea,'k');
            title([animalID ' ' strrep(fileID,'_',' ')])
            xlabel('Time (sec)')
            ylabel('Pupil area (pixels)')
            axis tight
            keepDiameter = input('Is this an accurate pupil tracking diameter? (y/n): ','s'); disp(' ')
            close(diameterCheck)
            if strcmp(keepDiameter,'y') == true || strcmp(keepDiameter,'n') == true
                ProcData.data.Pupil.diameterCheck = keepDiameter;
                check = true;
            end
        end
        ProcData.data.Pupil.diameterCheckComplete = 'y';
        save(procDataFileID,'ProcData')
    end
elseif strcmp(ProcData.data.Pupil.frameCheck,'n') == true
    ProcData.data.Pupil.diameterCheck = 'n';
    ProcData.data.Pupil.diameterCheckComplete = 'y';
    save(procDataFileID,'ProcData')
end

end
