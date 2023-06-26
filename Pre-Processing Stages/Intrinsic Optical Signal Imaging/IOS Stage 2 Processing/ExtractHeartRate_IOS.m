function [] = ExtractHeartRate_IOS(procDataFileIDs)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
load(procDataFileIDs(1,:));
imagingType = ProcData.notes.imagingType;
imagingWavelengths = ProcData.notes.imagingWavelengths;
for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    load(procDataFileID)
    if isfield(ProcData.data,'heartRate') == false
        disp(['Extracting heart rate from file (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
        if ProcData.notes.CBVCamSamplingRate >= 30
            if strcmpi(imagingType,'Single ROI (SI)') == true
                [~,~,~,HR] = FindHeartRate_IOS(ProcData.data.green.barrels,ProcData.notes.CBVCamSamplingRate);
            elseif strcmpi(imagingType,'Single ROI (SSS)') == true
                [~,~,~,HR] = FindHeartRate_IOS(ProcData.data.greenRefl.SSS,ProcData.notes.CBVCamSamplingRate);
            elseif strcmpi(imagingType,'Bilateral ROI (SI)') == true || strcmpi(imagingType,'Bilateral ROI (SI,FC)') == true
                % pull out the left and right window heart rate. they should be essentiall6 identical
                [~,~,~,LH_HR] = FindHeartRate_IOS(ProcData.data.green.LH,ProcData.notes.CBVCamSamplingRate);
                [~,~,~,RH_HR] = FindHeartRate_IOS(ProcData.data.green.RH,ProcData.notes.CBVCamSamplingRate);
                % average the two signals from the left and right windows
                HR = (LH_HR + RH_HR)/2;
            end
            % patch the missing data at the beginning and end of the signal
            patchedHR = horzcat(HR(1),HR,HR(end),HR(end));
            % smooth the signal with a 2 Hz low pass third-order butterworth filter
            [B,A] = butter(3,2/(ProcData.notes.CBVCamSamplingRate/2),'low');
            heartRate = filtfilt(B,A,patchedHR); % filtered heart rate signal
            ProcData.data.heartRate.frequency = heartRate;
            save(procDataFileID,'ProcData');
        else
            ProcData.data.heartRate.frequency = zeros(1,ProcData.notes.trialDuration_sec);
            save(procDataFileID,'ProcData');
        end
    else
        disp(['Heart rate for file (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ') already analyzed.']); disp(' ')
    end
end