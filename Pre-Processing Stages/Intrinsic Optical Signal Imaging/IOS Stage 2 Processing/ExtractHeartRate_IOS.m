function [] = ExtractHeartRate_IOS(procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Use the spectral properties of the CBV data to extract the heart rate.
%________________________________________________________________________________________________________________________

load(procDataFileIDs(1,:));
imagingType = ProcData.notes.imagingType;
for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Extracting heart rate from ProcData file ' num2str(a) ' of ' num2str(size(procDataFileIDs,1)) '...']); disp(' ')
    load(procDataFileID)
    if ProcData.notes.CBVCamSamplingRate >= 30
        if strcmpi(imagingType,'Single ROI (SI)') == true
            [~,~,~,HR] = FindHeartRate_IOS(ProcData.data.CBV.Barrels,ProcData.notes.CBVCamSamplingRate);
        elseif strcmpi(imagingType,'Single ROI (SSS)') == true
            [~,~,~,HR] = FindHeartRate_IOS(ProcData.data.CBV.SSS,ProcData.notes.CBVCamSamplingRate);
        elseif strcmpi(imagingType,'Bilateral ROI (SI)') == true || strcmpi(imagingType,'Bilateral ROI (SI,FC)') == true
            % pull out the left and right window heart rate. they should be essentiall6 identical
            [~,~,~,LH_HR] = FindHeartRate_IOS(ProcData.data.CBV.LH,ProcData.notes.CBVCamSamplingRate);
            [~,~,~,RH_HR] = FindHeartRate_IOS(ProcData.data.CBV.RH,ProcData.notes.CBVCamSamplingRate);
            % average the two signals from the left and right windows
            HR = (LH_HR + RH_HR)/2;
        end
        % patch the missing data at the beginning and end of the signal
        patchedHR = horzcat(HR(1),HR,HR(end),HR(end));
        % smooth the signal with a 2 Hz low pass third-order butterworth filter
        [B,A] = butter(3,2/(ProcData.notes.CBVCamSamplingRate/2),'low');
        heartRate = filtfilt(B,A,patchedHR); % filtered heart rate signal
        ProcData.data.heartRate = heartRate;
        save(procDataFileID,'ProcData');
    else
        ProcData.data.heartRate = zeros(1,ProcData.notes.trialDuration_sec);
        save(procDataFileID,'ProcData');
    end
end
