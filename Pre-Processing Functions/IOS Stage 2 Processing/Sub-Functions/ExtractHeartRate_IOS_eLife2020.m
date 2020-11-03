function [] = ExtractHeartRate_IOS_eLife2020(procDataFiles,imagingType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Qingguang Zhang
%________________________________________________________________________________________________________________________
%
%   Purpose: Use the spectral properties of the CBV data to extract the heart rate.
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFiles,1)
    procDataFile = procDataFiles(a,:);
    disp(['Extracting heart rate from ProcData file ' num2str(a) ' of ' num2str(size(procDataFiles, 1)) '...']); disp(' ')
    load(procDataFile)
    % control for imaging type
    if strcmp(imagingType,'bilateral') == true
        % Pull out the left and right window heart rate. They should be essentiall6 identical
        [~,~,~,LH_HR] = FindHeartRate_IOS_eLife2020(ProcData.data.CBV.LH,ProcData.notes.CBVCamSamplingRate);
        [~,~,~,RH_HR] = FindHeartRate_IOS_eLife2020(ProcData.data.CBV.RH,ProcData.notes.CBVCamSamplingRate);
        % Average the two signals from the left and right windows
        HR = (LH_HR + RH_HR)/2;
    elseif strcmp(imagingType,'single') == true
        [~,~,~,HR] = FindHeartRate_IOS_eLife2020(ProcData.data.CBV.Barrels,ProcData.notes.CBVCamSamplingRate);
    end
    % patch the missing data at the beginning and end of the signal
    patchedHR = horzcat(HR(1),HR,HR(end),HR(end));
    % Smooth the signal with a 2 Hz low pass third-order butterworth filter
    [B,A] = butter(3,2/(ProcData.notes.CBVCamSamplingRate/2),'low');
    heartRate = filtfilt(B,A,patchedHR);   % Filtered heart rate signal
    ProcData.data.heartRate = heartRate;
    save(procDataFile,'ProcData');
end
