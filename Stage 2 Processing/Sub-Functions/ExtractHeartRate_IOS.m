function [] = ExtractHeartRate_IOS(procDataFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Qingguang Zhang
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs: 
%
%   Outputs:
%
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFiles,1)
    procDataFile = procDataFiles(a, :);
    disp(['Extracting heart rate from ProcData file ' num2str(a) ' of ' num2str(size(procDataFiles, 1)) '...']); disp(' ')
    load(procDataFile)
    
    % Pull out the left and right window heart rate. They should be essentiall6 identical
    [~, ~, ~, LH_HR] = FindHeartRate_IOS(ProcData.data.CBV.LH, ProcData.notes.CBVCamSamplingRate);
    [~, ~, ~, RH_HR] = FindHeartRate_IOS(ProcData.data.CBV.RH, ProcData.notes.CBVCamSamplingRate);
    HR = (LH_HR + RH_HR)/2;   % Average the two signals from the left and right windows
    patchedHR = horzcat(HR(1), HR, HR(end), HR(end));
    
    % Smooth the signal with a 4 Hz low pass 4th-order butterworth filter
    [B, A] = butter(4, 2/(ProcData.notes.CBVCamSamplingRate/2), 'low');
    heartRate = filtfilt(B, A, patchedHR);    % Filtered heart rate signal
    ProcData.data.heartRate = heartRate;
    save(procDataFile, 'ProcData');
    
end