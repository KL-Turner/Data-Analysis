function AddHeartRate(fileName)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Add the heart rate from IOS data via FindHeartRate.m to the ProcData.mat structure.
%________________________________________________________________________________________________________________________
%
%   Inputs: ProcData.mat file name.
%
%   Outputs: None - but saves the updated ProcData.mat file to the current folder.
%
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

load(fileName)   % ProcData file

% Pull out the left and right window heart rate. They should be essentiall6 identical
[~, tr, ~, LH_HR] = FindHeartRate(ProcData.Data.CBV.LH, ProcData.Notes.CBVCamSamplingRate);
[~, ~, ~, RH_HR] = FindHeartRate(ProcData.Data.CBV.RH, ProcData.Notes.CBVCamSamplingRate);
HR = (LH_HR + RH_HR)/ 2;   % Average the two signals from the left and right windows

% Smooth the signal with a 4 Hz low pass 4th-order butterworth filter
[B, A] = butter(4, 2 / (ProcData.Notes.CBVCamSamplingRate / 2), 'low');
HeartRate = filtfilt(B, A, HR);    % Filtered heart rate signal

ProcData.Data.HeartRate.HR = HeartRate;
ProcData.Data.HeartRate.tr = tr;

save(fileName, 'ProcData');

end
