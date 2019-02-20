function AddHeartRate(fileName)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs: 
%
%   Outputs: 
%________________________________________________________________________________________________________________________

load(fileName)

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
