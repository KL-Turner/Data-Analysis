function [SleepData] = SleepScore(animal, hem, rawDataFiles, procDataFiles, electrodeInput, RestingBaselines, SpectrogramData)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: This function serves as the focal point for the initial stages of sleep-specific data analysis. Each block of
%            code contains a function or set of functions with a commmon purpose.
%________________________________________________________________________________________________________________________
%
%   Inputs: The List of the Raw and Proc data files.
%
%   Outputs: //
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [1] Add Sleep parameters and logicals for sleep scoring to each ProcData file
disp('Analyzing Block [1] Add sleep parameters and logicals.'); disp(' ')
% Loop through all ProcData files in the current folder
% Add the parameters needed for sleep scoring and those that pertain to post-scoring analysis to each ProcData file
for fileNumber = 1:size(procDataFiles, 1)
    procDataFile = procDataFiles(fileNumber, :);   % Obtain the ProcData filename from the list
    load(procDataFile);   % Load that ProcData file
    rawDataFile = rawDataFiles(fileNumber, :);   % Obtain the associated RawData filename from the list
    load(rawDataFile);   % Load that RawData file
    [animal, hem, fileDate, fileID] = GetFileInfo(procDataFile);   % Find the fileDate of the current loaded file
    strDay = ConvertDate(fileDate);   % Convert the fileDate to a word string
    disp(['Adding sleep parameters and logicals to ProcData file number ' num2str(fileNumber) ' of ' num2str(size(procDataFiles, 1)) '...']); disp(' ')
    [ProcData] = AddSleepParameters(animal, hem, fileID, strDay, RawData, ProcData, RestingBaselines);   % Add sleep parameters to each ProcData file
    [ProcData] = AddSleepLogicals(animal, hem, fileID, ProcData, electrodeInput);  % Add sleep logical to each ProcData file
end

%% BLOCK PURPOSE: [2] Create SleepData.mat struct
disp('Analyzing Block [2] Creating SleepData structure.'); disp(' ')
sleepTime = 30;   % seconds
[SleepData] = CreateSleepData(procDataFiles, sleepTime);   % Create struct containing sleep epochs

%% BLOCK PURPOSE: [3] Create Single Trial Checks for sleeping and first 30 minutes of rest
disp('Analyzing Block [3] Manually add REM sleep events to the SleepData struct.'); disp(' ')
if ~isempty(SleepData)
    % Manually add REM sleep events
    [SleepData] = AddREMSleepEvents(animal, hem, RestingBaselines, SleepData, SpectrogramData);
    DisplaySleepTimes(animal, SleepData)
    disp('SleepScore.mat - complete'); disp(' ')   
else
    disp('SleepData appears to be empty. Ending Sleep analysis.'); disp(' ')
end

end
