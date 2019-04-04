function [RestData] = ExtractRestingData_2P(mergedDataFiles, dataTypes)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Evalute the flags from each MergedData file and pull out the corresponding rest data associated with each
%            dataType.
%________________________________________________________________________________________________________________________
%
%   Inputs: List of all MergedData.mat files and the dataTypes that should be analyzed.
%
%   Outputs: A RestData.mat structure containing all periods of rest and their associated file IDs, event times, durations,
%            and corresponding data. This structure can later be filtered to extract the events greater than a certain duration.
%
%   Last Revised: March 21st, 2019
%________________________________________________________________________________________________________________________

RestData = [];

if not(iscell(dataTypes))
    dataTypes = {dataTypes};
end

for dT = 1:length(dataTypes)
    dataType = dataTypes(dT);
    restVals = cell(size(mergedDataFiles, 1), 1);
    eventTimes = cell(size(mergedDataFiles, 1), 1);
    durations = cell(size(mergedDataFiles, 1), 1);
    puffDistances = cell(size(mergedDataFiles, 1), 1);
    fileIDs = cell(size(mergedDataFiles, 1), 1);
    fileDates = cell(size(mergedDataFiles, 1), 1);
    vesselIDs = cell(size(mergedDataFiles, 1), 1);
    
    for f = 1:size(mergedDataFiles, 1)
        disp(['Gathering rest ' char(dataType) ' data from file ' num2str(f) ' of ' num2str(size(mergedDataFiles, 1)) '...']); disp(' ')
        filename = mergedDataFiles(f, :);
        load(filename);
        
        % Get the date and file identifier for the data to be saved with each resting event
        [animalID, fileDate, fileID, vesselID, ~] = GetFileInfo2_2P(filename);
        
        % Sampling frequency for element of dataTypes
        if strcmp(dataType, 'vesselDiameter')
            Fs = floor(MergedData.notes.p2Fs);
        else
            Fs = floor(MergedData.notes.dsFs);
        end
        
        % Expected number of samples for element of dataType
        expectedLength = MergedData.notes.trialDuration_Sec*Fs;
        
        % Get information about periods of rest from the loaded file
        trialEventTimes = MergedData.flags.rest.eventTime';
        trialPuffDistances = MergedData.flags.rest.puffDistance;
        trialDurations = MergedData.flags.rest.duration';
        
        % Initialize cell array for all periods of rest from the loaded file
        trialRestVals = cell(size(trialEventTimes'));
        for tET = 1:length(trialEventTimes)
            % Extract the whole duration of the resting event. Coerce the
            % start index to values above 1 to preclude rounding to 0.
            startInd = max(floor(trialEventTimes(tET)*Fs), 1);
            
            % Convert the duration from seconds to samples.
            dur = round(trialDurations(tET)*Fs);
            
            % Get ending index for data chunk. If event occurs at the end of
            % the trial, assume animal whisks as soon as the trial ends and
            % give a 200ms buffer.
            stopInd = min(startInd + dur, expectedLength - round(0.2*Fs));
            
            % Extract data from the trial and add to the cell array for the current loaded file
            trialRestVals{tET} = MergedData.data.(dataTypes{dT})(:, startInd:stopInd);
        end
        % Add all periods of rest to a cell array for all files
        restVals{f} = trialRestVals';
        
        % Transfer information about resting periods to the new structure
        eventTimes{f} = trialEventTimes';
        durations{f} = trialDurations';
        puffDistances{f} = trialPuffDistances';
        fileIDs{f} = repmat({fileID}, 1, length(trialEventTimes));
        fileDates{f} = repmat({fileDate}, 1, length(trialEventTimes));
        vesselIDs{f} = repmat({vesselID}, 1, length(trialEventTimes));
    end
    
    RestData.(dataTypes{dT}).data = [restVals{:}]';
    RestData.(dataTypes{dT}).eventTimes = cell2mat(eventTimes);
    RestData.(dataTypes{dT}).durations = cell2mat(durations);
    RestData.(dataTypes{dT}).puffDistances = [puffDistances{:}]';
    RestData.(dataTypes{dT}).fileIDs = [fileIDs{:}]';
    RestData.(dataTypes{dT}).fileDates = [fileDates{:}]';
    RestData.(dataTypes{dT}).vesselIDs = [vesselIDs{:}]';
    RestData.(dataTypes{dT}).samplingRate = Fs;
end

save([animalID '_RestData.mat'], 'RestData'); 

end
