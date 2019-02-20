function [RestData] = ExtractRestingData2(combDataFiles, dataTypes)
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
%   Outputs: RestData.mat
%________________________________________________________________________________________________________________________

if not(iscell(dataTypes))
    dataTypes = {dataTypes};
end

for dT = 1:length(dataTypes)
    dataType = dataTypes(dT);
    
    for f = 1:size(combDataFiles, 1)
        disp(['Gathering rest ' char(dataType) ' data from file ' num2str(f) ' of ' num2str(size(combDataFiles, 1)) '...']); disp(' ')
        filename = combDataFiles(f, :);
        load(filename);
        restVals = cell(size(combDataFiles, 1), 1);
        eventTimes = cell(size(combDataFiles, 1), 1);
        durations = cell(size(combDataFiles, 1), 1);
        puffDistances = cell(size(combDataFiles, 1), 1);
        fileIDs = cell(size(combDataFiles, 1), 1);
        fileDates = cell(size(combDataFiles, 1), 1);
        
        % Get the date and file identifier for the data to be saved with each resting event
        [animalID, fileDate, fileID, ~] = GetFileInfo2(filename);
        
        % Sampling frequency for element of dataTypes
        if strcmp(dataType, 'Vessel_Diameter')
            Fs = floor(CombData.Notes.MScan.frameRate);
        else
            Fs = 30;
        end
        
        % Expected number of samples for element of dataType
        expectedLength = (CombData.Notes.LabVIEW.trialDuration_Seconds-10)*Fs;
        
        % Get information about periods of rest from the loaded file
        trialEventTimes = CombData.Flags.rest.eventTime';
        trialPuffDistances = CombData.Flags.rest.puffDistance;
        trialDurations = CombData.Flags.rest.duration';
        
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
            trialRestVals{tET} = CombData.Data.(dataTypes{dT})(:, startInd:stopInd);
        end
        % Add all periods of rest to a cell array for all files
        restVals{f} = trialRestVals';
        
        % Transfer information about resting periods to the new structure
        eventTimes{f} = trialEventTimes';
        durations{f} = trialDurations';
        puffDistances{f} = trialPuffDistances';
        fileIDs{f} = repmat({fileID}, 1, length(trialEventTimes));
        fileDates{f} = repmat({fileDate}, 1, length(trialEventTimes));
    end
    % Combine the cells from separate files into a single cell array of all resting periods
    RestData.(dataTypes{dT}).data = [restVals{:}]';
    RestData.(dataTypes{dT}).eventTimes = cell2mat(eventTimes);
    RestData.(dataTypes{dT}).durations = cell2mat(durations);
    RestData.(dataTypes{dT}).puffDistances = [puffDistances{:}]';
    RestData.(dataTypes{dT}).fileIDs = [fileIDs{:}]';
    RestData.(dataTypes{dT}).fileDates = [fileDates{:}]';
    RestData.(dataTypes{dT}).samplingRate = Fs;
    
end

save([animalID '_RestData.mat'], 'RestData'); 

end
