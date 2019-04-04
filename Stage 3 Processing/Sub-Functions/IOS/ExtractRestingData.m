function [RestData] = ExtractRestingData(procDataFiles, dataTypes)
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
    if strcmp(dataType, 'CBV')
        hemisphereDataTypes = {'LH', 'LH_Electrode', 'RH', 'RH_Electrode'};
    else
        hemisphereDataTypes = {'LH', 'RH'};
    end
    
    for hDT = 1:length(hemisphereDataTypes)
        % Initialize cell arrays for resting data and other information.
        restVals = cell(size(procDataFiles, 1), 1);
        eventTimes = cell(size(procDataFiles, 1), 1);
        durations = cell(size(procDataFiles, 1), 1);
        puffDistances = cell(size(procDataFiles, 1), 1);
        fileIDs = cell(size(procDataFiles, 1), 1);
        fileDates = cell(size(procDataFiles, 1), 1);

        for f = 1:size(procDataFiles, 1)
            disp(['Gathering '  char(hemisphereDataTypes(hDT)) ' rest ' char(dataType) ' data from file ' num2str(f) ' of ' num2str(size(procDataFiles, 1)) '...']); disp(' ')
            filename = procDataFiles(f, :);
            load(filename);

            % Get the date and file identifier for the data to be saved with each resting event
            [animal, ~, fileDate, fileID] = GetFileInfo(filename);

            % Sampling frequency for element of dataTypes
            Fs = ProcData.Notes.CBVCamSamplingRate;                                     

            % Expected number of samples for element of dataType                      
            expectedLength = ProcData.Notes.trialDuration_Seconds*Fs;

            % Get information about periods of rest from the loaded file
            trialEventTimes = ProcData.Flags.rest.eventTime';
            trialPuffDistances = ProcData.Flags.rest.puffDistance;
            trialDurations = ProcData.Flags.rest.duration';

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
                trialRestVals{tET} = ProcData.Data.(dataTypes{dT}).(hemisphereDataTypes{hDT})(:, startInd:stopInd);            
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
        RestData.(dataTypes{dT}).(hemisphereDataTypes{hDT}).data = [restVals{:}]';
        RestData.(dataTypes{dT}).(hemisphereDataTypes{hDT}).eventTimes = cell2mat(eventTimes);
        RestData.(dataTypes{dT}).(hemisphereDataTypes{hDT}).durations = cell2mat(durations);
        RestData.(dataTypes{dT}).(hemisphereDataTypes{hDT}).puffDistances = [puffDistances{:}]';
        RestData.(dataTypes{dT}).(hemisphereDataTypes{hDT}).fileIDs = [fileIDs{:}]';
        RestData.(dataTypes{dT}).(hemisphereDataTypes{hDT}).fileDates = [fileDates{:}]';
        RestData.(dataTypes{dT}).(hemisphereDataTypes{hDT}).CBVCamSamplingRate = Fs;
    end
end

save([animal '_RestData.mat'], 'RestData'); 

end
