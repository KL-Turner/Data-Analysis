function [RestData] = ExtractRestingData_IOS(procdataFiles,dataTypes,imagingType)
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

if not(iscell(dataTypes))
    dataTypes = {dataTypes};
end

for a = 1:length(dataTypes)
    dataType = dataTypes(a);
    if strcmp(dataType,'CBV') == true || strcmp(dataType,'CBV_HbT') == true
        if strcmp(imagingType,'bilateral') == true
            subDataTypes = {'LH','adjLH','RH','adjRH'};
        elseif strcmp(imagingType,'single') == true
            subDataTypes = {'Barrels','adjBarrels'};
        end
    elseif strcmp(dataType,'EMG') == true
        subDataTypes = {'emg'};
    elseif strcmp(dataType,'flow') == true
        subDataTypes = {'data'};
    else
        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','muaPower'};
    end
    
    for b = 1:length(subDataTypes)
        % Initialize cell arrays for resting data and other information.
        restVals = cell(size(procdataFiles, 1), 1);
        eventTimes = cell(size(procdataFiles, 1), 1);
        durations = cell(size(procdataFiles, 1), 1);
        puffDistances = cell(size(procdataFiles, 1), 1);
        fileIDs = cell(size(procdataFiles, 1), 1);
        fileDates = cell(size(procdataFiles, 1), 1);

        for c = 1:size(procdataFiles, 1)
            disp(['Extracting ' char(dataType) ' ' char(subDataTypes(b)) ' rest data from file ' num2str(c) ' of ' num2str(size(procdataFiles, 1)) '...']); disp(' ')
            procdataFile = procdataFiles(c, :);
            load(procdataFile);

            % Get the date and file identifier for the data to be saved with each resting event
            [animal, fileDate, fileID] = GetFileInfo_IOS(procdataFile);

            % Sampling frequency for element of dataTypes
            Fs = ProcData.notes.CBVCamSamplingRate;                                     
            
            % Expected number of samples for element of dataType
            trialDuration_sec = ProcData.notes.trialDuration_sec;
            expectedLength = trialDuration_sec*Fs;
            % Get information about periods of rest from the loaded file
            trialEventTimes = ProcData.flags.rest.eventTime';
            trialPuffDistances = ProcData.flags.rest.puffDistance;
            trialDurations = ProcData.flags.rest.duration';

            % Initialize cell array for all periods of rest from the loaded file
            trialRestVals = cell(size(trialEventTimes'));
            for d = 1:length(trialEventTimes)
                % Extract the whole duration of the resting event. Coerce the 
                % start index to values above 1 to preclude rounding to 0.
                startInd = max(floor(trialEventTimes(d)*Fs), 1);

                % Convert the duration from seconds to samples.
                dur = round(trialDurations(d)*Fs); 

                % Get ending index for data chunk. If event occurs at the end of
                % the trial, assume animal whisks as soon as the trial ends and
                % give a 200ms buffer.
                stopInd = min(startInd + dur, expectedLength - round(0.2*Fs));
                
                try
                    % Extract data from the trial and add to the cell array for the current loaded file
                    trialRestVals{d} = ProcData.data.(dataTypes{a}).(subDataTypes{b})(:, startInd:stopInd);
                catch % some files don't have certain fields. Skip those
                    trialRestVals{d} = [];
                end
            end
                
            % Add all periods of rest to a cell array for all files
            restVals{c} = trialRestVals';

            % Transfer information about resting periods to the new structure
            eventTimes{c} = trialEventTimes';
            durations{c} = trialDurations';
            puffDistances{c} = trialPuffDistances';
            fileIDs{c} = repmat({fileID}, 1, length(trialEventTimes));
            fileDates{c} = repmat({fileDate}, 1, length(trialEventTimes));
        end
        % Combine the cells from separate files into a single cell array of all resting periods
        RestData.(dataTypes{a}).(subDataTypes{b}).data = [restVals{:}]';
        RestData.(dataTypes{a}).(subDataTypes{b}).eventTimes = cell2mat(eventTimes);
        RestData.(dataTypes{a}).(subDataTypes{b}).durations = cell2mat(durations);
        RestData.(dataTypes{a}).(subDataTypes{b}).puffDistances = [puffDistances{:}]';
        RestData.(dataTypes{a}).(subDataTypes{b}).fileIDs = [fileIDs{:}]';
        RestData.(dataTypes{a}).(subDataTypes{b}).fileDates = [fileDates{:}]';
        RestData.(dataTypes{a}).(subDataTypes{b}).CBVCamSamplingRate = Fs;
        RestData.(dataTypes{a}).(subDataTypes{b}).trialDuration_sec = trialDuration_sec;
    end
end

save([animal '_RestData.mat'], 'RestData','-v7.3'); 

end
