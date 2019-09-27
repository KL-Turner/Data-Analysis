function [EventData] = ExtractEventTriggeredData_IOS(procdataFiles, dataTypes)
%___________________________________________________________________________________________________
% Edited by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%
%   Last Revised: October 3rd, 2018
%___________________________________________________________________________________________________
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Separates data corresponding to various behaviors into
%   structures
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   ProcDataFiles - [matrix] names of files organized as rows.
%                   Files should already be processed using the script
%                   "ProcessRawdataFile.m" or "CalculatePredictedCBV.m".
%
%                   EventData - [struct] structure of existing data chunked
%                   around behavioral events
%
%                   dataTypes - [dataTypes] the measurements to be
%                   chunked into data epochs
%
%                   epoch - [struct] contains fields:
%                               duration - [double] the total time of the
%                               data epoch
%
%                               offset - [double] the amount of pre-event
%                               time (s) to include in the epoch.
%_______________________________________________________________
%   RETURN:                     
%                   EventData - [struct] data chunked around behavioral
%                   events
%_______________________________________________________________

EventData = [];
epoch.duration = 12;
epoch.offset = 2;

% Control for dataTypes as string
if not(iscell(dataTypes))
    dataTypes = {dataTypes};
end

for a = 1:length(dataTypes)
    dataType = char(dataTypes(a));
    if strcmp(dataType, 'CBV') == true || strcmp(dataType, 'CBV_HbT') == true
        subdataTypes = {'LH', 'LH_Electrode', 'RH', 'RH_Electrode'};
    elseif strcmp(dataType, 'EMG')
        subdataTypes = {'emg'};
    else
        subdataTypes = {'deltaBandPower', 'thetaBandPower', 'alphaBandPower', 'betaBandPower', 'gammaBandPower', 'muaPower'};
    end

    temp = struct();
    
    for b = 1:size(procdataFiles, 1)    
        % Load ProcData File
        filename = procdataFiles(b, :);
        load(filename);

        % Get the date and file ID to include in the EventData structure
        [animal, fileDate, fileID] = GetFileInfo_IOS(procdataFiles(b,:));

        % Get the types of behaviors present in the file (stim,whisk,rest)
        holddata = fieldnames(ProcData.flags);
        behaviorFields = holddata([1 2],1);

        for c = 1:length(subdataTypes)
            sDT = char(subdataTypes(c));

            % Set the sampling frequency for the dataType
            samplingRate = ProcData.notes.dsFs;
            trialDuration_sec = ProcData.notes.trialDuration_sec;

            % Loop over the behaviors present in the file
            for d = 1:length(behaviorFields)

                %Preallocate space for unknown number of events using a
                %'temporary' structure of cells
                if not(isfield(temp, sDT))
                    temp.(sDT) = [];
                end

                % Create behavioral subfields for the temp structure, if needed
                if not(isfield(temp.(sDT), behaviorFields{d}))
                    subFields = fieldnames(ProcData.flags.(behaviorFields{d}));
                    blankCell = cell(1, size(procdataFiles, 1));
                    structVals = cell(size(subFields));
                    structVals(:) = {blankCell};
                    temp.(sDT).(behaviorFields{d}) = cell2struct(structVals, subFields, 1)';
                    temp.(sDT).(behaviorFields{d}).fileIDs = blankCell;
                    temp.(sDT).(behaviorFields{d}).fileDates = blankCell;
                    temp.(sDT).(behaviorFields{d}).data = blankCell;
                end

                % Assemble a structure to send to the sub-functions
                fieldName2 = dataType;

                if isempty(fieldName2)
                    data = ProcData;
                else
                    data = ProcData.data.(fieldName2);
                end

                data.Flags = ProcData.flags;
                data.notes = ProcData.notes;

                % Extract the data from the epoch surrounding the event
                disp(['Extracting ' dataType ' ' sDT ' event-triggered ' behaviorFields{d} ' data from file ' num2str(b) ' of ' num2str(size(procdataFiles, 1)) '...']); disp(' ');
                [chunkdata, evFilter] = ExtractBehavioraldata(data, epoch, sDT, behaviorFields{d});

                % Add epoch details to temp struct
                [temp] = AddEpochInfo(data, sDT, behaviorFields{d}, temp, fileID, fileDate, evFilter, b);
                temp.(sDT).(behaviorFields{d}).data{b} = chunkdata;

                % Add the sampling frequency, assume all Fs are the same for given
                % dataType
                temp.(sDT).(behaviorFields{d}).samplingRate = {samplingRate};
                temp.(sDT).(behaviorFields{d}).trialDuration_sec = {trialDuration_sec};
            end 
        end
    end
    % Convert the temporary stuct into a final structure
    [EventData] = ProcessTempStruct(EventData, dataType, temp, epoch);
end
save([animal '_EventData.mat'], 'EventData');

function [chunkdata, evFilter] = ExtractBehavioraldata(data, epoch, dataType, behavior)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: 
%   
%_______________________________________________________________
%   PARAMETERS:             
%                               
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

% Setup variables
eventTimes = data.Flags.(behavior).eventTime;
trialDuration = data.notes.trialDuration_sec;
samplingRate = data.notes.dsFs;

% Get the content from data.(dataType)
data = getfield(data, {}, dataType, {});

% Calculate start/stop times (seconds) for the events
allEpochStarts = eventTimes - epoch.offset*ones(size(eventTimes));
allEpochEnds = allEpochStarts + epoch.duration*ones(size(eventTimes));

% Filter out events which are too close to the beginning or end of trials
startFilter = allEpochStarts > 0;
stopFilter = round(allEpochEnds) < trialDuration; % Apply "round" to give an 
                                              % extra half second buffer 
                                              % and prevent indexing errors
evFilter = logical(startFilter.*stopFilter);
% disp(['ExtractEventTriggereddata > ExtractBehavioraldata:' upper(behavior) ': Events at times: ' num2str(eventTimes(not(evFilter))') ' seconds omitted. Too close to beginning/end of trial.']);
% disp(' ');

% Convert the starts from seconds to samples, round down to the nearest
% sample, coerce the value above 1.

epochStarts = max(floor(allEpochStarts(evFilter)*samplingRate),1);

% Calculate stops indices using the duration of the epoch, this avoids
% potential matrix dimension erros caused by differences in rounding when
% converting from seconds to samples.
sampleDur = round(epoch.duration*samplingRate);
epochStops = epochStarts + sampleDur*ones(size(epochStarts));

% Extract the chunk of data from the trial
chunkdata = zeros(sampleDur + 1, length(epochStarts), size(data, 1));

for a = 1:length(epochStarts)
    chunkInds = epochStarts(a):epochStops(a);
    chunkdata(:, a, :) = data(:, chunkInds)';
end


function [temp] = AddEpochInfo(data, dataType, behavior, temp, fileID, fileDate, evFilter, f)
%   function [temp] = AddEpochInfo(data,dataType,Beh,temp,fileID,FileDate,EvFilter,f)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: 
%   
%_______________________________________________________________
%   PARAMETERS:             
%                               
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

% Get the field names for each behavior
fields = fieldnames(data.Flags.(behavior));

% Filter out the events which are too close to the trial edge
for a = 1:length(fields)
    field = fields{a};
    temp.(dataType).(behavior).(field){f} = data.Flags.(behavior).(field)(evFilter,:)';
end

% Tag each event with the file ID, arrange cell array horizontally for
% later processing.
temp.(dataType).(behavior).fileIDs{f} = repmat({fileID}, 1, sum(evFilter));
temp.(dataType).(behavior).fileDates{f} = repmat({fileDate}, 1, sum(evFilter));


function [EventData] = ProcessTempStruct(EventData, dataType, temp, epoch)
%   [EventData] = ProcessTempStruct(EventData,temp,epoch,freqs)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: 
%   
%_______________________________________________________________
%   PARAMETERS:             
%                               
%_______________________________________________________________
%   RETURN:                     
%                               
%______________________________________________________________

% Get the dataTypes from temp
dTs = fieldnames(temp);

for a = 1:length(dTs)
    dT = dTs{a};
    
    % Get dataType names
    behaviorFields = fieldnames(temp.(dT));
    
    % Intialize Behavior fields of the dataType sub-structure
    structArray2 = cell(size(behaviorFields));
    EventData.(dataType).(dT) = cell2struct(structArray2, behaviorFields, 1);
    
    for b = 1:length(behaviorFields)
        behavior = behaviorFields{b};
        
        % Get Behavior names
        eventFields = fieldnames(temp.(dT).(behavior));
        
        % Initialize Event fields for the Behavior sub-structure
        structArray3 = cell(size(eventFields));
        EventData.(dataType).(dT).(behavior) = cell2struct(structArray3, eventFields, 1);
        
        for c = 1:length(eventFields)
            evField = eventFields{c};
            transferArray = [temp.(dT).(behavior).(evField){:}];
            EventData.(dataType).(dT).(behavior).(evField) = permute(transferArray, unique([2, 1, ndims(transferArray)], 'stable'));
        end
        EventData.(dataType).(dT).(behavior).epoch = epoch;
    end
end