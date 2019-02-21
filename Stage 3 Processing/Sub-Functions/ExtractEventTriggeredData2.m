function [EventData] = ExtractEventTriggeredData2(MergedDataFiles, dataTypes)
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
%                   MergedDataFiles - [matrix] names of files organized as rows.
%                   Files should already be processed using the script
%                   "ProcessRawDataFile.m" or "CalculatePredictedCBV.m".
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

for dT = 1:length(dataTypes)
    temp = struct();
    dataType = dataTypes{dT};
    
    for f = 1:size(MergedDataFiles, 1)    
        % Load MergedData File
        filename = MergedDataFiles(f, :);
        load(filename);

        % Get the date and file ID to include in the EventData structure
        [animal, fileDate, fileID, ~] = GetFileInfo2(MergedDataFiles(f,:));

        % Get the types of behaviors present in the file (stim,whisk,rest)
        holdData = fieldnames(MergedData.Flags);
        behaviorFields = holdData([1 2],1);
        
        % Sampling frequency for element of dataTypes
        if strcmp(dataType, 'Vessel_Diameter')
            Fs = floor(MergedData.Notes.MScan.frameRate);
        else
            Fs = 30;
        end
            % Loop over the behaviors present in the file
            for bF = 1:length(behaviorFields)
                % Create behavioral subfields for the temp structure, if needed
                if not(isfield(temp, behaviorFields{bF}))
                    subFields = fieldnames(MergedData.Flags.(behaviorFields{bF}));
                    blankCell = cell(1, size(MergedDataFiles, 1));
                    structVals = cell(size(subFields));
                    structVals(:) = {blankCell};
                    temp.(behaviorFields{bF}) = cell2struct(structVals, subFields, 1)';
                    temp.(behaviorFields{bF}).fileIDs = blankCell;
                    temp.(behaviorFields{bF}).fileDates = blankCell;
                    temp.(behaviorFields{bF}).data = blankCell;
                end

                % Assemble a structure to send to the sub-functions
                data = MergedData.Data;
                data.Flags = MergedData.Flags;
                data.Notes = MergedData.Notes;

                % Extract the data from the epoch surrounding the event
                disp(['Extracting event-triggered ' dataType ' ' behaviorFields{bF} ' data from file ' num2str(f) ' of ' num2str(size(MergedDataFiles, 1)) '...']); disp(' ');
                [chunkData, evFilter] = ExtractBehavioralData(data, epoch, dataType, Fs, behaviorFields{bF});

                % Add epoch details to temp struct
                [temp] = AddEpochInfo(data, behaviorFields{bF}, temp, fileID, fileDate, evFilter, f);
                temp.(behaviorFields{bF}).data{f} = chunkData;
            end 
    end
    % Convert the temporary stuct into a final structure
    [EventData] = ProcessTempStruct(EventData, dataType, temp, epoch);
end
save([animal '_EventData.mat'], 'EventData');

function [chunkData, evFilter] = ExtractBehavioralData(data, epoch, dataType, Fs, behavior)
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
trialDuration = (data.Notes.LabVIEW.trialDuration_Seconds - 10);

% Get the content from Data.(dataType)
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
% disp(['ExtractEventTriggeredData > ExtractBehavioralData:' upper(behavior) ': Events at times: ' num2str(eventTimes(not(evFilter))') ' seconds omitted. Too close to beginning/end of trial.']);
% disp(' ');

% Convert the starts from seconds to samples, round down to the nearest
% sample, coerce the value above 1.

epochStarts = max(floor(allEpochStarts(evFilter)*Fs),1);

% Calculate stops indices using the duration of the epoch, this avoids
% potential matrix dimension erros caused by differences in rounding when
% converting from seconds to samples.
sampleDur = round(epoch.duration*Fs);
epochStops = epochStarts + sampleDur*ones(size(epochStarts));

% Extract the chunk of data from the trial
chunkData = zeros(sampleDur + 1, length(epochStarts), size(data, 1));

for eS = 1:length(epochStarts)
    chunkInds = epochStarts(eS):epochStops(eS);
    chunkData(:, eS, :) = data(:, chunkInds)';
end


function [temp] = AddEpochInfo(data, behavior, temp, fileID, fileDate, evFilter, f)
%   function [temp] = AddEpochInfo(Data,dataType,Beh,temp,fileID,FileDate,EvFilter,f)
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
for flds = 1:length(fields)
    field = fields{flds};
    temp.(behavior).(field){f} = data.Flags.(behavior).(field)(evFilter,:)';
end

% Tag each event with the file ID, arrange cell array horizontally for
% later processing.
temp.(behavior).fileIDs{f} = repmat({fileID}, 1, sum(evFilter));
temp.(behavior).fileDates{f} = repmat({fileDate}, 1, sum(evFilter));


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


% Get dataType names
behaviorFields = fieldnames(temp);

% Intialize Behavior fields of the dataType sub-structure
structArray2 = cell(size(behaviorFields));
EventData.(dataType) = cell2struct(structArray2, behaviorFields, 1);

for bF = 1:length(behaviorFields)
    behavior = behaviorFields{bF};
    
    % Get Behavior names
    eventFields = fieldnames(temp.(behavior));
    
    % Initialize Event fields for the Behavior sub-structure
    structArray3 = cell(size(eventFields));
    EventData.(dataType).(behavior) = cell2struct(structArray3, eventFields, 1);
    
    for eF = 1:length(eventFields)
        evField = eventFields{eF};
        transferArray = [temp.(behavior).(evField){:}];
        EventData.(dataType).(behavior).(evField) = permute(transferArray, unique([2, 1, ndims(transferArray)], 'stable'));
    end
    EventData.(dataType).(behavior).epoch = epoch;
end
