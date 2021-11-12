function [EventData] = ExtractPupilEventTriggeredData_IOS(procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
% Adapted from code written by Aaron T. Winder
%________________________________________________________________________________________________________________________
%
% Purpose: Extracts event-triggered pupil data using behavioral flags, add
%________________________________________________________________________________________________________________________

% load rest data file
eventDataFileID = ls('*_EventData.mat');
load(eventDataFileID)
EventData.Pupil = [];
epoch.duration = 12;
epoch.offset = 2;
% Control for dataTypes as string
dataTypes = {'Pupil'};
for a = 1:length(dataTypes)
    dataType = char(dataTypes(a));
    subDataTypes = {'pupilArea','diameter','mmArea','mmDiameter'};
    temp = struct();
    for b = 1:size(procDataFileIDs,1)
        % Load ProcData File
        procDataFileID = procDataFileIDs(b,:);
        load(procDataFileID);
        if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
            % Get the date and file ID to include in the EventData structure
            [animal,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
            % Get the types of behaviors present in the file (stim,whisk,rest)
            holdData = fieldnames(ProcData.flags);
            behaviorFields = holdData([1,2],1);
            for c = 1:length(subDataTypes)
                sDT = char(subDataTypes(c));
                % Set the sampling frequency for the dataType
                samplingRate = ProcData.notes.dsFs;
                trialDuration_sec = ProcData.notes.trialDuration_sec;
                % Loop over the behaviors present in the file
                for d = 1:length(behaviorFields)
                    %Preallocate space for unknown number of events using a
                    %'temporary' structure of cells
                    if not(isfield(temp,sDT))
                        temp.(sDT) = [];
                    end
                    % Assemble a structure to send to the sub-functions
                    fieldName2 = dataType;
                    try
                        data = ProcData.data.(fieldName2);
                    catch % some files don't have certain fields. Skip those
                        data = [];
                    end
                    data.Flags = ProcData.flags;
                    data.notes = ProcData.notes;
                    % Extract the data from the epoch surrounding the event
                    disp(['Extracting ' dataType ' ' sDT ' event-triggered ' behaviorFields{d} ' data from file ' num2str(b) ' of ' num2str(size(procDataFileIDs,1)) '...']); disp(' ');
                    try
                        [chunkdata,evFilter] = ExtractBehavioralData(data,epoch,sDT,behaviorFields{d});
                    catch
                        chunkdata = [];
                        evFilter = [];
                    end
                    % Add epoch details to temp struct
                    [temp] = AddEpochInfo(data,sDT,behaviorFields{d},temp,fileID,fileDate,evFilter,b);
                    temp.(sDT).(behaviorFields{d}).data{b} = chunkdata;
                    % Add the sampling frequency, assume all Fs are the same for given
                    % dataType
                    temp.(sDT).(behaviorFields{d}).samplingRate = {samplingRate};
                    temp.(sDT).(behaviorFields{d}).trialDuration_sec = {trialDuration_sec};
                end
            end
        end
    end
    % Convert the temporary stuct into a final structure
    [EventData] = ProcessTempStruct(EventData,dataType,temp,epoch);
end
save([animal '_EventData.mat'],'EventData','-v7.3');

end

function [chunkdata,evFilter] = ExtractBehavioralData(data,epoch,dataType,behavior)
% Setup variables
eventTime = data.Flags.(behavior).eventTime;
trialDuration = data.notes.trialDuration_sec;
samplingRate = data.notes.dsFs;
% Get the content from data.(dataType)
data = getfield(data,{},dataType,{});
% Calculate start/stop times (seconds) for the events
allEpochStarts = eventTime - epoch.offset*ones(size(eventTime));
allEpochEnds = allEpochStarts + epoch.duration*ones(size(eventTime));
% Filter out events which are too close to the beginning or end of trials
startFilter = allEpochStarts > 0;
stopFilter = round(allEpochEnds) < trialDuration; % Apply "round" to give an extra half second buffer and prevent indexing errors
evFilter = logical(startFilter.*stopFilter);
% Convert the starts from seconds to samples, round down to the nearest sample, coerce the value above 1.
epochStarts = max(floor(allEpochStarts(evFilter)*samplingRate),1);
% Calculate stops indices using the duration of the epoch, this avoids potential matrix dimension erros caused by differences in rounding when converting from seconds to samples.
sampleDur = round(epoch.duration*samplingRate);
epochStops = epochStarts + sampleDur*ones(size(epochStarts));
% Extract the chunk of data from the trial
chunkdata = zeros(sampleDur + 1,length(epochStarts),size(data,1));
for a = 1:length(epochStarts)
    chunkInds = epochStarts(a):epochStops(a);
    chunkdata(:,a,:) = data(:,chunkInds)';
end

end

function [temp] = AddEpochInfo(data,dataType,behavior,temp,fileID,fileDate,evFilter,f)
% Get the field names for each behavior
fields = fieldnames(data.Flags.(behavior));
% Filter out the events which are too close to the trial edge
for a = 1:length(fields)
    field = fields{a};
    temp.(dataType).(behavior).(field){f} = data.Flags.(behavior).(field)(evFilter,:)';
end
% Tag each event with the file ID, arrange cell array horizontally for
% later processing.
temp.(dataType).(behavior).fileIDs{f} = repmat({fileID},1,sum(evFilter));
temp.(dataType).(behavior).fileDates{f} = repmat({fileDate},1,sum(evFilter));

end

function [EventData] = ProcessTempStruct(EventData,dataType,temp,epoch)
% Get the dataTypes from temp
dTs = fieldnames(temp);
for a = 1:length(dTs)
    dT = dTs{a};
    % Get dataType names
    behaviorFields = fieldnames(temp.(dT));
    % Intialize Behavior fields of the dataType sub-structure
    structArray2 = cell(size(behaviorFields));
    EventData.(dataType).(dT) = cell2struct(structArray2,behaviorFields,1);
    for b = 1:length(behaviorFields)
        behavior = behaviorFields{b};
        % Get Behavior names
        eventFields = fieldnames(temp.(dT).(behavior));
        % Initialize Event fields for the Behavior sub-structure
        structArray3 = cell(size(eventFields));
        EventData.(dataType).(dT).(behavior) = cell2struct(structArray3,eventFields,1);
        for c = 1:length(eventFields)
            evField = eventFields{c};
            transferArray = [temp.(dT).(behavior).(evField){:}];
            EventData.(dataType).(dT).(behavior).(evField) = permute(transferArray,unique([2,1,ndims(transferArray)],'stable'));
        end
        EventData.(dataType).(dT).(behavior).epoch = epoch;
    end
end

end
