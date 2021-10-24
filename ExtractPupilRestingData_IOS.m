function [RestData] = ExtractPupilRestingData_IOS(procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
% Purpose: Extracts all resting data periods from the data using behavioral flags
%________________________________________________________________________________________________________________________

% load rest data file
restDataFileID = ls('*_RestData.mat');
load(restDataFileID)
% analyze each proc data file
zz = 1;
for c = 1:size(procDataFileIDs,1)
    disp(['Extracting resting pupil area from ProcData file ' num2str(c) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    procDataFileID = procDataFileIDs(c,:);
    load(procDataFileID);
    if strcmp(ProcData.data.Pupil.frameCheck,'y') == true
        % get the date and file identifier for the data to be saved with each resting event
        [animal,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
        % sampling frequency for element of dataTypes
        Fs = ProcData.notes.pupilCamSamplingRate;
        % expected number of samples for element of dataType
        trialDuration_sec = ProcData.notes.trialDuration_sec;
        expectedLength = trialDuration_sec*Fs;
        % get information about periods of rest from the loaded file
        trialEventTimes = ProcData.flags.rest.eventTime';
        trialPuffDistances = ProcData.flags.rest.puffDistance;
        trialDurations = ProcData.flags.rest.duration';
        % initialize cell array for all periods of rest from the loaded file
        trialRestVals = cell(size(trialEventTimes'));
        for d = 1:length(trialEventTimes)
            % extract the whole duration of the resting event. Coerce the
            % start index to values above 1 to preclude rounding to 0.
            startInd = max(floor(trialEventTimes(d)*Fs),1);
            % convert the duration from seconds to samples.
            dur = round(trialDurations(d)*Fs);
            % get ending index for data chunk. If event occurs at the end of
            % the trial, assume animal whisks as soon as the trial ends and
            % give a 200ms buffer.
            stopInd = min(startInd + dur,expectedLength - round(0.2*Fs));
            try
                % extract data from the trial and add to the cell array for the current loaded file
                trialRestVals{d} = ProcData.data.Pupil.pupilArea(:,startInd:stopInd);
            catch % some files don't have certain fields. Skip those
                trialRestVals{d} = [];
            end
        end
        % add all periods of rest to a cell array for all files
        restVals{zz,1} = trialRestVals'; %#ok<*AGROW>
        % transfer information about resting periods to the new structure
        eventTimes{zz,1} = trialEventTimes';
        durations{zz,1} = trialDurations';
        puffDistances{zz,1} = trialPuffDistances';
        fileIDs{zz,1} = repmat({fileID},1,length(trialEventTimes));
        fileDates{zz,1} = repmat({fileDate},1,length(trialEventTimes));
        zz = zz + 1;
    end
end
% combine the cells from separate files into a single cell array of all resting periods
RestData.Pupil.pupilArea.data = [restVals{:}]';
RestData.Pupil.pupilArea.eventTimes = cell2mat(eventTimes);
RestData.Pupil.pupilArea.durations = cell2mat(durations);
RestData.Pupil.pupilArea.puffDistances = [puffDistances{:}]';
RestData.Pupil.pupilArea.fileIDs = [fileIDs{:}]';
RestData.Pupil.pupilArea.fileDates = [fileDates{:}]';
RestData.Pupil.pupilArea.CBVCamSamplingRate = Fs;
RestData.Pupil.pupilArea.trialDuration_sec = trialDuration_sec;
% save updated structure
save([animal '_RestData.mat'],'RestData','-v7.3');

end
