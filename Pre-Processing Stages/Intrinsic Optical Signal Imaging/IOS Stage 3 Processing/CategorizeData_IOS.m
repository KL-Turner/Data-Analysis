function [] = CategorizeData_IOS(procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%
% Purpose: Catagorizes data based on behavioral flags from whisking/movement events
%________________________________________________________________________________________________________________________

for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    % load and Setup
    disp(['Categorizing data for: ' procDataFileID]); disp(' ')
    load(procDataFileID)
    % condense pulsed stimulations
    stimulationFields = fieldnames(ProcData.data.stimulations);
    for vv = 1:length(stimulationFields)
        stimulationTimes = ProcData.data.stimulations.(stimulationFields{vv,1});
        ProcData.data.stimulationsOriginal.(stimulationFields{vv,1}) = stimulationTimes;
        condensedStimulationTimes = [];
        cc = 1;
        if isempty(stimulationTimes) == false
            for bb = 1:length(stimulationTimes)
                if bb == 1
                    condensedStimulationTimes(1,bb) = stimulationTimes(1,bb);
                    cc = cc + 1;
                else
                    timeDifference = stimulationTimes(1,bb) - stimulationTimes(1,bb - 1);
                    if timeDifference > 1 % remove stimulations that are closer than 1 second to the previous
                        condensedStimulationTimes(1,cc) = stimulationTimes(1,bb);
                        cc = cc + 1;
                    end
                end
            end
            ProcData.data.stimulations.(stimulationFields{vv,1}) = condensedStimulationTimes;
        end
    end
    whiskerSamplingRate = ProcData.notes.dsFs;
    linkThresh = 0.5; % link events < 0.5 seconds apart
    breakThresh = 0;
    modBinWhiskers = ProcData.data.binWhiskerAngle;
    % link the binarized whisking
    binWhiskers = LinkBinaryEvents_IOS(gt(modBinWhiskers,0),[linkThresh breakThresh]*whiskerSamplingRate);
    % handle edge conditions
    if binWhiskers(1) == 0 && binWhiskers(2) == 1
        binWhiskers(1) = 1;
    elseif binWhiskers(1) == 1 && binWhiskers(2) == 0
        binWhiskers(1) = 0;
    end
    % handle edge conditions
    if binWhiskers(end) == 0 && binWhiskers(end - 1) == 1
        binWhiskers(end) = 1;
    elseif binWhiskers(end) == 1 && binWhiskers(end - 1) == 0
        binWhiskers(end) = 0;
    end
    % retrieve details on whisking events
    [ProcData.flags.whisk] = GetWhiskingdata_IOS(ProcData,binWhiskers);
    % retrieve details on stiming events
    [ProcData.flags.stim] = GetStimdata_IOS(ProcData);
    % identify and separate resting data
    [ProcData.flags.rest] = GetRestdata(ProcData);
    % Save ProcData structure
    save(procDataFileID,'ProcData');
end

end

function [stimTimes] = GetStimTimes_IOS(ProcData)
% get solenoid times
stimNames = fieldnames(ProcData.data.stimulations);
stimList = cell(1,length(stimNames));
for sN = 1:length(stimNames)
    stimList{sN} = ProcData.data.stimulations.(stimNames{sN});
end
stimTimes = cell2mat(stimList);

end

function [Stim] = GetStimdata_IOS(ProcData)
% get stimulation times
whiskerSamplingRate = ProcData.notes.dsFs;
forceSensorSamplingRate = ProcData.notes.dsFs;
stimTimes = GetStimTimes_IOS(ProcData);
trialDuration = ProcData.notes.trialDuration_sec;
% set time intervals for calculation of the whisk scores
preTime = 1;
postTime = 1;
% get stimer IDs
stimNames = fieldnames(ProcData.data.stimulations);
Stim.solenoidName = cell(length(stimTimes),1);
Stim.eventTime = zeros(length(stimTimes),1);
Stim.whiskScore_Pre = zeros(length(stimTimes),1);
Stim.whiskScore_Post = zeros(length(stimTimes),1);
Stim.movementScore_Pre = zeros(length(stimTimes),1);
Stim.movementScore_Post = zeros(length(stimTimes),1);
j = 1;
for sN = 1:length(stimNames)
    solStimTimes = ProcData.data.stimulations.(stimNames{sN});
    for spT = 1:length(solStimTimes)
        if trialDuration - solStimTimes(spT) <= postTime
            disp(['Stim at time: ' solStimTimes(spT) ' is too close to trial end'])
            continue;
        end
        % set indexes for pre and post periods
        wStimInd = round(solStimTimes(spT)*whiskerSamplingRate);
        mStimInd = round(solStimTimes(spT)*forceSensorSamplingRate);
        wPreStart = max(round((solStimTimes(spT) - preTime)*whiskerSamplingRate),1);
        mPreStart = max(round((solStimTimes(spT) - preTime)*forceSensorSamplingRate),1);
        wPostEnd = round((solStimTimes(spT) + postTime)*whiskerSamplingRate);
        mPostEnd = round((solStimTimes(spT) + postTime)*forceSensorSamplingRate);
        % calculate the percent of the pre-stim time that the animal moved or whisked
        whiskScorePre = sum(ProcData.data.binWhiskerAngle(wPreStart:wStimInd))/(preTime*whiskerSamplingRate);
        whiskScorePost = sum(ProcData.data.binWhiskerAngle(wStimInd:wPostEnd))/(postTime*whiskerSamplingRate);
        moveScorePre = sum(ProcData.data.binForceSensor(mPreStart:mStimInd))/(preTime*forceSensorSamplingRate);
        moveScorePost = sum(ProcData.data.binForceSensor(mStimInd:mPostEnd))/(postTime*forceSensorSamplingRate);
        % add to Stim structure
        Stim.solenoidName{j} = stimNames{sN};
        Stim.eventTime(j) = solStimTimes(spT)';
        Stim.whiskScore_Pre(j) = whiskScorePre';
        Stim.whiskScore_Post(j) = whiskScorePost';
        Stim.movementScore_Pre(j) = moveScorePre';
        Stim.movementScore_Post(j) = moveScorePost';
        j = j + 1;
    end
end
% calculate the time to the closest stim, omit comparison of stim to itself
stimMat = ones(length(stimTimes),1)*stimTimes;
timeElapsed = abs(nonzeros(stimMat - stimMat'));
% if no other stim occurred during the trial, store 0 as a place holder.
if isempty(timeElapsed)
    stimTimeElapsed = 0;
else
    % if not empty, Reshape the array to compensate for nonzeros command
    stimTimeElapsed = reshape(timeElapsed,numel(stimTimes) - 1,numel(stimTimes));
end
% convert to cell and add to struct, if length of Stim_Times = 0, coerce to
% 1 to accommodate the NaN entry.
stimTimeCell = mat2cell(stimTimeElapsed',ones(max(length(stimTimes),1),1));
Stim.StimDistance = stimTimeCell;
end

function [Whisk] = GetWhiskingdata_IOS(ProcData,binWhiskerAngle)
% setup
whiskerSamplingRate = ProcData.notes.dsFs;
forceSensorSamplingRate = ProcData.notes.dsFs;
% get Stim Times
[stimTimes] = GetStimTimes_IOS(ProcData);
% find the starts of whisking
whiskEdge = diff(binWhiskerAngle);
whiskSamples = find(whiskEdge > 0);
whiskStarts = whiskSamples/whiskerSamplingRate;
% classify each whisking event by duration, whisking intensity, rest durations
sampleVec = 1:length(binWhiskerAngle);
% identify periods of whisking/resting, include beginning and end of trial
% if needed (hence unique command) for correct interval calculation
highSamples = unique([1,sampleVec(binWhiskerAngle),sampleVec(end)]);
lowSamples = unique([1,sampleVec(not(binWhiskerAngle)),sampleVec(end)]);
% calculate the number of samples between consecutive high/low samples.
dHigh = diff(highSamples);
dLow = diff(lowSamples);
% identify skips in sample numbers which correspond to rests/whisks,
% convert from samples to seconds.
restLength = dHigh(dHigh > 1);
whiskLength = dLow(dLow > 1);
restDur = restLength/whiskerSamplingRate;
whiskDur = whiskLength/whiskerSamplingRate;
% control for the beginning/end of the trial to correctly map rests/whisks
% onto the whisk_starts.
if binWhiskerAngle(1)
    whiskDur(1) = [];
    whiskLength(1) = [];
end
if not(binWhiskerAngle(end))
    restDur(end) = [];
end
% calculate the whisking intensity -> sum(ProcData.Bin_wwf)/sum(Bin_wwf)
% over the duration of the whisk. Calculate the movement intensity over the same interval.
whiskInt = zeros(size(whiskStarts));
movementInt = zeros(size(whiskStarts));
for wS = 1:length(whiskSamples)
    % whisking intensity
    whiskInds = whiskSamples(wS):whiskSamples(wS) + whiskLength(wS);
    whiskInt(wS) = sum(ProcData.data.binWhiskerAngle(whiskInds))/numel(whiskInds);
    % movement intensity
    movementStart = round(whiskStarts(wS)*forceSensorSamplingRate);
    movementDur = round(whiskDur(wS)*forceSensorSamplingRate);
    movementInds = max(movementStart,1):min(movementStart + movementDur,length(ProcData.data.binForceSensor));
    movementInt(wS) = sum(ProcData.data.binForceSensor(movementInds))/numel(movementInds);
end
% calculate the time to the closest stim
% if no stim occurred during the trial, store 0 as a place holder.
if isempty(stimTimes)
    stimTimes = 0;
end
stimMat = ones(length(whiskSamples),1)*stimTimes;
whiskMat = whiskSamples'*ones(1,length(stimTimes))/whiskerSamplingRate;
stimTimeElapsed = abs(whiskMat - stimMat);
% convert to cell
stimTimeCell = mat2cell(stimTimeElapsed,ones(length(whiskStarts),1));
% error handle
if length(restDur) ~= length(whiskDur)
    disp('Error in GetWhiskdata! The number of whisks does not equal the number of rests...'); disp(' ')
    keyboard;
end
% compile into final structure
Whisk.eventTime = whiskStarts';
Whisk.duration = whiskDur';
Whisk.restTime = restDur';
Whisk.whiskScore = whiskInt';
Whisk.movementScore = movementInt';
Whisk.stimDistance = stimTimeCell;
end

function [Rest] = GetRestdata(ProcData)
% setup
whiskerSamplingRate = ProcData.notes.dsFs;
forceSensorSamplingRate = ProcData.notes.dsFs;
% get stimulation times
[stimTimes] = GetStimTimes_IOS(ProcData);
% recalculate linked binarized wwf without omitting any possible whisks
modBinarizedWhiskers = ProcData.data.binWhiskerAngle;
modBinarizedWhiskers([1,end]) = 1;
modBinarizedForceSensor = ProcData.data.binForceSensor;
modBinarizedForceSensor([1,end]) = 1;
linkThresh = 0.5; % seconds
breakThresh = 0; % seconds
binWhiskerAngle = LinkBinaryEvents_IOS(gt(modBinarizedWhiskers,0),[linkThresh breakThresh]*whiskerSamplingRate);
binForceSensor = LinkBinaryEvents_IOS(modBinarizedForceSensor,[linkThresh breakThresh]*forceSensorSamplingRate);
% combine binWhiskerAngle, binForceSensor, and stimTimes, to find periods of rest
sampleVec = 1:length(binWhiskerAngle);
whiskHigh = sampleVec(binWhiskerAngle)/whiskerSamplingRate;
dsBinarizedWhiskers = zeros(size(binForceSensor));
% find Bin_wwf == 1. Convert indexes into pswf time. Coerce converted indexes
% between 1 and length(Bin_pswf). Take only unique values.
dsInds = min(max(round(whiskHigh*forceSensorSamplingRate),1),length(binForceSensor));
dsBinarizedWhiskers(unique(dsInds)) = 1;
% combine binarized whisking and body movement
wfBin = logical(min(dsBinarizedWhiskers + binForceSensor,1));
samplingRate = forceSensorSamplingRate;
% add stim times into the Bin_wf
stimInds = round(stimTimes*samplingRate);
wfBin(stimInds) = 1;
% find index for end of whisking event
edge = diff(wfBin);
samples = find([not(wfBin(1)),edge < 0]);
stops = samples/samplingRate;
% identify periods of whisking/resting, include beginning and end of trial
% if needed (hence unique command) for correct interval calculation
sampleVec = 1:length(logical(wfBin));
highSamples = unique([1,sampleVec(wfBin),sampleVec(end)]);
lowSamples = unique([1,sampleVec(not(wfBin)),sampleVec(end)]);
% calculate the number of samples between consecutive high/low samples.
dHigh = diff(highSamples);
dLow = diff(lowSamples);
% identify skips in sample numbers which correspond to rests/whisks,
% convert from samples to seconds.
restLength = dHigh(dHigh > 1);
restDur = restLength/samplingRate;
whiskLength = dLow(dLow > 1);
whiskDur = whiskLength/samplingRate;
% control for the beginning/end of the trial to correctly map rests/whisks
if not(wfBin(2))
    whiskDur = [NaN,whiskDur];
end
% control for the beginning/end of the trial to correctly map rests/whisks
if wfBin(end - 1)
    whiskDur(end) = [];
end
% calculate the time to the closest stim
if isempty(stimTimes)
    stimTimes = 0;
end
stimMat = ones(length(samples),1)*stimTimes;
restMat = samples'*ones(1,length(stimTimes))/samplingRate;
stimTimeElapsed = abs(restMat - stimMat);
% convert to cell
stimTimeCell = mat2cell(stimTimeElapsed,ones(length(samples),1));
% compile into a structure
Rest.eventTime = stops';
Rest.duration = restDur';
Rest.stimDistance = stimTimeCell;
Rest.whiskDuration = whiskDur';

end
