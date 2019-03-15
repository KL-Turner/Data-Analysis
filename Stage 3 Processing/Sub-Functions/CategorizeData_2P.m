function [] = CategorizeData_2P(fileName)
%___________________________________________________________________________________________________
% Edited by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%
%   Last Revised: August 8th, 2018
%___________________________________________________________________________________________________
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Identifies periods of sensory stimulation, volitional
%   movement, and rest. Calculates relevant details for each behavioral
%   period:
%           Stimulation:    Whisk Score - measure of the intensity of 
%                           whisking before and after onset of a puff. 
%                           A 0 indicates no whisking, a 1 indicates 
%                           maximum whisking over a 1 second period.
%                           Movement Score - Same as whisk score except
%                           uses the force sensor beneath the animal to
%                           detect body movment.
%
%           Whisking:       Duration - the time, in seconds, from onset to
%                           cessation of a whisking event.
%                           Rest Time - the duration of resting behavior
%                           prior to onset of the volitional whisk
%                           Whisk Score - a measure of the intensity of
%                           whisking for the duration of the whisk a
%                           maximum whisk for the whole duration will give
%                           a score of 1. No whisking will give a score of
%                           0.
%                           Movement Score - Same as whisk score exept uses
%                           the force sensor beneath the animal to detect
%                           body movement over the duration of the
%                           volitional whisk.
%                           Puff Distance - The time, in seconds, between
%                           the onset of each whisk an every puff
%                           administered during the trial.
%
%
%          Rest:            Duration - the time, in seconds without any 
%                           detected whisking or body movement.
%                           Start Time - the trial time corresponding to
%                           the cessation of all volitional movement.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                           filename - [string] file identifier                      
%_______________________________________________________________
%   RETURN:                     
%                           None, output of the script is additions to the
%                           MergedData structure.
%_______________________________________________________________

%% Load and Setup
disp(['Categorizing data for: ' fileName]); disp(' ')
load(fileName)
whiskerSamplingRate = MergedData.Notes.dsFs;

%% Process binary whisking waveform to detect whisking events
% Setup parameters for link_binary_events
linkThresh = 0.5;   % seconds, Link events < 0.5 seconds apart
breakThresh = 0;   % seconds changed by atw on 2/6/18 from 0.07

% Assume that whisks at the beginning/end of trial continue outside of the
% trial time. This will link any event occurring within "link_thresh"
% seconds to the beginning/end of the trial rather than assuming that it is
% a new/isolated event.
modBinWhiskers = MergedData.Data.binWhisker_Angle;

% Link the binarized whisking for use in GetWhiskingData function
binWhiskers = LinkBinaryEvents(gt(modBinWhiskers,0), [linkThresh breakThresh]*whiskerSamplingRate);

% Added 2/6/18 with atw. Code throws errors if binWhiskers(1)=1 and binWhiskers(2) = 0, or if 
% binWhiskers(1) = 0 and binWhiskers(2) = 1. This happens in GetWhiskingData because starts of 
% whisks are detected by taking the derivative of binWhiskers. Purpose of following lines is to 
% handle trials where the above conditions occur and avoid difficult dimension errors.
if binWhiskers(1)==0 && binWhiskers(2)==1
    binWhiskers(1) = 1;
elseif binWhiskers(1)==1 && binWhiskers(2)==0
    binWhiskers(1) = 0;
end

if binWhiskers(end)==0 && binWhiskers(end-1)==1
    binWhiskers(end) = 1;
elseif binWhiskers(end)==1 && binWhiskers(end-1)==0
    binWhiskers(end) = 0;
end

%% Categorize data by behavior
% Retrieve details on whisking events
[MergedData.Flags.whisk] = GetWhiskingData(MergedData, binWhiskers);

% Identify and separate resting data
[MergedData.Flags.rest] = GetRestData(MergedData);

% Save MergedData structure
save(fileName, 'MergedData');

function [puffTimes] = GetPuffTimes(~)
%   function [Puff_Times] = GetPuffTimes(MergedData)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Gets the time in seconds of all puffs administered during
%   a trial.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                       MergedData - [struct] structure obtained using the 
%                       function ProcessRawDataFile.
%_______________________________________________________________
%   RETURN:                     
%                       Puff_Times - [array] time in seconds of all puffs              
%_______________________________________________________________

puffList = [];
puffTimes = cell2mat(puffList);

function [Whisk] = GetWhiskingData(MergedData, binarizedWhiskers)
%   function [Whisk] = GetWhiskingData(MergedData, Bin_wwf)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Returns details on whisks which occurred during a trial.
%   Including:
%                           Duration - the time, in seconds, from onset to
%                           cessation of a whisking event.
%                           Rest Time - the duration of resting behavior
%                           prior to onset of the volitional whisk
%                           Whisk Score - a measure of the intensity of
%                           whisking for the duration of the whisk a
%                           maximum whisk for the whole duration will give
%                           a score of 1. No whisking will give a score of
%                           0.
%                           Movement Score - Same as whisk score exept uses
%                           the force sensor beneath the animal to detect
%                           body movement over the duration of the
%                           volitional whisk.
%                           Puff Distance - The time, in seconds, between
%                           the onset of each whisk an every puff
%                           administered during the trial.
%_______________________________________________________________
%   PARAMETERS:             
%                       MergedData - [struct] structure obtained using the 
%                       function ProcessRawDataFile.    
%_______________________________________________________________
%   RETURN:                     
%                       Whisk - [struct] structure containing a nested 
%                       structure for each whisk performed.
%_______________________________________________________________

%% Setup
whiskerSamplingRate = MergedData.Notes.dsFs;
forceSensorSamplingRate = MergedData.Notes.dsFs;

%% Get Puff Times
[puffTimes] = GetPuffTimes(MergedData);

%% Find the starts of whisking
whiskEdge = diff(binarizedWhiskers);
whiskSamples = find(whiskEdge > 0);
whiskStarts = whiskSamples/whiskerSamplingRate;

%% Classify each whisking event by duration, whisking intensity, rest durations
sampleVec = 1:length(binarizedWhiskers); 

% Identify periods of whisking/resting, include beginning and end of trial
% if needed (hence unique command) for correct interval calculation
highSamples = unique([1, sampleVec(binarizedWhiskers), sampleVec(end)]); 
lowSamples = unique([1, sampleVec(not(binarizedWhiskers)), sampleVec(end)]);

% Calculate the number of samples between consecutive high/low samples.
dHigh = diff(highSamples);
dLow = diff(lowSamples);

% Identify skips in sample numbers which correspond to rests/whisks,
% convert from samples to seconds.
restLength = dHigh(dHigh > 1);
whiskLength = dLow(dLow > 1);
restDur = restLength/whiskerSamplingRate;
whiskDur = whiskLength/whiskerSamplingRate;

% Control for the beginning/end of the trial to correctly map rests/whisks
% onto the whisk_starts.
if binarizedWhiskers(1)
    whiskDur(1) = [];
    whiskLength(1) = [];
end

if not(binarizedWhiskers(end))
    restDur(end) = [];
end


% Calculate the whisking intensity -> sum(MergedData.Bin_wwf)/sum(Bin_wwf)
% over the duration of the whisk. Calculate the movement intensity over the
% same interval.
whiskInt = zeros(size(whiskStarts));
movementInt = zeros(size(whiskStarts));

for wS = 1:length(whiskSamples)
    % Whisking intensity
    whiskInds = whiskSamples(wS):whiskSamples(wS) + whiskLength(wS);
    whiskInt(wS) = sum(MergedData.Data.binWhisker_Angle(whiskInds))/numel(whiskInds);
    
    % Movement intensity
    movementStart = round(whiskStarts(wS)*forceSensorSamplingRate);
    movementDur = round(whiskDur(wS)*forceSensorSamplingRate);
    movementInds = max(movementStart, 1):min(movementStart + movementDur, length(MergedData.Data.binForce_Sensor_M));
    movementInt(wS) = sum(MergedData.Data.binForce_Sensor_M(movementInds))/numel(movementInds);
end

% Calculate the time to the closest puff
% If no puff occurred during the trial, store 0 as a place holder.
if isempty(puffTimes)
    puffTimes = 0;
end

puffMat = ones(length(whiskSamples), 1)*puffTimes;
whiskMat = whiskSamples'*ones(1, length(puffTimes)) / whiskerSamplingRate;
puffTimeElapsed = abs(whiskMat - puffMat);

% Convert to cell
puffTimeCell = mat2cell(puffTimeElapsed, ones(length(whiskStarts), 1));

%% Error handle
if length(restDur) ~= length(whiskDur)
    disp('Error in GetWhiskData! The number of whisks does not equal the number of rests...')
    disp(' ')
    keyboard;
end
%% Compile into final structure

Whisk.eventTime = whiskStarts';
Whisk.duration = whiskDur';
Whisk.restTime = restDur';
Whisk.whiskScore = whiskInt';
Whisk.movementScore = movementInt';
Whisk.puffDistance = puffTimeCell;

function [Rest] = GetRestData(MergedData)
%   function [Rest] = GetRestData(MergedData)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Returns details on periods of rest during a trial.
%   Including:
%          Rest:            Duration - the time, in seconds without any 
%                           detected whisking or body movement.
%                           Start Time - the trial time corresponding to
%                           the cessation of all volitional movement.   
%_______________________________________________________________
%   PARAMETERS:             
%                       MergedData - [struct] structure obtained using the 
%                       function ProcessRawDataFile.    
%_______________________________________________________________
%   RETURN:                     
%                       Rest - [struct] structure containing a nested 
%                       structure for each period of rest.
%_______________________________________________________________

% Setup
whiskerSamplingRate = MergedData.Notes.dsFs;
forceSensorSamplingRate = MergedData.Notes.dsFs;

%% Get stimulation times
[puffTimes] = GetPuffTimes(MergedData);

%% Recalculate linked binarized wwf without omitting any possible whisks,
% this avoids inclusion of brief whisker movements in periods of rest.

% Assume that whisks at the beginning/end of trial continue outside of the
% trial time. This will link any event occurring within "link_thresh"
% seconds to the beginning/end of the trial rather than assuming that it is
% a new/isolated event.
modBinarizedWhiskers = MergedData.Data.binWhisker_Angle;
modBinarizedWhiskers([1, end]) = 1;

modBinarizedForceSensor = MergedData.Data.binForce_Sensor_M;
modBinarizedForceSensor([1, end]) = 1;

linkThresh = 0.5; % seconds
breakThresh = 0;% seconds
binarizedWhiskers = LinkBinaryEvents(gt(modBinarizedWhiskers, 0),...
    [linkThresh breakThresh]*whiskerSamplingRate);
binarizedForceSensor = LinkBinaryEvents(modBinarizedForceSensor,...
    [linkThresh breakThresh]*forceSensorSamplingRate);

%% Combine binarizedWhiskers, binarizedForceSensor, and puffTimes, to find periods of rest. 

% Downsample bin_wwf to match length of bin_pswf
sampleVec = 1:length(binarizedWhiskers); 
whiskHigh = sampleVec(binarizedWhiskers)/whiskerSamplingRate;
dsBinarizedWhiskers = zeros(size(binarizedForceSensor));

% Find Bin_wwf == 1. Convert indexes into pswf time. Coerce converted indexes
% between 1 and length(Bin_pswf). Take only unique values.
dsInds = min(max(round(whiskHigh*forceSensorSamplingRate), 1), length(binarizedForceSensor));
dsBinarizedWhiskers(unique(dsInds)) = 1;

% Combine binarized whisking and body movement
wfBin = logical(min(dsBinarizedWhiskers + binarizedForceSensor, 1));
Fs = forceSensorSamplingRate;

% Add puff times into the Bin_wf
puffInds = round(puffTimes*Fs);
wfBin(puffInds) = 1;

% Find index for end of whisking event
edge = diff(wfBin);
samples = find([not(wfBin(1)) edge < 0]);
stops = samples / Fs;

% Identify periods of whisking/resting, include beginning and end of trial
% if needed (hence unique command) for correct interval calculation
sampleVec = 1:length(logical(wfBin));
highSamples = unique([1, sampleVec(wfBin), sampleVec(end)]); 
lowSamples = unique([1, sampleVec(not(wfBin)), sampleVec(end)]); 

% Calculate the number of samples between consecutive high/low samples.
dHigh = diff(highSamples);
dLow = diff(lowSamples);

% Identify skips in sample numbers which correspond to rests/whisks,
% convert from samples to seconds.
restLength = dHigh(dHigh > 1);
restDur = restLength / Fs;
whiskLength = dLow(dLow > 1);
whiskDur = whiskLength / Fs;

% Control for the beginning/end of the trial to correctly map rests/whisks
% onto the whisk_starts. Use index 2 and end-1 since it is assumed that the
% first and last indexes of a trial are the end/beginning of a volitional
% movement.
if not(wfBin(2)) 
    whiskDur = [NaN whiskDur];
end

if wfBin(end - 1)
    whiskDur(end) = [];
end

% Calculate the time to the closest puff
% If no puff occurred during the trial, store 0 as a place holder.
if isempty(puffTimes)
    puffTimes = 0;
end

puffMat = ones(length(samples), 1)*puffTimes;
restMat = samples'*ones(1, length(puffTimes)) / Fs;
puffTimeElapsed = abs(restMat - puffMat);

% Convert to cell
puffTimeCell = mat2cell(puffTimeElapsed, ones(length(samples), 1));

% Compile into a structure
Rest.eventTime = stops';
Rest.duration = restDur';
Rest.puffDistance = puffTimeCell;
Rest.whiskDuration = whiskDur';
