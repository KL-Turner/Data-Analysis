function ProcessRawDataFile(fileName)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%________________________________________________________________________________________________________________________
%
%   Purpose: Saves processed forms of measured data to the _RawData.mat file in order to reduce
%            run times in later scripts and increase consistency. No data in the appended data 
%            will be normalized. Run CreateRawDataStructure_parallel.m to create a RawData
%            structure. Note: That code was renamed to StageOneProcessing.m
%________________________________________________________________________________________________________________________
%
%   Inputs: filename: [str] name of _RawData.mat file to be processed.
%
%   Outputs: Function output is a saved structure called 'ProcData.mat'
%
%   Last Revised: August 8th, 2018
%________________________________________________________________________________________________________________________

disp(['Processing RawData File: ' fileName]); disp(' ')
load(fileName);
underscoreIndexes = strfind(fileName,'_');
[animal, ~, fileDate, fileID] = GetFileInfo(fileName);
strDay = ConvertDate(fileDate);

%% Load in ProcData notes from RawData structure:
ProcData.Notes.experimenter = RawData.Notes.experimenter;
ProcData.Notes.animalID = RawData.Notes.animalID;
ProcData.Notes.imagedHemisphere = RawData.Notes.imagedHemisphere;
ProcData.Notes.solenoidPressure_PSI = RawData.Notes.solenoidPressure_PSI;
ProcData.Notes.isofluraneTime_Military = RawData.Notes.isofluraneTime_Military;
ProcData.Notes.sessionID = RawData.Notes.sessionID;
ProcData.Notes.amplifierGain = RawData.Notes.amplifierGain;
ProcData.Notes.CBVCamSamplingRate = RawData.Notes.CBVCamSamplingRate;
ProcData.Notes.whiskerCamSamplingRate = RawData.Notes.whiskerCamSamplingRate;
ProcData.Notes.webCamSamplingRate = RawData.Notes.webCamSamplingRate;
ProcData.Notes.analogSamplingRate = RawData.Notes.analogSamplingRate;
ProcData.Notes.trialDuration_Seconds = RawData.Notes.trialDuration_Seconds;
ProcData.Notes.CBVCamPixelHeight = RawData.Notes.CBVCamPixelHeight;
ProcData.Notes.CBVCamPixelWidth = RawData.Notes.CBVCamPixelWidth;
ProcData.Notes.CBVCamBitDepth = RawData.Notes.CBVCamBitDepth;
ProcData.Notes.whiskerCamPixelHeight = RawData.Notes.whiskerCamPixelHeight;
ProcData.Notes.whiskerCamPixelWidth = RawData.Notes.whiskerCamPixelWidth;
ProcData.Notes.CBVCamExposureTime_Microseconds = RawData.Notes.CBVCamExposureTime_Microseconds;
ProcData.Notes.CBVCamBinning = RawData.Notes.CBVCamBinning;
ProcData.Notes.numberDroppedWhiskerCamFrames = RawData.Notes.numberDroppedWhiskerCamFrames;
ProcData.Notes.droppedWhiskerCamFrameIndex = RawData.Notes.droppedWhiskerCamFrameIndex;

%% Save solenoid times (in seconds).
% Identify the solenoids by amplitude
ProcData.Data.Sol.solenoidLeftPad = find(diff(RawData.Data.Solenoids) == 1) / RawData.Notes.analogSamplingRate;
ProcData.Data.Sol.solenoidRightPad = find(diff(RawData.Data.Solenoids) == 2) / RawData.Notes.analogSamplingRate;
ProcData.Data.Sol.solenoidTail = find(diff(RawData.Data.Solenoids) == 3) / RawData.Notes.analogSamplingRate;
ProcData.Data.Sol.solenoidAuditory = find(diff(RawData.Data.Solenoids) == 4) / RawData.Notes.analogSamplingRate;

%% CBV from ROIs.
CBVfields = fieldnames(RawData.Data.CBV);

for field = 1:length(CBVfields)
    ProcData.Data.CBV.([CBVfields{field}]) = RawData.Data.CBV.(CBVfields{field})(1:end - 1);
end

%% Process neural data into its various forms.
% Wide-Band LFP
[ProcData.Data.WideBandLFP_Power.LH, ProcData.Notes.wideBandLFPSamplingRate] = ...
    ProcessNeuro(RawData, 'Wideband_LFP', 'Neural_LH');

[ProcData.Data.WideBandLFP_Power.RH, ProcData.Notes.wideBandLFPSamplingRate] = ...
    ProcessNeuro(RawData, 'Wideband_LFP', 'Neural_RH');

% MUA Band [300 - 3000]
[ProcData.Data.MUA_Power.LH, ProcData.Notes.multiUnitSamplingRate] = ...
    ProcessNeuro(RawData, 'MUpower', 'Neural_LH');

[ProcData.Data.MUA_Power.RH, ProcData.Notes.multiUnitSamplingRate] = ...
    ProcessNeuro(RawData, 'MUpower', 'Neural_RH');

% Gamma Band [40 - 100]
[ProcData.Data.GammaBand_Power.LH, ProcData.Notes.gammaBandSamplingRate] = ...
    ProcessNeuro(RawData, 'Gam', 'Neural_LH');

[ProcData.Data.GammaBand_Power.RH, ProcData.Notes.gammaBandSamplingRate] = ...
    ProcessNeuro(RawData, 'Gam', 'Neural_RH');

% Theta [4 - 8 Hz]
[ProcData.Data.ThetaBand_Power.LH, ProcData.Notes.thetaBandSamplingRate] = ...
    ProcessNeuro(RawData, 'Theta', 'Neural_LH');

[ProcData.Data.ThetaBand_Power.RH, ProcData.Notes.thetaBandSamplingRate] = ...
    ProcessNeuro(RawData, 'Theta', 'Neural_RH');

% Delta [1 - 4 Hz]
[ProcData.Data.DeltaBand_Power.LH, ProcData.Notes.deltaBandSamplingRate] = ...
    ProcessNeuro(RawData, 'Delta', 'Neural_LH');

[ProcData.Data.DeltaBand_Power.RH, ProcData.Notes.deltaBandSamplingRate] = ...
    ProcessNeuro(RawData, 'Delta', 'Neural_RH');

%% Binarize the whisker angle and set the resting angle to zero degrees.
% Track dropped whisker camera frames
expectedLength = RawData.Notes.trialDuration_Seconds*RawData.Notes.analogSamplingRate;
expectedFrames = RawData.Notes.trialDuration_Seconds*RawData.Notes.whiskerCamSamplingRate;
ProcData.Notes.totalDroppedWhiskerFrames = expectedFrames - length(RawData.Data.WhiskerAngle);

% Trim any additional frames for resample
whiskerAngle = RawData.Data.WhiskerAngle(1:min(expectedFrames, length(RawData.Data.WhiskerAngle)));

% Create filter for whisking/movement
whiskerDownsampledSamplingRate = RawData.Notes.CBVCamSamplingRate;   % Downsample to CBV Camera Fs
whiskerFilterThreshold = 20;
whiskerFilterOrder = 2;
[z, p, k] = butter(whiskerFilterOrder, whiskerFilterThreshold / (RawData.Notes.whiskerCamSamplingRate / 2), 'low');
[sos, g] = zp2sos(z, p, k);
filteredWhiskers = filtfilt(sos, g, whiskerAngle - mean(whiskerAngle));
resampledWhiskers = resample(filteredWhiskers, whiskerDownsampledSamplingRate, RawData.Notes.whiskerCamSamplingRate);

% Binarize the whisker waveform (wwf)
threshfile = dir('*_Thresholds.mat');
if ~isempty(threshfile)
    load(threshfile.name)
end

[ok] = CheckForThreshold(['binarizedWhiskersLower_' strDay], animal);

if ok == 0
    [whiskersThresh1, whiskersThresh2] = CreateWhiskThreshold(resampledWhiskers, whiskerDownsampledSamplingRate);
    Thresholds.(['binarizedWhiskersLower_' strDay]) = whiskersThresh1;
    Thresholds.(['binarizedWhiskersUpper_' strDay]) = whiskersThresh2;
    save([animal '_Thresholds.mat'], 'Thresholds');
end

load([animal '_Thresholds.mat']);
binarizedWhiskers = BinarizeWhiskers(resampledWhiskers, whiskerDownsampledSamplingRate, Thresholds.(['binarizedWhiskersLower_' strDay]), Thresholds.(['binarizedWhiskersUpper_' strDay]));
[linkedBinarizedWhiskers] = LinkBinaryEvents(gt(binarizedWhiskers,0), [round(whiskerDownsampledSamplingRate/3), 0]);

inds = linkedBinarizedWhiskers == 0;
restAngle = mean(resampledWhiskers(inds));

ProcData.Data.Behavior.whiskers = resampledWhiskers - restAngle;
ProcData.Notes.downsampledWhiskerSamplingRate = whiskerDownsampledSamplingRate;
ProcData.Data.Behavior.binarizedWhiskers = binarizedWhiskers;
ProcData.Notes.binarizedWhiskerSamplingRate = whiskerDownsampledSamplingRate;

%% Downsample and binarize the force sensor.
% Trim any additional data points for resample
trimmedForceSensor = RawData.Data.Force_Sensor(1:min(expectedLength, length(RawData.Data.Force_Sensor)));

% Filter then downsample the Force Sensor waveform to desired frequency
forceSensorDownSampledSamplingRate = RawData.Notes.CBVCamSamplingRate;   % Downsample to CBV Camera Fs
forceSensorFilterThreshold = 20;
forceSensorFilterOrder = 2;
[z, p, k] = butter(forceSensorFilterOrder, forceSensorFilterThreshold / (RawData.Notes.analogSamplingRate / 2), 'low');
[sos, g] = zp2sos(z, p, k);
filteredForceSensor = filtfilt(sos, g, trimmedForceSensor);

ProcData.Data.Behavior.forceSensor = resample(filteredForceSensor, forceSensorDownSampledSamplingRate, RawData.Notes.analogSamplingRate);
ProcData.Notes.downsampledForceSensorSamplingRate = forceSensorDownSampledSamplingRate;

% Binarize the force sensor waveform
[ok] = CheckForThreshold(['binarizedForceSensor_' strDay], animal);

if ok == 0
    [forceSensorThreshold] = CreateForceSensorThreshold(ProcData.Data.Behavior.forceSensor);
    Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
    save([animal '_Thresholds.mat'], 'Thresholds');
end

ProcData.Data.Behavior.binarizedForceSensor = BinarizeForceSensor(ProcData.Data.Behavior.forceSensor, Thresholds.(['binarizedForceSensor_' strDay]));
ProcData.Notes.binarizedForceSensorSamplingRate = forceSensorDownSampledSamplingRate;

%% EMG.
if isfield(RawData.Data, 'EMG')
    %% Downsample and binarize the EMG.
    % Trim any additional data points for resample
    trimmedEMG = RawData.Data.EMG(1:min(expectedLength, length(RawData.Data.EMG)));
    
    % Filter then downsample the Force Sensor waveform to desired frequency
    emgDownSampledSamplingRate = RawData.Notes.CBVCamSamplingRate;   % Downsample to CBV Camera Fs
    emgFilterThreshold = 20;
    emgFilterOrder = 2;
    [z, p, k] = butter(emgFilterOrder, emgFilterThreshold / (RawData.Notes.analogSamplingRate / 2), 'low');
    [sos, g] = zp2sos(z, p, k);
    filteredEMG = filtfilt(sos, g, trimmedEMG);
    
    ProcData.Data.Behavior.EMG = resample(filteredEMG, emgDownSampledSamplingRate, RawData.Notes.analogSamplingRate);
    ProcData.Notes.downsampledEMGSamplingRate = emgDownSampledSamplingRate;
end

%% Save the processed data structure.
disp(['Saving ProcData file ' fileID '...']); disp(' ')
save([fileName(1:(underscoreIndexes(end) - 1)) '_ProcData.mat'], 'ProcData');

end