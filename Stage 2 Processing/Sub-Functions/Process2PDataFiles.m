function Process2PDataFiles(labviewDataFiles, mscanDataFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs: 
%
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

%% MScan data file analysis
for f = 1:size(mscanDataFiles, 1)
    %% Find offset between the two force sensor signals using the cross correlation
    disp(['Analyzing MScan neural bands and analog signals for file number ' num2str(f) ' of ' num2str(size(mscanDataFiles, 1)) '...']); disp(' '); 
    MScanDataFile = mscanDataFiles(f, :);
    load(MScanDataFile);
    animalID = MScanData.Notes.animalID;
    imageID = MScanData.Notes.imageID;
    date = MScanData.Notes.date;
    strDay = ConvertDate(date);

    %% Process neural data into its various forms.
    % MUA Band [300 - 3000]
    [MScanData.Data.MUA_Power, MScanData.Notes.MScan.multiUnitSamplingRate] = ...
        ProcessNeuro_2P(MScanData, 'MUApower', 'MScan_Neural_Data');

    % Gamma Band [40 - 100]
    [MScanData.Data.GammaBand_Power, MScanData.Notes.MScan.gammaBandSamplingRate] = ...
        ProcessNeuro_2P(MScanData, 'Gam', 'MScan_Neural_Data');

    % Beta [13 - 30 Hz]
    [MScanData.Data.BetaBand_Power, MScanData.Notes.MScan.betaBandSamplingRate] = ...
        ProcessNeuro_2P(MScanData, 'Beta', 'MScan_Neural_Data');

    % Alpha [8 - 12 Hz]
    [MScanData.Data.AlphaBand_Power, MScanData.Notes.MScan.alphaBandSamplingRate] = ...
        ProcessNeuro_2P(MScanData, 'Alpha', 'MScan_Neural_Data');

    % Theta [4 - 8 Hz]
    [MScanData.Data.ThetaBand_Power, MScanData.Notes.MScan.thetaBandSamplingRate] = ...
        ProcessNeuro_2P(MScanData, 'Theta', 'MScan_Neural_Data');

    % Delta [1 - 4 Hz]
    [MScanData.Data.DeltaBand_Power, MScanData.Notes.MScan.deltaBandSamplingRate] = ...
        ProcessNeuro_2P(MScanData, 'Delta', 'MScan_Neural_Data');
    %% Downsample and binarize the force sensor.
    % Trim any additional data points for resample
    mscanForce = MScanData.Data.MScan_Force_Sensor;
    
    % Filter then downsample the Force Sensor waveform to desired frequency
    forceSensorDownSampledSamplingRate = 30;   % Downsample to CBV Camera Fs
    forceSensorFilterThreshold = 20;
    forceSensorFilterOrder = 2;
    [z, p, k] = butter(forceSensorFilterOrder, forceSensorFilterThreshold / (MScanData.Notes.MScan_analogSamplingRate / 2), 'low');
    [sos, g] = zp2sos(z, p, k);
    filteredForceSensor_M = filtfilt(sos, g, mscanForce);
    
    MScanData.Data.dsForce_Sensor_M = resample(filteredForceSensor_M, forceSensorDownSampledSamplingRate, MScanData.Notes.MScan_analogSamplingRate);
    MScanData.Notes.downsampledForceSensorSamplingRate = forceSensorDownSampledSamplingRate;
    
    % Binarize the force sensor waveform
    threshfile = dir('*_Thresholds.mat');
    if ~isempty(threshfile)
        load(threshfile.name)
    end
    
    [ok] = CheckForThreshold(['binarizedForceSensor_' strDay], animalID);
    
    if ok == 0
        [forceSensorThreshold] = CreateForceSensorThreshold(MScanData.Data.dsForce_Sensor_M);
        Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
        save([animalID '_Thresholds.mat'], 'Thresholds');
    end
    
    MScanData.Data.binForce_Sensor_M = BinarizeForceSensor(MScanData.Data.dsForce_Sensor_M, Thresholds.(['binarizedForceSensor_' strDay]));

    % %% EMG.
    % if isfield(MScanData.Data, 'EMG')
    %     %% Downsample and binarize the EMG.
    %     % Trim any additional data points for resample
    %     trimmedEMG = MScanData.Data.EMG(1:min(expectedLength, length(MScanData.Data.EMG)));
    %
    %     % Filter then downsample the Force Sensor waveform to desired frequency
    %     emgDownSampledSamplingRate = MScanData.Notes.CBVCamSamplingRate;   % Downsample to CBV Camera Fs
    %     emgFilterThreshold = 20;
    %     emgFilterOrder = 2;
    %     [z, p, k] = butter(emgFilterOrder, emgFilterThreshold / (MScanData.Notes.analogSamplingRate / 2), 'low');
    %     [sos, g] = zp2sos(z, p, k);
    %     filteredEMG = filtfilt(sos, g, trimmedEMG);
    %
    %     MScanData.Data.Behavior.EMG = resample(filteredEMG, emgDownSampledSamplingRate, MScanData.Notes.analogSamplingRate);
    %     MScanData.Notes.downsampledEMGSamplingRate = emgDownSampledSamplingRate;
    % end
    
    save([animalID '_' date '_' imageID '_MScanData'], 'MScanData') 
end

%% LabVIEW data file analysis
for f = 1:size(labviewDataFiles, 1)
    %% Find offset between the two force sensor signals using the cross correlation
    disp(['Analyzing LabVIEW analog signals and whisker angle for file number ' num2str(f) ' of ' num2str(size(labviewDataFiles, 1)) '...']); disp(' '); 
    labviewDataFile = labviewDataFiles(f, :);
    load(labviewDataFile);
    [animalID, hem, fileDate, fileID] = GetFileInfo(labviewDataFile);
    strDay = ConvertDate(fileDate);

    %% Binarize the whisker angle and set the resting angle to zero degrees.
    % Trim any additional frames for resample
    whiskerAngle = LabVIEWData.Data.WhiskerAngle;
    
    % Create filter for whisking/movement
    whiskerDownsampledSamplingRate = 30;   % Downsample to CBV Camera Fs
    whiskerFilterThreshold = 20;
    whiskerFilterOrder = 2;
    [z, p, k] = butter(whiskerFilterOrder, whiskerFilterThreshold/(LabVIEWData.Notes.whiskerCamSamplingRate/2), 'low');
    [sos, g] = zp2sos(z, p, k);
    filteredWhiskers = filtfilt(sos, g, whiskerAngle - mean(whiskerAngle));
    resampledWhiskers = resample(filteredWhiskers, whiskerDownsampledSamplingRate, LabVIEWData.Notes.whiskerCamSamplingRate);
    
    % Binarize the whisker waveform (wwf)
    threshfile = dir('*_Thresholds.mat');
    if ~isempty(threshfile)
        load(threshfile.name)
    end
    
    [ok] = CheckForThreshold(['binarizedWhiskersLower_' strDay], animalID);
    
    if ok == 0
        [whiskersThresh1, whiskersThresh2] = CreateWhiskThreshold(resampledWhiskers, whiskerDownsampledSamplingRate);
        Thresholds.(['binarizedWhiskersLower_' strDay]) = whiskersThresh1;
        Thresholds.(['binarizedWhiskersUpper_' strDay]) = whiskersThresh2;
        save([animalID '_Thresholds.mat'], 'Thresholds');
    end
    
    load([animalID '_Thresholds.mat']);
    binarizedWhiskers = BinarizeWhiskers(resampledWhiskers, whiskerDownsampledSamplingRate, Thresholds.(['binarizedWhiskersLower_' strDay]), Thresholds.(['binarizedWhiskersUpper_' strDay]));
    [linkedBinarizedWhiskers] = LinkBinaryEvents(gt(binarizedWhiskers,0), [round(whiskerDownsampledSamplingRate/3), 0]);
    
    inds = linkedBinarizedWhiskers == 0;
    restAngle = mean(resampledWhiskers(inds));
    
    LabVIEWData.Data.dsWhisker_Angle = resampledWhiskers - restAngle;
    LabVIEWData.Data.binWhisker_Angle = binarizedWhiskers;
    LabVIEWData.Notes.downsampledWhiskerSamplingRate = whiskerDownsampledSamplingRate;
    
    %% Downsample and binarize the force sensor.
    % Trim any additional data points for resample
    labviewForce = LabVIEWData.Data.Force_Sensor;
    
    % Filter then downsample the Force Sensor waveform to desired frequency
    forceSensorDownSampledSamplingRate = 30;   % Downsample to CBV Camera Fs
    forceSensorFilterThreshold = 20;
    forceSensorFilterOrder = 2;
    [z, p, k] = butter(forceSensorFilterOrder, forceSensorFilterThreshold/(LabVIEWData.Notes.analogSamplingRate/2), 'low');
    [sos, g] = zp2sos(z, p, k);
    filteredForceSensor_L = filtfilt(sos, g, labviewForce);
    
    LabVIEWData.Data.dsForce_Sensor_L = resample(filteredForceSensor_L, forceSensorDownSampledSamplingRate, LabVIEWData.Notes.analogSamplingRate);
    LabVIEWData.Notes.downsampledForceSensorSamplingRate = forceSensorDownSampledSamplingRate;
    
    % Binarize the force sensor waveform
    [ok] = CheckForThreshold(['binarizedForceSensor_' strDay], animalID);
    
    if ok == 0
        [forceSensorThreshold] = CreateForceSensorThreshold(LabVIEWData.Data.dsForce_Sensor_L);
        Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
        save([animalID '_Thresholds.mat'], 'Thresholds');
    end
    
    LabVIEWData.Data.binForce_Sensor_L = BinarizeForceSensor(LabVIEWData.Data.dsForce_Sensor_L, Thresholds.(['binarizedForceSensor_' strDay]));
    
    save([animalID '_' hem '_' fileID '_LabVIEWData'], 'LabVIEWData')

end

end

