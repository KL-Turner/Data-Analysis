function ProcessCombDataFile(MergedDataFiles)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner 
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
%
%   Last Revised: February 20th, 2019    
%________________________________________________________________________________________________________________________

for f = 1:size(MergedDataFiles, 1)
    %% Find offset between the two force sensor signals using the cross correlation
    disp(['Analyzing neural bands, analog signals, and whisker angle ' num2str(f) ' of ' num2str(size(MergedDataFiles, 1)) '...']); disp(' '); 
    MergedDataFile = MergedDataFiles(f, :);
    load(MergedDataFile);
    [animalID, fileDate, fileID, vesselID] = GetFileInfo2(MergedDataFile);
    strDay = ConvertDate(fileDate);

    %% Process neural data into its various forms.
    % MUA Band [300 - 3000]
    [MergedData.Data.MUA_Power, MergedData.Notes.MScan.multiUnitSamplingRate] = ...
        ProcessNeuro2(MergedData, 'MUApower', 'Neural_Data');

    % Gamma Band [40 - 100]
    [MergedData.Data.GammaBand_Power, MergedData.Notes.MScan.gammaBandSamplingRate] = ...
        ProcessNeuro2(MergedData, 'Gam', 'Neural_Data');

    % Beta [13 - 30 Hz]
    [MergedData.Data.BetaBand_Power, MergedData.Notes.MScan.betaBandSamplingRate] = ...
        ProcessNeuro2(MergedData, 'Beta', 'Neural_Data');

    % Alpha [8 - 12 Hz]
    [MergedData.Data.AlphaBand_Power, MergedData.Notes.MScan.alphaBandSamplingRate] = ...
        ProcessNeuro2(MergedData, 'Alpha', 'Neural_Data');

    % Theta [4 - 8 Hz]
    [MergedData.Data.ThetaBand_Power, MergedData.Notes.MScan.thetaBandSamplingRate] = ...
        ProcessNeuro2(MergedData, 'Theta', 'Neural_Data');

    % Delta [1 - 4 Hz]
    [MergedData.Data.DeltaBand_Power, MergedData.Notes.MScan.deltaBandSamplingRate] = ...
        ProcessNeuro2(MergedData, 'Delta', 'Neural_Data');

    %% Binarize the whisker angle and set the resting angle to zero degrees.
    % Trim any additional frames for resample
    whiskerAngle = MergedData.Data.Whisker_Angle;
    
    % Create filter for whisking/movement
    whiskerDownsampledSamplingRate = 30;   % Downsample to CBV Camera Fs
    whiskerFilterThreshold = 20;
    whiskerFilterOrder = 2;
    [z, p, k] = butter(whiskerFilterOrder, whiskerFilterThreshold / (MergedData.Notes.LabVIEW.whiskerCamSamplingRate / 2), 'low');
    [sos, g] = zp2sos(z, p, k);
    filteredWhiskers = filtfilt(sos, g, whiskerAngle - mean(whiskerAngle));
    resampledWhiskers = resample(filteredWhiskers, whiskerDownsampledSamplingRate, MergedData.Notes.LabVIEW.whiskerCamSamplingRate);
    
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
    
    MergedData.Data.dsWhisker_Angle = resampledWhiskers - restAngle;
    MergedData.Notes.LabVIEW.downsampledWhiskerSamplingRate = whiskerDownsampledSamplingRate;
    MergedData.Data.binWhisker_Angle = binarizedWhiskers;
    
    %% Downsample and binarize the force sensor.
    % Trim any additional data points for resample
    mscanForce = MergedData.Data.Force_Sensor_M;
    labviewForce = MergedData.Data.Force_Sensor_L;
    
    % Filter then downsample the Force Sensor waveform to desired frequency
    forceSensorDownSampledSamplingRate = 30;   % Downsample to CBV Camera Fs
    forceSensorFilterThreshold = 20;
    forceSensorFilterOrder = 2;
    [z, p, k] = butter(forceSensorFilterOrder, forceSensorFilterThreshold / (MergedData.Notes.LabVIEW.analogSamplingRate / 2), 'low');
    [sos, g] = zp2sos(z, p, k);
    filteredForceSensor_M = filtfilt(sos, g, mscanForce);
    filteredForceSensor_L = filtfilt(sos, g, labviewForce);
    
    MergedData.Data.dsForce_Sensor_M = resample(filteredForceSensor_M, forceSensorDownSampledSamplingRate, MergedData.Notes.LabVIEW.analogSamplingRate);
    MergedData.Data.dsForce_Sensor_L = resample(filteredForceSensor_L, forceSensorDownSampledSamplingRate, MergedData.Notes.LabVIEW.analogSamplingRate);
    MergedData.Notes.LabVIEW.downsampledForceSensorSamplingRate = forceSensorDownSampledSamplingRate;
    
    % Binarize the force sensor waveform
    [ok] = CheckForThreshold(['binarizedForceSensor_' strDay], animalID);
    
    if ok == 0
        [forceSensorThreshold] = CreateForceSensorThreshold(MergedData.Data.dsForce_Sensor_M);
        Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
        save([animalID '_Thresholds.mat'], 'Thresholds');
    end
    
    MergedData.Data.binForce_Sensor_M = BinarizeForceSensor(MergedData.Data.dsForce_Sensor_M, Thresholds.(['binarizedForceSensor_' strDay]));
    MergedData.Data.binForce_Sensor_L = BinarizeForceSensor(MergedData.Data.dsForce_Sensor_L, Thresholds.(['binarizedForceSensor_' strDay]));

    % %% EMG.
    % if isfield(MergedData.Data, 'EMG')
    %     %% Downsample and binarize the EMG.
    %     % Trim any additional data points for resample
    %     trimmedEMG = MergedData.Data.EMG(1:min(expectedLength, length(MergedData.Data.EMG)));
    %
    %     % Filter then downsample the Force Sensor waveform to desired frequency
    %     emgDownSampledSamplingRate = MergedData.Notes.CBVCamSamplingRate;   % Downsample to CBV Camera Fs
    %     emgFilterThreshold = 20;
    %     emgFilterOrder = 2;
    %     [z, p, k] = butter(emgFilterOrder, emgFilterThreshold / (MergedData.Notes.analogSamplingRate / 2), 'low');
    %     [sos, g] = zp2sos(z, p, k);
    %     filteredEMG = filtfilt(sos, g, trimmedEMG);
    %
    %     MergedData.Data.Behavior.EMG = resample(filteredEMG, emgDownSampledSamplingRate, MergedData.Notes.analogSamplingRate);
    %     MergedData.Notes.downsampledEMGSamplingRate = emgDownSampledSamplingRate;
    % end
    
    save([animalID '_' fileID '_' vesselID '_MergedData'], 'MergedData')

end

end

