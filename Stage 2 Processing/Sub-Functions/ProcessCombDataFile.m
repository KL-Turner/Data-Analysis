function ProcessCombDataFile(combDataFiles)
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

for f = 1:size(combDataFiles, 1)
    %% Find offset between the two force sensor signals using the cross correlation
    disp(['Analyzing neural bands, analog signals, and whisker angle ' num2str(f) ' of ' num2str(size(combDataFiles, 1)) '...']); disp(' '); 
    combDataFile = combDataFiles(f, :);
    load(combDataFile);
    [animalID, fileDate, fileID, imageID] = GetFileInfo2(combDataFile);
    strDay = ConvertDate(fileDate);

    %% Process neural data into its various forms.
    % MUA Band [300 - 3000]
    [CombData.Data.MUA_Power, CombData.Notes.MScan.multiUnitSamplingRate] = ...
        ProcessNeuro2(CombData, 'MUApower', 'Neural_Data');

    % Gamma Band [40 - 100]
    [CombData.Data.GammaBand_Power, CombData.Notes.MScan.gammaBandSamplingRate] = ...
        ProcessNeuro2(CombData, 'Gam', 'Neural_Data');

    % Beta [13 - 30 Hz]
    [CombData.Data.BetaBand_Power, CombData.Notes.MScan.betaBandSamplingRate] = ...
        ProcessNeuro2(CombData, 'Beta', 'Neural_Data');

    % Alpha [8 - 12 Hz]
    [CombData.Data.AlphaBand_Power, CombData.Notes.MScan.alphaBandSamplingRate] = ...
        ProcessNeuro2(CombData, 'Alpha', 'Neural_Data');

    % Theta [4 - 8 Hz]
    [CombData.Data.ThetaBand_Power, CombData.Notes.MScan.thetaBandSamplingRate] = ...
        ProcessNeuro2(CombData, 'Theta', 'Neural_Data');

    % Delta [1 - 4 Hz]
    [CombData.Data.DeltaBand_Power, CombData.Notes.MScan.deltaBandSamplingRate] = ...
        ProcessNeuro2(CombData, 'Delta', 'Neural_Data');

    %% Binarize the whisker angle and set the resting angle to zero degrees.
    % Trim any additional frames for resample
    whiskerAngle = CombData.Data.Whisker_Angle;
    
    % Create filter for whisking/movement
    whiskerDownsampledSamplingRate = 30;   % Downsample to CBV Camera Fs
    whiskerFilterThreshold = 20;
    whiskerFilterOrder = 2;
    [z, p, k] = butter(whiskerFilterOrder, whiskerFilterThreshold / (CombData.Notes.LabVIEW.whiskerCamSamplingRate / 2), 'low');
    [sos, g] = zp2sos(z, p, k);
    filteredWhiskers = filtfilt(sos, g, whiskerAngle - mean(whiskerAngle));
    resampledWhiskers = resample(filteredWhiskers, whiskerDownsampledSamplingRate, CombData.Notes.LabVIEW.whiskerCamSamplingRate);
    
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
    
    CombData.Data.dsWhisker_Angle = resampledWhiskers - restAngle;
    CombData.Notes.LabVIEW.downsampledWhiskerSamplingRate = whiskerDownsampledSamplingRate;
    CombData.Data.binWhisker_Angle = binarizedWhiskers;
    
    %% Downsample and binarize the force sensor.
    % Trim any additional data points for resample
    mscanForce = CombData.Data.Force_Sensor_M;
    labviewForce = CombData.Data.Force_Sensor_L;
    
    % Filter then downsample the Force Sensor waveform to desired frequency
    forceSensorDownSampledSamplingRate = 30;   % Downsample to CBV Camera Fs
    forceSensorFilterThreshold = 20;
    forceSensorFilterOrder = 2;
    [z, p, k] = butter(forceSensorFilterOrder, forceSensorFilterThreshold / (CombData.Notes.LabVIEW.analogSamplingRate / 2), 'low');
    [sos, g] = zp2sos(z, p, k);
    filteredForceSensor_M = filtfilt(sos, g, mscanForce);
    filteredForceSensor_L = filtfilt(sos, g, labviewForce);
    
    CombData.Data.dsForce_Sensor_M = resample(filteredForceSensor_M, forceSensorDownSampledSamplingRate, CombData.Notes.LabVIEW.analogSamplingRate);
    CombData.Data.dsForce_Sensor_L = resample(filteredForceSensor_L, forceSensorDownSampledSamplingRate, CombData.Notes.LabVIEW.analogSamplingRate);
    CombData.Notes.LabVIEW.downsampledForceSensorSamplingRate = forceSensorDownSampledSamplingRate;
    
    % Binarize the force sensor waveform
    [ok] = CheckForThreshold(['binarizedForceSensor_' strDay], animalID);
    
    if ok == 0
        [forceSensorThreshold] = CreateForceSensorThreshold(CombData.Data.dsForce_Sensor_M);
        Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
        save([animalID '_Thresholds.mat'], 'Thresholds');
    end
    
    CombData.Data.binForce_Sensor_M = BinarizeForceSensor(CombData.Data.dsForce_Sensor_M, Thresholds.(['binarizedForceSensor_' strDay]));
    CombData.Data.binForce_Sensor_L = BinarizeForceSensor(CombData.Data.dsForce_Sensor_L, Thresholds.(['binarizedForceSensor_' strDay]));

    % %% EMG.
    % if isfield(CombData.Data, 'EMG')
    %     %% Downsample and binarize the EMG.
    %     % Trim any additional data points for resample
    %     trimmedEMG = CombData.Data.EMG(1:min(expectedLength, length(CombData.Data.EMG)));
    %
    %     % Filter then downsample the Force Sensor waveform to desired frequency
    %     emgDownSampledSamplingRate = CombData.Notes.CBVCamSamplingRate;   % Downsample to CBV Camera Fs
    %     emgFilterThreshold = 20;
    %     emgFilterOrder = 2;
    %     [z, p, k] = butter(emgFilterOrder, emgFilterThreshold / (CombData.Notes.analogSamplingRate / 2), 'low');
    %     [sos, g] = zp2sos(z, p, k);
    %     filteredEMG = filtfilt(sos, g, trimmedEMG);
    %
    %     CombData.Data.Behavior.EMG = resample(filteredEMG, emgDownSampledSamplingRate, CombData.Notes.analogSamplingRate);
    %     CombData.Notes.downsampledEMGSamplingRate = emgDownSampledSamplingRate;
    % end
    
    save([animalID '_' fileID '_' imageID '_CombData'], 'CombData')

end

end

