function ProcessRawDataFiles_IOS(rawDataFiles)
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
%   Last Revised: June 27th, 2019
%________________________________________________________________________________________________________________________

%% Raw data file analysis
for a = 1:size(rawDataFiles,1)
    rawDataFile = rawDataFiles(a,:);
    disp(['Analyzing RawData file ' num2str(a) ' of ' num2str(size(rawDataFiles, 1)) '...']); disp(' ')
    [animalID, fileDate, fileID] = GetFileInfo_IOS(rawDataFile);
    strDay = ConvertDate_IOS(fileDate);
    procDataFile = ([animalID '_' fileID '_ProcData.mat']);
    disp(['Generating ' procDataFile '...']); disp(' ')
    load(rawDataFile);
    
    %% Transfer RawData notes to ProcData structure.
    ProcData.notes = RawData.notes;
    
    %% Expected durations
    analogExpectedLength = ProcData.notes.trialDuration_sec*ProcData.notes.analogSamplingRate;
    
    %% Save solenoid times (in seconds). Identify the solenoids by amplitude.
    ProcData.data.solenoids.LPadSol = find(diff(RawData.data.solenoids) == 1)/RawData.notes.analogSamplingRate;
    ProcData.data.solenoids.RPadSol = find(diff(RawData.data.solenoids) == 2)/RawData.notes.analogSamplingRate;
    ProcData.data.solenoids.AudSol = find(diff(RawData.data.solenoids) == 3)/RawData.notes.analogSamplingRate;
    
    %% CBV from ROIs.
    CBVfields = fieldnames(RawData.data.CBV);
    for b = 1:length(CBVfields)
        ProcData.data.CBV.(CBVfields{b}(1:end-6)) = RawData.data.CBV.(CBVfields{b})(1:end - 1);
    end
    CheckForNaNs_IOS(ProcData);
    
    %% Process neural data into its various forms.
    ProcData.notes.dsFs = 30;   % downsampled Fs
    
    neuralDataTypes = {'cortical_LH', 'cortical_RH', 'hippocampus'};
    for c = 1:length(neuralDataTypes)
        neuralDataType = neuralDataTypes{1,c};
        % MUA Band [300 - 3000]
        [muaPower] = ProcessNeuro_IOS(RawData, analogExpectedLength, 'MUA', neuralDataType, ProcData.notes.dsFs);
        ProcData.data.(neuralDataType).muaPower = muaPower;
        
        % Gamma Band [40 - 100]
        [gammaBandPower] = ProcessNeuro_IOS(RawData, analogExpectedLength, 'Gamma', neuralDataType, ProcData.notes.dsFs);
        ProcData.data.(neuralDataType).gammaBandPower = gammaBandPower;
        
        % Beta [13 - 30 Hz]
        [betaBandPower] = ProcessNeuro_IOS(RawData, analogExpectedLength, 'Beta', neuralDataType, ProcData.notes.dsFs);
        ProcData.data.(neuralDataType).betaBandPower = betaBandPower;
        
        % Alpha [8 - 12 Hz]
        [alphaBandPower] = ProcessNeuro_IOS(RawData, analogExpectedLength, 'Alpha', neuralDataType, ProcData.notes.dsFs);
        ProcData.data.(neuralDataType).alphaBandPower = alphaBandPower;
        
        % Theta [4 - 8 Hz]
        [thetaBandPower] = ProcessNeuro_IOS(RawData, analogExpectedLength, 'Theta', neuralDataType, ProcData.notes.dsFs);
        ProcData.data.(neuralDataType).thetaBandPower = thetaBandPower;
        
        % Delta [1 - 4 Hz]
        [deltaBandPower] = ProcessNeuro_IOS(RawData, analogExpectedLength, 'Delta', neuralDataType, ProcData.notes.dsFs);
        ProcData.data.(neuralDataType).deltaBandPower = deltaBandPower;
    end
    
    %% Patch and binarize the whisker angle and set the resting angle to zero degrees.
    [patchedWhisk, droppedFrames] = PatchWhiskerAngle_IOS(RawData.data.whiskerAngle, RawData.notes.whiskCamSamplingRate, RawData.notes.trialDuration_sec, RawData.notes.droppedWhiskCamFrameIndex);
    RawData.data.patchedWhiskerAngle = patchedWhisk;
    if droppedFrames >= 200
        disp(['WARNING - ' num2str(droppedFrames) ' dropped whisker camera frames from file ID ' rawDataFile '.']); disp(' ')
    end
    
    % Create filter for whisking/movement
    filtThreshold = 20;
    filtOrder = 2;
    [z, p, k] = butter(filtOrder, filtThreshold/(RawData.notes.whiskCamSamplingRate/2), 'low');
    [sos, g] = zp2sos(z, p, k);
    filteredWhiskers = filtfilt(sos, g, patchedWhisk - mean(patchedWhisk));
    resampledWhisk = resample(filteredWhiskers, ProcData.notes.dsFs, RawData.notes.whiskCamSamplingRate);
    
    % Binarize the whisker waveform (wwf)
    threshfile = dir('*_Thresholds.mat');
    if ~isempty(threshfile)
        load(threshfile.name)
    end
    
    [ok] = CheckForThreshold_IOS(['binarizedWhiskersLower_' strDay], animalID);
    
    if ok == 0
        [whiskersThresh1, whiskersThresh2] = CreateWhiskThreshold_IOS(resampledWhisk, ProcData.notes.dsFs);
        Thresholds.(['binarizedWhiskersLower_' strDay]) = whiskersThresh1;
        Thresholds.(['binarizedWhiskersUpper_' strDay]) = whiskersThresh2;
        save([animalID '_Thresholds.mat'], 'Thresholds');
    end
    
    load([animalID '_Thresholds.mat']);
    binWhisk = BinarizeWhiskers_IOS(resampledWhisk, ProcData.notes.dsFs, Thresholds.(['binarizedWhiskersLower_' strDay]), Thresholds.(['binarizedWhiskersUpper_' strDay]));
    [linkedBinarizedWhiskers] = LinkBinaryEvents_IOS(gt(binWhisk,0), [round(ProcData.notes.dsFs/3), 0]);
    inds = linkedBinarizedWhiskers == 0;
    restAngle = mean(resampledWhisk(inds));
    
    ProcData.data.whiskerAngle = resampledWhisk - restAngle;
    ProcData.data.binWhiskerAngle = binWhisk;
    
    
    %% Downsample and binarize the force sensor.
    trimmedForce = RawData.data.forceSensor(1:min(analogExpectedLength, length(RawData.data.forceSensor)));
    
    % Filter then downsample the Force Sensor waveform to desired frequency
    filtThreshold = 20;
    filtOrder = 2;
    [z, p, k] = butter(filtOrder, filtThreshold/(ProcData.notes.analogSamplingRate/2), 'low');
    [sos, g] = zp2sos(z, p, k);
    filtForceSensor = filtfilt(sos, g, trimmedForce);
    ProcData.data.forceSensor = resample(filtForceSensor, ProcData.notes.dsFs, ProcData.notes.analogSamplingRate);
    
    % Binarize the force sensor waveform
    threshfile = dir('*_Thresholds.mat');
    if ~isempty(threshfile)
        load(threshfile.name)
    end
    
    [ok] = CheckForThreshold_IOS(['binarizedForceSensor_' strDay], animalID);
    
    if ok == 0
        [forceSensorThreshold] = CreateForceSensorThreshold_IOS(ProcData.data.forceSensor);
        Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
        save([animalID '_Thresholds.mat'], 'Thresholds');
    end
    
    ProcData.data.binForceSensor = BinarizeForceSensor_IOS(ProcData.data.forceSensor, Thresholds.(['binarizedForceSensor_' strDay]));
    
    %% EMG
    fpass = [300 3000];
    trimmedEMG = RawData.data.EMG(1:min(analogExpectedLength, length(RawData.data.EMG)));
    [z1, p1, k1] = butter(3, fpass/(ProcData.notes.analogSamplingRate/2));
    [sos1, g1] = zp2sos(z1, p1, k1);
    filtEMG = filtfilt(sos1, g1, trimmedEMG - mean(trimmedEMG));
    kernelWidth = 0.5;
    smoothingKernel = gausswin(kernelWidth*ProcData.notes.analogSamplingRate)/sum(gausswin(kernelWidth*ProcData.notes.analogSamplingRate));
    EMGPwr = log10(conv(filtEMG.^2, smoothingKernel, 'same'));
    resampEMG = resample(EMGPwr, ProcData.notes.dsFs, ProcData.notes.analogSamplingRate);
    ProcData.data.EMG.emg = resampEMG;
    
    % nancheck
    if sum(isnan(ProcData.data.EMG.emg)) > 0
        keyboard
    end
    
    if sum(isinf(ProcData.data.EMG.emg)) > 0
        keyboard
    end
    
    %% Laser Doppler
    if isfield(RawData.data, 'backScatter') == true
        trimmedBackScatter = RawData.data.backScatter(1:min(analogExpectedLength, length(RawData.data.backScatter)));
        trimmedFlow = RawData.data.flow(1:min(analogExpectedLength, length(RawData.data.flow)));
        ProcData.data.backScatter = max(resample(trimmedBackScatter, ProcData.notes.dsFs, ProcData.notes.analogSamplingRate), 0);
        ProcData.data.flow = max(resample(trimmedFlow, ProcData.notes.dsFs, ProcData.notes.analogSamplingRate), 0);
    end
    
    %% Save the processed data
    save(rawDataFile, 'RawData')
    save(procDataFile, 'ProcData')
end

end