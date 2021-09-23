function [] = ProcessRawDataFiles_FP(rawDataFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the force sensor and neural bands. Create a threshold for binarized movement/whisking if
%            one does not already exist.
%________________________________________________________________________________________________________________________

% Raw data file analysis
for a = 1:size(rawDataFiles,1)
    rawDataFile = rawDataFiles(a,:);
    disp(['Analyzing RawData file ' num2str(a) ' of ' num2str(size(rawDataFiles, 1)) '...']); disp(' ')
    [animalID, fileDate, fileID] = GetFileInfo_FP(rawDataFile);
    strDay = ConvertDate_FP(fileDate);
    procDataFile = ([animalID '_' fileID '_ProcData.mat']);
    disp(['Generating ' procDataFile '...']); disp(' ')
    load(rawDataFile);
    % Transfer RawData notes to ProcData structure.
    ProcData.notes = RawData.notes;
    ProcData.notes.dsFs = 30; % downsampled Fs
    analogExpectedLength = ProcData.notes.trialDuration_sec_trim*ProcData.notes.analogSamplingRate;
    dsExpectedLength = ProcData.notes.trialDuration_sec_trim*ProcData.notes.dsFs;
    %% identify the solenoids by amplitude.
    ProcData.data.stimulations.LPadSol = find(diff(RawData.data.stimulations) == 1)/RawData.notes.analogSamplingRate;
    ProcData.data.stimulations.RPadSol = find(diff(RawData.data.stimulations) == 2)/RawData.notes.analogSamplingRate;
    ProcData.data.stimulations.AudSol = find(diff(RawData.data.stimulations) == 3)/RawData.notes.analogSamplingRate;
    %% resample doric data to 30 Hz
    ProcData.data.HbT.LH = (resample(RawData.data.LH_560.zScored_trim,ProcData.notes.dsFs,ProcData.notes.doric.samplingRate));
    ProcData.data.HbT.LH = ProcData.data.HbT.LH(1:dsExpectedLength);
    ProcData.data.HbT.RH = (resample(RawData.data.RH_560.zScored_trim,ProcData.notes.dsFs,ProcData.notes.doric.samplingRate));
    ProcData.data.HbT.RH = ProcData.data.HbT.RH(1:dsExpectedLength);
    ProcData.data.GCaMP7s.LH = (resample(RawData.data.LH_465.zScored_trim,ProcData.notes.dsFs,ProcData.notes.doric.samplingRate));
    ProcData.data.GCaMP7s.LH = ProcData.data.GCaMP7s.LH(1:dsExpectedLength);
    ProcData.data.GCaMP7s.RH = (resample(RawData.data.RH_465.zScored_trim,ProcData.notes.dsFs,ProcData.notes.doric.samplingRate));
    ProcData.data.GCaMP7s.RH = ProcData.data.GCaMP7s.RH(1:dsExpectedLength);
    %% process neural data into its various forms.
    neuralDataTypes = {'cortical_LH','cortical_RH','hippocampus'};
    for c = 1:length(neuralDataTypes)
        neuralDataType = neuralDataTypes{1,c};
        % MUA Band [300 - 3000]
        [muaPower,~] = ProcessNeuro_FP(RawData,analogExpectedLength,'MUA',neuralDataType);
        ProcData.data.(neuralDataType).muaPower = muaPower;
        % Gamma Band [40 - 100]
        [gammaBandPower,~] = ProcessNeuro_FP(RawData,analogExpectedLength,'Gam',neuralDataType);
        ProcData.data.(neuralDataType).gammaBandPower = gammaBandPower;
        % Beta [13 - 30 Hz]
        [betaBandPower,~] = ProcessNeuro_FP(RawData,analogExpectedLength,'Beta',neuralDataType);
        ProcData.data.(neuralDataType).betaBandPower = betaBandPower;
        % Alpha [8 - 12 Hz]
        [alphaBandPower,~] = ProcessNeuro_FP(RawData,analogExpectedLength,'Alpha',neuralDataType);
        ProcData.data.(neuralDataType).alphaBandPower = alphaBandPower;
        % Theta [4 - 8 Hz]
        [thetaBandPower,~] = ProcessNeuro_FP(RawData,analogExpectedLength,'Theta',neuralDataType);
        ProcData.data.(neuralDataType).thetaBandPower = thetaBandPower;
        % Delta [1 - 4 Hz]
        [deltaBandPower,~] = ProcessNeuro_FP(RawData,analogExpectedLength,'Delta',neuralDataType);
        ProcData.data.(neuralDataType).deltaBandPower = deltaBandPower;
    end
    %% Patch and binarize the whisker angle and set the resting angle to zero degrees.
    filtThreshold = 20;
    filtOrder = 2;
    [z,p,k] = butter(filtOrder,filtThreshold/(RawData.notes.whiskCamSamplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    filteredWhiskers = filtfilt(sos,g,RawData.data.whiskerAngle - mean(RawData.data.whiskerAngle));
    resampledWhisk = resample(filteredWhiskers,ProcData.notes.dsFs,RawData.notes.whiskCamSamplingRate);
    % Binarize the whisker waveform (wwf)
    threshfile = dir('*_Thresholds.mat');
    if ~isempty(threshfile)
        load(threshfile.name)
    end
    [ok] = CheckForThreshold_FP(['binarizedWhiskersLower_' strDay],animalID);
    if ok == 0
        [whiskersThresh1,whiskersThresh2] = CreateWhiskThreshold_FP(resampledWhisk,ProcData.notes.dsFs);
        Thresholds.(['binarizedWhiskersLower_' strDay]) = whiskersThresh1;
        Thresholds.(['binarizedWhiskersUpper_' strDay]) = whiskersThresh2;
        save([animalID '_Thresholds.mat'], 'Thresholds');
    end
    load([animalID '_Thresholds.mat']);
    binWhisk = BinarizeWhiskers_FP(resampledWhisk,ProcData.notes.dsFs,Thresholds.(['binarizedWhiskersLower_' strDay]),Thresholds.(['binarizedWhiskersUpper_' strDay]));
    [linkedBinarizedWhiskers] = LinkBinaryEvents_FP(gt(binWhisk,0),[round(ProcData.notes.dsFs/3),0]);
    inds = linkedBinarizedWhiskers == 0;
    restAngle = mean(resampledWhisk(inds));
    ProcData.data.whiskerAngle = resampledWhisk - restAngle;
    ProcData.data.binWhiskerAngle = binWhisk;
    %% Downsample and binarize the force sensor.
    trimmedForce = RawData.data.forceSensor(1:min(analogExpectedLength,length(RawData.data.forceSensor)));
    % Filter then downsample the Force Sensor waveform to desired frequency
    filtThreshold = 20;
    filtOrder = 2;
    [z,p,k] = butter(filtOrder,filtThreshold/(ProcData.notes.analogSamplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    filtForceSensor = filtfilt(sos,g,trimmedForce);
    ProcData.data.forceSensor = resample(filtForceSensor,ProcData.notes.dsFs,ProcData.notes.analogSamplingRate);
    % Binarize the force sensor waveform
    threshfile = dir('*_Thresholds.mat');
    if ~isempty(threshfile)
        load(threshfile.name)
    end
    [ok] = CheckForThreshold_FP(['binarizedForceSensor_' strDay],animalID);
    if ok == 0
        [forceSensorThreshold] = CreateForceSensorThreshold_FP(ProcData.data.forceSensor);
        Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
        save([animalID '_Thresholds.mat'],'Thresholds');
    end
    ProcData.data.binForceSensor = BinarizeForceSensor_FP(ProcData.data.forceSensor,Thresholds.(['binarizedForceSensor_' strDay]));
    %% EMG
    fpass = [300,3000];
    trimmedEMG = RawData.data.EMG(1:min(analogExpectedLength,length(RawData.data.EMG)));
    [z,p,k] = butter(3,fpass/(ProcData.notes.analogSamplingRate/2));
    [sos,g] = zp2sos(z,p,k);
    filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
    kernelWidth = 0.5;
    smoothingKernel = gausswin(kernelWidth*ProcData.notes.analogSamplingRate)/sum(gausswin(kernelWidth*ProcData.notes.analogSamplingRate));
    EMGPwr = log10(conv(filtEMG.^2,smoothingKernel,'same'));
    resampEMG = resample(EMGPwr,ProcData.notes.dsFs,ProcData.notes.analogSamplingRate);
    ProcData.data.EMG.emg = resampEMG;
    %% save the processed data
    save(procDataFile, 'ProcData')
end

end
