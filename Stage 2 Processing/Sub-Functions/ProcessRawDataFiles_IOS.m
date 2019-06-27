function ProcessRawdataFiles_IOS(rawdataFiles)
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
%
%   Inputs: List of LabVIEW and MScan data files.
%
%   Outputs: Saves updates to both files in the current directory.
%
%   Last Revised: March 21st, 2019
%________________________________________________________________________________________________________________________

%% Raw data file analysis
for a = 1:size(rawdataFiles,1)
    rawdataFile = rawdataFiles(a,:);
    disp(['Processing Rawdata File: ' rawdataFile]); disp(' ')
    [animal, hemisphere, fileDate, fileID] = GetFileInfo_IOS(fileName);
    strDay = ConvertDate(fileDate);
    procdataFile = ([animal '_' hemisphere '_' fileID '_Procdata.mat']);
    
    % Skip processing the file if it already exists
    if ~exist(procdataFile, 'file')
        load(rawdataFile);
        
        %% Transfer Rawdata notes to Procdata structure.
        Procdata.notes = Rawdata.notes;
        
        %% Save solenoid times (in seconds). Identify the solenoids by amplitude.
        Procdata.data.solenoids.LPadSol = find(diff(Rawdata.data.solenoids) == 1) / Rawdata.notes.analogSamplingRate;
        Procdata.data.solenoids.RPadSol = find(diff(Rawdata.data.solenoids) == 2) / Rawdata.notes.analogSamplingRate;
        Procdata.data.solenoids.AudSol = find(diff(Rawdata.data.solenoids) == 3) / Rawdata.notes.analogSamplingRate;
        
        %% CBV from ROIs.
        CBVfields = fieldnames(Rawdata.data.CBV);
        for b = 1:length(CBVfields)
            Procdata.data.CBV.([CBVfields{b}]) = Rawdata.data.CBV.(CBVfields{b})(1:end - 1);
        end
        
        %% Process neural data into its various forms.
        expectedLength = Procdata.notes.trialDuration_sec*Procdata.notes.analogSamplingRate;
        dsFs = 30;   % downsampled Fs

        neuralDataTypes = {'cortical_LH', 'cortical_RH', 'hippocampus'};
        for c = 1:length(neuralDataTypes)
            neuralDataType = neuralDataTypes{c,1};
        % MUA Band [300 - 3000]
        [muaPower] = ProcessNeuro_IOS(Procdata, expectedLength, 'MUA', 'rawNeuraldata');
        Procdata.data.(neuralDataType).muaPower = muaPower;
        
        % Gamma Band [40 - 100]
        [Procdata.data.gammaPower, ~] = ProcessNeuro_IOS(Procdata, expectedLength, 'Gam', 'rawNeuraldata');
        
        % Beta [13 - 30 Hz]
        [Procdata.data.betaPower, ~] = ProcessNeuro_IOS(Procdata, expectedLength, 'Beta', 'rawNeuraldata');
        
        % Alpha [8 - 12 Hz]
        [Procdata.data.alphaPower, ~] = ProcessNeuro_IOS(Procdata, expectedLength, 'Alpha', 'rawNeuraldata');
        
        % Theta [4 - 8 Hz]
        [Procdata.data.thetaPower, ~] = ProcessNeuro_IOS(Procdata, expectedLength, 'Theta', 'rawNeuraldata');
        
        % Delta [1 - 4 Hz]
        [Procdata.data.deltaPower, ~] = ProcessNeuro_IOS(Procdata, expectedLength, 'Delta', 'rawNeuraldata');
        
        %% Patch and binarize the whisker angle and set the resting angle to zero degrees.
        [patchedWhisk] = PatchWhiskerAngle_2P(LabVIEWdata.data.whiskerAngle, LabVIEWdata.notes.whiskerCamSamplingRate_Hz, LabVIEWdata.notes.trialDuration_Seconds, LabVIEWdata.notes.droppedWhiskerCamFrameIndex);
        
        % Create filter for whisking/movement
        dsFs = 30;
        filtThreshold = 20;
        filtOrder = 2;
        [z, p, k] = butter(filtOrder, filtThreshold/(LabVIEWdata.notes.whiskerCamSamplingRate_Hz/2), 'low');
        [sos, g] = zp2sos(z, p, k);
        filteredWhiskers = filtfilt(sos, g, patchedWhisk - mean(patchedWhisk));
        resampledWhisk = resample(filteredWhiskers, dsFs, LabVIEWdata.notes.whiskerCamSamplingRate_Hz);
        
        % Binarize the whisker waveform (wwf)
        threshfile = dir('*_Thresholds.mat');
        if ~isempty(threshfile)
            load(threshfile.name)
        end
        
        [ok] = CheckForThreshold_2P(['binarizedWhiskersLower_' strDay], animalID);
        
        if ok == 0
            [whiskersThresh1, whiskersThresh2] = CreateWhiskThreshold_2P(resampledWhisk, dsFs);
            Thresholds.(['binarizedWhiskersLower_' strDay]) = whiskersThresh1;
            Thresholds.(['binarizedWhiskersUpper_' strDay]) = whiskersThresh2;
            save([animalID '_Thresholds.mat'], 'Thresholds');
        end
        
        load([animalID '_Thresholds.mat']);
        binWhisk = BinarizeWhiskers_2P(resampledWhisk, dsFs, Thresholds.(['binarizedWhiskersLower_' strDay]), Thresholds.(['binarizedWhiskersUpper_' strDay]));
        [linkedBinarizedWhiskers] = LinkBinaryEvents_2P(gt(binWhisk,0), [round(dsFs/3), 0]);
        inds = linkedBinarizedWhiskers == 0;
        restAngle = mean(resampledWhisk(inds));
        
        LabVIEWdata.data.dsWhiskerAngle = resampledWhisk - restAngle;
        LabVIEWdata.data.binWhiskerAngle = binWhisk;
        
        
        %% Downsample and binarize the force sensor.
        trimmedForceM = Procdata.data.forceSensor(1:min(expectedLength, length(Procdata.data.forceSensor)));
        
        % Filter then downsample the Force Sensor waveform to desired frequency
        filtThreshold = 20;
        filtOrder = 2;
        [z, p, k] = butter(filtOrder, filtThreshold/(Procdata.notes.analogSamplingRate/2), 'low');
        [sos, g] = zp2sos(z, p, k);
        filtForceSensorM = filtfilt(sos, g, trimmedForceM);
        Procdata.data.dsForceSensorM = resample(filtForceSensorM, dsFs, Procdata.notes.analogSamplingRate);
        
        % Binarize the force sensor waveform
        threshfile = dir('*_Thresholds.mat');
        if ~isempty(threshfile)
            load(threshfile.name)
        end
        
        [ok] = CheckForThreshold_2P(['binarizedForceSensor_' strDay], animalID);
        
        if ok == 0
            [forceSensorThreshold] = CreateForceSensorThreshold_2P(Procdata.data.dsForceSensorM);
            Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
            save([animalID '_Thresholds.mat'], 'Thresholds');
        end
        
        Procdata.data.binForceSensorM = BinarizeForceSensor_2P(Procdata.data.dsForceSensorM, Thresholds.(['binarizedForceSensor_' strDay]));
        
        %% EMG
        fpass = [30 300];
        trimmedEMG = (1:min(expectedLength, length(Procdata.data.EMG)));
        [z1, p1, k1] = butter(4, fpass/(Procdata.notes.analogSamplingRate/2));
        [sos1, g1] = zp2sos(z1, p1, k1);
        filtEMG = filtfilt(sos1, g1, trimmedEMG - mean(trimmedEMG));
        [z2, p2, k2] = butter(4, 10/(Procdata.notes.analogSamplingRate/2), 'low');
        [sos2, g2] = zp2sos(z2, p2, k2);
        smoothEMGPower = filtfilt(sos2, g2, filtEMG.^2);
        Procdata.data.dsEMG = max(resample(smoothEMGPower, dsFs, Procdata.notes.analogSamplingRate), 0);
        
        %% Save the processed data
        save([animalID '_' hemisphere '_' fileID '_Procdata'], 'Procdata')
    else
        disp([rawdataFile ' has already been processed. Continuing...']); disp(' ');
    end
end

end