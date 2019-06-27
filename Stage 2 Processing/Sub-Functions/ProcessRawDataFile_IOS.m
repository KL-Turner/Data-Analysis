function ProcessRawDataFile_IOS(fileName)
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
for a = 1:size(rawDataFiles,1)
    rawDataFile = rawDataFiles(a,:);
    load(rawDataFile);
    
    disp(['Processing RawData File: ' fileName]); disp(' ')
load(fileName);
underscoreIndexes = strfind(fileName,'_');
[animal, hemisphere, fileDate, fileID] = GetFileInfo_IOS(fileName);
strDay = ConvertDate(fileDate);

%% Load in ProcData notes from RawData structure:
LabVIEWProcData.notes = RawData.notes;

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
    % Skip the file if it has already been processed
        disp(['Analyzing MScan neural bands and analog signals for file number ' num2str(a) ' of ' num2str(size(rawDataFiles, 1)) '...']); disp(' ');
        animalID = ProcData.notes.animalID;
        imageID = ProcData.notes.imageID;
        vesselID = ProcData.notes.vesselID;
        date = ProcData.notes.date;
        strDay = ConvertDate_2P(date);
        
        expectedLength = (ProcData.notes.numberOfFrames/ProcData.notes.frameRate)*ProcData.notes.analogSamplingRate;
        %% Process neural data into its various forms.
        % MUA Band [300 - 3000]
        [ProcData.data.muaPower, ProcData.notes.downSampledFs] = ProcessNeuro_2P(ProcData, expectedLength, 'MUA', 'rawNeuralData');
        downSampledFs = ProcData.notes.downSampledFs;

        % Gamma Band [40 - 100]
        [ProcData.data.gammaPower, ~] = ProcessNeuro_2P(ProcData, expectedLength, 'Gam', 'rawNeuralData');
        
        % Beta [13 - 30 Hz]
        [ProcData.data.betaPower, ~] = ProcessNeuro_2P(ProcData, expectedLength, 'Beta', 'rawNeuralData');
        
        % Alpha [8 - 12 Hz]
        [ProcData.data.alphaPower, ~] = ProcessNeuro_2P(ProcData, expectedLength, 'Alpha', 'rawNeuralData');
        
        % Theta [4 - 8 Hz]
        [ProcData.data.thetaPower, ~] = ProcessNeuro_2P(ProcData, expectedLength, 'Theta', 'rawNeuralData');
        
        % Delta [1 - 4 Hz]
        [ProcData.data.deltaPower, ~] = ProcessNeuro_2P(ProcData, expectedLength, 'Delta', 'rawNeuralData');
        
        %% Patch and binarize the whisker angle and set the resting angle to zero degrees.
        [patchedWhisk] = PatchWhiskerAngle_2P(LabVIEWData.data.whiskerAngle, LabVIEWData.notes.whiskerCamSamplingRate_Hz, LabVIEWData.notes.trialDuration_Seconds, LabVIEWData.notes.droppedWhiskerCamFrameIndex);
        
        % Create filter for whisking/movement
        downSampledFs = 30;
        filtThreshold = 20;
        filtOrder = 2;
        [z, p, k] = butter(filtOrder, filtThreshold/(LabVIEWData.notes.whiskerCamSamplingRate_Hz/2), 'low');
        [sos, g] = zp2sos(z, p, k);
        filteredWhiskers = filtfilt(sos, g, patchedWhisk - mean(patchedWhisk));
        resampledWhisk = resample(filteredWhiskers, downSampledFs, LabVIEWData.notes.whiskerCamSamplingRate_Hz);
        
        % Binarize the whisker waveform (wwf)
        threshfile = dir('*_Thresholds.mat');
        if ~isempty(threshfile)
            load(threshfile.name)
        end
        
        [ok] = CheckForThreshold_2P(['binarizedWhiskersLower_' strDay], animalID);
        
        if ok == 0
            [whiskersThresh1, whiskersThresh2] = CreateWhiskThreshold_2P(resampledWhisk, downSampledFs);
            Thresholds.(['binarizedWhiskersLower_' strDay]) = whiskersThresh1;
            Thresholds.(['binarizedWhiskersUpper_' strDay]) = whiskersThresh2;
            save([animalID '_Thresholds.mat'], 'Thresholds');
        end
        
        load([animalID '_Thresholds.mat']);
        binWhisk = BinarizeWhiskers_2P(resampledWhisk, downSampledFs, Thresholds.(['binarizedWhiskersLower_' strDay]), Thresholds.(['binarizedWhiskersUpper_' strDay]));
        [linkedBinarizedWhiskers] = LinkBinaryEvents_2P(gt(binWhisk,0), [round(downSampledFs/3), 0]);
        inds = linkedBinarizedWhiskers == 0;
        restAngle = mean(resampledWhisk(inds));
        
        LabVIEWData.data.dsWhiskerAngle = resampledWhisk - restAngle;
        LabVIEWData.data.binWhiskerAngle = binWhisk;
        
        
        %% Downsample and binarize the force sensor.
        trimmedForceM = ProcData.data.forceSensor(1:min(expectedLength, length(ProcData.data.forceSensor)));
        
        % Filter then downsample the Force Sensor waveform to desired frequency
        filtThreshold = 20;
        filtOrder = 2;
        [z, p, k] = butter(filtOrder, filtThreshold/(ProcData.notes.analogSamplingRate/2), 'low');
        [sos, g] = zp2sos(z, p, k);
        filtForceSensorM = filtfilt(sos, g, trimmedForceM);
        ProcData.data.dsForceSensorM = resample(filtForceSensorM, downSampledFs, ProcData.notes.analogSamplingRate);
        
        % Binarize the force sensor waveform
        threshfile = dir('*_Thresholds.mat');
        if ~isempty(threshfile)
            load(threshfile.name)
        end
        
        [ok] = CheckForThreshold_2P(['binarizedForceSensor_' strDay], animalID);
        
        if ok == 0
            [forceSensorThreshold] = CreateForceSensorThreshold_2P(ProcData.data.dsForceSensorM);
            Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
            save([animalID '_Thresholds.mat'], 'Thresholds');
        end
        
        ProcData.data.binForceSensorM = BinarizeForceSensor_2P(ProcData.data.dsForceSensorM, Thresholds.(['binarizedForceSensor_' strDay]));
        
        %% EMG
        fpass = [30 300];
        trimmedEMG = (1:min(expectedLength, length(ProcData.data.EMG)));
        [z1, p1, k1] = butter(4, fpass/(ProcData.notes.analogSamplingRate/2));
        [sos1, g1] = zp2sos(z1, p1, k1);
        filtEMG = filtfilt(sos1, g1, trimmedEMG - mean(trimmedEMG));
        [z2, p2, k2] = butter(4, 10/(ProcData.notes.analogSamplingRate/2), 'low');
        [sos2, g2] = zp2sos(z2, p2, k2);
        smoothEMGPower = filtfilt(sos2, g2, filtEMG.^2);
        ProcData.data.dsEMG = max(resample(smoothEMGPower, downSampledFs, ProcData.notes.analogSamplingRate), 0);

        %% Save the processed data
        save([animalID '_' hemisphere '_' fileID '_ProcData'], 'ProcData')
    else
        disp([rawDataFile ' has already been processed. Continuing...']); disp(' ');
    end
end

end