function [] = AddSleepParameters_SVM(procDataFileIDs, RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs: List of processed data file IDs (char) and the loaded RestingBaselines structure
%
%   Outputs: 
%
%   Last Revised: July 26th, 2019
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs, 1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Adding sleep scoring parameters to ' procDataFileID '... (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
    [~, fileDate, ~] = GetFileInfo_IOS(procDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    load(procDataFileID)
    
    %% BLOCK PURPOSE: Create folder for the Neural data of each electrode
    % delta
    LH_Delta = ProcData.data.cortical_LH.deltaBandPower;
    RH_Delta = ProcData.data.cortical_RH.deltaBandPower;
    LH_baselineDelta = RestingBaselines.cortical_LH.deltaBandPower.(strDay);
    RH_baselineDelta = RestingBaselines.cortical_RH.deltaBandPower.(strDay);
    LH_NormDelta = (LH_Delta-LH_baselineDelta)/LH_baselineDelta;
    RH_NormDelta = (RH_Delta-RH_baselineDelta)/RH_baselineDelta;
    % theta
    LH_Theta = ProcData.data.cortical_LH.thetaBandPower;
    RH_Theta = ProcData.data.cortical_RH.thetaBandPower;
    LH_baselineTheta = RestingBaselines.cortical_LH.thetaBandPower.(strDay);
    RH_baselineTheta = RestingBaselines.cortical_RH.thetaBandPower.(strDay);
    LH_NormTheta = (LH_Theta-LH_baselineTheta)/LH_baselineTheta;
    RH_NormTheta = (RH_Theta-RH_baselineTheta)/RH_baselineTheta;
    % gamma
    LH_Gamma = ProcData.data.cortical_LH.gammaBandPower;
    RH_Gamma = ProcData.data.cortical_RH.gammaBandPower;
    LH_baselineGamma = RestingBaselines.cortical_LH.gammaBandPower.(strDay);
    RH_baselineGamma = RestingBaselines.cortical_LH.gammaBandPower.(strDay);
    LH_NormGamma = (LH_Gamma-LH_baselineGamma)/LH_baselineGamma;
    RH_NormGamma = (RH_Gamma-RH_baselineGamma)/RH_baselineGamma;
    
    % Smooth the signal with a 1 Hz low pass 4th-order butterworth filter
    [B, A] = butter(4, 1/(30/2), 'low');
    LH_DeltaNeuro = filtfilt(B, A, LH_NormDelta);
    RH_DeltaNeuro = filtfilt(B, A, RH_NormDelta);
    LH_ThetaNeuro = filtfilt(B, A, LH_NormTheta);
    RH_ThetaNeuro = filtfilt(B, A, RH_NormTheta);  
    LH_GammaNeuro = filtfilt(B, A, LH_NormGamma);
    RH_GammaNeuro = filtfilt(B, A, RH_NormGamma);
    
    % Divide the neural signals into five second bins and put them in a cell array
    LH_tempDeltaStruct = cell(180,1);   
    RH_tempDeltaStruct = cell(180,1);
    LH_tempThetaStruct = cell(180,1);
    RH_tempThetaStruct = cell(180,1);
    LH_tempGammaStruct = cell(180,1);
    RH_tempGammaStruct = cell(180,1);
    
    for neuralBins = 1:180   % loop through all 9000 samples across 5 minutes in 5 second bins (180 total)
        if neuralBins == 1
            LH_tempDeltaStruct(neuralBins,1) = {LH_DeltaNeuro(neuralBins:150)};
            RH_tempDeltaStruct(neuralBins,1) = {RH_DeltaNeuro(neuralBins:150)};    
            LH_tempThetaStruct(neuralBins,1) = {LH_ThetaNeuro(neuralBins:150)};
            RH_tempThetaStruct(neuralBins,1) = {RH_ThetaNeuro(neuralBins:150)};          
            LH_tempGammaStruct(neuralBins,1) = {LH_GammaNeuro(neuralBins:150)};
            RH_tempGammaStruct(neuralBins,1) = {RH_GammaNeuro(neuralBins:150)};
        elseif neuralBins == 180
            LH_tempDeltaStruct(neuralBins,1) = {LH_NormDelta((((150*(neuralBins-1))+1)):end)};
            RH_tempDeltaStruct(neuralBins,1) = {RH_NormDelta((((150*(neuralBins-1))+1)):end)};  
            LH_tempThetaStruct(neuralBins,1) = {LH_ThetaNeuro((((150*(neuralBins-1))+1)):end)};
            RH_tempThetaStruct(neuralBins,1) = {RH_ThetaNeuro((((150*(neuralBins-1))+1)):end)};     
            LH_tempGammaStruct(neuralBins,1) = {LH_GammaNeuro((((150*(neuralBins-1))+1)):end)};
            RH_tempGammaStruct(neuralBins,1) = {RH_GammaNeuro((((150*(neuralBins-1))+1)):end)};
        else
            LH_tempDeltaStruct(neuralBins,1) = {LH_NormDelta((((150*(neuralBins-1))+1)):(150*neuralBins))};
            RH_tempDeltaStruct(neuralBins,1) = {RH_NormDelta((((150*(neuralBins-1))+1)):(150*neuralBins))};
            LH_tempThetaStruct(neuralBins,1) = {LH_ThetaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            RH_tempThetaStruct(neuralBins,1) = {RH_ThetaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};            
            LH_tempGammaStruct(neuralBins,1) = {LH_GammaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            RH_tempGammaStruct(neuralBins,1) = {RH_GammaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};   
        end
    end
    
    ProcData.sleep.parameters.deltaBandPower.LH = LH_tempDeltaStruct;
    ProcData.sleep.parameters.deltaBandPower.RH = RH_tempDeltaStruct;    
    ProcData.sleep.parameters.thetaBandPower.LH = LH_tempThetaStruct;
    ProcData.sleep.parameters.thetaBandPower.RH = RH_tempThetaStruct;   
    ProcData.sleep.parameters.gammaBandPower.LH = LH_tempGammaStruct;
    ProcData.sleep.parameters.gammaBandPower.RH = RH_tempGammaStruct;

    %% BLOCK PURPOSE: Create folder for binarized whisking and binarized force sensor
    binWhiskerAngle = ProcData.data.binWhiskerAngle;
    binForceSensor = ProcData.data.binForceSensor;
    
    % Find the number of whiskerBins due to frame drops.
    whiskerBinNumber = ceil(length(binWhiskerAngle)/150);

    % Divide the signal into five second bins and put them in a cell array
    tempWhiskerStruct = cell(whiskerBinNumber, 1);   % Pre-allocate cell array
    tempForceStruct = cell(whiskerBinNumber, 1);   % Pre-allocate cell array
    
    for whiskerBins = 1:whiskerBinNumber
        if whiskerBins == 1
            tempWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle(whiskerBins:150)};
            tempForceStruct(whiskerBins, 1) = {binForceSensor(whiskerBins:150)};
        elseif whiskerBins == whiskerBinNumber
            tempWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle((((150*(whiskerBins-1))+1)):end)};
            tempForceStruct(whiskerBins, 1) = {binForceSensor((((150*(whiskerBins-1))+1)):end)};
        else
            tempWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle((((150*(whiskerBins-1))+1)):(150*whiskerBins))};
            tempForceStruct(whiskerBins, 1) = {binForceSensor((((150*(whiskerBins-1))+1)):(150*whiskerBins))};
        end
    end
    
    ProcData.sleep.parameters.binWhiskerAngle = tempWhiskerStruct;
    ProcData.sleep.parameters.binForceSensor = tempForceStruct;
    
    %% Create folder for the EMG
    EMG = log(ProcData.data.EMG);
    tempEMGStruct = cell(180, 1);
    
    for EMGBins = 1:180
        if EMGBins == 1
            tempEMGStruct(EMGBins,1) = {EMG(EMGBins:150)};
        else
            tempEMGStruct(EMGBins,1) = {EMG((((150*(EMGBins-1))+1)):(150*EMGBins))};
        end
    end
    
    ProcData.sleep.parameters.EMG = tempEMGStruct;
    
    %% BLOCK PURPOSE: Create folder for the Heart Rate
    % Find the heart rate from the current ProcData file
    HeartRate = ProcData.data.heartRate;
    
    % Divide the signal into five second bins and put them in a cell array
    tempHRStruct = cell(180, 1);   % Pre-allocate cell array
    
    for HRBins = 1:180  % loop through all 297 samples across 5 minutes in 5 second bins (180 total)
        if HRBins == 1
            tempHRStruct(HRBins, 1) = {HeartRate(HRBins:5)};  % Samples 1 to 5
        else
            tempHRStruct(HRBins, 1) = {HeartRate((((5*(HRBins-1))+1)):(5*HRBins))};  % Samples 6 to 10, etc...
        end
    end
    
    ProcData.sleep.parameters.heartRate = tempHRStruct;   % Place the data in the ProcData struct to later be saved
    
    %% BLOCK PURPOSE: Create folder for the left and right CBV data
    LH_CBV = ProcData.data.CBV.LH;
    RH_CBV = ProcData.data.CBV.RH;
    LH_ElectrodeCBV = ProcData.data.CBV.LH_Electrode;
    RH_ElectrodeCBV = ProcData.data.CBV.RH_Electrode;
    
    LH_NormCBV = (LH_CBV-RestingBaselines.CBV.LH.(strDay))/RestingBaselines.CBV.LH.(strDay);
    RH_NormCBV = (RH_CBV-RestingBaselines.CBV.RH.(strDay))/RestingBaselines.CBV.RH.(strDay);
    LH_NormElectrodeCBV = (LH_ElectrodeCBV-RestingBaselines.CBV.LH_Electrode.(strDay))/RestingBaselines.CBV.LH_Electrode.(strDay);
    RH_NormElectrodeCBV = (RH_ElectrodeCBV-RestingBaselines.CBV.RH_Electrode.(strDay))/RestingBaselines.CBV.RH_Electrode.(strDay);
    
    [D, C] = butter(4, 1/(30/2), 'low');
    LH_FiltCBV = filtfilt(D, C, LH_NormCBV);
    RH_FiltCBV = filtfilt(D, C, RH_NormCBV);
    LH_ElectrodeFiltCBV = filtfilt(D, C, LH_NormElectrodeCBV);
    RH_ElectrodeFiltCBV = filtfilt(D, C, RH_NormElectrodeCBV);
    
    LH_tempCBVStruct = cell(180,1);   % Pre-allocate cell array
    RH_tempCBVStruct = cell(180,1);
    LH_tempElectrodeCBVStruct = cell(180,1);
    RH_tempElectrodeCBVStruct = cell(180,1);
    
    for CBVBins = 1:180   % loop through all 9000 samples across 5 minutes in 5 second bins (180 total)
        if CBVBins == 1
            LH_tempCBVStruct(CBVBins,1) = {LH_FiltCBV(CBVBins:150)};  % Samples 1 to 150
            RH_tempCBVStruct(CBVBins,1) = {RH_FiltCBV(CBVBins:150)};
            LH_tempElectrodeCBVStruct(CBVBins,1) = {LH_ElectrodeFiltCBV(CBVBins:150)};
            RH_tempElectrodeCBVStruct(CBVBins,1) = {RH_ElectrodeFiltCBV(CBVBins:150)};
        else
            LH_tempCBVStruct(CBVBins, 1) = {LH_FiltCBV((((150*(CBVBins-1))+1)):(150*CBVBins))};  % Samples 151 to 300, etc...
            RH_tempCBVStruct(CBVBins, 1) = {RH_FiltCBV((((150*(CBVBins-1))+1)):(150*CBVBins))};
            LH_tempElectrodeCBVStruct(CBVBins, 1) = {LH_ElectrodeFiltCBV((((150*(CBVBins - 1)) + 1)):(150*CBVBins))};
            RH_tempElectrodeCBVStruct(CBVBins, 1) = {RH_ElectrodeFiltCBV((((150*(CBVBins - 1)) + 1)):(150*CBVBins))};
        end
    end
    
    ProcData.sleep.parameters.CBV.LH = LH_tempCBVStruct;   % Place the data in the ProcData struct to later be saved
    ProcData.sleep.parameters.CBV.RH = RH_tempCBVStruct;
    ProcData.sleep.parameters.CBV.LH_Electrode = LH_tempElectrodeCBVStruct;
    ProcData.sleep.parameters.CBV.RH_Electrode = RH_tempElectrodeCBVStruct;
    
    save(procDataFileID, 'ProcData');
end

end
