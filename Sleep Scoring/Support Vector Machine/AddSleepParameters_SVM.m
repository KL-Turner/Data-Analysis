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
    specDataFileID = [procDataFileID(1:end-12) 'SpecData.mat'];
    load(specDataFileID)
    
    %% BLOCK PURPOSE: Create folder for the Neural data of each electrode
    % delta
    LH_Delta = ProcData.data.cortical_LH.deltaBandPower;
    RH_Delta = ProcData.data.cortical_RH.deltaBandPower;
    LH_baselineDelta = RestingBaselines.cortical_LH.deltaBandPower.(strDay);
    RH_baselineDelta = RestingBaselines.cortical_RH.deltaBandPower.(strDay);
    LH_NormDelta = (LH_Delta-LH_baselineDelta)/LH_baselineDelta;
    RH_NormDelta = (RH_Delta-RH_baselineDelta)/RH_baselineDelta;
    % theta
    hippTheta = ProcData.data.hippocampus.thetaBandPower;
    hippBaselineTheta = RestingBaselines.hippocampus.thetaBandPower.(strDay);
    hippNormTheta = (hippTheta-hippBaselineTheta)/hippBaselineTheta;
    % beta
    LH_Beta = ProcData.data.cortical_LH.betaBandPower;
    RH_Beta = ProcData.data.cortical_RH.betaBandPower;
    LH_baselineBeta = RestingBaselines.cortical_LH.betaBandPower.(strDay);
    RH_baselineBeta = RestingBaselines.cortical_LH.betaBandPower.(strDay);
    LH_NormBeta = (LH_Beta-LH_baselineBeta)/LH_baselineBeta;
    RH_NormBeta = (RH_Beta-RH_baselineBeta)/RH_baselineBeta;
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
    hippThetaNeuro = filtfilt(B, A, hippNormTheta);
    LH_BetaNeuro = filtfilt(B, A, LH_NormBeta);
    RH_BetaNeuro = filtfilt(B, A, RH_NormBeta);
    LH_GammaNeuro = filtfilt(B, A, LH_NormGamma);
    RH_GammaNeuro = filtfilt(B, A, RH_NormGamma);
    
    % Divide the neural signals into five second bins and put them in a cell array
    LH_tempDeltaStruct = cell(180,1);
    RH_tempDeltaStruct = cell(180,1);
    hipptempThetaStruct = cell(180,1);
    LH_tempBetaStruct = cell(180,1);
    RH_tempBetaStruct = cell(180,1);
    LH_tempGammaStruct = cell(180,1);
    RH_tempGammaStruct = cell(180,1);
    
    for neuralBins = 1:180   % loop through all 9000 samples across 5 minutes in 5 second bins (180 total)
        if neuralBins == 1
            LH_tempDeltaStruct(neuralBins,1) = {LH_DeltaNeuro(neuralBins:150)};
            RH_tempDeltaStruct(neuralBins,1) = {RH_DeltaNeuro(neuralBins:150)};
            hipptempThetaStruct(neuralBins,1) = {hippThetaNeuro(neuralBins:150)};
            LH_tempBetaStruct(neuralBins,1) = {LH_BetaNeuro(neuralBins:150)};
            RH_tempBetaStruct(neuralBins,1) = {RH_BetaNeuro(neuralBins:150)};
            LH_tempGammaStruct(neuralBins,1) = {LH_GammaNeuro(neuralBins:150)};
            RH_tempGammaStruct(neuralBins,1) = {RH_GammaNeuro(neuralBins:150)};
        elseif neuralBins == 180
            LH_tempDeltaStruct(neuralBins,1) = {LH_DeltaNeuro((((150*(neuralBins-1))+1)):end)};
            RH_tempDeltaStruct(neuralBins,1) = {RH_DeltaNeuro((((150*(neuralBins-1))+1)):end)};
            hipptempThetaStruct(neuralBins,1) = {hippThetaNeuro((((150*(neuralBins-1))+1)):end)};
            LH_tempBetaStruct(neuralBins,1) = {LH_BetaNeuro((((150*(neuralBins-1))+1)):end)};
            RH_tempBetaStruct(neuralBins,1) = {RH_BetaNeuro((((150*(neuralBins-1))+1)):end)};
            LH_tempGammaStruct(neuralBins,1) = {LH_GammaNeuro((((150*(neuralBins-1))+1)):end)};
            RH_tempGammaStruct(neuralBins,1) = {RH_GammaNeuro((((150*(neuralBins-1))+1)):end)};
        else
            LH_tempDeltaStruct(neuralBins,1) = {LH_DeltaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            RH_tempDeltaStruct(neuralBins,1) = {RH_DeltaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            hipptempThetaStruct(neuralBins,1) = {hippThetaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            LH_tempBetaStruct(neuralBins,1) = {LH_BetaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            RH_tempBetaStruct(neuralBins,1) = {RH_BetaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            LH_tempGammaStruct(neuralBins,1) = {LH_GammaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            RH_tempGammaStruct(neuralBins,1) = {RH_GammaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
        end
    end
    
    ProcData.sleep.parameters.cortical_LH.deltaBandPower = LH_tempDeltaStruct;
    ProcData.sleep.parameters.cortical_RH.deltaBandPower = RH_tempDeltaStruct;
    ProcData.sleep.parameters.hippocampus.thetaBandPower = hipptempThetaStruct;
    ProcData.sleep.parameters.cortical_LH.betaBandPower = LH_tempBetaStruct;
    ProcData.sleep.parameters.cortical_RH.betaBandPower = RH_tempBetaStruct;
    ProcData.sleep.parameters.cortical_LH.gammaBandPower = LH_tempGammaStruct;
    ProcData.sleep.parameters.cortical_RH.gammaBandPower = RH_tempGammaStruct;
    
    %% BLOCK PURPOSE: Create folder for the Neural spectrogram data of each electrode
    T = SpecData.cortical_LH.fiveSec.T;
    F = SpecData.cortical_LH.fiveSec.F;
    specLH = SpecData.cortical_LH.fiveSec.normS;
    specRH = SpecData.cortical_RH.fiveSec.normS;
    specHip = SpecData.hippocampus.fiveSec.normS;
    binWidth = ceil((length(T)/900))*5;
    freqFloor = floor(F);
    
    % delta
    deltaLow = freqFloor == 1;
    deltaHigh = freqFloor == 4;
    deltaLowStart = find(deltaLow, 1, 'first');
    deltaLowEnd = find(deltaHigh, 1, 'last');
    deltaSpecLH = specLH(deltaLowStart:deltaLowEnd, :);
    meanDeltaSpecLH = mean(deltaSpecLH,1);
    deltaSpecRH = specRH(deltaLowStart:deltaLowEnd, :);
    meanDeltaSpecRH = mean(deltaSpecRH,1);

    % theta
    thetaLow = freqFloor == 4;
    thetaHigh = freqFloor == 10;
    thetaLowStart = find(thetaLow, 1, 'first');
    thetaLowEnd = find(thetaHigh, 1, 'last');
    thetaSpecHip = specHip(thetaLowStart:thetaLowEnd, :);
    meanThetaSpecHip = mean(thetaSpecHip,1);
    
    % beta
    betaLow = freqFloor == 13;
    betaHigh = freqFloor == 30;
    betaLowStart = find(betaLow, 1, 'first');
    betaLowEnd = find(betaHigh, 1, 'last');
    betaSpecLH = specLH(betaLowStart:betaLowEnd, :);
    meanBetaSpecLH = mean(betaSpecLH,1);
    betaSpecRH = specRH(betaLowStart:betaLowEnd, :);
    meanBetaSpecRH = mean(betaSpecRH,1);

    % gamma
    gammaLow = freqFloor == 30;
    gammaHigh = freqFloor == 99;
    gammaLowStart = find(gammaLow, 1, 'first');
    gammaLowEnd = find(gammaHigh, 1, 'last');
    gammaSpecLH = specLH(gammaLowStart:gammaLowEnd, :);
    meanGammaSpecLH = mean(gammaSpecLH,1);
    gammaSpecRH = specRH(gammaLowStart:gammaLowEnd, :);
    meanGammaSpecRH = mean(gammaSpecRH,1);

    % Divide the neural signals into five second bins and put them in a cell array
    LH_tempDeltaSpecStruct = cell(180,1);
    RH_tempDeltaSpecStruct = cell(180,1);
    hipptempThetaSpecStruct = cell(180,1);
    LH_tempBetaSpecStruct = cell(180,1);
    RH_tempBetaSpecStruct = cell(180,1);
    LH_tempGammaSpecStruct = cell(180,1);
    RH_tempGammaSpecStruct = cell(180,1);
    
    for neuralBins = 1:180   % loop through all 9000 samples across 5 minutes in 5 second bins (180 total)
        if neuralBins == 1
            LH_tempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecLH(neuralBins:binWidth)};
            RH_tempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecRH(neuralBins:binWidth)};
            hipptempThetaSpecStruct(neuralBins,1) = {meanThetaSpecHip(neuralBins:binWidth)};
            LH_tempBetaSpecStruct(neuralBins,1) = {meanBetaSpecLH(neuralBins:binWidth)};
            RH_tempBetaSpecStruct(neuralBins,1) = {meanBetaSpecRH(neuralBins:binWidth)};
            LH_tempGammaSpecStruct(neuralBins,1) = {meanGammaSpecLH(neuralBins:binWidth)};
            RH_tempGammaSpecStruct(neuralBins,1) = {meanGammaSpecRH(neuralBins:binWidth)};
        elseif neuralBins == 180
            LH_tempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecLH((((binWidth*(neuralBins-1))+1)):end)};
            RH_tempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecRH((((binWidth*(neuralBins-1))+1)):end)};
            hipptempThetaSpecStruct(neuralBins,1) = {meanThetaSpecHip((((binWidth*(neuralBins-1))+1)):end)};
            LH_tempBetaSpecStruct(neuralBins,1) = {meanBetaSpecLH((((binWidth*(neuralBins-1))+1)):end)};
            RH_tempBetaSpecStruct(neuralBins,1) = {meanBetaSpecRH((((binWidth*(neuralBins-1))+1)):end)};
            LH_tempGammaSpecStruct(neuralBins,1) = {meanGammaSpecLH((((binWidth*(neuralBins-1))+1)):end)};
            RH_tempGammaSpecStruct(neuralBins,1) = {meanGammaSpecRH((((binWidth*(neuralBins-1))+1)):end)};
        else
            LH_tempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecLH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            RH_tempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecRH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            hipptempThetaSpecStruct(neuralBins,1) = {meanThetaSpecHip((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            LH_tempBetaSpecStruct(neuralBins,1) = {meanBetaSpecLH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            RH_tempBetaSpecStruct(neuralBins,1) = {meanBetaSpecRH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            LH_tempGammaSpecStruct(neuralBins,1) = {meanGammaSpecLH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            RH_tempGammaSpecStruct(neuralBins,1) = {meanGammaSpecRH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
        end
    end
    
    ProcData.sleep.parameters.cortical_LH.specDeltaBandPower = LH_tempDeltaSpecStruct;
    ProcData.sleep.parameters.cortical_RH.specDeltaBandPower = RH_tempDeltaSpecStruct;
    ProcData.sleep.parameters.hippocampus.specThetaBandPower = hipptempThetaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specBetaBandPower = LH_tempBetaSpecStruct;
    ProcData.sleep.parameters.cortical_RH.specBetaBandPower = RH_tempBetaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specGammaBandPower = LH_tempGammaSpecStruct;
    ProcData.sleep.parameters.cortical_RH.specGammaBandPower = RH_tempGammaSpecStruct;
    
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
    EMG = ProcData.data.EMG.emg;
    normEMG = (EMG-RestingBaselines.EMG.emg.(strDay))/RestingBaselines.EMG.emg.(strDay);  
    tempEMGStruct = cell(180, 1);
    
    for EMGBins = 1:180
        if EMGBins == 1
            tempEMGStruct(EMGBins,1) = {normEMG(EMGBins:150)};
        else
            tempEMGStruct(EMGBins,1) = {normEMG((((150*(EMGBins-1))+1)):(150*EMGBins))};
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
    LH_FiltCBV = detrend(filtfilt(D, C, LH_NormCBV), 'constant');
    RH_FiltCBV = detrend(filtfilt(D, C, RH_NormCBV), 'constant');
    LH_ElectrodeFiltCBV = detrend(filtfilt(D, C, LH_NormElectrodeCBV), 'constant');
    RH_ElectrodeFiltCBV = detrend(filtfilt(D, C, RH_NormElectrodeCBV), 'constant');
    
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
