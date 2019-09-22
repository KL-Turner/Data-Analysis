function [] = AddSleepParameters_SVM(procDataFileIDs, RestingBaselines, baselineType)
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
    hippDelta = ProcData.data.hippocampus.deltaBandPower;
    hippBaselineDelta = RestingBaselines.(baselineType).hippocampus.deltaBandPower.(strDay);
    hippNormDelta = (hippDelta-hippBaselineDelta)/hippBaselineDelta;
    
    LH_Delta = ProcData.data.cortical_LH.deltaBandPower;
    RH_Delta = ProcData.data.cortical_RH.deltaBandPower;
    LH_baselineDelta = RestingBaselines.(baselineType).cortical_LH.deltaBandPower.(strDay);
    RH_baselineDelta = RestingBaselines.(baselineType).cortical_RH.deltaBandPower.(strDay);
    LH_NormDelta = (LH_Delta-LH_baselineDelta)/LH_baselineDelta;
    RH_NormDelta = (RH_Delta-RH_baselineDelta)/RH_baselineDelta;
    
    % theta
    hippTheta = ProcData.data.hippocampus.thetaBandPower;
    hippBaselineTheta = RestingBaselines.(baselineType).hippocampus.thetaBandPower.(strDay);
    hippNormTheta = (hippTheta-hippBaselineTheta)/hippBaselineTheta;
    
    LH_Theta = ProcData.data.cortical_LH.thetaBandPower;
    RH_Theta = ProcData.data.cortical_RH.thetaBandPower;
    LH_baselineTheta = RestingBaselines.(baselineType).cortical_LH.thetaBandPower.(strDay);
    RH_baselineTheta = RestingBaselines.(baselineType).cortical_LH.thetaBandPower.(strDay);
    LH_NormTheta = (LH_Theta-LH_baselineTheta)/LH_baselineTheta;
    RH_NormTheta = (RH_Theta-RH_baselineTheta)/RH_baselineTheta;
    
    % alpha
    hippAlpha = ProcData.data.hippocampus.alphaBandPower;
    hippBaselineAlpha = RestingBaselines.(baselineType).hippocampus.alphaBandPower.(strDay);
    hippNormAlpha = (hippAlpha-hippBaselineAlpha)/hippBaselineAlpha;
    
    LH_Alpha = ProcData.data.cortical_LH.alphaBandPower;
    RH_Alpha = ProcData.data.cortical_RH.alphaBandPower;
    LH_baselineAlpha = RestingBaselines.(baselineType).cortical_LH.alphaBandPower.(strDay);
    RH_baselineAlpha = RestingBaselines.(baselineType).cortical_LH.alphaBandPower.(strDay);
    LH_NormAlpha = (LH_Alpha-LH_baselineAlpha)/LH_baselineAlpha;
    RH_NormAlpha = (RH_Alpha-RH_baselineAlpha)/RH_baselineAlpha;
    
    % beta
    hippBeta = ProcData.data.hippocampus.betaBandPower;
    hippBaselineBeta = RestingBaselines.(baselineType).hippocampus.betaBandPower.(strDay);
    hippNormBeta = (hippBeta-hippBaselineBeta)/hippBaselineBeta;
    
    LH_Beta = ProcData.data.cortical_LH.betaBandPower;
    RH_Beta = ProcData.data.cortical_RH.betaBandPower;
    LH_baselineBeta = RestingBaselines.(baselineType).cortical_LH.betaBandPower.(strDay);
    RH_baselineBeta = RestingBaselines.(baselineType).cortical_LH.betaBandPower.(strDay);
    LH_NormBeta = (LH_Beta-LH_baselineBeta)/LH_baselineBeta;
    RH_NormBeta = (RH_Beta-RH_baselineBeta)/RH_baselineBeta;
    
    % gamma
    hippGamma = ProcData.data.hippocampus.gammaBandPower;
    hippBaselineGamma = RestingBaselines.(baselineType).hippocampus.gammaBandPower.(strDay);
    hippNormGamma = (hippGamma-hippBaselineGamma)/hippBaselineGamma;
    
    LH_Gamma = ProcData.data.cortical_LH.gammaBandPower;
    RH_Gamma = ProcData.data.cortical_RH.gammaBandPower;
    LH_baselineGamma = RestingBaselines.(baselineType).cortical_LH.gammaBandPower.(strDay);
    RH_baselineGamma = RestingBaselines.(baselineType).cortical_LH.gammaBandPower.(strDay);
    LH_NormGamma = (LH_Gamma-LH_baselineGamma)/LH_baselineGamma;
    RH_NormGamma = (RH_Gamma-RH_baselineGamma)/RH_baselineGamma;
    
    % Smooth the signal with a 1 Hz low pass 4th-order butterworth filter
    [B, A] = butter(4, 1/(30/2), 'low');
    hippDeltaNeuro = filtfilt(B, A, hippNormDelta);
    hippThetaNeuro = filtfilt(B, A, hippNormTheta);
    hippAlphaNeuro = filtfilt(B, A, hippNormAlpha);
    hippBetaNeuro = filtfilt(B, A, hippNormBeta);
    hippGammaNeuro = filtfilt(B, A, hippNormGamma);
    
    LH_DeltaNeuro = filtfilt(B, A, LH_NormDelta);
    RH_DeltaNeuro = filtfilt(B, A, RH_NormDelta);
    LH_ThetaNeuro = filtfilt(B, A, LH_NormTheta);
    RH_ThetaNeuro = filtfilt(B, A, RH_NormTheta);
    LH_AlphaNeuro = filtfilt(B, A, LH_NormAlpha);
    RH_AlphaNeuro = filtfilt(B, A, RH_NormAlpha);
    LH_BetaNeuro = filtfilt(B, A, LH_NormBeta);
    RH_BetaNeuro = filtfilt(B, A, RH_NormBeta);
    LH_GammaNeuro = filtfilt(B, A, LH_NormGamma);
    RH_GammaNeuro = filtfilt(B, A, RH_NormGamma);
    
    % Divide the neural signals into five second bins and put them in a cell array
    hipptempDeltaStruct = cell(180,1);
    hipptempThetaStruct = cell(180,1);
    hipptempAlphaStruct = cell(180,1);
    hipptempBetaStruct = cell(180,1);
    hipptempGammaStruct = cell(180,1);
    
    LH_tempDeltaStruct = cell(180,1);
    RH_tempDeltaStruct = cell(180,1);
    LH_tempThetaStruct = cell(180,1);
    RH_tempThetaStruct = cell(180,1);
    LH_tempAlphaStruct = cell(180,1);
    RH_tempAlphaStruct = cell(180,1);
    LH_tempBetaStruct = cell(180,1);
    RH_tempBetaStruct = cell(180,1);
    LH_tempGammaStruct = cell(180,1);
    RH_tempGammaStruct = cell(180,1);
    
    for neuralBins = 1:180   % loop through all 9000 samples across 5 minutes in 5 second bins (180 total)
        if neuralBins == 1
            hipptempDeltaStruct(neuralBins,1) = {hippDeltaNeuro(neuralBins:150)};
            hipptempThetaStruct(neuralBins,1) = {hippThetaNeuro(neuralBins:150)};
            hipptempAlphaStruct(neuralBins,1) = {hippAlphaNeuro(neuralBins:150)};
            hipptempBetaStruct(neuralBins,1) = {hippBetaNeuro(neuralBins:150)};
            hipptempGammaStruct(neuralBins,1) = {hippGammaNeuro(neuralBins:150)};
            
            LH_tempDeltaStruct(neuralBins,1) = {LH_DeltaNeuro(neuralBins:150)};
            RH_tempDeltaStruct(neuralBins,1) = {RH_DeltaNeuro(neuralBins:150)};
            LH_tempThetaStruct(neuralBins,1) = {LH_ThetaNeuro(neuralBins:150)};
            RH_tempThetaStruct(neuralBins,1) = {RH_ThetaNeuro(neuralBins:150)};
            LH_tempAlphaStruct(neuralBins,1) = {LH_AlphaNeuro(neuralBins:150)};
            RH_tempAlphaStruct(neuralBins,1) = {RH_AlphaNeuro(neuralBins:150)};
            LH_tempBetaStruct(neuralBins,1) = {LH_BetaNeuro(neuralBins:150)};
            RH_tempBetaStruct(neuralBins,1) = {RH_BetaNeuro(neuralBins:150)};
            LH_tempGammaStruct(neuralBins,1) = {LH_GammaNeuro(neuralBins:150)};
            RH_tempGammaStruct(neuralBins,1) = {RH_GammaNeuro(neuralBins:150)};
        elseif neuralBins == 180
            hipptempDeltaStruct(neuralBins,1) = {hippDeltaNeuro((((150*(neuralBins-1))+1)):end)};
            hipptempThetaStruct(neuralBins,1) = {hippThetaNeuro((((150*(neuralBins-1))+1)):end)};
            hipptempAlphaStruct(neuralBins,1) = {hippAlphaNeuro((((150*(neuralBins-1))+1)):end)};
            hipptempBetaStruct(neuralBins,1) = {hippBetaNeuro((((150*(neuralBins-1))+1)):end)};
            hipptempGammaStruct(neuralBins,1) = {hippGammaNeuro((((150*(neuralBins-1))+1)):end)};
            
            LH_tempDeltaStruct(neuralBins,1) = {LH_DeltaNeuro((((150*(neuralBins-1))+1)):end)};
            RH_tempDeltaStruct(neuralBins,1) = {RH_DeltaNeuro((((150*(neuralBins-1))+1)):end)};
            LH_tempThetaStruct(neuralBins,1) = {LH_ThetaNeuro((((150*(neuralBins-1))+1)):end)};
            RH_tempThetaStruct(neuralBins,1) = {RH_ThetaNeuro((((150*(neuralBins-1))+1)):end)};
            LH_tempAlphaStruct(neuralBins,1) = {LH_AlphaNeuro((((150*(neuralBins-1))+1)):end)};
            RH_tempAlphaStruct(neuralBins,1) = {RH_AlphaNeuro((((150*(neuralBins-1))+1)):end)};
            LH_tempBetaStruct(neuralBins,1) = {LH_BetaNeuro((((150*(neuralBins-1))+1)):end)};
            RH_tempBetaStruct(neuralBins,1) = {RH_BetaNeuro((((150*(neuralBins-1))+1)):end)};
            LH_tempGammaStruct(neuralBins,1) = {LH_GammaNeuro((((150*(neuralBins-1))+1)):end)};
            RH_tempGammaStruct(neuralBins,1) = {RH_GammaNeuro((((150*(neuralBins-1))+1)):end)};
        else
            hipptempDeltaStruct(neuralBins,1) = {hippDeltaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            hipptempThetaStruct(neuralBins,1) = {hippThetaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            hipptempAlphaStruct(neuralBins,1) = {hippAlphaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            hipptempBetaStruct(neuralBins,1) = {hippBetaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            hipptempGammaStruct(neuralBins,1) = {hippGammaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            
            
            LH_tempDeltaStruct(neuralBins,1) = {LH_DeltaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            RH_tempDeltaStruct(neuralBins,1) = {RH_DeltaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            LH_tempThetaStruct(neuralBins,1) = {LH_ThetaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            RH_tempThetaStruct(neuralBins,1) = {RH_ThetaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            LH_tempAlphaStruct(neuralBins,1) = {LH_AlphaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            RH_tempAlphaStruct(neuralBins,1) = {RH_AlphaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            LH_tempBetaStruct(neuralBins,1) = {LH_BetaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            RH_tempBetaStruct(neuralBins,1) = {RH_BetaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            LH_tempGammaStruct(neuralBins,1) = {LH_GammaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
            RH_tempGammaStruct(neuralBins,1) = {RH_GammaNeuro((((150*(neuralBins-1))+1)):(150*neuralBins))};
        end
    end
    
    ProcData.sleep.parameters.hippocampus.deltaBandPower = hipptempDeltaStruct;
    ProcData.sleep.parameters.hippocampus.thetaBandPower = hipptempThetaStruct;
    ProcData.sleep.parameters.hippocampus.alphaBandPower = hipptempThetaStruct;
    ProcData.sleep.parameters.hippocampus.betaBandPower = hipptempBetaStruct;
    ProcData.sleep.parameters.hippocampus.gammaBandPower = hipptempGammaStruct;
    
    ProcData.sleep.parameters.cortical_LH.deltaBandPower = LH_tempDeltaStruct;
    ProcData.sleep.parameters.cortical_RH.deltaBandPower = RH_tempDeltaStruct;
    ProcData.sleep.parameters.cortical_LH.thetaBandPower = LH_tempThetaStruct;
    ProcData.sleep.parameters.cortical_RH.thetaBandPower = RH_tempThetaStruct;
    ProcData.sleep.parameters.cortical_LH.alphaBandPower = LH_tempAlphaStruct;
    ProcData.sleep.parameters.cortical_RH.alphaBandPower = RH_tempAlphaStruct;
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
    
    deltaSpecHip = specHip(deltaLowStart:deltaLowEnd, :);
    deltaSpecLH = specLH(deltaLowStart:deltaLowEnd, :);
    deltaSpecRH = specRH(deltaLowStart:deltaLowEnd, :);
    
    meanDeltaSpecHip = mean(deltaSpecHip,1);
    meanDeltaSpecLH = mean(deltaSpecLH,1);
    meanDeltaSpecRH = mean(deltaSpecRH,1);

    
    % theta
    thetaLow = freqFloor == 7;
    thetaHigh = freqFloor == 10;
    thetaLowStart = find(thetaLow, 1, 'first');
    thetaLowEnd = find(thetaHigh, 1, 'last');
    
    thetaSpecHip = specHip(thetaLowStart:thetaLowEnd, :);
    thetaSpecLH = specLH(thetaLowStart:thetaLowEnd, :);
    thetaSpecRH = specRH(thetaLowStart:thetaLowEnd, :);
    
    meanThetaSpecHip = mean(thetaSpecHip,1);
    meanThetaSpecLH = mean(thetaSpecLH,1);
    meanThetaSpecRH = mean(thetaSpecRH,1);
    
    % alpha
    alphaLow = freqFloor == 10;
    alphaHigh = freqFloor == 13;
    alphaLowStart = find(alphaLow, 1, 'first');
    alphaLowEnd = find(alphaHigh, 1, 'last');
    
    alphaSpecHip = specHip(alphaLowStart:alphaLowEnd, :);
    alphaSpecLH = specLH(alphaLowStart:alphaLowEnd, :);
    alphaSpecRH = specRH(alphaLowStart:alphaLowEnd, :);
    
    meanAlphaSpecHip = mean(alphaSpecHip,1);
    meanAlphaSpecLH = mean(alphaSpecLH,1);
    meanAlphaSpecRH = mean(alphaSpecRH,1);

    % beta
    betaLow = freqFloor == 13;
    betaHigh = freqFloor == 30;
    betaLowStart = find(betaLow, 1, 'first');
    betaLowEnd = find(betaHigh, 1, 'last');
    
    betaSpecHip = specHip(betaLowStart:betaLowEnd, :);
    betaSpecLH = specLH(betaLowStart:betaLowEnd, :);
    betaSpecRH = specRH(betaLowStart:betaLowEnd, :);
    
    meanBetaSpecHip = mean(betaSpecHip,1);
    meanBetaSpecLH = mean(betaSpecLH,1);
    meanBetaSpecRH = mean(betaSpecRH,1);

    
    % gamma
    gammaLow = freqFloor == 30;
    gammaHigh = freqFloor == 99;
    gammaLowStart = find(gammaLow, 1, 'first');
    gammaLowEnd = find(gammaHigh, 1, 'last');
    
    gammaSpecHip = specHip(gammaLowStart:gammaLowEnd, :);
    gammaSpecLH = specLH(gammaLowStart:gammaLowEnd, :);
    gammaSpecRH = specRH(gammaLowStart:gammaLowEnd, :);
    
    meanGammaSpecHip = mean(gammaSpecHip,1);
    meanGammaSpecRH = mean(gammaSpecRH,1);
    meanGammaSpecLH = mean(gammaSpecLH,1);
    
    % Divide the neural signals into five second bins and put them in a cell array
    hipptempDeltaSpecStruct = cell(180,1);
    hipptempThetaSpecStruct = cell(180,1);
    hipptempAlphaSpecStruct = cell(180,1);
    hipptempBetaSpecStruct = cell(180,1);
    hipptempGammaSpecStruct = cell(180,1);
    
    LH_tempDeltaSpecStruct = cell(180,1);
    RH_tempDeltaSpecStruct = cell(180,1);
    LH_tempThetaSpecStruct = cell(180,1);
    RH_tempThetaSpecStruct = cell(180,1);
    LH_tempAlphaSpecStruct = cell(180,1);
    RH_tempAlphaSpecStruct = cell(180,1);
    LH_tempBetaSpecStruct = cell(180,1);
    RH_tempBetaSpecStruct = cell(180,1);
    LH_tempGammaSpecStruct = cell(180,1);
    RH_tempGammaSpecStruct = cell(180,1);
    
    for neuralBins = 1:180   % loop through all 9000 samples across 5 minutes in 5 second bins (180 total)
        if neuralBins == 1
            hipptempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecHip(neuralBins:binWidth)};
            hipptempThetaSpecStruct(neuralBins,1) = {meanThetaSpecHip(neuralBins:binWidth)};
            hipptempAlphaSpecStruct(neuralBins,1) = {meanAlphaSpecHip(neuralBins:binWidth)};
            hipptempBetaSpecStruct(neuralBins,1) = {meanBetaSpecHip(neuralBins:binWidth)};
            hipptempGammaSpecStruct(neuralBins,1) = {meanGammaSpecHip(neuralBins:binWidth)};
            
            LH_tempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecLH(neuralBins:binWidth)};
            RH_tempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecRH(neuralBins:binWidth)};
            LH_tempThetaSpecStruct(neuralBins,1) = {meanThetaSpecLH(neuralBins:binWidth)};
            RH_tempThetaSpecStruct(neuralBins,1) = {meanThetaSpecRH(neuralBins:binWidth)};
            LH_tempAlphaSpecStruct(neuralBins,1) = {meanAlphaSpecLH(neuralBins:binWidth)};
            RH_tempAlphaSpecStruct(neuralBins,1) = {meanAlphaSpecRH(neuralBins:binWidth)};
            LH_tempBetaSpecStruct(neuralBins,1) = {meanBetaSpecLH(neuralBins:binWidth)};
            RH_tempBetaSpecStruct(neuralBins,1) = {meanBetaSpecRH(neuralBins:binWidth)};
            LH_tempGammaSpecStruct(neuralBins,1) = {meanGammaSpecLH(neuralBins:binWidth)};
            RH_tempGammaSpecStruct(neuralBins,1) = {meanGammaSpecRH(neuralBins:binWidth)};
        elseif neuralBins == 180
            hipptempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecHip((((binWidth*(neuralBins-1))+1)):end)};
            hipptempThetaSpecStruct(neuralBins,1) = {meanThetaSpecHip((((binWidth*(neuralBins-1))+1)):end)};
            hipptempAlphaSpecStruct(neuralBins,1) = {meanAlphaSpecHip((((binWidth*(neuralBins-1))+1)):end)};
            hipptempBetaSpecStruct(neuralBins,1) = {meanBetaSpecHip((((binWidth*(neuralBins-1))+1)):end)};
            hipptempGammaSpecStruct(neuralBins,1) = {meanGammaSpecHip((((binWidth*(neuralBins-1))+1)):end)};
            
            LH_tempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecLH((((binWidth*(neuralBins-1))+1)):end)};
            RH_tempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecRH((((binWidth*(neuralBins-1))+1)):end)};
            LH_tempThetaSpecStruct(neuralBins,1) = {meanThetaSpecLH((((binWidth*(neuralBins-1))+1)):end)};
            RH_tempThetaSpecStruct(neuralBins,1) = {meanThetaSpecRH((((binWidth*(neuralBins-1))+1)):end)};
            LH_tempAlphaSpecStruct(neuralBins,1) = {meanAlphaSpecLH((((binWidth*(neuralBins-1))+1)):end)};
            RH_tempAlphaSpecStruct(neuralBins,1) = {meanAlphaSpecRH((((binWidth*(neuralBins-1))+1)):end)};
            LH_tempBetaSpecStruct(neuralBins,1) = {meanBetaSpecLH((((binWidth*(neuralBins-1))+1)):end)};
            RH_tempBetaSpecStruct(neuralBins,1) = {meanBetaSpecRH((((binWidth*(neuralBins-1))+1)):end)};
            LH_tempGammaSpecStruct(neuralBins,1) = {meanGammaSpecLH((((binWidth*(neuralBins-1))+1)):end)};
            RH_tempGammaSpecStruct(neuralBins,1) = {meanGammaSpecRH((((binWidth*(neuralBins-1))+1)):end)};
        else
            hipptempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecHip((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            hipptempThetaSpecStruct(neuralBins,1) = {meanThetaSpecHip((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            hipptempAlphaSpecStruct(neuralBins,1) = {meanAlphaSpecHip((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            hipptempBetaSpecStruct(neuralBins,1) = {meanBetaSpecHip((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            hipptempGammaSpecStruct(neuralBins,1) = {meanGammaSpecHip((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            
            LH_tempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecLH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            RH_tempDeltaSpecStruct(neuralBins,1) = {meanDeltaSpecRH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            LH_tempThetaSpecStruct(neuralBins,1) = {meanThetaSpecLH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            RH_tempThetaSpecStruct(neuralBins,1) = {meanThetaSpecRH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            LH_tempAlphaSpecStruct(neuralBins,1) = {meanAlphaSpecLH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            RH_tempAlphaSpecStruct(neuralBins,1) = {meanAlphaSpecRH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            LH_tempBetaSpecStruct(neuralBins,1) = {meanBetaSpecLH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            RH_tempBetaSpecStruct(neuralBins,1) = {meanBetaSpecRH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            LH_tempGammaSpecStruct(neuralBins,1) = {meanGammaSpecLH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
            RH_tempGammaSpecStruct(neuralBins,1) = {meanGammaSpecRH((((binWidth*(neuralBins-1))+1)):(binWidth*neuralBins))};
        end
    end
    
    ProcData.sleep.parameters.hippocampus.specDeltaBandPower = hipptempDeltaSpecStruct;
    ProcData.sleep.parameters.hippocampus.specThetaBandPower = hipptempThetaSpecStruct;
    ProcData.sleep.parameters.hippocampus.specAlphaBandPower = hipptempAlphaSpecStruct;
    ProcData.sleep.parameters.hippocampus.specBetaBandPower = hipptempBetaSpecStruct;
    ProcData.sleep.parameters.hippocampus.specGammaBandPower = hipptempGammaSpecStruct;
    
    ProcData.sleep.parameters.cortical_LH.specDeltaBandPower = LH_tempDeltaSpecStruct;
    ProcData.sleep.parameters.cortical_RH.specDeltaBandPower = RH_tempDeltaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specThetaBandPower = LH_tempThetaSpecStruct;
    ProcData.sleep.parameters.cortical_RH.specThetaBandPower = RH_tempThetaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specAlphaBandPower = LH_tempAlphaSpecStruct;
    ProcData.sleep.parameters.cortical_RH.specAlphaBandPower = RH_tempAlphaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specBetaBandPower = LH_tempBetaSpecStruct;
    ProcData.sleep.parameters.cortical_RH.specBetaBandPower = RH_tempBetaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specGammaBandPower = LH_tempGammaSpecStruct;
    ProcData.sleep.parameters.cortical_RH.specGammaBandPower = RH_tempGammaSpecStruct;
    
    %% BLOCK PURPOSE: Create folder for binarized whisking and binarized force sensor
    binWhiskerAngle = ProcData.data.binWhiskerAngle;
    binForceSensor = ProcData.data.binForceSensor;
    
    whiskerAngle = ProcData.data.whiskerAngle;
    whiskerAcceleration = diff(whiskerAngle, 2);
    
    % Find the number of whiskerBins due to frame drops.
    whiskerBinNumber = ceil(length(binWhiskerAngle)/150);
    
    % Divide the signal into five second bins and put them in a cell array
    tempWhiskerStruct = cell(whiskerBinNumber, 1);   % Pre-allocate cell array
    tempWhiskerAccelStruct = cell(whiskerBinNumber, 1);   % Pre-allocate cell array
    tempBinWhiskerStruct = cell(whiskerBinNumber, 1);   % Pre-allocate cell array
    tempForceStruct = cell(whiskerBinNumber, 1);   % Pre-allocate cell array
    
    for whiskerBins = 1:whiskerBinNumber
        if whiskerBins == 1
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle(whiskerBins:150)};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration(whiskerBins:150)};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle(whiskerBins:150)};
            tempForceStruct(whiskerBins, 1) = {binForceSensor(whiskerBins:150)};
        elseif whiskerBins == whiskerBinNumber
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle((((150*(whiskerBins-1))+1)):end)};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration((((150*(whiskerBins-1))+1)):end)};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle((((150*(whiskerBins-1))+1)):end)};
            tempForceStruct(whiskerBins, 1) = {binForceSensor((((150*(whiskerBins-1))+1)):end)};
        else
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle((((150*(whiskerBins-1))+1)):(150*whiskerBins))};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration((((150*(whiskerBins-1))+1)):(150*whiskerBins))};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle((((150*(whiskerBins-1))+1)):(150*whiskerBins))};
            tempForceStruct(whiskerBins, 1) = {binForceSensor((((150*(whiskerBins-1))+1)):(150*whiskerBins))};
        end
    end
    
    ProcData.sleep.parameters.whiskerAngle = tempWhiskerStruct;
    ProcData.sleep.parameters.whiskerAcceleration = tempWhiskerAccelStruct;
    ProcData.sleep.parameters.binWhiskerAngle = tempBinWhiskerStruct;
    ProcData.sleep.parameters.binForceSensor = tempForceStruct;
    
    %% Create folder for the EMG
    EMG = ProcData.data.EMG.emg;
    normEMG = (EMG-RestingBaselines.(baselineType).EMG.emg.(strDay))/RestingBaselines.(baselineType).EMG.emg.(strDay);  
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
    
    LH_NormCBV = (LH_CBV-RestingBaselines.(baselineType).CBV.LH.(strDay))/RestingBaselines.(baselineType).CBV.LH.(strDay);
    RH_NormCBV = (RH_CBV-RestingBaselines.(baselineType).CBV.RH.(strDay))/RestingBaselines.(baselineType).CBV.RH.(strDay);
    LH_NormElectrodeCBV = (LH_ElectrodeCBV-RestingBaselines.(baselineType).CBV.LH_Electrode.(strDay))/RestingBaselines.(baselineType).CBV.LH_Electrode.(strDay);
    RH_NormElectrodeCBV = (RH_ElectrodeCBV-RestingBaselines.(baselineType).CBV.RH_Electrode.(strDay))/RestingBaselines.(baselineType).CBV.RH_Electrode.(strDay);
    
    [D, C] = butter(4, 1/(30/2), 'low');
    ledType = 'M530L3';
    bandfilterType = 'FB530-10';
    cutfilterType = 'EO46540';
    conv2um = 1e6;
    [~,~,weightedcoeffHbT] = getHbcoeffs_IOS(ledType,bandfilterType,cutfilterType);
    
    hbtLH_CBV = (log(LH_CBV/RestingBaselines.(baselineType).CBV.LH.(strDay)))*weightedcoeffHbT*conv2um;
    filtHbtLH_CBV = detrend(filtfilt(D, C, hbtLH_CBV), 'constant');
    hbtRH_CBV = (log(RH_CBV/ RestingBaselines.(baselineType).CBV.RH.(strDay)))*weightedcoeffHbT*conv2um;
    filtHbtRH_CBV = detrend(filtfilt(D, C, hbtRH_CBV), 'constant');
    
    hbtLH_ElectrodeCBV = (log(LH_CBV/RestingBaselines.(baselineType).CBV.LH_Electrode.(strDay)))*weightedcoeffHbT*conv2um;
    filtHbtLH_ElectrodeCBV = detrend(filtfilt(D, C, hbtLH_ElectrodeCBV), 'constant');
    hbtRH_ElectrodeCBV = (log(RH_ElectrodeCBV/ RestingBaselines.(baselineType).CBV.RH_Electrode.(strDay)))*weightedcoeffHbT*conv2um;
    filtHbtRH_ElectrodeCBV = detrend(filtfilt(D, C, hbtRH_ElectrodeCBV), 'constant');
    
    LH_FiltCBV = detrend(filtfilt(D, C, LH_NormCBV), 'constant');
    RH_FiltCBV = detrend(filtfilt(D, C, RH_NormCBV), 'constant');
    LH_ElectrodeFiltCBV = detrend(filtfilt(D, C, LH_NormElectrodeCBV), 'constant');
    RH_ElectrodeFiltCBV = detrend(filtfilt(D, C, RH_NormElectrodeCBV), 'constant');
    
    LH_tempCBVStruct = cell(180,1);   % Pre-allocate cell array
    RH_tempCBVStruct = cell(180,1);
    LH_tempElectrodeCBVStruct = cell(180,1);
    RH_tempElectrodeCBVStruct = cell(180,1);
    
    hbtLH_tempCBVStruct = cell(180,1);   % Pre-allocate cell array
    hbtRH_tempCBVStruct = cell(180,1);
    hbtLH_tempElectrodeCBVStruct = cell(180,1);
    hbtRH_tempElectrodeCBVStruct = cell(180,1);
    
    for CBVBins = 1:180   % loop through all 9000 samples across 5 minutes in 5 second bins (180 total)
        if CBVBins == 1
            LH_tempCBVStruct(CBVBins,1) = {LH_FiltCBV(CBVBins:150)};  % Samples 1 to 150
            RH_tempCBVStruct(CBVBins,1) = {RH_FiltCBV(CBVBins:150)};
            LH_tempElectrodeCBVStruct(CBVBins,1) = {LH_ElectrodeFiltCBV(CBVBins:150)};
            RH_tempElectrodeCBVStruct(CBVBins,1) = {RH_ElectrodeFiltCBV(CBVBins:150)};
            
            hbtLH_tempCBVStruct(CBVBins,1) = {filtHbtLH_CBV(CBVBins:150)};  % Samples 1 to 150
            hbtRH_tempCBVStruct(CBVBins,1) = {filtHbtRH_CBV(CBVBins:150)};
            hbtLH_tempElectrodeCBVStruct(CBVBins,1) = {filtHbtLH_ElectrodeCBV(CBVBins:150)};
            hbtRH_tempElectrodeCBVStruct(CBVBins,1) = {filtHbtRH_ElectrodeCBV(CBVBins:150)};
        else
            LH_tempCBVStruct(CBVBins, 1) = {LH_FiltCBV((((150*(CBVBins-1))+1)):(150*CBVBins))};  % Samples 151 to 300, etc...
            RH_tempCBVStruct(CBVBins, 1) = {RH_FiltCBV((((150*(CBVBins-1))+1)):(150*CBVBins))};
            LH_tempElectrodeCBVStruct(CBVBins, 1) = {LH_ElectrodeFiltCBV((((150*(CBVBins - 1)) + 1)):(150*CBVBins))};
            RH_tempElectrodeCBVStruct(CBVBins, 1) = {RH_ElectrodeFiltCBV((((150*(CBVBins - 1)) + 1)):(150*CBVBins))};
            
            hbtLH_tempCBVStruct(CBVBins, 1) = {filtHbtLH_CBV((((150*(CBVBins-1))+1)):(150*CBVBins))};  % Samples 151 to 300, etc...
            hbtRH_tempCBVStruct(CBVBins, 1) = {filtHbtRH_CBV((((150*(CBVBins-1))+1)):(150*CBVBins))};
            hbtLH_tempElectrodeCBVStruct(CBVBins, 1) = {filtHbtLH_ElectrodeCBV((((150*(CBVBins - 1)) + 1)):(150*CBVBins))};
            hbtRH_tempElectrodeCBVStruct(CBVBins, 1) = {filtHbtRH_ElectrodeCBV((((150*(CBVBins - 1)) + 1)):(150*CBVBins))};
        end
    end
    
    ProcData.sleep.parameters.CBV.LH = LH_tempCBVStruct;   % Place the data in the ProcData struct to later be saved
    ProcData.sleep.parameters.CBV.RH = RH_tempCBVStruct;
    ProcData.sleep.parameters.CBV.LH_Electrode = LH_tempElectrodeCBVStruct;
    ProcData.sleep.parameters.CBV.RH_Electrode = RH_tempElectrodeCBVStruct;
    
    ProcData.sleep.parameters.CBV.hbtLH = hbtLH_tempCBVStruct;   % Place the data in the ProcData struct to later be saved
    ProcData.sleep.parameters.CBV.hbtRH = hbtRH_tempCBVStruct;
    ProcData.sleep.parameters.CBV.hbtLH_Electrode = hbtLH_tempElectrodeCBVStruct;
    ProcData.sleep.parameters.CBV.hbtRH_Electrode = hbtRH_tempElectrodeCBVStruct;
    
    save(procDataFileID, 'ProcData');
end

end
