function [ProcData] = AddSleepParameters(animal, hem, fileID, strDay, RawData, ProcData, RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: This function serves to add the relevant sleep parameters to each ProcData file, the most notable being the
%            delta power that will be normalized by the current day's DeltaBandBaseline, which was obtained from the 
%            whisking events. The other parameters include the raw CBV from both hemispheres, the raw neural data from
%            both hemispheres (as well as the non-normalized Gamma band power), the Heart Rate, and whisker acceleration.
%________________________________________________________________________________________________________________________
%
%   Inputs: ProcData and RawData files from the same 5 minute imaging session. The RawData file contains the raw neural
%           data, while the ProcData file contains the Delta and Gamma bands, and the whisker position used to calculate
%           the acceleration. The raw CBV is in both files - either can be used. The Heart Rate is then calculated from 
%           the raw CBV.
%
%   Outputs: 5 second bins of each of the previously mentioned parameters that are either relevenat for sleep scoring, 
%            or for post-scoring analysis.
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: Create folder for the Neural Data of each electrode
LH_Delta = ProcData.Data.DeltaBand_Power.LH;   % Left electrode Delta power signal
RH_Delta = ProcData.Data.DeltaBand_Power.RH;   % Right electrode Delta power signal

LH_Theta = ProcData.Data.ThetaBand_Power.LH;   % Left electrode Theta power signal
RH_Theta = ProcData.Data.ThetaBand_Power.RH;   % Right electrode Theta power signal

LH_Gamma = ProcData.Data.GammaBand_Power.LH;   % Left electrode Gamma power signal
RH_Gamma = ProcData.Data.GammaBand_Power.RH;   % Right electrode Gamma power signal

LH_Neuro = RawData.Data.Neural_LH;   % Left electrode raw data (full signal)
RH_Neuro = RawData.Data.Neural_RH;   % Right electrode raw data (full signal)

% Normalize the filtered signal by the power during resting baseline
LH_baselineDelta = RestingBaselines.DeltaBand_Power.LH.(strDay);   % Take the average for the specific day of this particular file (left)
RH_baselineDelta = RestingBaselines.DeltaBand_Power.RH.(strDay);   % Take the average for the specific day of this particular file (right)

LH_NormDelta = (LH_Delta - LH_baselineDelta) / LH_baselineDelta;   % Normalize left electrode by baseline delta power
RH_NormDelta = (RH_Delta - RH_baselineDelta) / RH_baselineDelta;   % Normalize right electrode by baseline delta power

LH_baselineTheta = RestingBaselines.ThetaBand_Power.LH.(strDay);   % Take the average for the specific day of this particular file (left)
RH_baselineTheta = RestingBaselines.ThetaBand_Power.RH.(strDay);   % Take the average for the specific day of this particular file (right)

LH_NormTheta = (LH_Theta - LH_baselineTheta) / LH_baselineTheta;   % Normalize left electrode by baseline theta power
RH_NormTheta = (RH_Theta - RH_baselineTheta) / RH_baselineTheta;   % Normalize right electrode by baseline theta power

LH_baselineGamma = RestingBaselines.GammaBand_Power.LH.(strDay);   % Take the average for the specific day of this particular file (left)
RH_baselineGamma = RestingBaselines.GammaBand_Power.RH.(strDay);   % Take the average for the specific day of this particular file (right)

LH_NormGamma = (LH_Gamma - LH_baselineGamma) / LH_baselineGamma;   % Normalize left electrode by baseline gamma power
RH_NormGamma = (RH_Gamma - RH_baselineGamma) / RH_baselineGamma;   % Normalize right electrode by baseline gamma power

% Smooth the signal with a 1 Hz low pass 4th-order butterworth filter
% Sampling Rate is 30 Hz for Delta, Theta, and Gamma signals
[B, A] = butter(4, 1 / (30 / 2), 'low');  
LH_DeltaNeuro = filtfilt(B, A, LH_NormDelta);   % Filtered left electrode Delta power signal
RH_DeltaNeuro = filtfilt(B, A, RH_NormDelta);   % Filtered right electrode Delta power signal

LH_ThetaNeuro = filtfilt(B, A, LH_NormTheta);   % Filtered left electrode Delta power signal
RH_ThetaNeuro = filtfilt(B, A, RH_NormTheta);   % Filtered right electrode Delta power signal

LH_GammaNeuro = filtfilt(B, A, LH_NormGamma);   % Filtered left electrode Gamma power signal
RH_GammaNeuro = filtfilt(B, A, RH_NormGamma);   % Filtered right electrode Gamma power signal

% Divide the neural signals into five second bins and put them in a cell array
LH_tempDeltaStruct = cell(60, 1);   % Pre-allocate cell array
RH_tempDeltaStruct = cell(60, 1);

LH_tempThetaStruct = cell(60, 1);
RH_tempThetaStruct = cell(60, 1);

LH_tempGammaStruct = cell(60, 1);
RH_tempGammaStruct = cell(60, 1);

LH_tempRawNeuroStruct = cell(60, 1);
RH_tempRawNeuroStruct = cell(60, 1);

for neuralBins = 1:60   % loop through all 9000 samples across 5 minutes in 5 second bins (60 total)
    if neuralBins == 1
        LH_tempDeltaStruct(neuralBins, 1) = {LH_NormDelta(neuralBins:150)};   % Samples 1 to 150
        RH_tempDeltaStruct(neuralBins, 1) = {RH_NormDelta(neuralBins:150)};
        
        LH_tempThetaStruct(neuralBins, 1) = {LH_ThetaNeuro(neuralBins:150)};
        RH_tempThetaStruct(neuralBins, 1) = {RH_ThetaNeuro(neuralBins:150)};
        
        LH_tempGammaStruct(neuralBins, 1) = {LH_GammaNeuro(neuralBins:150)};
        RH_tempGammaStruct(neuralBins, 1) = {RH_GammaNeuro(neuralBins:150)};
        
        LH_tempRawNeuroStruct(neuralBins, 1) = {LH_Neuro(neuralBins:100000)};
        RH_tempRawNeuroStruct(neuralBins, 1) = {RH_Neuro(neuralBins:100000)};
        
    elseif neuralBins == 60
        LH_tempDeltaStruct(neuralBins, 1) = {LH_NormDelta((((150*(neuralBins - 1)) + 1)):end)};   % Samples 151 to 300, etc...
        RH_tempDeltaStruct(neuralBins, 1) = {RH_NormDelta((((150*(neuralBins - 1)) + 1)):end)};
        
        LH_tempThetaStruct(neuralBins, 1) = {LH_ThetaNeuro((((150*(neuralBins - 1)) + 1)):end)};
        RH_tempThetaStruct(neuralBins, 1) = {RH_ThetaNeuro((((150*(neuralBins - 1)) + 1)):end)};
        
        LH_tempGammaStruct(neuralBins, 1) = {LH_GammaNeuro((((150*(neuralBins - 1)) + 1)):end)};
        RH_tempGammaStruct(neuralBins, 1) = {RH_GammaNeuro((((150*(neuralBins - 1)) + 1)):end)};
        
        LH_tempRawNeuroStruct(neuralBins, 1) = {LH_Neuro((((100000*(neuralBins - 1)) + 1)):end)};
        RH_tempRawNeuroStruct(neuralBins, 1) = {RH_Neuro((((100000*(neuralBins - 1)) + 1)):end)};
        
    else
        LH_tempDeltaStruct(neuralBins, 1) = {LH_NormDelta((((150*(neuralBins - 1)) + 1)):(150*neuralBins))};  % Samples 151 to 300, etc...
        RH_tempDeltaStruct(neuralBins, 1) = {RH_NormDelta((((150*(neuralBins - 1)) + 1)):(150*neuralBins))};
        
        LH_tempThetaStruct(neuralBins, 1) = {LH_ThetaNeuro((((150*(neuralBins - 1)) + 1)):(150*neuralBins))};
        RH_tempThetaStruct(neuralBins, 1) = {RH_ThetaNeuro((((150*(neuralBins - 1)) + 1)):(150*neuralBins))};
        
        LH_tempGammaStruct(neuralBins, 1) = {LH_GammaNeuro((((150*(neuralBins - 1)) + 1)):(150*neuralBins))};
        RH_tempGammaStruct(neuralBins, 1) = {RH_GammaNeuro((((150*(neuralBins - 1)) + 1)):(150*neuralBins))};
        
        LH_tempRawNeuroStruct(neuralBins, 1) = {LH_Neuro((((100000*(neuralBins - 1)) + 1)):(100000*neuralBins))};
        RH_tempRawNeuroStruct(neuralBins, 1) = {RH_Neuro((((100000*(neuralBins - 1)) + 1)):(100000*neuralBins))};
        
    end
end

ProcData.Sleep.Parameters.DeltaBand_Power.LH = LH_tempDeltaStruct;   % Place the data in the ProcData struct to later be saved
ProcData.Sleep.Parameters.DeltaBand_Power.RH = RH_tempDeltaStruct;

ProcData.Sleep.Parameters.ThetaBand_Power.LH = LH_tempThetaStruct;
ProcData.Sleep.Parameters.ThetaBand_Power.RH = RH_tempThetaStruct;

ProcData.Sleep.Parameters.GammaBand_Power.LH = LH_tempGammaStruct;
ProcData.Sleep.Parameters.GammaBand_Power.RH = RH_tempGammaStruct;

ProcData.Sleep.Parameters.RawNeuro.LH = LH_tempRawNeuroStruct;
ProcData.Sleep.Parameters.RawNeuro.RH = RH_tempRawNeuroStruct;

%% BLOCK PURPOSE: Create folder for the Whisker Acceleration
whiskerAngle = ProcData.Data.Behavior.whiskers;     % Unfiltered whisker angle

%     % Smooth the signal with a 10 Hz low pass 4th-order butterworth filter
%     [D, C] = butter(4, 10 / (ProcData.Notes.whiskerCamSamplingRate / 2), 'low');
%     filteredWhiskerAngle = filtfilt(D, C, whiskerAngle);    % Filtered whisker angle
whiskerAcceleration = diff(whiskerAngle, 2);

% Find the number of whiskerBins due to frame drops.
whiskerBinNumber = ceil(length(whiskerAcceleration) / 150);

% Divide the signal into five second bins and put them in a cell array
tempWhiskerStruct = cell(whiskerBinNumber, 1);   % Pre-allocate cell array

for whiskerBins = 1:whiskerBinNumber  % loop through all 9000 samples across 5 minutes in 5 second bins (60 total)
    if whiskerBins == 1
        tempWhiskerStruct(whiskerBins, 1) = {abs(whiskerAcceleration(whiskerBins:150))};  % Samples 1 to 150
    elseif whiskerBins == whiskerBinNumber
        tempWhiskerStruct(whiskerBins, 1) = {abs(whiskerAcceleration((((150*(whiskerBins - 1)) + 1)):end))};  % Samples 8701 to end. which changes due to dropped frames
    else
        tempWhiskerStruct(whiskerBins, 1) = {abs(whiskerAcceleration((((150*(whiskerBins - 1)) + 1)):(150*whiskerBins)))};  % Samples 151 to 300, etc...
    end
end

ProcData.Sleep.Parameters.WhiskerAcceleration = tempWhiskerStruct;   % Place the data in the ProcData struct to later be saved

%% BLOCK PURPOSE: Create folder for the Heart Rate
% Find the heart rate from the current ProcData file
HeartRate = ProcData.Data.HeartRate.HR;

% Divide the signal into five second bins and put them in a cell array
tempHRStruct = cell(60, 1);   % Pre-allocate cell array

for HRBins = 1:60  % loop through all 297 samples across 5 minutes in 5 second bins (60 total)
    if HRBins == 1
        tempHRStruct(HRBins, 1) = {HeartRate(HRBins:5)};  % Samples 1 to 5
    elseif HRBins == 60
        tempHRStruct(HRBins, 1) = {HeartRate((((5*(HRBins - 1)) + 1)):end)};  % Samples 297 to end.
    else
        tempHRStruct(HRBins, 1) = {HeartRate((((5*(HRBins - 1)) + 1)):(5*HRBins))};  % Samples 6 to 10, etc...
    end
end

ProcData.Sleep.Parameters.HeartRate = tempHRStruct;   % Place the data in the ProcData struct to later be saved

%% BLOCK PURPOSE: Create folder for the left and right CBV data
LH_CBV = ProcData.Data.CBV.LH;
RH_CBV = ProcData.Data.CBV.RH;
LH_ElectrodeCBV = ProcData.Data.CBV.LH_Electrode;
RH_ElectrodeCBV = ProcData.Data.CBV.RH_Electrode;

LH_NormCBV = (LH_CBV - RestingBaselines.CBV.LH.(strDay)) / RestingBaselines.CBV.LH.(strDay);
RH_NormCBV = (RH_CBV - RestingBaselines.CBV.RH.(strDay)) / RestingBaselines.CBV.RH.(strDay);
LH_NormElectrodeCBV = (RH_ElectrodeCBV - RestingBaselines.CBV.LH_Electrode.(strDay)) / RestingBaselines.CBV.LH_Electrode.(strDay);
RH_NormElectrodeCBV = (RH_ElectrodeCBV - RestingBaselines.CBV.RH_Electrode.(strDay)) / RestingBaselines.CBV.RH_Electrode.(strDay);

[D, C] = butter(4, 1 / (30 / 2), 'low');  
LH_FiltCBV = filtfilt(D, C, LH_NormCBV);
RH_FiltCBV = filtfilt(D, C, RH_NormCBV);
LH_ElectrodeFiltCBV = filtfilt(D, C, LH_NormElectrodeCBV);
RH_ElectrodeFiltCBV = filtfilt(D, C, RH_NormElectrodeCBV);

LH_tempCBVStruct = cell(60, 1);   % Pre-allocate cell array
RH_tempCBVStruct = cell(60, 1);   
LH_tempElectrodeCBVStruct = cell(60, 1);
RH_tempElectrodeCBVStruct = cell(60, 1);   

for CBVBins = 1:60   % loop through all 9000 samples across 5 minutes in 5 second bins (60 total)
    if CBVBins == 1
        LH_tempCBVStruct(CBVBins, 1) = {LH_FiltCBV(CBVBins:150)};  % Samples 1 to 150
        RH_tempCBVStruct(CBVBins, 1) = {RH_FiltCBV(CBVBins:150)};
        LH_tempElectrodeCBVStruct(CBVBins, 1) = {LH_ElectrodeFiltCBV(CBVBins:150)};
        RH_tempElectrodeCBVStruct(CBVBins, 1) = {RH_ElectrodeFiltCBV(CBVBins:150)};
    else
        LH_tempCBVStruct(CBVBins, 1) = {LH_FiltCBV((((150*(CBVBins - 1)) + 1)):(150*CBVBins))};  % Samples 151 to 300, etc...
        RH_tempCBVStruct(CBVBins, 1) = {RH_FiltCBV((((150*(CBVBins - 1)) + 1)):(150*CBVBins))};
        LH_tempElectrodeCBVStruct(CBVBins, 1) = {LH_ElectrodeFiltCBV((((150*(CBVBins - 1)) + 1)):(150*CBVBins))};  
        RH_tempElectrodeCBVStruct(CBVBins, 1) = {RH_ElectrodeFiltCBV((((150*(CBVBins - 1)) + 1)):(150*CBVBins))};
    end
end

ProcData.Sleep.Parameters.CBV.LH = LH_tempCBVStruct;   % Place the data in the ProcData struct to later be saved
ProcData.Sleep.Parameters.CBV.RH = RH_tempCBVStruct;
ProcData.Sleep.Parameters.CBV.LH_Electrode = LH_tempElectrodeCBVStruct;
ProcData.Sleep.Parameters.CBV.RH_Electrode = RH_tempElectrodeCBVStruct;

save([animal '_' hem '_' fileID '_', 'ProcData']);

end
