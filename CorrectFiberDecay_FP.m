function [lowPassData,FiberData] = CorrectFiberDecay_FP(filtData,channelData,sosLow,gLow,FiberData)

%% READ ME
% This subfunction is used to remove slow trends from fiber photometry
% Slow decay is fitted using a 2 exponent regression on 0.05Hz low pass
% filtered data (filtData). This fit is then subtracted from the
% raw data (channelData) prior to low-pass filtering @1Hz.
%% INPUTS
%filtData: optical signals low-pass filtered below 0.05Hz to regression
%of metabolic decay/photobleaching
%
%channelData: raw optical channel data to have baseline drift
%subtracted from
%
%sosLow/gLow: filter parameters for low-pass filtering detrended
%signals
%
%FiberData: large data structure for loading fitting parameters and
%results for storage and analysis.
%
%trialNum: which trial number within imaging session was fit on
%function call
%
%sessionNum: which imaging session is being processed

%% Remove photobleaching/metabolism baseline drift RH
spacing = 1:1:length(filtData(:,3));
% RH correct TRITC blood volume
[fitVals] = fit(spacing',filtData(:,3),'exp2');
coeffVals = coeffvalues(fitVals);
RH_predictedCBV = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
RH_correctedCBV = channelData(:,3) - RH_predictedCBV';
FiberData.FitStruct.RH_CBV = fitVals;
% RH correct Ca2+ dependent GCaMP
[fitVals] = fit(spacing',filtData(:,2),'exp2');
coeffVals = coeffvalues(fitVals);
RH_predictedGCaMP = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
RH_correctedGCaMP = channelData(:,2) - RH_predictedGCaMP';
FiberData.FitStruct.RH_GCaMP = fitVals;
% RH correct Ca2+ independent GCaMP
[fitVals] = fit(spacing',filtData(:,1),'exp2');
coeffVals = coeffvalues(fitVals);
RH_predictedGFP = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
RH_correctedGFP = channelData(:,1) - RH_predictedGFP';
FiberData.FitStruct.RH_GFP = fitVals;
%% Remove photobleaching/metaloism baseline drift RH
% LH correct TRITC blood volume
[fitVals] = fit(spacing',filtData(:,6),'exp2');
coeffVals = coeffvalues(fitVals);
LH_predictedCBV = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
LH_correctedCBV = channelData(:,6) - LH_predictedCBV';
FiberData.FitStruct.LH_CBV = fitVals;
% LH correct Ca2+ dependent GCaMP
[fitVals] = fit(spacing',filtData(:,5),'exp2');
coeffVals = coeffvalues(fitVals);
LH_predictedGCaMP = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
LH_correctedGCaMP = channelData(:,5) - LH_predictedGCaMP';
FiberData.FitStruct.LH_GCaMP = fitVals;
% LH correct Ca2+ independent GCaMP
[fitVals] = fit(spacing',filtData(:,4),'exp2');
coeffVals = coeffvalues(fitVals);
LH_predictedGFP = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
LH_correctedGFP = channelData(:,4) - LH_predictedGFP';
FiberData.FitStruct.LH_GFP = fitVals;
%% detrend, low-pass filter, rescale data
detrendData = [RH_correctedGFP,RH_correctedGCaMP,RH_correctedCBV,LH_correctedGFP,LH_correctedGCaMP,LH_correctedCBV]; % This is the optical data with photobleaching and metabolism dependent decreases removed- KWG
lowPassData = filtfilt(sosLow,gLow,detrendData); % This is the data low-pass filtered below 1Hz for each of the sessions- KWG

end
