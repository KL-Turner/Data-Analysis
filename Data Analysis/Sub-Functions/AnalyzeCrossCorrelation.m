%___________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%___________________________________________________________________________________________________
%
%   Purpose: To analyze the CBV and correlation coefficient for bilateral hemispheres
%            in a single trial
%___________________________________________________________________________________________________
%
%   Inputs: 
%          
%
%   Outputs: 
%
%___________________________________________________________________________________________________

%% Close all open figures, clear the Workspace and command window of previous data
close all;
clc;
clear;

CBVEventData = ls('*_CBVEventData.mat');
load(CBVEventData);

GammaBandNeuralEventData = ls('*_GammaBandNeuralEventData.mat');
load(GammaBandNeuralEventData);

CBVRestData = ls('*_CBVRestData.mat');
load(CBVRestData);

GammaBandNeuralRestData = ls('*_GammaBandNeuralRestData.mat');
load(GammaBandNeuralRestData);

MultiUnitActivityNeuralRestData = ls('*_MultiUnitActivityNeuralRestData.mat');
load(MultiUnitActivityNeuralRestData);

%%
samplingRate = CBVEventData.RH.stim.samplingRate;     
timeVector = (0:(CBVEventData.RH.stim.epoch.duration*samplingRate)) / samplingRate - CBVEventData.RH.stim.epoch.offset; % 12 seconds

% Set up fields for FilterEvents
% RH_Criteria.Fieldname = {'solenoidName'};
% RH_Criteria.Comparison = {'equal'};
% RH_Criteria.Value = {'solenoidRightPad'};

RestCriteria.Fieldname = {'duration', 'duration'};
RestCriteria.Comparison = {'gt', 'lt'};
RestCriteria.Value = {15, 20};

% CBV stimulus data from right whisker puffs
% [RH_CBVFromRightStims] = FilterEvents(CBVEventData.LH.stim, RH_Criteria);  
% RH_CBVStimData = CBVEventData.LH.stim.NormData(RH_CBVFromRightStims, :);

% Gamma-band stimulus data from right whisker puffs
% [RH_GammaFromRightStims] = FilterEvents(GammaBandNeuralEventData.LH.stim, RH_Criteria);   
% RH_GammaStimData = GammaBandNeuralEventData.LH.stim.NormData(RH_GammaFromRightStims, :);

% Rest Data Logicals
[LH_RestCBVLogical] = FilterEvents(CBVRestData.LH, RestCriteria);
[LH_RestGammaLogical] = FilterEvents(GammaBandNeuralRestData.LH, RestCriteria);
[LH_RestMUALogical] = FilterEvents(MultiUnitActivityNeuralRestData.LH, RestCriteria);
[RH_RestCBVLogical] = FilterEvents(CBVRestData.RH, RestCriteria);
[RH_RestGammaLogical] = FilterEvents(GammaBandNeuralRestData.RH, RestCriteria);
[RH_RestMUALogical] = FilterEvents(MultiUnitActivityNeuralRestData.RH, RestCriteria);

LH_RestCBV = CBVRestData.LH.NormData(LH_RestCBVLogical, :);
LH_RestGamma = GammaBandNeuralRestData.LH.NormData(LH_RestGammaLogical, :);
LH_RestMUA = MultiUnitActivityNeuralRestData.LH.NormData(LH_RestMUALogical, :);
RH_RestCBV = CBVRestData.RH.NormData(RH_RestCBVLogical, :);
RH_RestGamma = GammaBandNeuralRestData.RH.NormData(RH_RestGammaLogical, :);
RH_RestMUA = MultiUnitActivityNeuralRestData.RH.NormData(RH_RestMUALogical, :);

%% Filter both signals below 2 Hz
[B, A] = butter(4, 2/(30/2), 'low');
% filteredRH_CBVStimData = filtfilt(B, A, RH_CBVStimData')';
% filteredRH_GammaStimData = filtfilt(B, A, RH_GammaStimData')';

for f = 1:length(LH_RestCBV)
    filteredLH_RestCBVData{f, :} = filtfilt(B, A, LH_RestCBV{f, :});
    filteredRH_RestCBVData{f, :} = filtfilt(B, A, RH_RestCBV{f, :});
    filteredLH_RestGammaData{f, :} = filtfilt(B, A, LH_RestGamma{f, :});
    filteredRH_RestGammaData{f, :} = filtfilt(B, A, RH_RestGamma{f, :});
    filteredLH_RestMUAData{f, :} = filtfilt(B, A, LH_RestMUA{f, :});
    filteredRH_RestMUAData{f, :} = filtfilt(B, A, RH_RestMUA{f, :});
end
 
for r = 1:length(filteredLH_RestCBVData)
    LH_restCBVData(r, :) = filteredLH_RestCBVData{r}(1:300);
    RH_restCBVData(r, :) = filteredRH_RestCBVData{r}(1:300); 
    LH_restGammaData(r, :) = filteredLH_RestGammaData{r}(1:300);
    RH_restGammaData(r, :) = filteredRH_RestGammaData{r}(1:300); 
    LH_restMUAData(r, :) = filteredLH_RestMUAData{r}(1:300);
    RH_restMUAData(r, :) = filteredRH_RestMUAData{r}(1:300); 
end

LH_restCBVMeans = mean(LH_restCBVData, 2);
RH_restCBVMeans = mean(RH_restCBVData, 2);
LH_restGammaMeans = mean(LH_restGammaData, 2);
RH_restGammaMeans = mean(RH_restGammaData, 2);
LH_restMUAMeans = mean(LH_restMUAData, 2);
RH_restMUAMean = mean(RH_restMUAData, 2);

for m = 1:length(LH_restCBVMeans)
    if LH_restCBVMeans(m, 1) >= 0
        msLH_restCBVData(m, :) = LH_restCBVData(m, :) - LH_restCBVMeans(m, 1);
    else
        msLH_restCBVData(m, :) = LH_restCBVData(m, :) + abs(LH_restCBVMeans(m, 1));
    end
    
    if RH_restCBVMeans(m, 1) >= 0
        msRH_restCBVData(m, :) = RH_restCBVData(m, :) - RH_restCBVMeans(m, 1); 
    else
        msRH_restCBVData(m, :) = RH_restCBVData(m, :) + abs(RH_restCBVMeans(m, 1));
    end
    
    if LH_restGammaMeans(m, 1) >= 0
        msLH_restGammaData(m, :) = LH_restGammaData(m, :) - LH_restGammaMeans(m, 1); 
    else
        msLH_restGammaData(m, :) = LH_restGammaData(m, :) + abs(LH_restGammaMeans(m, 1)); 
    end
    
    if RH_restGammaMeans(m, 1) >= 0
            msRH_restGammaData(m, :) = RH_restGammaData(m, :) - RH_restGammaMeans(m, 1); 
    else
            msRH_restGammaData(m, :) = RH_restGammaData(m, :) + abs(RH_restGammaMeans(m, 1)); 
    end
            
    if LH_restMUAMeans(m, 1) >= 0
            msLH_restMUAData(m, :) = LH_restMUAData(m, :) - LH_restMUAMeans(m, 1);
    else
            msLH_restMUAData(m, :) = LH_restMUAData(m, :) + abs(LH_restMUAMeans(m, 1));
    end
    
    if LH_restMUAMeans(m, 1) >= 0
            msRH_restMUAData(m, :) = RH_restMUAData(m, :) - RH_restMUAMean(m, 1);
    else
            msRH_restMUAData(m, :) = RH_restMUAData(m, :) + abs(RH_restMUAMean(m, 1));
    end
end

for xc = 1:size(msLH_restCBVData, 1)
    [LH_acor_CBV_Gam(xc, :), LH_lag_CBV_Gam(xc, :)] = xcorr(msLH_restCBVData(xc, :), msLH_restGammaData(xc, :), 150, 'coeff');
end

for xc = 1:size(msLH_restCBVData, 1)
    [LH_acor_CBV_MUA(xc, :), LH_lag_CBV_MUA(xc, :)] = xcorr(msLH_restCBVData(xc, :), msLH_restMUAData(xc, :), 150, 'coeff');
end

LH_acor_CBV_Gam_mean = mean(LH_acor_CBV_Gam, 1);
LH_lag_CBV_Gam_mean = mean(LH_lag_CBV_Gam, 1);

LH_acor_CBV_MUA_mean = mean(LH_acor_CBV_MUA, 1);
LH_lag_CBV_MUA_mean = mean(LH_lag_CBV_MUA, 1);

figure;
plot(LH_lag_CBV_Gam_mean, LH_acor_CBV_Gam_mean)
hold on;
plot(LH_lag_CBV_MUA_mean, LH_acor_CBV_MUA_mean)
