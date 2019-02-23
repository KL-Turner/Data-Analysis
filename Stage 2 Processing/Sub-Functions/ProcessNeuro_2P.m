function [neuro, neuroFs] = ProcessNeuro_2P(MScanData, NeurType, Neural_Hem)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%
%   Last Revised: August 8th, 2018
%________________________________________________________________________________________________________________________

%% Thresholds and Neurtype switch
neuralData = MScanData.Data.(Neural_Hem);
analogFs = MScanData.Notes.MScan_analogSamplingRate;

switch NeurType
    case 'MUApower'
        fpass = [300 3000];
    case 'Gam'
        fpass = [40 100];
    case 'Beta'
        fpass = [13 30];
    case 'Alpha'
        fpass = [8 12];
    case 'Theta'
        fpass = [4 8];
    case 'Delta'
        fpass = [1 4];
end

%% CALCULATE NEURAL POWER
if ismember(NeurType, [{'MUApower'}, {'Gam'}, {'Beta'}, {'Alpha'}, {'Theta'}, {'Delta'}])
    disp(['ProcessNeuro.m: Processing ' Neural_Hem ' ' NeurType]); disp(' ')
    neuroFs = 30;
    [z, p, k] = butter(4, fpass / (analogFs / 2));
    [sos, g] = zp2sos(z, p, k);
    filtNeuro = filtfilt(sos, g, neuralData - mean(neuralData));
    [z1, p1, k1] = butter(4, 10 / (analogFs / 2), 'low');
    [sos1, g1] = zp2sos(z1, p1, k1);
    Long_Neuro = filtfilt(sos1, g1, filtNeuro.^2);
    neuro = max(resample(Long_Neuro, neuroFs, analogFs), 0);
end

end