function [Neuro, NeurFs] = ProcessNeuro(RawData, NeurType, Neural_Hem)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%
%   Last Revised: August 8th, 2018
%________________________________________________________________________________________________________________________
%
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Nov 2013
%   Version 1
%
%   SUMMARY: Processes neural data according to variable NeurType
%_______________________________________________________________
%   INPUTS:
%                           RawData - [struct] a structure containing
%                           neural data (RawData.Neuro) and sampling
%                           frequency (RawData.an_fs)
%
%                           NeurType - [string] the type of processing
%                           desired. This variable references parameters
%                           which are stored in Global variables and
%                           ensures that all neural analysis is done in a
%                           conserved manner.
%                               Acceptable inputs are:
%                                   -HiGam
%_______________________________________________________________
%   OUTPUTS:
%                           Neuro - [array] the processed neural data
%
%                           NeurFs - [integer] the sampling frequency of
%                           the new neural signal
%_______________________________________________________________
%   REQUIRED SCRIPTS:
%       [res_data1, res_data2] = MatchDataLengths(data1,data2)
%_______________________________________________________________
%   CALLED BY:
%                           Process RawDataFile.m
%_______________________________________________________________
%   FUTURE VERSIONS:
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%_______________________________________________________________

%% Thresholds and Neurtype switch
ExpectedLength = RawData.Notes.trialDuration_Seconds*RawData.Notes.analogSamplingRate;
Trimmed_Neuro = RawData.Data.(Neural_Hem)(1:min(ExpectedLength, length(RawData.Data.(Neural_Hem))));

switch NeurType
    case 'MUpower'
        fpass = [300 3000];
    case 'Gam'
        fpass = [40 100];
    case 'HiGam'
        fpass = [60 100];
    case 'LoGam'
        fpass = [40 60];
    case 'Theta'
        fpass = [4 8];
    case 'Delta'
        fpass = [1 4];
end

%% CALCULATE NEURAL POWER
if ismember(NeurType, [{'MUpower'}, {'Gam'}, {'HiGam'}, {'LoGam'}, {'Theta'}, {'Delta'}])
    disp(['ProcessNeuro.m: Processing ' Neural_Hem ' ' NeurType]); disp(' ')
    Fs = RawData.Notes.analogSamplingRate;
    [z, p, k] = butter(4, fpass / (Fs / 2));
    [sos, g] = zp2sos(z, p, k);
    filtNeuro = filtfilt(sos, g, Trimmed_Neuro - mean(Trimmed_Neuro));
    [z1, p1, k1] = butter(4, 10 / (Fs / 2), 'low');
    [sos1, g1] = zp2sos(z1, p1, k1);
    Long_Neuro = filtfilt(sos1, g1, filtNeuro.^2);
    Neuro = max(resample(Long_Neuro, RawData.Notes.CBVCamSamplingRate, RawData.Notes.analogSamplingRate), 0);
    NeurFs = RawData.Notes.CBVCamSamplingRate;
    
% FILTER WIDEBAND LFP    
elseif strcmp(NeurType, 'Wideband_LFP')
    display(['ProcessNeuro.m: Processing ' Neural_Hem ' ' NeurType]); disp(' ')
    Fs = RawData.Notes.analogSamplingRate;
    [z,p,k] = butter(4,300/(Fs/2),'low');
    [sos,g] = zp2sos(z,p,k);
    filtNeuro = filtfilt(sos, g, Trimmed_Neuro-mean(Trimmed_Neuro));
    NeurFs = 1500;   % 300 Hz * 5 will give adequate oversampling
    Neuro = resample(filtNeuro, NeurFs, RawData.Notes.analogSamplingRate);
end

end
