function [diamPerc, S, f] = DiamPercPower(rawDiameter, baseDiameter, samplingRate)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%
% Originally written by Patrick J. Drew and Yurong Gao
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised: February 19th, 2019    
%________________________________________________________________________________________________________________________

L_Diam = length(rawDiameter);
duration = L_Diam/samplingRate;
the_fs = 1/samplingRate;

% change raw diameter to diameter percentage with a median filter and a low pass filter
d_Amp = (rawDiameter/baseDiameter - 1)*100;

diamMedf = medfilt1(d_Amp, 5);   % 5-point median filter
[B, A] = butter(3, 3/samplingRate, 'low');   % change from 6 to 3
diamPerc = filtfilt(B, A, diamMedf);

% power spectrum multitaper
NW = floor(max(1, the_fs*duration/2));

params.Fs = samplingRate;
params.tapers = [NW 2*NW - 1];
params.err = [1 .01];
params.fpass = [0.05 3];
[S,f,~] = mtspectrumc(diamPerc - mean(diamPerc), params);

end
