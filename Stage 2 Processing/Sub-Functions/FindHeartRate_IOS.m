function [Sr, tr, fr, HR] = FindHeartRate(r, Fr)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Qingguang Zhang
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs: 
%
%   Outputs:
%
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

r = r - mean(r); % mean subtract to remove slow drift

% Select taper parameters
tapers_r=[2 3]; % [time band width, number of tapers]
movingwin_r=[3.33,1];
params_r.Fs=Fr; % Frame rate
params_r.fpass=[5 15];
params_r.tapers=tapers_r;

[Sr, tr, fr] = mtspecgramc(r, movingwin_r, params_r);
% Sr: spectrum; tr: time; fr: frequency
[~, ridx] = max(Sr, [], 2); % largest elements along the frequency direction
HR = fr(ridx); % heart rate, in Hz
end
