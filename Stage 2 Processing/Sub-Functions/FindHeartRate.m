function [Sr, tr, fr, HR] = FindHeartRate(r, Fr)

% MultiTaperIOS: multi-taper analysis for IOS signal
% Written by Qingguang Zhang
% Adapted by K.L. Turner

% INPUTS:
%       r: dRR0 data for all pixels within each POI
%       Fr: frame rate
%
% OUTPUTS:
%       Sr: spectrum, a.u.
%       tr: time, in seconds
%       fr: frequency, in Hz
%       HR: heart rate, in Hz

% pre-process of IOS signal
% r = mean(r); % get average of all pixels within ROI
% r = r(:); % make sure is column vector
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
