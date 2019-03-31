function [binForceSensor] = BinarizeForceSensor_2P(forceSensor, thresh)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Binarize the Force sensor with a given threshold.
%________________________________________________________________________________________________________________________
%
%   Inputs: Force sensor data and threshold value.
%
%   Outputs: Binarized value for each array point.
%
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

y = hilbert(diff(forceSensor));
env = abs(y);
binForceSensor = gt(env, thresh);

end
