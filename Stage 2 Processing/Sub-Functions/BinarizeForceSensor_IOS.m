function [binPSWF] = BinarizeForceSensor(PSWF, thresh)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
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

y = hilbert(diff(PSWF));
env = abs(y);
binPSWF = gt(env, thresh);

end
