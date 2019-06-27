function  [bin_wwf] = BinarizeWhiskers(angl, fs, thresh1, thresh2)
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

% Differentiate, rectify, subtract off noise, rectify
dd_wwf = abs((diff(angl, 2)))*fs^2;
bin_wwf1 = gt(dd_wwf, thresh1);   % Acceleration exceeds lower threshold
bin_wwf2 = gt(dd_wwf, thresh2);   % Acceleration exceeds upper threshold

% Combine the two waveforms
bin_wwf = (bin_wwf1 + bin_wwf2) / 2;

end