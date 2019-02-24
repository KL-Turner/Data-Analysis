function [refl] = BinToIntensity(fileName, ROImask, frames)
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

% Import camera frames from dalsa file
if nargin < 3
    [frames] = ReadDalsaBinary(fileName, 256, 256);
end

nFrames = length(frames);


refl = zeros(1, nFrames);
for n = 1:nFrames
    mask = ROImask.*double(frames{n});
    refl(n) = mean(nonzeros(mask));
end

end
