function [angle] = WhiskerTrackerParallel(fileName)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%
%   Last Revised: August 4th, 2018
%________________________________________________________________________________________________________________________
%
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University
%
%   SUMMARY: Tracks the approximate movement of whiskers using the radon
%           transform of an image of the whiskers. Identifies the whisker
%           angles by analyzing the column variance of the transformed
%           image. Columns with the highest variance will correspond to the
%           correct angle where the highest value of the radon transform
%           will correspond to the correct angle and y-position in the
%           image, while the lowest value will correspond to the correct
%           angle but an incorrect image position.
%
%           This version of the code uses an onboard GPU to speed up the
%           calculation of the whisker angles.
%_______________________________________________________________
%   PARAMETER TYPE:             
%                               filename - [string] list of filames
%                               without the extension.
%_______________________________________________________________
%   RETURN:                     
%                               TDMSFile - [struct] contains measured
%                               analog data and trial notes from the
%                               LabVIEW acquisition program
%_______________________________________________________________

% Variable Setup
theta = -40:80; % Angles used for radon

% Import whisker movie
importStart = tic;
baslerFrames = ReadBinFileU8MatrixGradient([fileName '_WhiskerCam.bin'], 350, 30);
importTime = toc(importStart);
disp(['WhiskerTrackerParallel: Binary file import time was ' num2str(importTime) ' seconds.']); disp(' ')

% Transfer the images to the GPU
gpuTrans1 = tic;
gpuFrame = gpuArray(baslerFrames);
gpuTransfer = toc(gpuTrans1);
disp(['WhiskerTrackerParallel: GPU transfer time was ' num2str(gpuTransfer) ' seconds.']); disp(' ')

% PreAllocate array of whisker angles, use NaN as a place holder
angle = NaN*ones(1, length(baslerFrames));
radonTime1 = tic;
for f = 1:(length(baslerFrames) - 1)
    % Radon on individual frame
    [R, ~] = radon(gpuFrame(:, :, f), theta);
    % Get transformed image from GPU and calculate the variance
    colVar = var(gather(R));
    % Sort the columns according to variance
    ordVar = sort(colVar);
    % Choose the top 0.1*number columns which show the highest variance
    thresh = round(numel(ordVar)*0.9);
    sieve = gt(colVar, ordVar(thresh));
    % Associate the columns with the corresponding whisker angle
    angles = nonzeros(theta.*sieve);
    % Calculate the average of the whisker angles
    angle(f) = mean(angles);
end
radonTime = toc(radonTime1);
disp(['WhiskerTrackerParallel: Whisker Tracking time was ' num2str(radonTime) ' seconds.']); disp(' ')

inds = isnan(angle)==1;
angle(inds) = [];

end
