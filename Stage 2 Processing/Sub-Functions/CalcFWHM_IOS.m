function width = CalcFWHM_IOS(data,smoothing,threshold)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: Calculates the full width at half max to determine changes in vessel diameter
%________________________________________________________________________________________________________________________
%
%   Inputs: data, smoothing and threshold conditions
%
%   Outputs: width
%
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

% function which takes data and calculates the full-width, half max value
% half-max values are found looking in from the sides, i.e., the program will work
% even if the data dips to a lower value in the middle
% 2009-06-17 - no longer subracting the min from the entire data set
% 2009-06-17 - changed convolution to conv2, and added 'valid' parameter

data = double(data(:));     % make sure this is column, and cast to double

% smooth data, if appropriate
if nargin < 2
    % smoothing not passed in, set to default (none)
    smoothing = 1;
end

if smoothing > 1
    data = conv2(data,rectwin(smoothing) ./ smoothing,'valid');
end

% subtract out baseline
%data = data - min(data);   %jd - don't subract out min, in case threshold was set externally

if nargin < 3
    %threshold = max(data)/2;
    offset = min(data);                           % find the baseline 
    threshold = max(data - offset) / 2 + offset;  % threshold is half max, taking offset into account    
end

aboveI = find(data > threshold);    % all the indices where the data is above half max

if isempty(aboveI)
    % nothing was above threshold!
    width = 0;
    return
end

firstI = aboveI(1);                 % index of the first point above threshold
lastI = aboveI(end);                % index of the last point above threshold

if (firstI-1 < 1) | (lastI+1) > length(data)
    % interpolation would result in error, set width to zero and just return ...
    width = 0;
    return
end

% use linear interpolation to get a more accurate picture of where the max was
% find value difference between the point and the threshold value,
% and scale this by the difference between integer points ...
point1offset = (threshold-data(firstI-1)) / (data(firstI)-data(firstI-1));
point2offset = (threshold-data(lastI)) / (data(lastI+1)-data(lastI));

point1 = firstI-1 + point1offset;
point2 = lastI + point2offset;

width = point2-point1;
