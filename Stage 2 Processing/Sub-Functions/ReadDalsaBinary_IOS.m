function [image] = ReadDalsaBinary_IOS(file, imageHeight, imageWidth)
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

% imagebasics
pixelsPerFrame = imageWidth*imageHeight;
% open the file , get file size , back to the begining
fid = fopen(file);
fseek(fid, 0, 'eof');
fileSize = ftell(fid);
fseek(fid, 0, 'bof');

% identify the number of frames to read. Each frame has a previously
% defined width and height (as inputs), along with a grayscale "depth" of 2"

nFramesToRead = floor(fileSize / (2*pixelsPerFrame));
% preallocate memory
image = cell(1, nFramesToRead);
for n = 1:nFramesToRead
    z = fread(fid, pixelsPerFrame, '*int16', 'b');
    img = reshape(z(1:pixelsPerFrame), imageHeight, imageWidth);
    image{n} = rot90(img', 2);
end

fclose('all');

end
