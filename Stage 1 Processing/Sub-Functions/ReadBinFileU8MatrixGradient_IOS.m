function [imageGrad] = ReadBinFileU8MatrixGradient_IOS(fileName, height, width)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the pixel intensity of the whisker movie and output the the intensity vals as w x h x time.
%________________________________________________________________________________________________________________________
%
%   Inputs: File name ending in '.WhiskerCam.bin' that contains a movie of the whiskers. image height and width in pixels.
%
%   Outputs: imageout - [array, u8 integer] the measured intensity values organized as
%            [image width, image height, frame number].
%
%   Last Revised: February 23rd, 2019
%________________________________________________________________________________________________________________________
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   SUMMARY: Reads binary file with little endian ordering and 8-bit
%   grayscale formatting, calculates the gradient of the image to
%   emphasize the edges and outputs a cell array of images.
%_______________________________________________________________
%   PARAMETERS:
%               filename - [string] complete file name including extension
%
%               height - [integer] the height of the recorded image in
%               pixels
%
%               width - [integer]  the width of the recorded image in
%               pixels
%
%   RETURN:
%               imageout - [array, u8 integer] the measured intensity
%               values organized as [image width, image height, frame
%               number]
%_______________________________________________________________

% Calculate pixels per frame for fread
pixelsPerFrame = width*height;

% open the file , get file size , back to the begining
fid = fopen(fileName);
fseek(fid,0, 'eof');
fileSize = ftell(fid);
fseek(fid,0, 'bof');

% Identify the number of frames to read. Each frame has a previously
% defined width and height (as inputs), U8 has a depth of 1.
nFrameToRead = floor(fileSize / (pixelsPerFrame));
disp(['ReadBinFileU8MatrixGradient: ' num2str(nFrameToRead) ' frames to read.']); disp(' ')

% PreAllocate
imageGrad = int8(zeros(width, height, nFrameToRead));
for n = 1:nFrameToRead
    z = fread(fid, pixelsPerFrame, '*uint8', 0, 'l');
    indImg = reshape(z(1:pixelsPerFrame), width, height);
    imageGrad(:, :, n) = int8(gradient(double(indImg)));
end
fclose(fid);

end

