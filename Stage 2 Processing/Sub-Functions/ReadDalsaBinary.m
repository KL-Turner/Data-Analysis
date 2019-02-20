function [image] = ReadDalsaBinary(file, imageHeight, imageWidth)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%
%   Last Revised: August 8th, 2018
%________________________________________________________________________________________________________________________
%
%   Adapted from ReadDalsaBinary2 by Aaron Winder, Drew Lab, ESM, 
%   Penn State University, Nov 2013
%   Version 1
%
%   SUMMARY: Converts a binary file containing images from a Dalsa 1M60
%   pantera camera into an image with proper orientation.
%_______________________________________________________________
%   INPUTS:
%                       thefile - the name of the binary file including the
%                       extension
%
%                       image_height - the height of the image in pixels
%
%                       image_width - the width of the image in pixels
%_______________________________________________________________
%   OUTPUTS:
%                       theimage - a cell array containing a stack of
%                       images
%_______________________________________________________________
%   REQUIRED SCRIPTS:
%_______________________________________________________________
%   CALLED BY:
%_______________________________________________________________
%   FUTURE VERSIONS:
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%_______________________________________________________________

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
