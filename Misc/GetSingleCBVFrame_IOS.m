function [Frame] = GetSingleCBVFrame(firstFile, imageWidth, imageHeight, animal, strDay)
%_______________________________________________________________________________________________
% Edited by Kevin L. Turner, Jr. 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%_______________________________________________________________________________________________%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Reads in a single frame of a binary image
%   
%_______________________________________________________________
%   PARAMETERS:      
%               filename - [string] binary file name with extension
%
%               image_width - [double] number of pixels in width of image
%
%               image_height - [double] number of pixels in height of image
%                               
%_______________________________________________________________
%   RETURN:                     
%               Frame - [array] single image of the binary file           
%_______________________________________________________________

% Calculate the number of pixels in a single frame
pixels_per_frame=imageWidth*imageHeight;

% Open the Binary File
fid = fopen(firstFile);

% Read the image from the binary file
FramePix = fread(fid, pixels_per_frame,'*int16','b');

% Reshape the image into rows and columns
img=reshape(FramePix,imageHeight,imageWidth);

% Orient the frame so that rostral is up
Frame = rot90(img',2);

fig = figure;
imagesc(Frame);           
colormap('gray')                    
axis image;
title([animal ' ' strDay ' Windows']);
set(gca, 'Ticklength', [0 0]);
colorbar
caxis([0 4096])

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Windows/'];

if ~exist(dirpath, 'dir') 
    mkdir(dirpath); 
end

savefig(fig, [dirpath animal ' ' strDay ' Windows']);

end
