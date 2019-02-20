function [Frame] = UniqueDayROIs(firstFile, imageWidth, imageHeight, animal, strDay)
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

ROIFile = ls('*_ROIs.mat');
load(ROIFile);

% Calculate the number of pixels in a single frame
pixels_per_frame=imageWidth*imageHeight;

% Open the Binary File
fid=fopen(firstFile);

% Read the image from the binary file
FramePix=fread(fid, pixels_per_frame,'*int16','b');

% Reshape the image into rows and columns
img=reshape(FramePix,imageHeight,imageWidth);

% Orient the frame so that rostral is up
Frame = rot90(img',2);

ROIFig = figure; 
imagesc(Frame);  
colormap('gray')  
hold on;
plot(ROIs.(['LH_' strDay]).xi, ROIs.(['LH_' strDay]).yi);
plot(ROIs.(['RH_' strDay]).xi, ROIs.(['RH_' strDay]).yi);
plot(ROIs.(['LH_Electrode_' strDay]).xi, ROIs.(['LH_Electrode_' strDay]).yi);
plot(ROIs.(['RH_Electrode_' strDay]).xi, ROIs.(['RH_Electrode_' strDay]).yi);
title([animal ' ' strDay ' Windows with ROIs'])
xlabel('Caudal')
ylabel('Left')
axis image
set(gca, 'Ticklength', [0 0])

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Windows/'];

if ~exist(dirpath, 'dir') 
    mkdir(dirpath); 
end

savefig(ROIFig, [dirpath animal ' ' strDay ' Windows with ROIs'])

end