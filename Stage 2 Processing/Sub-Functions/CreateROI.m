function [mask] = CreateROI(img, ROIname, animal, hem)
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
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Prompts user to identify a region of interest on an image
%   and saves the selection to the shared variables.
%_______________________________________________________________
%   PARAMETERS:
%                   img - [matrix] the as output by ReadDalsaBinary.m
%
%                   ROIname - [string] a designation for the region of
%                   interest, usually should have some description and a
%                   date (i.g. 'Barrels_May20')
%
%                   animal - [string] ID for the animal
%
%                   hem - [string] hemisphere recorded
%_______________________________________________________________
%   RETURN:
%                   mask - [matrix] same size as image where pixels within
%                   the ROI are designated with a '1'. All other elements
%                   are '0'.
%_______________________________________________________________
ROIFile = ls('*ROIs.mat');
if not(isempty(ROIFile))
    load(ROIFile)
else
    ROIs = [];
end

disp(['Please select your region of interest for ' animal ' ' ROIname '.']); disp(' ')
figure;
imagesc(img);
colormap(gray);
axis image;
xlabel('Caudal');

if strcmp(hem,'LH')
    ylabel('Lateral');
elseif strcmp(hem,'RH')
    ylabel('Medial')
end

[mask, xi, yi] = roipoly;
ROIs.(ROIname).xi = xi;
ROIs.(ROIname).yi = yi;
save([animal '_ROIs.mat'], 'ROIs');
impoly(gca, [xi, yi]);

end
