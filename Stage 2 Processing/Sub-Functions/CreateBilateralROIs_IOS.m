function [ROIs] = CreateBilateralROIs_IOS(img, ROIname, animalID, ROIs)
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
%   Last Revised: June 27th, 2019
%________________________________________________________________________________________________________________________

roiFig = figure;
imagesc(img)
colormap(gray)
axis image
xlabel('Caudal')

disp(['Please select your region of interest for ' animalID ' ' ROIname '.']); disp(' ')
[~, xi, yi] = roipoly;
ROIs.(ROIname).xi = xi;
ROIs.(ROIname).yi = yi;
close(roiFig)

end
