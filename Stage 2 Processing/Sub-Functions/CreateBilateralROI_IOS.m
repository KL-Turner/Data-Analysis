function [mask] = CreateBilateralROIs_IOS(img, ROIname, animal)
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
