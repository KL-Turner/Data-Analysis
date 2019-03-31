function [mask] = GetROI(img, ROIname, animal, hem)
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

isok = 'n';
ROIfile = ls('*ROIs.mat');
if not(isempty(ROIfile))
    load(ROIfile)
else
    ROIs = [];
end

if isfield(ROIs,ROIname)
    xi = ROIs.(ROIname).xi;
    yi = ROIs.(ROIname).yi;
    disp([animal ' ROI for ' ROIname ' found, saving to RawData file.']); disp(' ')

    if strcmp(isok,'y')
        mask = roipoly(img,xi,yi);
    end
    while strcmp(isok,'n') == 1
        isok = 'y';
        if strcmp(isok,'y')
            mask = roipoly(img, xi, yi);
            continue;
        elseif strcmp(isok,'n')
            [mask] = CreateROI(img,animal,hem);
        end
    end
else
    [mask] = CreateROI(img, ROIname, animal, hem);
end

end
