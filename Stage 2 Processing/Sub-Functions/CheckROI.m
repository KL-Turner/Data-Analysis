function [isok] = CheckROI(ROIname)
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

ROIfile = ls('*ROIs.mat');
if not(isempty(ROIfile))
    load(ROIfile)
else
    ROIs = [];
end

if isfield(ROIs, ROIname)
    isok = 1;
else
    isok = 0;
end

end
