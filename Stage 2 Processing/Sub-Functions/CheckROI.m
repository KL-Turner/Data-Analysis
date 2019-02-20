function [isok] = CheckROI(ROIname)
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
%   DESCRIPTION: Checks whether the given ROI exists in shared
%   variables.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   ROIname - [string] a name designating the ROI
%
%                   animal - [string] identifier for the animal
%
%                   hem - [string] hemisphere recorded
%_______________________________________________________________
%   RETURN:                     
%                   isok - [int] designates whether ROI exists or not (1/0)        
%_______________________________________________________________

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
