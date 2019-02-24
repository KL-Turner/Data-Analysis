function [days] = ConvertDate(dateTag)
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

days = cell(size(dateTag, 1), 1);
for f = 1:size(dateTag, 1)
    days{f} = datestr([2000 + str2double(dateTag(f, 1:2)) str2double(dateTag(f, 3:4)) str2double(dateTag(f, 5:6)) 00 00 00], 'mmmdd');
end

if length(days) == 1
    days = days{1};
end

end
