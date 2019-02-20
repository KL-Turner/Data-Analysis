function [days] = ConvertDate(dateTag)
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
%   DESCRIPTION: Converts the date format output by the LabVIEW acquisition
%   program into a date string.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   DateTag - [string] date vector output by the LabVIEW
%                       acquisition
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

days = cell(size(dateTag, 1), 1);
for f = 1:size(dateTag, 1)
    days{f} = datestr([2000 + str2double(dateTag(f, 1:2)) str2double(dateTag(f, 3:4)) str2double(dateTag(f, 5:6)) 00 00 00], 'mmmdd');
end

if length(days) == 1
    days = days{1};
end

end
