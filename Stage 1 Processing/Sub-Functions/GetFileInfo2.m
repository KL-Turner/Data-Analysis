function [animalID, fileDate, fileID, vesselID] = GetFileInfo2(fileName)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%
%   Last Revised: February 20th, 2019
%________________________________________________________________________________________________________________________
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Uses file name in a standard format to get the animal ID
%   and hemisphere recorded.
%
%_______________________________________________________________
%   PARAMETERS:
%                  filename - [string] filename in a standard format.
%                       'SubjectID_Hemisphere_Date_HH_mm_ssdd'
%_______________________________________________________________
%   RETURN:
%                   Animal_ID - [string] the subject identifier
%
%                   hem - [string] the hemisphere recorded
%
%                   filedate - [string] the date the file was recorded
%
%                   fileID - [string] the timestamp of the file
%_______________________________________________________________

% Identify the extension
extInd = strfind(fileName(1, :), '.');
extension = fileName(1, extInd + 1:end);

% Identify the underscores
fileBreaks = strfind(fileName(1, :), '_');

switch extension
    case 'bin'
        animalID = [];
        vesselID = [];
        fileDate = fileName(:, 1:fileBreaks(1) - 1);
        fileID = fileName(:, 1:fileBreaks(4) - 1);
    case 'mat'
        % Use the known format to parse
        animalID = fileName(:, 1:fileBreaks(1) - 1);
        fileDate = fileName(:, fileBreaks(1) + 1:fileBreaks(2) - 1);
        if numel(fileBreaks) > 3
            fileID = fileName(:, fileBreaks(1) + 1:fileBreaks(5) - 1);
            vesselID = fileName(:, fileBreaks(5) + 1:fileBreaks(6) - 1);
        else
            fileDate = [];
            fileID = [];
        end
end

end
