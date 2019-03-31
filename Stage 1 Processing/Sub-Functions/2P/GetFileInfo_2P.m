function [animalID, hem, fileDate, fileID] = GetFileInfo_2P(fileName)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Identify important aspects of a file name and output each individually.
%________________________________________________________________________________________________________________________
%
%   Inputs: Any filename in the format AnimalID_Hemisphere_YYMMDD_HH_MM_SS independent of extension.
%
%   Outputs: animalID - typically Letter(s) followed by number(s)
%            hemisphere - typically 'LH' (left hem), RH, (right), or 'Both' (bilateral imaging)
%            fileDate - date in the form Year, Month, Day (YYMMDD) typically 6 numbers
%            fileID - date followed by the underscores of hour (military time) minutes and seconds. 'YYMMDD_HH_MM_SS'
%
%   Last Revised: March 21st, 2019
%________________________________________________________________________________________________________________________

% Identify the extension
extInd = strfind(fileName(1, :), '.');
extension = fileName(1, extInd + 1:end);

% Identify the underscores
fileBreaks = strfind(fileName(1, :), '_');

switch extension
    case 'bin'
        animalID = [];
        hem = [];
        fileDate = fileName(:, 1:fileBreaks(1) - 1);
        fileID = fileName(:, 1:fileBreaks(4) - 1);
    case 'mat'
        % Use the known format to parse
        animalID = fileName(:, 1:fileBreaks(1) - 1);
        hem = fileName(:, fileBreaks(1) + 1:fileBreaks(2) - 1);
        if numel(fileBreaks) > 3
            fileDate = fileName(:, fileBreaks(2) + 1:fileBreaks(3) - 1);
            fileID = fileName(:, fileBreaks(2) + 1:fileBreaks(6) - 1);
        else
            fileDate = [];
            fileID = [];
        end
end

end
