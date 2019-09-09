function [TDMSFile] = ReadInTDMSWhiskerTrials_LD_IOS(fileName)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Pull the data and notes from the LabVIEW '.tdms' files into a Matlab structure.
%________________________________________________________________________________________________________________________
%
%   Inputs: File name ending in '.tdms' that contains the LabVIEW aquired analog data and notes from the session.
%
%   Outputs: Structure containing the data (arranged into rows with corresponding labels in a different field)
%            and various descriptive variables/strings of the session notes.
%
%   Last Revised: February 23rd, 2019
%________________________________________________________________________________________________________________________

%% Convert the .tdms file into something that Matlab understands
[TempStruct, ~] = ConvertTDMS_IOS(0, fileName);

% Data is contained in .Vals folder in rows with corresponding labels in .Names
TDMSFile.data.vals = NaN*ones(length(TempStruct.Data.MeasuredData), length(TempStruct.Data.MeasuredData(1).Data));
TDMSFile.data.names = cell(length(TempStruct.Data.MeasuredData), 1);

for k = 1:length(TempStruct.Data.MeasuredData)
    TDMSFile.data.vals(k,:) = TempStruct.Data.MeasuredData(k).Data;
    TDMSFile.data.names{k} = strrep(TempStruct.Data.MeasuredData(k).Name, 'Analog_Data', '');
end

end

