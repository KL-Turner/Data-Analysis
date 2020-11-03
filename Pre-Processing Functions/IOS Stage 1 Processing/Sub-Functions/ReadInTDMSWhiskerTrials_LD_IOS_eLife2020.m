function [TDMSFile] = ReadInTDMSWhiskerTrials_LD_IOS_eLife2020(fileName)
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

% Convert the .tdms file into something that Matlab understands
[tempStruct,~] = ConvertTDMS_IOS_eLife2020(0,fileName);
% Pre-allocate - Data is contained in .vals folder in rows with corresponding labels in .names
TDMSFile.data.vals = NaN*ones(length(tempStruct.Data.MeasuredData),length(tempStruct.Data.MeasuredData(1).Data));
TDMSFile.data.names = cell(length(tempStruct.Data.MeasuredData),1);
% Pull data from tempStruct and allocate it in the proper areas 
for k = 1:length(tempStruct.Data.MeasuredData)
    TDMSFile.data.vals(k,:) = tempStruct.Data.MeasuredData(k).Data;
    TDMSFile.data.names{k} = strrep(tempStruct.Data.MeasuredData(k).Name,'Analog_Data','');
end

end
