function [TDMSFile] = ReadInTDMSWhiskerTrials_FP(fileName)
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
[tempStruct,~] = ConvertTDMS_FP(0,fileName);
% Extract Whisker Camera info and transfer from TempStruct
TDMSFile.experimenter = tempStruct.Data.Root.Experimenter;
TDMSFile.animalID = tempStruct.Data.Root.Animal_ID;
TDMSFile.hemisphere = tempStruct.Data.Root.Hemisphere;
TDMSFile.solenoidPSI = tempStruct.Data.Root.Solenoid_PSI;
TDMSFile.isofluraneTime = tempStruct.Data.Root.Isoflurane_time;
TDMSFile.sessionID = tempStruct.Data.Root.Session_ID;
TDMSFile.amplifierGain = tempStruct.Data.Root.Amplifier_Gain;
TDMSFile.whiskCamSamplingRate = tempStruct.Data.Root.WhiskerCam_Fs;
TDMSFile.analogSamplingRate = tempStruct.Data.Root.Analog_Fs;
TDMSFile.trialDuration_sec = tempStruct.Data.Root.TrialDuration_sec;
TDMSFile.whiskCamPixelWidth = tempStruct.Data.Root.Whisker_Cam_Width_pix;
TDMSFile.whiskCamPixelHeight = tempStruct.Data.Root.Whisker_Cam_Height_pix;
% Pre-allocate - Data is contained in .vals folder in rows with corresponding labels in .names
TDMSFile.data.vals = NaN*ones(length(tempStruct.Data.MeasuredData),length(tempStruct.Data.MeasuredData(1).Data));
TDMSFile.data.names = cell(length(tempStruct.Data.MeasuredData),1);
% Pull data from tempStruct and allocate it in the proper areas 
for k = 1:length(tempStruct.Data.MeasuredData)
    TDMSFile.data.vals(k,:) = tempStruct.Data.MeasuredData(k).Data;
    TDMSFile.data.names{k} = strrep(tempStruct.Data.MeasuredData(k).Name,'Analog_Data','');
end

end

