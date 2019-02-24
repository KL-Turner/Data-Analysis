function [TDMSFile] = ReadInTDMSWhiskerTrials(fileName)
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
[TempStruct, ~] = ConvertTDMS(0, fileName);

% Extract Whisker Camera info and transfer from TempStruct
TDMSFile.experimenter = TempStruct.Data.Root.Experimenter;
TDMSFile.animalID = TempStruct.Data.Root.Animal_ID;
TDMSFile.imagedHemisphere = TempStruct.Data.Root.Imaged_Hemisphere;
TDMSFile.solenoidPressure_PSI = TempStruct.Data.Root.Solenoid_Pressure;
TDMSFile.isofluraneTime_Military = TempStruct.Data.Root.Isoflurane_time;
TDMSFile.sessionID = TempStruct.Data.Root.Session_ID;
TDMSFile.amplifierGain = TempStruct.Data.Root.Amplifier_Gain;
TDMSFile.CBVCamSamplingRate = TempStruct.Data.Root.CBVCam_Rate;
TDMSFile.whiskerCamSamplingRate = TempStruct.Data.Root.WhiskerCam_Rate;
TDMSFile.webCamSamplingRate = TempStruct.Data.Root.WebCam_Rate;
TDMSFile.analogSamplingRate = TempStruct.Data.Root.Analog_Sampling_Rate;
TDMSFile.trialDuration_Seconds = TempStruct.Data.Root.TrialDuration_sec;
TDMSFile.CBVCamPixelHeight = TempStruct.Data.Root.CBV_Cam_Height_pix;
TDMSFile.CBVCamPixelWidth = TempStruct.Data.Root.CBV_Cam_Width_pix;
TDMSFile.CBVCamBitDepth = TempStruct.Data.Root.CBV_Cam_Bit_Depth;
TDMSFile.whiskerCamPixelHeight = TempStruct.Data.Root.Whisker_Cam_Height_pix;
TDMSFile.whiskerCamPixelWidth = TempStruct.Data.Root.Whisker_Cam_Width_pix;
TDMSFile.CBVCamExposureTime_Microseconds = TempStruct.Data.Root.CBVCam_Exposure_Time_microsec;
TDMSFile.CBVCamBinning = TempStruct.Data.Root.CBVCam_Binning;
TDMSFile.numberDroppedWhiskerCamFrames = TempStruct.Data.Root.WhiskerCam_NumberDropped;
TDMSFile.droppedWhiskerCamFrameIndex = TempStruct.Data.Root.WhiskerCam_DroppedFrameIndex;

% Data is contained in .Vals folder in rows with corresponding labels in .Names
TDMSFile.Data.Vals = NaN*ones(length(TempStruct.Data.MeasuredData), length(TempStruct.Data.MeasuredData(1).Data));
TDMSFile.Data.Names = cell(length(TempStruct.Data.MeasuredData), 1) ;

for k = 1:length(TempStruct.Data.MeasuredData)
    TDMSFile.Data.Vals(k,:) = TempStruct.Data.MeasuredData(k).Data;
    TDMSFile.Data.Names{k} = strrep(TempStruct.Data.MeasuredData(k).Name, 'Analog_Data', '');
end

end

