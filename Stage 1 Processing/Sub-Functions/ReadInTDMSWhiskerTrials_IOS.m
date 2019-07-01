function [TDMSFile] = ReadInTDMSWhiskerTrials_IOS(fileName)
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

% Extract Whisker Camera info and transfer from TempStruct
TDMSFile.experimenter = TempStruct.Data.Root.Experimenter;
TDMSFile.animalID = TempStruct.Data.Root.Animal_ID;
TDMSFile.hemisphere = TempStruct.Data.Root.Hemisphere;
TDMSFile.solenoidPSI = TempStruct.Data.Root.Solenoid_PSI;
TDMSFile.isofluraneTime = TempStruct.Data.Root.Isoflurane_time;
TDMSFile.sessionID = TempStruct.Data.Root.Session_ID;
TDMSFile.amplifierGain = TempStruct.Data.Root.Amplifier_Gain;
TDMSFile.CBVCamSamplingRate = TempStruct.Data.Root.CBV_Cam_Fs;
TDMSFile.whiskCamSamplingRate = TempStruct.Data.Root.Whisk_Cam_Fs;
TDMSFile.webCamSamplingRate = TempStruct.Data.Root.Web_Cam_Fs;
TDMSFile.pupilCamSamplingRate = TempStruct.Data.Root.Pupil_Cam_Fs;
TDMSFile.analogSamplingRate = TempStruct.Data.Root.Analog_Fs;
TDMSFile.trialDuration_sec = TempStruct.Data.Root.TrialDuration_sec;
TDMSFile.CBVCamPixelWidth = TempStruct.Data.Root.CBVCam_Width_pix;
TDMSFile.CBVCamPixelHeight = TempStruct.Data.Root.CBVCam_Height_pix;
TDMSFile.CBVCamBitDepth = TempStruct.Data.Root.CBVCam_Bit_Depth;
TDMSFile.pupilCamPixelWidth = TempStruct.Data.Root.PupilCam_Width_pix;
TDMSFile.pupilCamPixelHeight = TempStruct.Data.Root.PupilCam_Height_pix;
TDMSFile.whiskCamPixelWidth = TempStruct.Data.Root.WhiskCam_Width_pix;
TDMSFile.whiskCamPixelHeight = TempStruct.Data.Root.WhiskCam_Height_pix;
TDMSFile.CBVCamExposureTime_microsec = TempStruct.Data.Root.CBVCam_Exposure_Time_microsec;
TDMSFile.CBVCamBinning = TempStruct.Data.Root.CBVCam_Binning;
TDMSFile.droppedPupilCamFrameIndex = TempStruct.Data.Root.PupilCam_DroppedFrameIndex;
TDMSFile.droppedWhiskCamFrameIndex = TempStruct.Data.Root.WhiskCam_DroppedFrameIndex;

% Data is contained in .Vals folder in rows with corresponding labels in .Names
TDMSFile.data.vals = NaN*ones(length(TempStruct.Data.MeasuredData), length(TempStruct.Data.MeasuredData(1).Data));
TDMSFile.data.names = cell(length(TempStruct.Data.MeasuredData), 1) ;

for k = 1:length(TempStruct.Data.MeasuredData)
    TDMSFile.data.vals(k,:) = TempStruct.Data.MeasuredData(k).Data;
    TDMSFile.data.names{k} = strrep(TempStruct.Data.MeasuredData(k).Name, 'Analog_Data', '');
end

end

