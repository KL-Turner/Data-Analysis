function [TDMSFile] = ReadInTDMSWhiskerTrials(fileName)
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
%   Purpose: Reads in .TDMS files from a LABVIEW acquisition program using the convertTDMS.m 
%            script (acquired through Mathworks file exchange)and organizes the acquired data 
%            into a structure.
%
%            .bin - Cameras
%            .tdms - Digital and Analog Data
%
%            http://www.mathworks.com/matlabcentral/fileexchange/44206-converttdms--v10-
%___________________________________________________________________________________________________
%
%   Inputs: filenames - 
%
%   Outputs: TDMSFile - [struct] contains measured analog data and trial notes from the
%            LabVIEW acquisition program

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
        
TDMSFile.Data.Vals = NaN*ones(length(TempStruct.Data.MeasuredData), length(TempStruct.Data.MeasuredData(1).Data));
TDMSFile.Data.Names = cell(length(TempStruct.Data.MeasuredData), 1) ;

for k = 1:length(TempStruct.Data.MeasuredData)
    TDMSFile.Data.Vals(k,:) = TempStruct.Data.MeasuredData(k).Data;
    TDMSFile.Data.Names{k} = strrep(TempStruct.Data.MeasuredData(k).Name, 'Analog_Data', '');
end

end

