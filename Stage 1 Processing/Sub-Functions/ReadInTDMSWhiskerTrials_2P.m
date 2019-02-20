function [TDMSFile] = ReadInTDMSWhiskerTrials_2P(fileName)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%
%   Last Revised: January 18th, 2019
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
TDMSFile.imagedHemisphere = TempStruct.Data.Root.Hemisphere;
TDMSFile.isofluraneTime_Military = TempStruct.Data.Root.Isoflurane_time;
TDMSFile.sessionID = TempStruct.Data.Root.Session_ID;
TDMSFile.amplifierGain = TempStruct.Data.Root.Amplifier_Gain;
TDMSFile.whiskerCamSamplingRate = TempStruct.Data.Root.WhiskerCam_Fs;
TDMSFile.analogSamplingRate = TempStruct.Data.Root.Analog_Fs;
TDMSFile.trialDuration_Seconds = TempStruct.Data.Root.TrialDuration_sec;
TDMSFile.whiskerCamPixelHeight = TempStruct.Data.Root.Whisker_Cam_Height_pix;
TDMSFile.whiskerCamPixelWidth = TempStruct.Data.Root.Whisker_Cam_Width_pix;
TDMSFile.numberDroppedWhiskerCamFrames = TempStruct.Data.Root.WhiskerCam_NumberDropped;
TDMSFile.droppedWhiskerCamFrameIndex = TempStruct.Data.Root.WhiskerCam_DroppedFrameIndex;
        
TDMSFile.Data.Vals = NaN*ones(length(TempStruct.Data.MeasuredData), length(TempStruct.Data.MeasuredData(1).Data));
TDMSFile.Data.Names = cell(length(TempStruct.Data.MeasuredData), 1) ;

for k = 1:length(TempStruct.Data.MeasuredData)
    TDMSFile.Data.Vals(k,:) = TempStruct.Data.MeasuredData(k).Data;
    TDMSFile.Data.Names{k} = strrep(TempStruct.Data.MeasuredData(k).Name, 'Analog_Data', '');
end

end

