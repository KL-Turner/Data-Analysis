function CombineLabVIEWMScanFiles(labviewDataFiles, mscanDataFiles)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised: February 20th, 2019    
%________________________________________________________________________________________________________________________

for f = 1:size(labviewDataFiles, 1)
    disp(['Combining the data from LabVIEW and MScan file(s) number ' num2str(f) ' of ' num2str(size(labviewDataFiles, 1)) '...']); disp(' ');
    labviewDataFile = labviewDataFiles(f, :);
    mscanDataFile = mscanDataFiles(f, :);
    load(labviewDataFile);
    load(mscanDataFile);
    
    [animalID, ~, ~, fileID] = GetFileInfo(labviewDataFile);
    vesselID = MScanData.Notes.vesselID;
    
    MergedData.Notes.LabVIEW = LabVIEWData.Notes;
    MergedData.Data.Whisker_Angle = LabVIEWData.Data.dsWhisker_Angle;
    MergedData.Data.binWhisker_Angle = LabVIEWData.Data.binWhisker_Angle;
    MergedData.Data.Force_Sensor_L = LabVIEWData.Data.dsForce_Sensor_L;
    MergedData.Data.binForce_Sensor_L = LabVIEWData.Data.binForce_Sensor_L;
    
    MergedData.Notes.MScan = MScanData.Notes;
    MergedData.Data.Raw_NeuralData = MScanData.Data.MScan_Neural_Data;
    MergedData.Data.MUA_Power = MScanData.Data.MUA_Power;
    MergedData.Data.GammaBand_Power = MScanData.Data.GammaBand_Power;
    MergedData.Data.BetaBand_Power = MScanData.Data.BetaBand_Power;
    MergedData.Data.AlphaBand_Power = MScanData.Data.AlphaBand_Power;
    MergedData.Data.ThetaBand_Power = MScanData.Data.ThetaBand_Power;
    MergedData.Data.DeltaBand_Power = MScanData.Data.DeltaBand_Power;
    MergedData.Data.Force_Sensor_M = MScanData.Data.dsForce_Sensor_M;
    MergedData.Data.binForce_Sensor_M = MScanData.Data.binForce_Sensor_M;
    MergedData.Data.Vessel_Diameter = MScanData.Data.Vessel_Diameter;
    
    save([animalID '_' fileID '_' vesselID '_MergedData'], 'MergedData')
end
 
end
