function CombineLabVIEWMScanFiles(labviewDataFiles, mscanDataFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs: 
%
%   Last Revised: February 29th, 2019
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
    MergedData.Data.Whisker_Angle = LabVIEWData.Data.dsWhisker_Angle2;
    MergedData.Data.binWhisker_Angle = LabVIEWData.Data.binWhisker_Angle2;
    MergedData.Data.Force_Sensor_L = LabVIEWData.Data.dsForce_Sensor_L2;
    MergedData.Data.binForce_Sensor_L = LabVIEWData.Data.binForce_Sensor_L2;
    
    MergedData.Notes.MScan = MScanData.Notes;
    MergedData.Data.Raw_NeuralData = MScanData.Data.MScan_Neural_Data2;
    MergedData.Data.MUA_Power = MScanData.Data.MUA_Power2;
    MergedData.Data.GammaBand_Power = MScanData.Data.GammaBand_Power2;
    MergedData.Data.BetaBand_Power = MScanData.Data.BetaBand_Power2;
    MergedData.Data.AlphaBand_Power = MScanData.Data.AlphaBand_Power2;
    MergedData.Data.ThetaBand_Power = MScanData.Data.ThetaBand_Power2;
    MergedData.Data.DeltaBand_Power = MScanData.Data.DeltaBand_Power2;
    MergedData.Data.Force_Sensor_M = MScanData.Data.dsForce_Sensor_M2;
    MergedData.Data.binForce_Sensor_M = MScanData.Data.binForce_Sensor_M2;
    MergedData.Data.Vessel_Diameter = MScanData.Data.Vessel_Diameter2;
    MergedData.Data.EMG = MScanData.Data.filtEMG2;
    
    save([animalID '_' vesselID '_' fileID '_MergedData'], 'MergedData')
end

end
