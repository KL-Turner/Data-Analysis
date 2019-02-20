function Combine_LabVIEW_MScan_Files(labviewDataFiles, mscanDataFiles)
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
    imageID = MScanData.Notes.imageID;
    
    CombData.Notes.LabVIEW = LabVIEWData.Notes;
    CombData.Data.Whisker_Angle = LabVIEWData.Data.WhiskerAngle;
    CombData.Data.Force_Sensor_L = LabVIEWData.Data.Force_Sensor;
    
    CombData.Notes.MScan = MScanData.Notes;
    CombData.Data.Neural_Data = MScanData.Data.MScan_Neural_Data';
    CombData.Data.Force_Sensor_M = MScanData.Data.MScan_Force_Sensor';
    CombData.Data.Vessel_Diameter = MScanData.Data.Vessel_Diameter;
    
    save([animalID '_' fileID '_' imageID '_CombData'], 'CombData')
end
 
end
