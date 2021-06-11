function Analyze2PDiameter_2P(mscanDataFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Run the sub-functions necessary to analyze the TIFF stack's changes in vessel diameter.
%________________________________________________________________________________________________________________________

% Analyze every image (surface/penatrating/capillary)
for a = 1:size(mscanDataFiles,1)
    load(mscanDataFiles(a,:),'-mat');
    if ~isfield(MScanData,'notes') && isfield(MScanData, 'Header') %Because I'm too tired to fix this properly
        MScanData.notes=MScanData.Header.notes;
    end
    if MScanData.notes.checklist.analyzeDiam == false
        MScanData_backup = MScanData;
        if strcmp(MScanData.notes.movieType,'MS') == true || strcmp(MScanData.notes.movieType,'MD') == true
            [MScanData] = ExtractTiffAnalogData_2P(MScanData,[MScanData.notes.date '_' MScanData.notes.imageID]);
            MScanData.notes.checklist.analyzeDiam = true;
            save([MScanData.notes.animalID '_' MScanData_backup.notes.date '_' MScanData.notes.imageID '_' MScanData.notes.vesselID '_MScanData'],'MScanData')        
        elseif strcmp(MScanData.notes.movieType,'MP') == true
            [MScanData] = CalcPenVesselArea_2P(MScanData,[MScanData.notes.date '_' MScanData.notes.imageID]);
            MScanData.notes.checklist.analyzeDiam = true;         
            save([MScanData.notes.animalID '_' MScanData_backup.notes.date '_' MScanData.notes.imageID '_' MScanData.notes.vesselID '_MScanData'],'MScanData')        
        elseif strcmp(MScanData.notes.movieType,'MC') == true
            
            [MScanData] = GetVelocityLineScan_2P(MScanData,[MScanData.notes.date '_' MScanData.notes.imageID]);
            MScanData.Header.notes.checklist.analyzeDiam = true;        
            save([MScanData.Header.notes.animalID '_' MScanData_backup.notes.date '_' MScanData.Header.notes.imageID '_' MScanData.Header.notes.vesselID '_MScanData'],'MScanData')        
        end
    else
        disp([mscanDataFiles(a,:) ' vessel diameter already analyzed. Continuing...']); disp(' ')
    end
end

end
