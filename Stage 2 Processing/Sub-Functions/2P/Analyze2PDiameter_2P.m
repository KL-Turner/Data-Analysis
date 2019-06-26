function Analyze2PDiameter_2P(mscanDataFiles)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab and Yurong Gao
%________________________________________________________________________________________________________________________
%
%   Purpose: Run the sub-functions necessary to analyze the TIFF stack's changes in vessel diameter.
%________________________________________________________________________________________________________________________
%
%   Inputs: List of all MScan files.
%
%   Outputs: Saves the updated files to the current directory.
%
%   Last Revised: January 29th, 2019
%________________________________________________________________________________________________________________________

%% Analyze every image (surface/penatrating/capillary)
for a = 1:size(mscanDataFiles, 1)
    load(mscanDataFiles(a,:), '-mat');
    if MScanData.notes.checklist.analyzeDiam == false
        if strcmp(MScanData.notes.movieType, 'MS') || strcmp(MScanData.notes.movieType, 'MD')
            [MScanData] = ExtractTiffAnalogData_2P(MScanData, [MScanData.notes.date '_' MScanData.notes.imageID]);
            MScanData.notes.checklist.analyzeDiam = true;
            save([MScanData.notes.animalID '_' MScanData.notes.date '_' MScanData.notes.imageID '_' MScanData.notes.vesselID '_MScanData'], 'MScanData')
            
        % elseif strcmp(MScanData.notes.movieType, 'MP')   % Not used in this analysis
            % [MScanData] = GetAreaPA(MScanData,[MScanData(1).ImageID '.TIF']);
            % MScanData.notes.checklist.analyzeDiam = true;         
            % save([MScanData.notes.animalID '_' MScanData.notes.date '_' MScanData.notes.imageID '_MScanData'], 'MScanData')

        % elseif strcmp(MScanData.notes.movieType, 'C')   % Not used in this analysis
            % [MScanData] = GetVelocityLineScan(MScanData,[MScanData(1).ImageID '.TIF']);
            % MScanData.notes.checklist.analyzeDiam = true;        
            % save([MScanData.notes.animalID '_' MScanData.notes.date '_' MScanData.notes.imageID '_MScanData'], 'MScanData')
        end
    else
        disp([mscanDataFiles(a,:) ' vessel diameter already analyzed. Continuing...']); disp(' ')
    end
end

end