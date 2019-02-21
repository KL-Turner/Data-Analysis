function Analyze2PDiameter(directoryInfo)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%
% Originally written by Patrick J. Drew and Yurong Gao
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised: February 19th, 2019    
%________________________________________________________________________________________________________________________

%% Analyze every image (pial surface/dural/intracortical)
for f = 1:length(directoryInfo)
    load(directoryInfo(f).name);
    
    if strcmp(MScanData.Notes.movieType, 'MS') || strcmp(MScanData.Notes.movieType, 'MD')
        [MScanData] = ExtractTiffAnalogData(MScanData, [MScanData.Notes.date '_' MScanData.Notes.imageID]);
        save([MScanData.Notes.animalID '_' MScanData.Notes.date '_' MScanData.Notes.imageID '_MScanData'], 'MScanData')
        
    elseif strcmp(MScanData.Notes.movieType, 'MP')
        [MScanData] = GetArea_PA_Tiff_auto02_new_soft_nowhisk(MScanData,[MScanData(1).ImageID '.TIF']);
        save([MScanData.Notes.animalID '_' MScanData.Notes.date '_' MScanData.Notes.imageID '_MScanData'], 'MScanData')

    elseif strcmp(MScanData.Notes.movieType, 'C')
        [MScanData] = GetVelocity_LineScan_new_soft_nowhisk_new(MScanData,[MScanData(1).ImageID '.TIF']);
        save([MScanData.Notes.animalID '_' MScanData.Notes.date '_' MScanData.Notes.imageID '_MScanData'], 'MScanData')
    end
end

% if strcmp(MScanData.Data.movieType, 'C') ~= 1
%     [MScanData] = post_process_pial_single3(MScanData);
% end

end