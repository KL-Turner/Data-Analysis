function Analyze2PDiameter(directoryInfo)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab and Yurong Gao
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised: January 29th, 2019
%________________________________________________________________________________________________________________________

%% Analyze every image (pial surface/dural/intracortical)
for f = 1:length(directoryInfo)
    load(directoryInfo(f).name);
    if MScanData.Notes.checklist.analyzeDiam == false
        if strcmp(MScanData.Notes.movieType, 'MS') || strcmp(MScanData.Notes.movieType, 'MD')
            [MScanData] = ExtractTiffAnalogData(MScanData, [MScanData.Notes.date '_' MScanData.Notes.imageID]);
            MScanData.Notes.checklist.analizeDiam = true;
            save([MScanData.Notes.animalID '_' MScanData.Notes.date '_' MScanData.Notes.imageID '_MScanData'], 'MScanData')
            
        elseif strcmp(MScanData.Notes.movieType, 'MP')
            [MScanData] = GetArea_PA_Tiff_auto02_new_soft_nowhisk(MScanData,[MScanData(1).ImageID '.TIF']);
            MScanData.Notes.checklist.analizeDiam = true;         
            save([MScanData.Notes.animalID '_' MScanData.Notes.date '_' MScanData.Notes.imageID '_MScanData'], 'MScanData')

        elseif strcmp(MScanData.Notes.movieType, 'C')
            [MScanData] = GetVelocity_LineScan_new_soft_nowhisk_new(MScanData,[MScanData(1).ImageID '.TIF']);
            MScanData.Notes.checklist.analizeDiam = true;        
            save([MScanData.Notes.animalID '_' MScanData.Notes.date '_' MScanData.Notes.imageID '_MScanData'], 'MScanData')
        end
    end
end
% if strcmp(MScanData.Data.movieType, 'C') ~= 1
%     [MScanData] = post_process_pial_single3(MScanData);
% end

end