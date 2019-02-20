function Analyze_2P_Diameter(directoryInfo)
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
    
    if strcmp(OrgData.Notes.movieType, 'MS') || strcmp(OrgData.Notes.movieType, 'MD')
        [OrgData] = Extract_Tiff_Analog_Data(OrgData, [OrgData.Notes.date '_' OrgData.Notes.imageID]);
        save([OrgData.Notes.animalID '_' OrgData.Notes.date '_' OrgData.Notes.imageID '_OrgData'], 'OrgData')
        
    elseif strcmp(OrgData.Notes.movieType, 'MP')
        [OrgData] = GetArea_PA_Tiff_auto02_new_soft_nowhisk(OrgData,[OrgData(1).ImageID '.TIF']);
        save([OrgData.Notes.animalID '_' OrgData.Notes.date '_' OrgData.Notes.imageID '_OrgData'], 'OrgData')

    elseif strcmp(OrgData.Notes.movieType, 'C')
        [OrgData] = GetVelocity_LineScan_new_soft_nowhisk_new(OrgData,[OrgData(1).ImageID '.TIF']);
        save([OrgData.Notes.animalID '_' OrgData.Notes.date '_' OrgData.Notes.imageID '_OrgData'], 'OrgData')
    end
end

% if strcmp(OrgData.Data.movieType, 'C') ~= 1
%     [OrgData] = post_process_pial_single3(OrgData);
% end

end