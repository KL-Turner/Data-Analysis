function [] = CreateTrainingDataSet_SVM
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
%   Last Revised: July 26th, 2019
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs, 1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Adding sleep scoring parameters to ' procDataFileID '... (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
    [~, fileDate, ~] = GetFileInfo_IOS(procDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    load(procDataFileID)
    
    variableNames = {'maxLH_CBV', 'maxRH_CBV', 'maxLH_Delta', 'maxRH_Delta', 'maxLH_Theta', 'maxRH_Theta',...
        'numWhiskEvents', 'numForceEvents', 'numEMGEvents'};
%     T = table(

end

end
