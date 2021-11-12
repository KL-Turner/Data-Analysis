function [] = AnalyzePupilSleepModelAccuracy_Pupil_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_PupilSleepModel = [];
elseif runFromStart == false
    % load existing results structure, if it exists
    if exist('Results_PupilSleepModel.mat','file') == 2
        load('Results_PupilSleepModel.mat','-mat')
    else
        Results_PupilSleepModel = [];
    end
end
folderList = dir('Data');
folderList = folderList(~startsWith({folderList.name}, '.'));
animalIDs = {folderList.name};
% run analysis for each animal in the group
if isempty(Results_PupilSleepModel) == true
    AnalzyePupilSleepModelAccuracy2_Pupil(animalIDs,rootFolder,delim);
end

end
