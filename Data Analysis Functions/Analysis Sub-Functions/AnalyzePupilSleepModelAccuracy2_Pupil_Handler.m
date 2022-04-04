function [] = AnalyzePupilSleepModelAccuracy2_Pupil_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: 
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_SleepModel = [];
elseif runFromStart == false
    % load existing results structure, if it exists
    if exist('Results_SleepModel.mat','file') == 2
        load('Results_SleepModel.mat','-mat')
    else
        Results_SleepModel = [];
    end
end
% determine waitbar length
waitBarLength = 0;
folderList = dir('Data');
folderList = folderList(~startsWith({folderList.name}, '.'));
animalIDs = {folderList.name};
waitBarLength = waitBarLength + length(animalIDs);
% run analysis for each animal in the group
aa = 1;
multiWaitbar('Analyzing pupil sleep model accuracy',0,'Color','P'); pause(0.25);
for bb = 1:length(animalIDs)
    if isfield(Results_SleepModel,(animalIDs{1,bb})) == false
        [Results_SleepModel] = AnalzyePupilSleepModelAccuracy_Pupil(animalIDs{1,bb},rootFolder,delim,Results_SleepModel);
    end
    multiWaitbar('Analyzing pupil sleep model accuracy','Value',aa/waitBarLength);
    aa = aa + 1;
end

end
