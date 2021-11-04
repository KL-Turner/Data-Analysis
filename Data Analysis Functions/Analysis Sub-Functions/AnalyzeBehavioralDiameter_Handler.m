function [] = AnalyzeBehavioralDiameter_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the hemodynamic signal [Diameter] during different behavioral states (IOS)
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_BehavDiameter = [];
elseif runFromStart == false
    % load existing results structure, if it exists
    if exist('Results_BehavDiameter.mat','file') == 2
        load('Results_BehavDiameter.mat','-mat')
    else
        Results_BehavDiameter = [];
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
multiWaitbar('Analyzing behavioral pupil diameter',0,'Color','P'); pause(0.25);
for bb = 1:length(animalIDs)
    if isfield(Results_BehavDiameter,(animalIDs{1,bb})) == false
        [Results_BehavDiameter] = AnalyzeBehavioralDiameter(animalIDs{1,bb},rootFolder,delim,Results_BehavDiameter);
    end
    multiWaitbar('Analyzing behavioral pupil diameter','Value',aa/waitBarLength);
    aa = aa + 1;
end
multiWaitbar('CloseAll');

end
