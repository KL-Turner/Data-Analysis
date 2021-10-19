function [] = AnalyzeArousalTransitions_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the transitions between different arousal-states (IOS)
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_Transitions = [];
elseif runFromStart == false
    % load existing results structure, if it exists
    if exist('Results_Transitions.mat','file') == 2
        load('Results_Transitions.mat','-mat')
    else
        Results_Transitions = [];
    end
end
setName = 'IOS Set A';
expGroups = {'Naive','SSP-SAP','Blank-SAP'};
% determine waitbar length
waitBarLength = 0;
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name}, '.'));
    folderAnimalIDs = {folderList.name};
    waitBarLength = waitBarLength + length(folderAnimalIDs);
end
% run analysis for each animal in the group
cc = 1;
for aa = 1:length(expGroups)  
    folderList = dir([expGroups{1,aa} delim setName]); 
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    multiWaitbar('Analyzing arousal-state transitions',0,'Color','P'); pause(0.25);
    for bb = 1:length(animalIDs)
        if isfield(Results_Transitions,(animalIDs{1,bb})) == false
            [Results_Transitions] = AnalyzeArousalTransitions(animalIDs{1,bb},[expGroups{1,aa} delim 'IOS Set A'],rootFolder,delim,Results_Transitions);
        end
        multiWaitbar('Analyzing arousal-state transitions','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
multiWaitbar('CloseAll');

end
