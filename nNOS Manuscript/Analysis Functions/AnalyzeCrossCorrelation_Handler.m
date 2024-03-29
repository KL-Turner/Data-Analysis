function [] = AnalyzeCrossCorrelation_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the cross-correlation between neural activity and hemodynamics [HbT] (IOS)
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_CrossCorrEphys = [];
elseif runFromStart == false
    cd([rootFolder delim 'Results_Turner']);
    % load existing results structure, if it exists
    if exist('Results_CrossCorrEphys.mat','file') == 2
        load('Results_CrossCorrEphys.mat','-mat')
    else
        Results_CrossCorrEphys = [];
    end
end
cd([rootFolder delim 'Data']);
expGroups = {'Naive','SSP_SAP','Blank_SAP'};
setName = 'IOS_Ephys';
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
multiWaitbar('Analyzing cross correlation for IOS_Ephys',0,'Color','G');
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(Results_CrossCorrEphys,(animalIDs{1,bb})) == false
            [Results_CrossCorrEphys] = AnalyzeCrossCorrelation_Ephys(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,delim,Results_CrossCorrEphys);
        end
        multiWaitbar('Analyzing cross correlation for IOS_Ephys','Value',cc/waitBarLength); pause(0.5);
        cc = cc + 1;
    end
end
