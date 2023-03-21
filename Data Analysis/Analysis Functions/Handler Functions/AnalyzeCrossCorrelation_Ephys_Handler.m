function [] = AnalyzeCrossCorrelation_Ephys_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the cross-correlation between neural activity and hemodynamics [HbT] (IOS)
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_CrossCorr_Ephys = [];
elseif runFromStart == false
    cd([rootFolder delim 'Results_Turner']);
    % load existing results structure, if it exists
    if exist('Results_CrossCorr_Ephys.mat','file') == 2
        load('Results_CrossCorr_Ephys.mat','-mat')
    else
        Results_CrossCorr_Ephys = [];
    end
end
cd([rootFolder delim 'Data']);
groups = {'Naive','SSP_SAP','Blank_SAP'};
set = 'IOS_Ephys';
% determine waitbar length
waitBarLength = 0;
for aa = 1:length(groups)
    folderList = dir([groups{1,aa} delim set]);
    folderList = folderList(~startsWith({folderList.name}, '.'));
    folderAnimalIDs = {folderList.name};
    waitBarLength = waitBarLength + length(folderAnimalIDs);
end
% run analysis for each animal in the group
cc = 1;
multiWaitbar('Analyzing cross correlation for IOS_Ephys',0,'Color','G');
for aa = 1:length(groups)
    folderList = dir([groups{1,aa} delim set]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(Results_CrossCorr_Ephys,(animalIDs{1,bb})) == false
            [Results_CrossCorr_Ephys] = AnalyzeCrossCorrelation_Ephys(animalIDs{1,bb},groups{1,aa},set,rootFolder,delim,Results_CrossCorr_Ephys);
        end
        multiWaitbar('Analyzing cross correlation for IOS_Ephys','Value',cc/waitBarLength); pause(0.5);
        cc = cc + 1;
    end
end
