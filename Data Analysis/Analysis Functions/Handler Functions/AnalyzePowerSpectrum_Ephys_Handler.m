function [] = AnalyzePowerSpectrum_Ephys_Handler(rootFolder,delim,runFromStart)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

% create or load results structure
if runFromStart == true
    Results_PowerSpec_Ephys = [];
elseif runFromStart == false
    cd([rootFolder delim 'Results_Turner']);
    % load existing results structure, if it exists
    if exist('Results_PowerSpec_Ephys.mat','file') == 2
        load('Results_PowerSpec_Ephys.mat','-mat')
    else
        Results_PowerSpec_Ephys.Naive = [];
        Results_PowerSpec_Ephys.SSP_SAP = [];
        Results_PowerSpec_Ephys.Blank_SAP = [];
    end
end
cd([rootFolder delim 'Data']);
groups = {'Naive','SSP_SAP','Blank_SAP'};
set = 'Ephys';
% determine waitbar length
waitBarLength = 0;
for aa = 1:length(groups)
    group = groups{1,aa};
    folderList = dir([group delim set]);
    folderList = folderList(~startsWith({folderList.name}, '.'));
    folderAnimalIDs = {folderList.name};
    waitBarLength = waitBarLength + length(folderAnimalIDs);
end
% run analysis for each animal in the group
cc = 1;
multiWaitbar('Analyzing power spectrum for IOS_Ephys',0,'Color','P');
for aa = 1:length(groups)
    group = groups{1,aa};
    folderList = dir([expGroup delim set]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        animalID = animalIDs{1,bb};
        if isfield(Results_PowerSpec_Ephys.(group),animalID) == false
            [Results_PowerSpec_Ephys] = AnalyzePowerSpectrum_Ephys(animalID,group,set,rootFolder,delim,Results_PowerSpec_Ephys);
        end
        multiWaitbar('Analyzing power spectrum for IOS_Ephys','Value',cc/waitBarLength); pause(0.5);
        cc = cc + 1;
    end
end