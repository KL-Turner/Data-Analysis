function [] = AnalyzeBilateralCoherence_GCaMP_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the spectral coherence between bilateral hemodynamic [HbT] and neural signals (IOS)
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_BilatCoherGCaMP = [];
elseif runFromStart == false
    cd([rootFolder delim 'Results_Turner']);
    % load existing results structure, if it exists
    if exist('Results_BilatCoherGCaMP.mat','file') == 2
        load('Results_BilatCoherGCaMP.mat','-mat')
    else
        Results_BilatCoherGCaMP = [];
    end
end
cd([rootFolder delim 'Data']);
expGroups = {'SSP_SAP','Blank_SAP'};
setName = 'IOS_GCaMP7s';
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
multiWaitbar('Analyzing bilateral coherence for IOS_GCaMP',0,'Color','B');
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(Results_BilatCoherGCaMP,(animalIDs{1,bb})) == false
            [Results_BilatCoherGCaMP] = AnalyzeBilateralCoherence_GCaMP(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,delim,Results_BilatCoherGCaMP);
        end
        multiWaitbar('Analyzing bilateral coherence for IOS_GCaMP','Value',cc/waitBarLength); pause(0.5);
        cc = cc + 1;
    end
end
