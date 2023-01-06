function [] = AnalyzeBilateralCoherence_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the spectral coherence between bilateral hemodynamic [HbT] and neural signals (IOS)
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_BilatCoher = [];
elseif runFromStart == false
    % load existing results structure, if it exists
    if exist('Results_BilatCoher.mat','file') == 2
        load('Results_BilatCoher.mat','-mat')
    else
        Results_BilatCoher = [];
    end
end
expGroups = {'SSP_SAP','Blank_SAP','NPY_Gi_DREADDs','NPY_Gq_DREADDs'};
sets = {'IOS_Ephys','IOS_GCaMP7s';'IOS_Ephys','IOS_GCaMP7s';'CNO','DMSO';'CNO','DMSO'};
% determine waitbar length
waitBarLength = 0;
for aa = 1:length(expGroups)
    setNames = sets(aa,:);
    for bb = 1:length(setNames)
        folderList = dir([expGroups{1,aa} delim setNames{1,bb}]);
        folderList = folderList(~startsWith({folderList.name}, '.'));
        folderAnimalIDs = {folderList.name};
        waitBarLength = waitBarLength + length(folderAnimalIDs);
    end
end
% run analysis for each animal in the group
dd = 1;
multiWaitbar('Analyzing bilateral coherence',0,'Color','P');
for aa = 1:length(expGroups)
    setNames = sets(aa,:);
    for bb = 1:length(setNames)
        folderList = dir([expGroups{1,aa} delim setNames{1,bb}]);
        folderList = folderList(~startsWith({folderList.name},'.'));
        animalIDs = {folderList.name};
        for cc = 1:length(animalIDs)
            if isfield(Results_BilatCoher,(animalIDs{1,cc})) == false
                [Results_BilatCoher] = AnalyzeBilateralCoherence(animalIDs{1,cc},[expGroups{1,aa} delim setNames{1,bb}],rootFolder,delim,Results_BilatCoher);
            end
            multiWaitbar('Analyzing bilateral coherence','Value',dd/waitBarLength); pause(0.5);
            dd = dd + 1;
        end
    end
end