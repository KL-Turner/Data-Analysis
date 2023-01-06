function [] = AnalyzePowerSpectrum_LFP_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the LFP power spectrum (IOS)
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_PowerSpecLFP = [];
elseif runFromStart == false
    % load existing results structure, if it exists
    if exist('Results_PowerSpecLFP.mat','file') == 2
        load('Results_PowerSpecLFP.mat','-mat')
    else
        Results_PowerSpecLFP = [];
    end
end
expGroups = {'SSP_SAP','Blank_SAP'};
sets = {'IOS_Ephys';'IOS_Ephys'};
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
multiWaitbar('Analyzing LFP power spectrum',0,'Color','P');
for aa = 1:length(expGroups)
    setNames = sets(aa,:);
    for bb = 1:length(setNames)
        folderList = dir([expGroups{1,aa} delim setNames{1,bb}]);
        folderList = folderList(~startsWith({folderList.name},'.'));
        animalIDs = {folderList.name};
        for cc = 1:length(animalIDs)
            if isfield(Results_PowerSpecLFP,(animalIDs{1,cc})) == false
                [Results_PowerSpecLFP] = AnalyzePowerSpectrum_LFP(animalIDs{1,cc},[expGroups{1,aa} delim setNames{1,bb}],rootFolder,delim,Results_PowerSpecLFP);
            end
            multiWaitbar('Analyzing LFP power spectrum','Value',dd/waitBarLength); pause(0.5);
            dd = dd + 1;
        end
    end
end