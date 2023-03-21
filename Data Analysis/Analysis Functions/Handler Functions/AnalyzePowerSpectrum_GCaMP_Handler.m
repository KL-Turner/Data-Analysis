function [] = AnalyzePowerSpectrum_GCaMP_Handler(rootFolder,delim,runFromStart)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

% create or load results structure
if runFromStart == true
    Results_PowerSpec_GCaMP = [];
elseif runFromStart == false
    cd([rootFolder delim 'Results_Turner']);
    % load existing results structure, if it exists
    if exist('Results_PowerSpec_GCaMP.mat','file') == 2
        load('Results_PowerSpec_GCaMP.mat','-mat')
    else
        Results_PowerSpec_GCaMP.Blank_SAP = [];
        Results_PowerSpec_GCaMP.SSP_SAP = [];
    end
end
cd([rootFolder delim 'Data']);
expGroups = {'Blank_SAP','SSP_SAP'};
setName = 'GCaMP';
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
multiWaitbar('Analyzing power spectrum for IOS_GCaMP',0,'Color','B');
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(Results_PowerSpec_GCaMP,(animalIDs{1,bb})) == false
            [Results_PowerSpec_GCaMP] = AnalyzePowerSpectrum_GCaMP(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,delim,Results_PowerSpec_GCaMP);
        end
        multiWaitbar('Analyzing power spectrum for IOS_GCaMP','Value',cc/waitBarLength); pause(0.5);
        cc = cc + 1;
    end
end