function [] = AnalyzeVesselEvokedResponses_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the whisking-evoked and stimulus-evoked arteriole D/D responses (2PLSM)
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_VesselEvoked = [];
elseif runFromStart == false
    % load existing results structure, if it exists
    if exist('Results_VesselEvoked.mat','file') == 2
        load('Results_VesselEvoked.mat','-mat')
    else
        Results_VesselEvoked = [];
    end
end
setName = '2PLSM Set B';
expGroups = {'SSP-SAP','Blank-SAP'};
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
    multiWaitbar('Analyzing vessel stimulus and whisking evoked responses (A)',0,'Color','P'); pause(0.25);
    for bb = 1:length(animalIDs)
        if isfield(Results_VesselEvoked,(animalIDs{1,bb})) == false
            [Results_VesselEvoked] = AnalyzeVesselEvokedResponses(animalIDs{1,bb},[expGroups{1,aa} delim '2PLSM Set B'],rootFolder,delim,Results_VesselEvoked);
        end
        multiWaitbar('Analyzing vessel stimulus and whisking evoked responses (A)','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
multiWaitbar('CloseAll');

end
