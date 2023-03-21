function [] = AnalyzeNeuralHemoCoherence_GCaMP_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the spectral coherence between neural-hemodynamic [HbT] signals (IOS)
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_NeuralHemoCoherGCaMP = [];
elseif runFromStart == false
    cd([rootFolder delim 'Results_Turner']);
    % load existing results structure, if it exists
    if exist('Results_NeuralHemoCoherGCaMP.mat','file') == 2
        load('Results_NeuralHemoCoherGCaMP.mat','-mat')
    else
        Results_NeuralHemoCoherGCaMP = [];
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
multiWaitbar('Analyzing neural-hemo coherence for IOS_GCaMP',0,'Color','B');
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(Results_NeuralHemoCoherGCaMP,(animalIDs{1,bb})) == false
            [Results_NeuralHemoCoherGCaMP] = AnalyzeNeuralHemoCoherence_GCaMP(animalIDs{1,bb},[expGroups{1,aa} delim setName],rootFolder,delim,Results_NeuralHemoCoherGCaMP);
        end
        multiWaitbar('Analyzing neural-hemo coherence for IOS_GCaMP','Value',cc/waitBarLength); pause(0.5);
        cc = cc + 1;
    end
end
