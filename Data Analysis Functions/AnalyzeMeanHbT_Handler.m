function [] = AnalyzeMeanHbT_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the hemodynamic signal [HbT] during different arousal states (IOS)
%________________________________________________________________________________________________________________________

% create or load results structure
if strcmp(runFromStart,'y') == true
    AnalysisResults = [];
elseif strcmp(runFromStart,'n') == true
    % load existing results structure, if it exists
    if exist('AnalysisResults_MeanHbT.mat','file') == 2
        load('AnalysisResults_MeanHbT.mat','-mat')
    else
        AnalysisResults = [];
    end
end
% analyze the hemodynamic signal [HbT] during different arousal states (IOS)
setName = 'IOS Set A';
expGroups = {'C57BL6J','SSP-SAP','Blank-SAP'};
% determine waitbar length
waitBarLength = 0;
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name}, '.'));
    animalIDs = {folderList.name};
    waitBarLength = waitBarLength + length(animalIDs);
end
% run analysis for each animal in the group
cc = 1;
for aa = 1:length(expGroups)
    multiWaitbar('Analyzing behavioral hemodynamics',0,'Color','P'); pause(0.25);
    for bb = 1:length(animalIDs)
        if isfield(AnalysisResults,(animalIDs{1,bb})) == false
            [AnalysisResults] = AnalyzeMeanHbT(animalIDs{1,bb},[expGroups{1,aa} delim 'IOS Set A'],rootFolder,AnalysisResults);
        end
        multiWaitbar('Analyzing behavioral hemodynamics','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
multiWaitbar('CloseAll');
end
