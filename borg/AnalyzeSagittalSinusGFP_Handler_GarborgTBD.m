function [] = AnalyzeSagittalSinusGFP_Handler_GarborgTBD(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Handler function for AnalyzeEvokedResponses.mat
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_ImageStim = [];
elseif runFromStart == false
    % load existing results structure, if it exists
    if exist('Results_ImageStim.mat','file') == 2
        load('Results_ImageStim.mat','-mat')
    else
        Results_ImageStim = [];
    end
end
cd(rootFolder)
% determine waitbar length
waitBarLength = 0;
folderList = dir('Data');
folderList = folderList(~startsWith({folderList.name},'.'));
animalIDs = {folderList.name};
waitBarLength = waitBarLength + length(animalIDs);
% run analysis for each animal in the group
aa = 1;
multiWaitbar('Analyzing image-wise sagittal sinus fluorescence',0,'Color','P'); pause(0.25);
for bb = 1:length(animalIDs)
    if isfield(Results_ImageStim,(animalIDs{1,bb})) == false
        [Results_ImageStim] = AnalyzeSagittalSinusGFP_GarborgTBD(animalIDs{1,bb},rootFolder,delim,Results_ImageStim);
    end
    multiWaitbar('Analyzing image-wise sagittal sinus fluorescence','Value',aa/waitBarLength);
    aa = aa + 1;
end

end
