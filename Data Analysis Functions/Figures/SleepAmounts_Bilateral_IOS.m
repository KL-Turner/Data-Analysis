function [] = SleepAmounts_Bilateral_IOS(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

colorRfcAwake = [(0/256),(64/256),(64/256)];
colorRfcNREM = [(0/256),(174/256),(239/256)];
colorRfcREM = [(190/256),(30/256),(45/256)];
%% set-up and process data
resultsStruct = 'Results_ArousalStateProb';
load(resultsStruct);
expGroups = {'Naive','SSP-SAP','Blank-SAP'};
setName = 'IOS Set A';
animalIDs.all = {};
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs.all = horzcat(animalIDs.all,{folderList.name});
    animalIDs.(strrep(expGroups{1,aa},'-','_')) = {folderList.name};
end
treatments = {'Naive','SSP_SAP','Blank_SAP'};
%% Pearson's correlations during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.Naive = []; data.SSP_SAP = []; data.Blank_SAP = [];
for aa = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,aa};
    % recognize treatment based on animal group
    if ismember(animalID,animalIDs.Naive) == true
        treatment = 'Naive';
    elseif ismember(animalID,animalIDs.SSP_SAP) == true
        treatment = 'SSP_SAP';
    elseif ismember(animalID,animalIDs.Blank_SAP) == true
        treatment = 'Blank_SAP';
    end
    if isfield(data.(treatment),'awakePercent') == false
        data.(treatment).awakePercent = [];
        data.(treatment).nremPercent = [];
        data.(treatment).remPercent = [];
    end
    % concatenate mean R and the animalID/behavior
    data.(treatment).awakePercent = cat(1,data.(treatment).awakePercent,Results_ArousalStateProb.(animalID).awakePercent);
    data.(treatment).nremPercent = cat(1,data.(treatment).nremPercent,Results_ArousalStateProb.(animalID).nremPercent);
    data.(treatment).remPercent = cat(1,data.(treatment).remPercent,Results_ArousalStateProb.(animalID).remPercent);
end
%% take mean/STD of R
for qq = 1:length(treatments)
    treatment = treatments{1,qq};
    data.(treatment).meanAwakePercent = mean(data.(treatment).awakePercent,1);
    data.(treatment).stdAwakePercent = std(data.(treatment).awakePercent,0,1);
    data.(treatment).meanNremPercent = mean(data.(treatment).nremPercent,1);
    data.(treatment).stdNremPercent = std(data.(treatment).nremPercent,0,1);
    data.(treatment).meanRemPercent = mean(data.(treatment).remPercent,1);
    data.(treatment).stdRemPercent = std(data.(treatment).remPercent,0,1);
end
meanPercsNaive = horzcat(data.Naive.meanAwakePercent,data.Naive.meanNremPercent,data.Naive.meanRemPercent);
meanPercsSAP = horzcat(data.SSP_SAP.meanAwakePercent,data.SSP_SAP.meanNremPercent,data.SSP_SAP.meanRemPercent);
meanPercsBlank = horzcat(data.Blank_SAP.meanAwakePercent,data.Blank_SAP.meanNremPercent,data.Blank_SAP.meanRemPercent);
%% arousal state percentages
summaryFigure = figure;
%% Naive
subplot(1,3,1);
p1 = pie(meanPercsNaive);
pText = findobj(p1,'Type','text');
percentValues = get(pText,'String');
txt = {'rfc-Awake: ';'rfc-NREM: ';'rfc-REM: '};
combinedtxt = strcat(txt,percentValues);
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
title({'Naive'})
%% SSP-SAP
subplot(1,3,2);
p1 = pie(meanPercsSAP);
pText = findobj(p1,'Type','text');
percentValues = get(pText,'String');
txt = {'rfc-Awake: ';'rfc-NREM: ';'rfc-REM: '};
combinedtxt = strcat(txt,percentValues);
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
title({'SSP-SAP'})
%% Blank-SAP
subplot(1,3,3);
p1 = pie(meanPercsBlank);
pText = findobj(p1,'Type','text');
percentValues = get(pText,'String');
txt = {'rfc-Awake: ';'rfc-NREM: ';'rfc-REM: '};
combinedtxt = strcat(txt,percentValues);
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
title({'Blank-SAP'})
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Arousal State Probability - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'ArousalProbability']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'ArousalProbability'])
end

end
