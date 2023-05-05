function [] = ArousalStateProb_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_ArousalStateProb_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Naive','Blank_SAP','SSP_SAP'};
dataTypes = {'awakePercent','nremPercent','remPercent'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_ArousalStateProb_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        data.(group).dummCheck = 1;
        data.(group).group = {};
        data.(group).animalID = {};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            if isfield(data.(group),dataType) == false
                data.(group).(dataType) = [];
            end
            data.(group).(dataType) = cat(1,data.(group).(dataType),Results_ArousalStateProb_Ephys.(group).(animalID).(dataType));
        end
        data.(group).group = cat(1,data.(group).group,group);
        data.(group).animalID = cat(1,data.(group).animalID,animalID);
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.(group).(['mean_' dataType]) = mean(data.(group).(dataType),1);
        data.(group).(['std_' dataType]) = std(data.(group).(dataType),0,1);
    end
    data.(group).meanPercs = cat(1,data.(group).mean_awakePercent,data.(group).mean_nremPercent,data.(group).mean_remPercent);
end
% figure
summaryFigure = figure;
sgtitle('Arousal state probability [Ephys]')
for aa = 1:length(groups)
    group = groups{1,aa};
    subplot(1,3,aa);
    p1 = pie(data.(group).meanPercs);
    pText = findobj(p1,'Type','text');
    percentValues = get(pText,'String');
    txt = {'Awake: ';'NREM: ';'REM: '};
    combinedtxt = strcat(txt,percentValues);
    pText(1).String = combinedtxt(1);
    pText(2).String = combinedtxt(2);
    pText(3).String = combinedtxt(3);
    title(strrep(group,'_',' '))
end
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures' delim 'Behavior' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'ArousalStateProb_Ephys']);
end
% statistics - unpaired ttest
[awakeStats1.h,awakeStats1.p,awakeStats1.ci,awakeStats1.stats] = ttest2(data.Naive.awakePercent,data.Blank_SAP.awakePercent);
[awakeStats2.h,awakeStats2.p,awakeStats2.ci,awakeStats2.stats] = ttest2(data.Blank_SAP.awakePercent,data.SSP_SAP.awakePercent);
[nremStats1.h,nremStats1.p,nremStats1.ci,nremStats1.stats] = ttest2(data.Naive.nremPercent,data.Blank_SAP.nremPercent);
[nremStats2.h,nremStats2.p,nremStats2.ci,nremStats2.stats] = ttest2(data.Blank_SAP.nremPercent,data.SSP_SAP.nremPercent);
[remStats1.h,remStats1.p,remStats1.ci,remStats1.stats] = ttest2(data.Naive.nremPercent,data.Blank_SAP.nremPercent);
[remStats2.h,remStats2.p,remStats2.ci,remStats2.stats] = ttest2(data.Blank_SAP.remPercent,data.SSP_SAP.remPercent);
% statistical diary
if saveFigs == true
    % statistical diary
    diaryFile = [dirpath 'StimEvoked_Ephys_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    disp('======================================================================================================================')
    disp('ttest2 statistics:')
    disp('======================================================================================================================')
    disp(['Naive vs. Blank (awake) p < ' num2str(awakeStats1.p)]); disp(' ')
    disp(['Blank vs. SAP (awake) p < ' num2str(awakeStats2.p)]); disp(' ')
    disp(['Naive vs. Blank (nrem) p < ' num2str(nremStats1.p)]); disp(' ')
    disp(['Blank vs. SAP (nrem) p < ' num2str(nremStats2.p)]); disp(' ')
    disp(['Naive vs. Blank (rem) p < ' num2str(remStats1.p)]); disp(' ')
    disp(['Blank vs. SAP (rem) p < ' num2str(remStats2.p)]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
end