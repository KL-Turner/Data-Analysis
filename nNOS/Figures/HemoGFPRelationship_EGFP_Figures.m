function [] = HemoGFPRelationship_EGFP_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_GFP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'EGFP'};
hemispheres = {'LH','RH'};
data.blue = [];
data.green = [];
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_GFP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            data.blue = cat(2,data.blue,Results_GFP.(group).(animalID).(hemisphere).blue);
            data.green = cat(2,data.green,Results_GFP.(group).(animalID).(hemisphere).green);
        end
    end
end
% figure
summaryFigure = figure;
histogram2(data.blue,data.green,'DisplayStyle','tile','ShowEmptyBins','on')
axis square
set(gca,'box','off')
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures' delim 'Intrinsic Signals' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'IntrinsicSignals_GFP']);
end