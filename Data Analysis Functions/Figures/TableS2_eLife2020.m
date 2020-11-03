function [AnalysisResults] = TableS2_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table S2 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
TwoP_animalIDs = AnalysisResults.ArterioleDurations.TwoP_animalIDs;
TwoP_baselineDiams = AnalysisResults.ArterioleDurations.TwoP_baselineDiams;
TwoP_totalTimeAwake = AnalysisResults.ArterioleDurations.TwoP_totalTimeAwake;
TwoP_totalTimeNREM = AnalysisResults.ArterioleDurations.TwoP_totalTimeNREM;
TwoP_totalTimeREM = AnalysisResults.ArterioleDurations.TwoP_totalTimeREM ;
TwoP_totalTimeMins = AnalysisResults.ArterioleDurations.TwoP_totalTimeMins;
TwoP_allTimeHours = AnalysisResults.ArterioleDurations.TwoP_allTimeHours;
TwoP_meanTimeHours = AnalysisResults.ArterioleDurations.TwoP_meanTimeHours;
TwoP_stdTimeHours = AnalysisResults.ArterioleDurations.TwoP_stdTimeHours;
%% Table S2
summaryTable = figure('Name','TableS2');
sgtitle('Table S2 - Turner et al. 2020')
variableNames = {'TotalTimeMins','BaseDiamUm','AwakeTimeMins','NREMTimeMins','REMTimeMins'};
T = table(TwoP_totalTimeMins,TwoP_baselineDiams,TwoP_totalTimeAwake,TwoP_totalTimeNREM,TwoP_totalTimeREM,'RowNames',TwoP_animalIDs,'VariableNames',variableNames);
T2 = sortrows(T,'RowNames');
uitable('Data',T2{:,:},'ColumnName',T2.Properties.VariableNames,'RowName',T2.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
uicontrol('Style','text','Position',[700,600,100,150],'String',{'Total Time (Hrs): ' num2str(TwoP_allTimeHours),'Mean time per animal (Hrs): ' num2str(TwoP_meanTimeHours) ' +/- ' num2str(TwoP_stdTimeHours)});
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryTable,[dirpath 'TableS2']);
end

end
