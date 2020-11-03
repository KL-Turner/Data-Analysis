function [AnalysisResults] = TableS5_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table S5 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
% columnNames = AnalysisResults.PSD.columnNames;
columnNames = {'Rest','NREM','REM','Alert','Asleep','All'};
rowNames = {'Delta_S01_meanStD','Delta_S01_pVal','Theta_S01_meanStD','Theta_S01_pVal'...
    'Alpha_S01_meanStD','Alpha_S01_pVal','Beta_S01_meanStD','Beta_S01_pVal'};
T(1,:) = cell2table(AnalysisResults.PSD.deltaBandPower.meanStD01);
T(2,:) = cell2table(AnalysisResults.PSD.deltaBandPower.p01);
T(3,:) = cell2table(AnalysisResults.PSD.thetaBandPower.meanStD01);
T(4,:) = cell2table(AnalysisResults.PSD.thetaBandPower.p01);
T(5,:) = cell2table(AnalysisResults.PSD.alphaBandPower.meanStD01);
T(6,:) = cell2table(AnalysisResults.PSD.alphaBandPower.p01);
T(7,:) = cell2table(AnalysisResults.PSD.betaBandPower.meanStD01);
T(8,:) = cell2table(AnalysisResults.PSD.betaBandPower.p01);
T.Properties.RowNames = rowNames;
T.Properties.VariableNames = columnNames;
%% Table S5
summaryTable = figure('Name','TableS5');
sgtitle('Table S5 - Turner et al. 2020')
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryTable,[dirpath 'TableS5']);
end

end
