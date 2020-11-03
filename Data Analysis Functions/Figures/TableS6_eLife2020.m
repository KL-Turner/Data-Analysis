function [AnalysisResults] = TableS6_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table S6 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
% columnNames = AnalysisResults.PSD.columnNames;
columnNames = {'Rest','NREM','REM','Alert','Asleep','All'};
rowNames = {'Delta_S001_meanStD','Delta_S001_pVal','Theta_S001_meanStD','Theta_S001_pVal'...
    'Alpha_S001_meanStD','Alpha_S001_pVal','Beta_S001_meanStD','Beta_S001_pVal'};
T(1,:) = cell2table(AnalysisResults.PSD.deltaBandPower.meanStD001);
T(2,:) = cell2table(AnalysisResults.PSD.deltaBandPower.p001);
T(3,:) = cell2table(AnalysisResults.PSD.thetaBandPower.meanStD001);
T(4,:) = cell2table(AnalysisResults.PSD.thetaBandPower.p001);
T(5,:) = cell2table(AnalysisResults.PSD.alphaBandPower.meanStD001);
T(6,:) = cell2table(AnalysisResults.PSD.alphaBandPower.p001);
T(7,:) = cell2table(AnalysisResults.PSD.betaBandPower.meanStD001);
T(8,:) = cell2table(AnalysisResults.PSD.betaBandPower.p001);
T.Properties.RowNames = rowNames;
T.Properties.VariableNames = columnNames;
%% Table S6
summaryTable = figure('Name','TableS6');
sgtitle('Table S6 - Turner et al. 2020')
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryTable,[dirpath 'TableS6']);
end

end
