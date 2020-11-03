function [AnalysisResults] = TableS11_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table S11 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
% columnNames = AnalysisResults.NeuralHemoCoherence.columnNames;
columnNames = {'Rest','NREM','REM','Alert','Asleep','All'};
rowNames = {'Delta_C001_meanStD','Delta_C001_pVal','Theta_C001_meanStD','Theta_C001_pVal'...
    'Alpha_C001_meanStD','Alpha_C001_pVal','Beta_C001_meanStD','Beta_C001_pVal','Gamma_C001_meanStD','Gamma_C001_pVal'};
T(1,:) = cell2table(AnalysisResults.NeuralHemoCoherence.deltaBandPower.meanStD001);
T(2,:) = cell2table(AnalysisResults.NeuralHemoCoherence.deltaBandPower.p001);
T(3,:) = cell2table(AnalysisResults.NeuralHemoCoherence.thetaBandPower.meanStD001);
T(4,:) = cell2table(AnalysisResults.NeuralHemoCoherence.thetaBandPower.p001);
T(5,:) = cell2table(AnalysisResults.NeuralHemoCoherence.alphaBandPower.meanStD001);
T(6,:) = cell2table(AnalysisResults.NeuralHemoCoherence.alphaBandPower.p001);
T(7,:) = cell2table(AnalysisResults.NeuralHemoCoherence.betaBandPower.meanStD001);
T(8,:) = cell2table(AnalysisResults.NeuralHemoCoherence.betaBandPower.p001);
T(9,:) = cell2table(AnalysisResults.NeuralHemoCoherence.gammaBandPower.meanStD001);
T(10,:) = cell2table(AnalysisResults.NeuralHemoCoherence.gammaBandPower.p001);
T.Properties.RowNames = rowNames;
T.Properties.VariableNames = columnNames;
%% Table S11
summaryTable = figure('Name','TableS11');
sgtitle('Table S11 - Turner et al. 2020')
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryTable,[dirpath 'TableS11']);
end

end
