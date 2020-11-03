function [AnalysisResults] = TableS10_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table S10 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
% columnNames = AnalysisResults.NeuralHemoCoherence.columnNames;
columnNames = {'Rest','NREM','REM','Alert','Asleep','All'};
rowNames = {'Delta_C01_meanStD','Delta_C01_pVal','Theta_C01_meanStD','Theta_C01_pVal'...
    'Alpha_C01_meanStD','Alpha_C01_pVal','Beta_C01_meanStD','Beta_C01_pVal','Gamma_C01_meanStD','Gamma_C01_pVal'};
T(1,:) = cell2table(AnalysisResults.NeuralHemoCoherence.deltaBandPower.meanStD01);
T(2,:) = cell2table(AnalysisResults.NeuralHemoCoherence.deltaBandPower.p01);
T(3,:) = cell2table(AnalysisResults.NeuralHemoCoherence.thetaBandPower.meanStD01);
T(4,:) = cell2table(AnalysisResults.NeuralHemoCoherence.thetaBandPower.p01);
T(5,:) = cell2table(AnalysisResults.NeuralHemoCoherence.alphaBandPower.meanStD01);
T(6,:) = cell2table(AnalysisResults.NeuralHemoCoherence.alphaBandPower.p01);
T(7,:) = cell2table(AnalysisResults.NeuralHemoCoherence.betaBandPower.meanStD01);
T(8,:) = cell2table(AnalysisResults.NeuralHemoCoherence.betaBandPower.p01);
T(9,:) = cell2table(AnalysisResults.NeuralHemoCoherence.gammaBandPower.meanStD01);
T(10,:) = cell2table(AnalysisResults.NeuralHemoCoherence.gammaBandPower.p01);
T.Properties.RowNames = rowNames;
T.Properties.VariableNames = columnNames;
%% Table S10
summaryTable = figure('Name','TableS10');
sgtitle('Table S10 - Turner et al. 2020')
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryTable,[dirpath 'TableS10']);
end

end
