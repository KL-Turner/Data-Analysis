function [AnalysisResults] = Table1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table 1 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
% columnNames = AnalysisResults.PSD.columnNames;
columnNames = {'Rest','NREM','REM','Alert','Asleep','All'};
rowNames = {'Gamma_S01_meanStD','Gamma_S01_pVal','HbT_S01_meanStD','HbT_S01_pVal','TwoP_S01_meanStD','TwoP_S01_pVal'};
T(1,:) = cell2table(AnalysisResults.PSD.gammaBandPower.meanStD01);
T(2,:) = cell2table(AnalysisResults.PSD.gammaBandPower.p01);
T(3,:) = cell2table(AnalysisResults.PSD.CBV_HbT.meanStD01);
T(4,:) = cell2table(AnalysisResults.PSD.CBV_HbT.p01);
T(5,:) = cell2table(AnalysisResults.PSD.TwoP.meanStD01);
T(6,:) = cell2table(AnalysisResults.PSD.TwoP.p01);
T.Properties.RowNames = rowNames;
T.Properties.VariableNames = columnNames;
%% Table 1
summaryTable = figure('Name','Table1');
sgtitle('Table 1 - Turner et al. 2020')
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryTable,[dirpath 'Table1']);
end

end
