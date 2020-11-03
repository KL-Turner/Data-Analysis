function [AnalysisResults] = Table2_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table 2 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
% columnNames = AnalysisResults.PSD.columnNames;
columnNames = {'Rest','NREM','REM','Alert','Asleep','All'};
rowNames = {'Gamma_S001_meanStD','Gamma_S001_pVal','HbT_S001_meanStD','HbT_S001_pVal','TwoP_S001_meanStD','TwoP_S001_pVal'};
T(1,:) = cell2table(AnalysisResults.PSD.gammaBandPower.meanStD001);
T(2,:) = cell2table(AnalysisResults.PSD.gammaBandPower.p001);
T(3,:) = cell2table(AnalysisResults.PSD.CBV_HbT.meanStD001);
T(4,:) = cell2table(AnalysisResults.PSD.CBV_HbT.p001);
T(5,:) = cell2table(AnalysisResults.PSD.TwoP.meanStD001);
T(6,:) = cell2table(AnalysisResults.PSD.TwoP.p001);
T.Properties.RowNames = rowNames;
T.Properties.VariableNames = columnNames;
%% Table 2
summaryTable = figure('Name','Table2');
sgtitle('Table 2 - Turner et al. 2020')
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryTable,[dirpath 'Table2']);
end

end
