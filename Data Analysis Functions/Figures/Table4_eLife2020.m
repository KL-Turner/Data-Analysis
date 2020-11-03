function [AnalysisResults] = Table4_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table 4 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
% columnNames = AnalysisResults.Coherr.columnNames;
columnNames = {'Rest','NREM','REM','Alert','Asleep','All'};
rowNames = {'Gamma_C001_meanStD','Gamma_C001_pVal','HbT_C001_meanStD','HbT_C001_pVal'};
T(1,:) = cell2table(AnalysisResults.Coherr.gammaBandPower.meanStD001);
T(2,:) = cell2table(AnalysisResults.Coherr.gammaBandPower.p001);
T(3,:) = cell2table(AnalysisResults.Coherr.CBV_HbT.meanStD001);
T(4,:) = cell2table(AnalysisResults.Coherr.CBV_HbT.p001);
T.Properties.RowNames = rowNames;
T.Properties.VariableNames = columnNames;
%% Table 4
summaryTable = figure('Name','Table4');
sgtitle('Table 4 - Turner et al. 2020')
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryTable,[dirpath 'Table4']);
end

end
