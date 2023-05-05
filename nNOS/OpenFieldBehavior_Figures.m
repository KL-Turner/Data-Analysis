function [] = OpenFieldBehavior_Figures(rootFolder,~,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Hossain'];
cd(path)
plotStatistics()
bodyparts_heatmap_plot()