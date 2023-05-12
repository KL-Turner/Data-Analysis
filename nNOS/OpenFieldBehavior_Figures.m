function [] = OpenFieldBehavior_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Hossain'];
cd(path)
plotStatistics_first5min(rootFolder,saveFigs,delim); % plot the stats for the first 5 minutes
vectormap_plot(rootFolder,saveFigs,delim) % plot the vector map
