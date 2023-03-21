function [] = OpenFieldBehavior_Figures(rootFolder,~,delim)
path = [rootFolder delim 'Results_Hossain'];
cd(path)
plotStatistics() % plot the statistics
bodyparts_heatmap_plot() % plot the heatmap
cd(rootFolder)
