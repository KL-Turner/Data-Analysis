function [] = PupilSleepProbability_Pupil(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

resultsStruct = 'Results_SleepProbability';
load(resultsStruct);
diameterAllCatMeans = Results_SleepProbability.diameterCatMeans;
awakeProbPerc = Results_SleepProbability.awakeProbPerc;
nremProbPerc = Results_SleepProbability.nremProbPerc;
remProbPerc = Results_SleepProbability.remProbPerc;
edges = 100:10:3500;
summaryFigure = figure;
ax1 = subplot(1,1,1);
yyaxis right
h1 = histogram(diameterAllCatMeans,edges,'Normalization','probability','EdgeColor','k','FaceColor',colors('dark candy apple red'));
ylabel('Probability','rotation',-90,'VerticalAlignment','bottom')
yyaxis left
p1 = plot(edges,sgolayfilt(medfilt1(awakeProbPerc,10,'truncate'),3,17),'-','color','b','LineWidth',2);
hold on
p2 = plot(edges,sgolayfilt(medfilt1(nremProbPerc,10,'truncate'),3,17),'-','color','r','LineWidth',2);
p3 = plot(edges,sgolayfilt(medfilt1(remProbPerc,10,'truncate'),3,17),'-','color','k','LineWidth',2);
ylabel({'Arousal-state probability (%)'})
xlim([100,3500])
ylim([0,100])
legend([p1,p2,p3,h1],'Awake','NREM','REM','\DeltaArea','Location','NorthEast')
title('Diameter (Z Units) vs. arousal state probability')
xlabel('Diameter (Z Units)')
axis square
set(gca,'box','off')
set(gca,'TickLength',[0.03,0.03]);
set(h1,'facealpha',0.2);
ax1.TickLength = [0.03,0.03];
ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = colors_eLife2020('dark candy apple red');
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Diameter vs. Sleep Probability' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Diameter_SleepProbability']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Diameter_SleepProbability'])
end

end
