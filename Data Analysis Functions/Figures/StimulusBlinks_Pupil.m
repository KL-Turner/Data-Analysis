function [] = StimulusBlinks_Pupil(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% variables for loops
resultsStruct = 'Results_StimulusBlinks';
load(resultsStruct);
animalIDs = fieldnames(Results_StimulusBlinks);
%% pre-allocate data structure
data.stimPerc = []; data.binProb = []; data.indBinProb = []; data.duration = [];
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    data.stimPerc = cat(1,data.stimPerc,Results_StimulusBlinks.(animalID).stimPercentage);
    data.duration = cat(1,data.duration,Results_StimulusBlinks.(animalID).stimPercentageDuration);
    data.binProb = cat(1,data.binProb,Results_StimulusBlinks.(animalID).binProbability);
    data.indBinProb = cat(1,data.indBinProb,Results_StimulusBlinks.(animalID).indBinProbability);
end
data.meanStimPerc = mean(data.stimPerc,1);
data.stdStimPerc = std(data.stimPerc,0,1);
data.meanDuration = mean(data.duration,1);
data.meanBinProb = mean(data.binProb,1);
data.stdBinProb = std(data.binProb,0,1);
data.meanIndBinProb = mean(data.indBinProb,1);
data.stdIndBinProb = std(data.indBinProb,0,1);
%% figures
summaryFigure = figure;
subplot(1,2,1)
scatter(data.duration,data.stimPerc,75,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold on
e1 = errorbar(data.meanDuration,data.meanStimPerc,data.stdStimPerc,'d','MarkerEdgeColor','k','MarkerFaceColor','g');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
title('Probability of blinking post-whisker stimulus')
ylabel('Probability (%)')
set(gca,'box','off')
% xlim([0.5,1.5]);
axis square
subplot(1,2,2)
plot(0.5:0.5:5,data.meanBinProb)
hold on; 
plot(0.5:0.5:5,data.meanBinProb + data.stdBinProb,'r')
plot(0.5:0.5:5,data.meanBinProb - data.stdBinProb,'b')
title('Probability of first blink post-stimulus')
ylabel('Time (s)')
xlabel('Probability (%)')
set(gca,'box','off')
xlim([0,5]);
axis square
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Post-Stimulus Blink Probability' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'PostStimulusBlinkProbability']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'PostStimulusBlinkProbability'])
end

end
