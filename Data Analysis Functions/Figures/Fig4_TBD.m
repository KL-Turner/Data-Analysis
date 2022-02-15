function [] = Fig4_TBD(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% variables for loops
resultsStruct = 'Results_BlinkPeriodogram';
load(resultsStruct);
animalIDs = fieldnames(Results_BlinkPeriodogram); %#ok<NODEF>
animalIDs = animalIDs(1:end - 1);
%% pre-allocate data structure
data.f1 = []; data.S = []; data.blinkArray = [];
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    data.blinkArray = cat(2,data.blinkArray,Results_BlinkPeriodogram.(animalID).blinkArray);
    data.S = cat(2,data.S,Results_BlinkPeriodogram.(animalID).S);
    data.f1 = cat(1,data.f1,Results_BlinkPeriodogram.(animalID).f);
end
data.meanS = mean(data.S,2);
data.meanF1 = mean(data.f1,1);
if isfield(Results_BlinkPeriodogram,'results') == false 
    %% mean/std
    [data.pxx,data.f2] = plomb(data.blinkArray,2);
    bb = 1; pxx2 = [];
    for aa = 1:length(animalIDs)
        animalID = animalIDs{aa,1};
        avgLen = size(Results_BlinkPeriodogram.(animalID).blinkArray,2);
        if bb == 1
            pxx2(:,aa) = mean(data.pxx(:,bb:bb + avgLen),2);
        else
            pxx2(:,aa) = mean(data.pxx(:,bb + 1:bb + avgLen - 1),2);
        end
        bb = bb + avgLen;
    end
    Results_BlinkPeriodogram.results.f = data.f2;
    Results_BlinkPeriodogram.results.pxx = pxx2;
    save('Results_BlinkPeriodogram.mat','Results_BlinkPeriodogram')
else
    data.f2 = Results_BlinkPeriodogram.results.f;
    data.pxx = Results_BlinkPeriodogram.results.pxx;
end
data.meanPxx = mean(data.pxx,2,'omitnan');
data.meanF2 = data.f2;
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
sgtitle('Figure 4')
subplot(2,2,1)
semilogx(data.meanF1,data.meanS)
title('Power Spectrum')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
xlim([0.003,1]);
axis square
subplot(2,2,2)
semilogx(data.meanF2,data.meanPxx)
title('Lomb-Scargle Periodogram')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
xlim([0.003,1]);
axis square
subplot(2,2,3)
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
subplot(2,2,4)
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
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig4_TBD']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig4_TBD'])
end

end
