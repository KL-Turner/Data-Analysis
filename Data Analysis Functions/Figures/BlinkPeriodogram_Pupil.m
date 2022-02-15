function [] = BlinkPeriodogram_Pupil(rootFolder,saveFigs,delim)
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
%% figures
summaryFigure = figure;
subplot(1,2,1)
semilogx(data.meanF1,data.meanS)
title('Power Spectrum')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
xlim([0.003,1]);
axis square
subplot(1,2,2)
semilogx(data.meanF2,data.meanPxx)
title('Lomb-Scargle Periodogram')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
xlim([0.003,1]);
axis square
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Blink Lomb-Scargle Periodogram' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Blink_Periodogram']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Blink_Periodogram'])
end

end
