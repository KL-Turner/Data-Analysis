function [] = BlinkResponses_Pupil(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

resultsStruct = 'Results_BlinkResponses';
load(resultsStruct);
animalIDs = fieldnames(Results_BlinkResponses);
data.HbT = []; data.cort = []; data.hip = []; data.EMG = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    data.pupilArea  = cat(1,data.pupilArea,Results_BlinkResponses.(animalID).pupilArea);
    data.HbT  = cat(1,data.HbT,Results_BlinkResponses.(animalID).LH_HbT,Results_BlinkResponses.(animalID).RH_HbT);
    data.cort = cat(3,data.cort,Results_BlinkResponses.(animalID).LH_cort,Results_BlinkResponses.(animalID).RH_cort);
    data.hip = cat(3,data.hip,Results_BlinkResponses.(animalID).hip);
    data.EMG = cat(1,data.EMG,Results_BlinkResponses.(animalID).EMG);
    T = Results_BlinkResponses.(animalID).T;
    F = Results_BlinkResponses.(animalID).F;
end
%
data.meanArea = mean(data.pupilArea,1);
data.stdArea = std(data.pupilArea,0,1)./sqrt(size(data.pupilArea,1));
data.meanHbT = mean(data.HbT,1);
data.stdHbT = std(data.HbT,0,1)./sqrt(size(data.HbT,1));
data.meanCort = mean(data.cort,3);
data.meanHip = mean(data.hip,3);
data.meanEMG = mean(data.EMG,1);
data.stdEMG = std(data.EMG,0,1)./sqrt(size(data.EMG,1));
% 
summaryFigure = figure;
subplot(5,1,1);
plot(timeVector,data.meanArea,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.meanArea + data.stdArea,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.meanArea - data.stdArea,'color',colors('smoky black'),'LineWidth',0.5)
title('Pupil area')
ylabel('\DeltaArea (pixels)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
subplot(5,1,2);
plot(timeVector,data.meanEMG,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.meanEMG + data.stdEMG,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.meanEMG - data.stdEMG,'color',colors('smoky black'),'LineWidth',0.5)
title('EMG')
ylabel('Power (a.u.)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
subplot(5,1,3);
plot(timeVector,data.meanHbT,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.meanHbT + data.stdHbT,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.meanHbT - data.stdHbT,'color',colors('smoky black'),'LineWidth',0.5)
title('HbT')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
subplot(5,1,4);
imagesc(T,F,data.meanCort)
title('Cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
axis xy
subplot(5,1,5);
imagesc(T,F,data.meanHip)
title('Hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square
axis xy
% % save figure(s)
% if saveFigs == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Stimulus-evoked Pupil Area' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure2,[dirpath 'Stimulus_PupilArea']);
%     set(summaryFigure2,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-bestfit',[dirpath 'Stimulus_PupilArea'])
% end

end
