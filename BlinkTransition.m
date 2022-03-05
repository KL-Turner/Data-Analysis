function [] = BlinkTransition(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generate figure panel 8 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% variables for loops
resultsStruct = 'Results_BlinkTransition';
load(resultsStruct);
animalIDs = fieldnames(Results_BlinkTransition);
%% take data from each animal corresponding to the CBV-gamma relationship
catAwakeMat = [];
catNremMat = [];
catRemMat = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    catAwakeMat = cat(1,catAwakeMat,Results_BlinkTransition.(animalID).awakeProbabilityMatrix);
    catNremMat = cat(1,catNremMat,Results_BlinkTransition.(animalID).nremProbabilityMatrix);
    catRemMat = cat(1,catRemMat,Results_BlinkTransition.(animalID).remProbabilityMatrix);
end
% average probability
awakeProbability = smooth(mean(catAwakeMat,1))*100;
nremProbability = smooth(mean(catNremMat,1))*100;
remProbability = smooth(mean(catRemMat,1))*100;
% figure
figure;
p1 = plot(awakeProbability);
hold on
p2 = plot(nremProbability);
p3 = plot(remProbability);
x1 = xline(7);
title('Arousal state probability adjacent to blinking')
xlabel('Time (sec)')
ylabel('Probability (%)')
legend([p1,p2,p3,x1],'Awake','NREM','REM','Blink')
xticks([1,3,5,7,9,11,13])
xticklabels({'-30','-20','-10','0','10','20','30'})
xlim([1,13])
ylim([0,100])
%% save figure(s)
% if saveFigs == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Pupil-Gamma Relationship' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     imwrite(h1Img,[dirpath 'Fig8d_Awake.png'])
%     imwrite(h2Img,[dirpath 'Fig8d_NREM.png'])
%     imwrite(h3Img,[dirpath 'Fig8d_REM.png'])
%     savefig(summaryFigure,[dirpath 'PupilGamma_Relationship']);
%     set(summaryFigure,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'PupilGamma_Relationship'])
% end

end
